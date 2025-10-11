#' Stan Utility Functions
#'
#' Functions for preparing data, fitting Stan models, and extracting results
#' for radiocarbon dating contemporaneity analysis.
#'
#' @author Project Team
#' @date 2025-10-11

library(cmdstanr)
library(posterior)
library(loo)

#' Prepare Stan data from radiocarbon dataset
#'
#' Converts calibrated radiocarbon dates and deposit assignments into
#' Stan-compatible data format.
#'
#' @param c14_data Data frame with radiocarbon dates
#' @param cal_result Calibrated dates from calibrate_dates()
#' @param prior_params List with prior parameters (optional)
#' @return List formatted for Stan models
#' @export
prepare_stan_data <- function(c14_data, cal_result, prior_params = NULL) {

  # Validate inputs
  validate_c14_data(c14_data)

  if (!inherits(cal_result, "calibrated_c14")) {
    stop("cal_result must be output from calibrate_dates()")
  }

  # Get calibration data
  calib_data <- prepare_stan_calibration_data(cal_result, method = "mixture", n_components = 3)

  # Count deposits
  n_deposits <- length(unique(c14_data$deposit))
  deposit_ids <- as.integer(factor(c14_data$deposit))

  # Set default priors if not provided
  if (is.null(prior_params)) {
    # Use data-driven priors
    age_range <- range(c14_data$age)

    # Approximate calendar range (rough conversion)
    cal_range_approx <- age_range

    prior_params <- list(
      prior_start_mean = max(cal_range_approx) + 200,
      prior_start_sd = 500,
      prior_end_mean = min(cal_range_approx) - 200,
      prior_end_sd = 500,
      prior_duration_mean = 300,
      prior_duration_sd = 200,
      prior_gap_mean = 100,
      prior_gap_sd = 100
    )
  }

  # Combine into Stan data list
  stan_data <- c(
    list(
      N = nrow(c14_data),
      n_deposits = n_deposits,
      c14_age = c14_data$age,
      c14_error = c14_data$error,
      deposit_id = deposit_ids
    ),
    calib_data,
    prior_params
  )

  return(stan_data)
}

#' Compile Stan model
#'
#' Compiles a Stan model file using cmdstanr.
#'
#' @param model_file Path to .stan file
#' @param force_recompile Logical, force recompilation (default FALSE)
#' @param cpp_options List of C++ compiler options (optional)
#' @return Compiled Stan model object
#' @export
compile_stan_model <- function(model_file,
                                 force_recompile = FALSE,
                                 cpp_options = list()) {

  if (!file.exists(model_file)) {
    stop("Stan model file not found: ", model_file)
  }

  message("Compiling Stan model: ", basename(model_file))

  model <- cmdstan_model(
    stan_file = model_file,
    compile = TRUE,
    force_recompile = force_recompile,
    cpp_options = cpp_options
  )

  message("Compilation successful")

  return(model)
}

#' Fit Stan model
#'
#' Fits a Stan model to radiocarbon data with automated diagnostics.
#'
#' @param model Compiled Stan model (from compile_stan_model)
#' @param stan_data Data list for Stan
#' @param chains Number of MCMC chains (default 4)
#' @param parallel_chains Number of chains to run in parallel (default 4)
#' @param iter_warmup Number of warmup iterations (default 1000)
#' @param iter_sampling Number of sampling iterations (default 1000)
#' @param adapt_delta Target acceptance rate (default 0.95)
#' @param max_treedepth Maximum tree depth (default 12)
#' @param seed Random seed for reproducibility
#' @param init Initial values (optional)
#' @param verbose Print sampling progress (default TRUE)
#' @return Fitted Stan model object with diagnostics
#' @export
fit_stan_model <- function(model,
                             stan_data,
                             chains = 4,
                             parallel_chains = 4,
                             iter_warmup = 1000,
                             iter_sampling = 1000,
                             adapt_delta = 0.95,
                             max_treedepth = 12,
                             seed = NULL,
                             init = NULL,
                             verbose = TRUE) {

  if (verbose) {
    message("Fitting Stan model...")
    message("  Chains: ", chains)
    message("  Iterations: ", iter_warmup, " warmup + ", iter_sampling, " sampling")
  }

  # Fit model
  fit <- model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    init = init,
    refresh = ifelse(verbose, 200, 0)
  )

  # Check diagnostics
  diagnostics <- check_stan_diagnostics(fit, verbose = verbose)

  # Store diagnostics as an attribute (not in locked environment)
  attr(fit, "diagnostics_summary") <- diagnostics

  if (verbose) {
    message("Fitting complete")
  }

  return(fit)
}

#' Check Stan model diagnostics
#'
#' Checks convergence diagnostics and reports warnings.
#'
#' @param fit Fitted Stan model object
#' @param rhat_threshold Threshold for Rhat (default 1.01)
#' @param ess_threshold Threshold for ESS (default 400)
#' @param verbose Print diagnostic summary (default TRUE)
#' @return List with diagnostic results and flags
#' @export
check_stan_diagnostics <- function(fit,
                                     rhat_threshold = 1.01,
                                     ess_threshold = 400,
                                     verbose = TRUE) {

  # Get diagnostics
  diag <- fit$diagnostic_summary()
  summ <- fit$summary()

  # Check for divergences
  n_divergent <- sum(diag$num_divergent)

  # Check for max treedepth
  n_max_treedepth <- sum(diag$num_max_treedepth)

  # Check Rhat
  max_rhat <- max(summ$rhat, na.rm = TRUE)
  n_high_rhat <- sum(summ$rhat > rhat_threshold, na.rm = TRUE)

  # Check ESS
  min_ess_bulk <- min(summ$ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(summ$ess_tail, na.rm = TRUE)
  n_low_ess <- sum(summ$ess_bulk < ess_threshold | summ$ess_tail < ess_threshold, na.rm = TRUE)

  # Determine if diagnostics passed
  all_good <- (n_divergent == 0 &&
               n_max_treedepth == 0 &&
               max_rhat < rhat_threshold &&
               min_ess_bulk >= ess_threshold &&
               min_ess_tail >= ess_threshold)

  diagnostics <- list(
    n_divergent = n_divergent,
    n_max_treedepth = n_max_treedepth,
    max_rhat = max_rhat,
    n_high_rhat = n_high_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    n_low_ess = n_low_ess,
    all_diagnostics_passed = all_good
  )

  if (verbose) {
    message("\n=== Stan Diagnostics ===")

    if (all_good) {
      message("All diagnostics PASSED")
    } else {
      message("WARNING: Some diagnostics FAILED")
    }

    message("  Divergent transitions: ", n_divergent)
    message("  Max treedepth hits: ", n_max_treedepth)
    message("  Max Rhat: ", round(max_rhat, 4))
    message("  Parameters with Rhat > ", rhat_threshold, ": ", n_high_rhat)
    message("  Min ESS bulk: ", round(min_ess_bulk, 0))
    message("  Min ESS tail: ", round(min_ess_tail, 0))
    message("  Parameters with low ESS: ", n_low_ess)
    message("")

    if (!all_good) {
      message("Consider:")
      if (n_divergent > 0) {
        message("  - Increase adapt_delta (e.g., 0.99)")
      }
      if (n_max_treedepth > 0) {
        message("  - Increase max_treedepth (e.g., 15)")
      }
      if (n_high_rhat > 0 || n_low_ess > 0) {
        message("  - Run more iterations")
        message("  - Check for label switching or multimodality")
      }
    }
  }

  return(diagnostics)
}

#' Compute LOO for model comparison
#'
#' Calculates leave-one-out cross-validation using PSIS-LOO.
#'
#' @param fit Fitted Stan model object
#' @param cores Number of cores for parallel computation (default 2)
#' @return LOO object from loo package
#' @export
compute_loo <- function(fit, cores = 2) {

  # Extract log likelihood
  log_lik <- fit$draws("log_lik", format = "matrix")

  # Compute LOO
  loo_result <- loo::loo(log_lik, cores = cores)

  # Check for high Pareto k values
  high_k <- sum(loo_result$diagnostics$pareto_k > 0.7)

  if (high_k > 0) {
    warning("Found ", high_k, " observations with high Pareto k (>0.7). ",
            "LOO estimates may be unreliable for these observations.")
  }

  return(loo_result)
}

#' Compare multiple models using LOO
#'
#' Performs model comparison using LOO and computes model weights.
#'
#' @param ... Named fitted Stan models
#' @param cores Number of cores for parallel computation (default 2)
#' @return List with comparison results and model weights
#' @export
compare_models_loo <- function(..., cores = 2) {

  models <- list(...)

  if (length(models) < 2) {
    stop("Provide at least two models for comparison")
  }

  model_names <- names(models)
  if (is.null(model_names)) {
    model_names <- paste0("Model", 1:length(models))
  }

  message("Computing LOO for ", length(models), " models...")

  # Compute LOO for each model
  loo_list <- vector("list", length(models))
  names(loo_list) <- model_names

  for (i in seq_along(models)) {
    message("  ", model_names[i], "...")
    loo_list[[i]] <- compute_loo(models[[i]], cores = cores)
  }

  # Compare models
  comparison <- loo::loo_compare(loo_list)

  # Compute model weights
  weights <- loo::loo_model_weights(loo_list)

  message("\n=== Model Comparison (LOO) ===")
  print(comparison)
  message("\nModel weights:")
  print(round(weights, 3))
  message("")

  return(list(
    comparison = comparison,
    weights = weights,
    loo_list = loo_list
  ))
}

#' Extract posterior samples
#'
#' Extracts posterior samples for specific parameters from fitted model.
#'
#' @param fit Fitted Stan model object
#' @param pars Character vector of parameter names (optional, default all)
#' @param format Return format: "matrix", "data.frame", or "draws" (default "data.frame")
#' @return Posterior samples in requested format
#' @export
extract_posterior <- function(fit, pars = NULL, format = "data.frame") {

  if (is.null(pars)) {
    draws <- fit$draws(format = format)
  } else {
    draws <- fit$draws(variables = pars, format = format)
  }

  return(draws)
}

#' Extract boundary estimates
#'
#' Extracts and summarizes occupation boundary parameters.
#'
#' @param fit Fitted Stan model object
#' @param probs Quantiles to compute (default c(0.025, 0.16, 0.5, 0.84, 0.975))
#' @return Data frame with boundary summaries
#' @export
extract_boundaries <- function(fit, probs = c(0.025, 0.16, 0.5, 0.84, 0.975)) {

  summ <- fit$summary(
    variables = c("theta_start", "theta_end", "durations"),
    ~quantile(.x, probs = probs),
    mean,
    sd
  )

  return(summ)
}

#' Generate initial values for Stan model
#'
#' Creates reasonable initial values to aid convergence.
#'
#' @param stan_data Data list for Stan
#' @param model_type Type of model: "contemporaneous", "sequential", or "partial"
#' @return Function that returns initial values for each chain
#' @export
generate_inits <- function(stan_data, model_type = "partial") {

  init_fun <- function() {
    # Estimate reasonable starting values from data
    cal_means <- apply(matrix(unlist(stan_data$mix_means), nrow = stan_data$N), 1, mean)

    if (model_type == "contemporaneous") {
      list(
        theta_start = max(cal_means) + 100,
        theta_end = min(cal_means) - 100,
        calendar_dates = cal_means
      )
    } else {
      # For deposit-specific models
      n_deps <- stan_data$n_deposits
      dep_means <- numeric(n_deps)

      for (k in 1:n_deps) {
        idx <- which(stan_data$deposit_id == k)
        dep_means[k] <- mean(cal_means[idx])
      }

      list(
        theta_start = dep_means + 100,
        theta_end = dep_means - 100,
        calendar_dates_raw = rnorm(stan_data$N, 0, 0.1)
      )
    }
  }

  return(init_fun)
}

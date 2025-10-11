#' Model Validation Framework
#'
#' Functions for systematic validation of contemporaneity models using
#' simulation studies with known ground truth.
#'
#' @author Project Team
#' @date 2025-10-11

library(dplyr)
library(tidyr)
library(parallel)

#' Run single validation trial
#'
#' Generates simulated data, fits all models, and evaluates performance.
#'
#' @param scenario Character: "contemporaneous", "sequential", or "partial"
#' @param n_deposits Number of deposits
#' @param n_dates_per_deposit Sample size per deposit
#' @param measurement_error Measurement error (14C years)
#' @param seed Random seed for this trial
#' @param compile_models Logical, whether to compile models (default FALSE for pre-compiled)
#' @param model_contemp Pre-compiled contemporaneous model (optional)
#' @param model_seq Pre-compiled sequential model (optional)
#' @param model_partial Pre-compiled partial overlap model (optional)
#' @return List with trial results
#' @export
run_validation_trial <- function(scenario,
                                  n_deposits,
                                  n_dates_per_deposit,
                                  measurement_error,
                                  seed,
                                  compile_models = FALSE,
                                  model_contemp = NULL,
                                  model_seq = NULL,
                                  model_partial = NULL) {

  set.seed(seed)

  # Generate data based on scenario
  if (scenario == "contemporaneous") {
    sim_data <- simulate_contemporaneous_deposits(
      n_deposits = n_deposits,
      shared_window = c(2700, 2400),
      n_dates_per_deposit = n_dates_per_deposit,
      measurement_error = measurement_error,
      seed = seed
    )
  } else if (scenario == "sequential") {
    sim_data <- simulate_sequential_deposits(
      n_deposits = n_deposits,
      deposit_duration = 200,
      gap_between_deposits = 100,
      earliest_start = 3000,
      n_dates_per_deposit = n_dates_per_deposit,
      measurement_error = measurement_error,
      seed = seed
    )
  } else if (scenario == "partial") {
    sim_data <- simulate_partial_overlap(
      n_deposits = n_deposits,
      overlap_percentage = 50,
      deposit_duration = 300,
      earliest_start = 2800,
      n_dates_per_deposit = n_dates_per_deposit,
      measurement_error = measurement_error,
      seed = seed
    )
  } else {
    stop("Invalid scenario: ", scenario)
  }

  c14_data <- sim_data$dates[, c("lab_code", "age", "error", "deposit")]

  # Calibrate
  cal_result <- calibrate_dates(
    ages = c14_data$age,
    errors = c14_data$error,
    lab_codes = c14_data$lab_code,
    verbose = FALSE
  )

  # Prepare Stan data
  stan_data <- prepare_stan_data(c14_data, cal_result)

  # Compile models if needed
  if (compile_models) {
    model_contemp <- compile_stan_model("stan/contemporaneous_model.stan")
    model_seq <- compile_stan_model("stan/sequential_model.stan")
    model_partial <- compile_stan_model("stan/partial_overlap_model.stan")
  }

  # Fit models with reduced iterations for validation
  results <- list()

  tryCatch({
    fit_contemp <- fit_stan_model(
      model_contemp,
      stan_data,
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      verbose = FALSE
    )
    results$fit_contemp <- fit_contemp
    results$contemp_converged <- fit_contemp$diagnostics_summary$all_diagnostics_passed
  }, error = function(e) {
    results$fit_contemp <- NULL
    results$contemp_converged <- FALSE
  })

  tryCatch({
    fit_seq <- fit_stan_model(
      model_seq,
      stan_data,
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      verbose = FALSE
    )
    results$fit_seq <- fit_seq
    results$seq_converged <- fit_seq$diagnostics_summary$all_diagnostics_passed
  }, error = function(e) {
    results$fit_seq <- NULL
    results$seq_converged <- FALSE
  })

  tryCatch({
    fit_partial <- fit_stan_model(
      model_partial,
      stan_data,
      chains = 2,
      iter_warmup = 500,
      iter_sampling = 500,
      verbose = FALSE
    )
    results$fit_partial <- fit_partial
    results$partial_converged <- fit_partial$diagnostics_summary$all_diagnostics_passed
  }, error = function(e) {
    results$fit_partial <- NULL
    results$partial_converged <- FALSE
  })

  # Compare models if all converged
  if (results$contemp_converged && results$seq_converged && results$partial_converged) {
    comparison <- compare_models_loo(
      contemporaneous = results$fit_contemp,
      sequential = results$fit_seq,
      partial = results$fit_partial,
      cores = 1
    )

    best_model <- names(comparison$weights)[which.max(comparison$weights)]
    results$selected_model <- best_model
    results$model_weights <- comparison$weights
  } else {
    results$selected_model <- NA
    results$model_weights <- c(NA, NA, NA)
  }

  # Calculate metrics if partial model converged
  if (results$partial_converged) {
    # Estimate occupation spans
    spans <- calculate_occupation_spans(results$fit_partial, n_deposits)

    # Compare to ground truth
    true_windows <- sim_data$ground_truth$deposit_windows

    results$boundary_bias <- list()
    results$boundary_rmse <- list()
    results$coverage_95 <- list()

    for (k in 1:n_deposits) {
      # Start boundary
      start_bias <- mean(extract_posterior(results$fit_partial, "theta_start", "matrix")[, k]) -
                    true_windows[k, 1]
      start_covered <- (true_windows[k, 1] >= spans$start_q2.5[k] &
                        true_windows[k, 1] <= spans$start_q97.5[k])

      # End boundary
      end_bias <- mean(extract_posterior(results$fit_partial, "theta_end", "matrix")[, k]) -
                  true_windows[k, 2]
      end_covered <- (true_windows[k, 2] >= spans$end_q2.5[k] &
                      true_windows[k, 2] <= spans$end_q97.5[k])

      results$boundary_bias[[k]] <- list(start = start_bias, end = end_bias)
      results$coverage_95[[k]] <- list(start = start_covered, end = end_covered)
    }

    # Overall coverage rate
    all_covered <- unlist(lapply(results$coverage_95, function(x) c(x$start, x$end)))
    results$overall_coverage <- mean(all_covered)

  } else {
    results$boundary_bias <- NA
    results$boundary_rmse <- NA
    results$coverage_95 <- NA
    results$overall_coverage <- NA
  }

  # Ground truth for validation
  results$true_scenario <- scenario
  results$true_windows <- sim_data$ground_truth$deposit_windows

  return(results)
}

#' Run comprehensive validation study
#'
#' Conducts systematic validation across multiple scenarios and parameters.
#'
#' @param n_trials Number of trials per condition (default 100)
#' @param scenarios Character vector of scenarios to test
#' @param sample_sizes Vector of sample sizes to test
#' @param n_deposits_range Vector of deposit numbers to test
#' @param measurement_errors Vector of measurement errors to test
#' @param parallel Logical, use parallel processing (default TRUE)
#' @param n_cores Number of cores for parallel processing (default 4)
#' @return Data frame with validation results
#' @export
run_validation_study <- function(n_trials = 100,
                                  scenarios = c("contemporaneous", "sequential", "partial"),
                                  sample_sizes = c(5, 10, 20),
                                  n_deposits_range = c(2, 3),
                                  measurement_errors = c(30),
                                  parallel = TRUE,
                                  n_cores = 4) {

  message("=== Starting Validation Study ===")
  message("Trials per condition: ", n_trials)
  message("Scenarios: ", paste(scenarios, collapse = ", "))
  message("Sample sizes: ", paste(sample_sizes, collapse = ", "))

  # Pre-compile models once
  message("\nCompiling Stan models...")
  model_contemp <- compile_stan_model("stan/contemporaneous_model.stan")
  model_seq <- compile_stan_model("stan/sequential_model.stan")
  model_partial <- compile_stan_model("stan/partial_overlap_model.stan")

  # Create parameter grid
  param_grid <- expand.grid(
    scenario = scenarios,
    n_deposits = n_deposits_range,
    sample_size = sample_sizes,
    measurement_error = measurement_errors,
    trial = 1:n_trials,
    stringsAsFactors = FALSE
  )

  message("\nTotal conditions: ", nrow(param_grid))

  # Run trials
  if (parallel) {
    message("Running in parallel with ", n_cores, " cores...")
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("run_validation_trial", "simulate_contemporaneous_deposits",
                        "simulate_sequential_deposits", "simulate_partial_overlap",
                        "calibrate_dates", "prepare_stan_data", "fit_stan_model",
                        "compare_models_loo", "calculate_occupation_spans",
                        "extract_posterior"))

    results_list <- parLapply(cl, 1:nrow(param_grid), function(i) {
      row <- param_grid[i, ]
      seed <- i + 1000

      result <- run_validation_trial(
        scenario = row$scenario,
        n_deposits = row$n_deposits,
        n_dates_per_deposit = row$sample_size,
        measurement_error = row$measurement_error,
        seed = seed,
        compile_models = FALSE,
        model_contemp = model_contemp,
        model_seq = model_seq,
        model_partial = model_partial
      )

      return(result)
    })

    stopCluster(cl)

  } else {
    message("Running sequentially...")
    results_list <- lapply(1:nrow(param_grid), function(i) {
      row <- param_grid[i, ]
      seed <- i + 1000

      if (i %% 10 == 0) {
        message("  Trial ", i, " of ", nrow(param_grid))
      }

      result <- run_validation_trial(
        scenario = row$scenario,
        n_deposits = row$n_deposits,
        n_dates_per_deposit = row$sample_size,
        measurement_error = row$measurement_error,
        seed = seed,
        compile_models = FALSE,
        model_contemp = model_contemp,
        model_seq = model_seq,
        model_partial = model_partial
      )

      return(result)
    })
  }

  # Compile results
  message("\nCompiling results...")
  results_df <- param_grid
  results_df$selected_model <- sapply(results_list, function(x) x$selected_model)
  results_df$overall_coverage <- sapply(results_list, function(x) x$overall_coverage)
  results_df$contemp_converged <- sapply(results_list, function(x) x$contemp_converged)
  results_df$seq_converged <- sapply(results_list, function(x) x$seq_converged)
  results_df$partial_converged <- sapply(results_list, function(x) x$partial_converged)

  # Calculate accuracy metrics
  results_df$correct_model <- with(results_df, {
    correct <- rep(FALSE, nrow(results_df))
    correct[scenario == "contemporaneous" & selected_model == "contemporaneous"] <- TRUE
    correct[scenario == "sequential" & selected_model == "sequential"] <- TRUE
    correct[scenario == "partial" & selected_model == "partial"] <- TRUE
    correct
  })

  message("Validation study complete!")

  return(results_df)
}

#' Summarize validation results
#'
#' Computes summary statistics from validation study.
#'
#' @param results_df Output from run_validation_study()
#' @return Data frame with summary statistics
#' @export
summarize_validation_results <- function(results_df) {

  summary_df <- results_df %>%
    group_by(scenario, n_deposits, sample_size, measurement_error) %>%
    summarize(
      n_trials = n(),
      convergence_rate = mean(partial_converged, na.rm = TRUE),
      selection_accuracy = mean(correct_model, na.rm = TRUE),
      coverage_rate = mean(overall_coverage, na.rm = TRUE),
      .groups = "drop"
    )

  return(summary_df)
}

#' Print validation summary
#'
#' @param summary_df Output from summarize_validation_results()
#' @export
print_validation_summary <- function(summary_df) {

  cat("=== Validation Study Summary ===\n\n")

  for (sc in unique(summary_df$scenario)) {
    cat("Scenario:", sc, "\n")

    sub_df <- summary_df[summary_df$scenario == sc, ]

    cat("  Sample Size | Convergence | Accuracy | Coverage\n")
    cat("  ", rep("-", 50), "\n", sep = "")

    for (i in 1:nrow(sub_df)) {
      cat(sprintf("  %11d | %10.1f%% | %7.1f%% | %7.1f%%\n",
                  sub_df$sample_size[i],
                  sub_df$convergence_rate[i] * 100,
                  sub_df$selection_accuracy[i] * 100,
                  sub_df$coverage_rate[i] * 100))
    }
    cat("\n")
  }
}

#' Plot validation results
#'
#' Creates visualization of validation study results.
#'
#' @param results_df Output from run_validation_study()
#' @return List of ggplot objects
#' @export
plot_validation_study <- function(results_df) {

  summary_df <- summarize_validation_results(results_df)

  plots <- list()

  # Accuracy by sample size
  plots$accuracy <- ggplot(summary_df, aes(x = sample_size, y = selection_accuracy,
                                            color = scenario, group = scenario)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
    ylim(0, 1) +
    labs(
      x = "Sample size per deposit",
      y = "Model selection accuracy",
      title = "Validation: Model Selection Accuracy",
      subtitle = "Dashed line: 70% threshold",
      color = "Scenario"
    ) +
    theme_minimal()

  # Coverage by sample size
  plots$coverage <- ggplot(summary_df, aes(x = sample_size, y = coverage_rate,
                                            color = scenario, group = scenario)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue") +
    ylim(0.7, 1) +
    labs(
      x = "Sample size per deposit",
      y = "95% credible interval coverage",
      title = "Validation: Interval Coverage",
      subtitle = "Dashed line: nominal 95% coverage",
      color = "Scenario"
    ) +
    theme_minimal()

  # Convergence rate
  plots$convergence <- ggplot(summary_df, aes(x = sample_size, y = convergence_rate,
                                               color = scenario, group = scenario)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    ylim(0, 1) +
    labs(
      x = "Sample size per deposit",
      y = "Convergence rate",
      title = "Validation: Stan Convergence Rate",
      color = "Scenario"
    ) +
    theme_minimal()

  return(plots)
}

#' Radiocarbon Calibration Functions
#'
#' Wrapper functions around rcarbon for calibrating radiocarbon dates
#' and preparing calibrated distributions for Stan models.
#'
#' References:
#'   Reimer et al. (2020). The IntCal20 Northern Hemisphere radiocarbon
#'   age calibration curve (0-55 cal kBP). Radiocarbon, 62(4), 725-757.
#'
#' @author Project Team
#' @date 2025-10-11

library(rcarbon)

#' Calibrate radiocarbon dates
#'
#' Wrapper around rcarbon::calibrate() with metadata tracking and
#' standardized output format.
#'
#' @param ages Numeric vector of radiocarbon ages in BP
#' @param errors Numeric vector of standard errors in radiocarbon years
#' @param lab_codes Character vector of laboratory codes (optional)
#' @param calCurve Character, calibration curve name (default "intcal20")
#' @param normalize Logical, normalize calibrated distributions (default TRUE)
#' @param verbose Logical, print calibration messages (default FALSE)
#' @return List containing:
#'   - calDates: calibrated dates object from rcarbon
#'   - metadata: list with calibration settings
#'   - summary: data frame with summary statistics
#' @examples
#' ages <- c(2450, 2380, 2420)
#' errors <- c(30, 25, 35)
#' cal_result <- calibrate_dates(ages, errors)
#' @export
calibrate_dates <- function(ages,
                             errors,
                             lab_codes = NULL,
                             calCurve = "intcal20",
                             normalize = TRUE,
                             verbose = FALSE) {

  # Input validation
  if (length(ages) != length(errors)) {
    stop("ages and errors must have the same length")
  }

  if (any(is.na(ages)) || any(is.na(errors))) {
    stop("ages and errors cannot contain NA values")
  }

  if (any(errors <= 0)) {
    stop("errors must be positive")
  }

  # Set default lab codes if not provided
  if (is.null(lab_codes)) {
    lab_codes <- paste0("Date_", seq_along(ages))
  }

  if (length(lab_codes) != length(ages)) {
    stop("lab_codes must have same length as ages")
  }

  # Record metadata
  metadata <- list(
    calibration_curve = calCurve,
    date_of_calibration = Sys.Date(),
    rcarbon_version = as.character(packageVersion("rcarbon")),
    normalize = normalize,
    n_dates = length(ages)
  )

  if (verbose) {
    message("Calibrating ", length(ages), " radiocarbon dates...")
    message("Using calibration curve: ", calCurve)
  }

  # Calibrate using rcarbon
  calDates <- rcarbon::calibrate(
    x = ages,
    errors = errors,
    calCurves = calCurve,
    ids = lab_codes,
    normalised = normalize,
    verbose = verbose
  )

  # Extract summary statistics
  summary_stats <- extract_calibration_summary(calDates, lab_codes)

  result <- list(
    calDates = calDates,
    metadata = metadata,
    summary = summary_stats
  )

  class(result) <- c("calibrated_c14", "list")

  if (verbose) {
    message("Calibration complete")
  }

  return(result)
}

#' Extract summary statistics from calibrated dates
#'
#' Calculates median, mode, and credible intervals from calibrated
#' probability distributions.
#'
#' @param calDates Calibrated dates object from rcarbon
#' @param lab_codes Character vector of laboratory codes
#' @return Data frame with summary statistics for each date
#' @export
extract_calibration_summary <- function(calDates, lab_codes = NULL) {

  n_dates <- length(calDates$grids)

  if (is.null(lab_codes)) {
    lab_codes <- paste0("Date_", 1:n_dates)
  }

  summary_df <- data.frame(
    lab_code = lab_codes,
    median_cal_bp = NA_real_,
    mean_cal_bp = NA_real_,
    hdi68_lower = NA_real_,
    hdi68_upper = NA_real_,
    hdi95_lower = NA_real_,
    hdi95_upper = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in 1:n_dates) {
    grid <- calDates$grids[[i]]
    cal_bp <- grid$calBP
    prob_dens <- grid$PrDens

    # Calculate median (50th percentile)
    cumsum_prob <- cumsum(prob_dens) / sum(prob_dens)
    summary_df$median_cal_bp[i] <- cal_bp[which.min(abs(cumsum_prob - 0.5))]

    # Calculate mean
    summary_df$mean_cal_bp[i] <- sum(cal_bp * prob_dens) / sum(prob_dens)

    # Calculate HDI (Highest Density Intervals)
    hdi68 <- calculate_hdi(cal_bp, prob_dens, prob = 0.68)
    hdi95 <- calculate_hdi(cal_bp, prob_dens, prob = 0.95)

    summary_df$hdi68_lower[i] <- hdi68[1]
    summary_df$hdi68_upper[i] <- hdi68[2]
    summary_df$hdi95_lower[i] <- hdi95[1]
    summary_df$hdi95_upper[i] <- hdi95[2]
  }

  return(summary_df)
}

#' Calculate Highest Density Interval (HDI)
#'
#' Finds the narrowest interval containing the specified probability mass.
#'
#' @param x Numeric vector of values (e.g., calendar years)
#' @param density Numeric vector of probability densities
#' @param prob Probability mass to include (default 0.95)
#' @return Numeric vector of length 2 with lower and upper bounds
#' @export
calculate_hdi <- function(x, density, prob = 0.95) {

  # Normalize density
  density <- density / sum(density)

  # Sort by density (highest first)
  ord <- order(density, decreasing = TRUE)
  x_sorted <- x[ord]
  dens_sorted <- density[ord]

  # Find cumulative probability
  cumsum_dens <- cumsum(dens_sorted)

  # Find indices that contain target probability
  idx <- which(cumsum_dens <= prob)

  if (length(idx) == 0) {
    idx <- 1
  }

  # Get the range of these values
  hdi_values <- x_sorted[idx]
  hdi <- c(max(hdi_values), min(hdi_values))  # Note: reversed due to BP convention

  return(hdi)
}

#' Approximate calibrated distribution with mixture of normals
#'
#' Fits a mixture of normal distributions to a calibrated radiocarbon
#' probability distribution for use in Stan models. This is more efficient
#' than including the full calibration curve in Stan.
#'
#' @param calDate Single calibrated date object from rcarbon
#' @param n_components Number of mixture components (default 3)
#' @return List with mixture parameters (weights, means, sds)
#' @export
approximate_calibration_mixture <- function(calDate, n_components = 3) {

  grid <- calDate$grids[[1]]
  cal_bp <- grid$calBP
  prob_dens <- grid$PrDens

  # Normalize
  prob_dens <- prob_dens / sum(prob_dens)

  # Fit mixture model using EM algorithm (simple implementation)
  # For production, consider using mixtools package
  mixture <- fit_gaussian_mixture_simple(cal_bp, prob_dens, n_components)

  return(mixture)
}

#' Simple Gaussian mixture fitting via EM
#'
#' Fits a mixture of Gaussians to discrete probability distribution
#' using Expectation-Maximization.
#'
#' @param x Numeric vector of values
#' @param weights Numeric vector of weights/probabilities
#' @param n_components Number of mixture components
#' @param max_iter Maximum EM iterations (default 100)
#' @param tol Convergence tolerance (default 1e-6)
#' @return List with mixture parameters
#' @keywords internal
fit_gaussian_mixture_simple <- function(x, weights, n_components,
                                         max_iter = 100, tol = 1e-6) {

  n <- length(x)
  weights <- weights / sum(weights)

  # Initialize parameters
  # Use quantiles to initialize means
  init_probs <- seq(0, 1, length.out = n_components + 2)[2:(n_components + 1)]
  means <- quantile(x, probs = init_probs)
  sds <- rep(sd(x) / sqrt(n_components), n_components)
  mix_weights <- rep(1 / n_components, n_components)

  log_lik_old <- -Inf

  for (iter in 1:max_iter) {

    # E-step: Calculate responsibilities
    resp_matrix <- matrix(0, n, n_components)
    for (k in 1:n_components) {
      resp_matrix[, k] <- mix_weights[k] * dnorm(x, means[k], sds[k])
    }

    # Normalize responsibilities
    resp_sum <- rowSums(resp_matrix)
    resp_matrix <- resp_matrix / (resp_sum + 1e-10)

    # Weight by original probabilities
    weighted_resp <- resp_matrix * weights

    # M-step: Update parameters
    for (k in 1:n_components) {
      nk <- sum(weighted_resp[, k])

      mix_weights[k] <- nk / sum(weights)
      means[k] <- sum(weighted_resp[, k] * x) / nk
      sds[k] <- sqrt(sum(weighted_resp[, k] * (x - means[k])^2) / nk)
    }

    # Calculate log likelihood
    log_lik <- sum(weights * log(resp_sum + 1e-10))

    # Check convergence
    if (abs(log_lik - log_lik_old) < tol) {
      break
    }

    log_lik_old <- log_lik
  }

  # Ensure weights sum to exactly 1.0 (avoid numerical precision issues)
  mix_weights <- mix_weights / sum(mix_weights)

  # Return mixture parameters
  list(
    weights = mix_weights,
    means = means,
    sds = sds,
    n_components = n_components,
    log_likelihood = log_lik,
    n_iterations = iter
  )
}

#' Convert calibrated dates to Stan-friendly format
#'
#' Prepares calibrated radiocarbon distributions for use in Stan models.
#' Two approaches:
#' 1. Mixture approximation (faster, recommended)
#' 2. Full calibration curve (more accurate, slower)
#'
#' @param cal_result Calibrated dates from calibrate_dates()
#' @param method Character: "mixture" or "full" (default "mixture")
#' @param n_components Number of mixture components if method="mixture" (default 3)
#' @return List formatted for Stan
#' @export
prepare_stan_calibration_data <- function(cal_result,
                                           method = "mixture",
                                           n_components = 3) {

  if (!inherits(cal_result, "calibrated_c14")) {
    stop("cal_result must be output from calibrate_dates()")
  }

  n_dates <- length(cal_result$calDates$grids)

  if (method == "mixture") {
    # Approximate each calibrated distribution with mixture of normals
    mixture_params <- vector("list", n_dates)

    for (i in 1:n_dates) {
      mixture_params[[i]] <- approximate_calibration_mixture(
        list(grids = list(cal_result$calDates$grids[[i]])),
        n_components = n_components
      )
    }

    # Format for Stan
    stan_data <- list(
      K = n_components,
      mix_weights = do.call(rbind, lapply(mixture_params, function(x) x$weights)),
      mix_means = do.call(rbind, lapply(mixture_params, function(x) x$means)),
      mix_sds = do.call(rbind, lapply(mixture_params, function(x) x$sds))
    )

  } else if (method == "full") {
    # Load full calibration curve
    cal_curve <- get_calibration_curve(cal_result$metadata$calibration_curve)

    stan_data <- list(
      N_cal_curve = nrow(cal_curve),
      cal_curve_bp = cal_curve$cal_bp,
      cal_curve_c14 = cal_curve$c14_age,
      cal_curve_error = cal_curve$error
    )

  } else {
    stop("method must be 'mixture' or 'full'")
  }

  return(stan_data)
}

#' Load calibration curve data
#'
#' Retrieves calibration curve data from rcarbon package.
#'
#' @param curve_name Character, name of calibration curve
#' @return Data frame with calibration curve
#' @export
get_calibration_curve <- function(curve_name = "intcal20") {

  # Get calibration curve from rcarbon
  if (curve_name == "intcal20") {
    data("intcal20", package = "rcarbon", envir = environment())
    curve <- intcal20
  } else if (curve_name == "intcal13") {
    data("intcal13", package = "rcarbon", envir = environment())
    curve <- intcal13
  } else if (curve_name == "shcal20") {
    data("shcal20", package = "rcarbon", envir = environment())
    curve <- shcal20
  } else {
    stop("Unsupported calibration curve: ", curve_name)
  }

  # Standardize column names
  cal_curve_df <- data.frame(
    cal_bp = curve[, 1],
    c14_age = curve[, 2],
    error = curve[, 3]
  )

  return(cal_curve_df)
}

#' Print calibration results
#'
#' @param cal_result Output from calibrate_dates()
#' @export
print.calibrated_c14 <- function(cal_result) {
  cat("=== Calibrated Radiocarbon Dates ===\n\n")
  cat("Number of dates:", cal_result$metadata$n_dates, "\n")
  cat("Calibration curve:", cal_result$metadata$calibration_curve, "\n")
  cat("Calibration date:", as.character(cal_result$metadata$date_of_calibration), "\n")
  cat("rcarbon version:", cal_result$metadata$rcarbon_version, "\n\n")

  cat("Summary statistics (cal BP):\n")
  print(cal_result$summary, row.names = FALSE)
  cat("\n")
}

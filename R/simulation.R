#' Simulation Framework for Radiocarbon Dating Analysis
#'
#' Functions to generate synthetic radiocarbon datasets with known
#' temporal structure for model validation and testing.
#'
#' @author Project Team
#' @date 2025-10-11

#' Simulate contemporaneous deposits
#'
#' Generates radiocarbon dates from deposits that share a common
#' occupation window (fully or partially overlapping).
#'
#' @param n_deposits Number of deposits to simulate
#' @param shared_window Numeric vector of length 2: c(start, end) in cal BP
#' @param n_dates_per_deposit Integer or vector of sample sizes per deposit
#' @param deposit_offsets Numeric vector of offsets for each deposit (optional)
#' @param measurement_error Standard deviation of measurement error in 14C years (default 30)
#' @param calCurve Calibration curve to use (default "intcal20")
#' @param seed Random seed for reproducibility
#' @return List with simulated data and ground truth
#' @examples
#' sim <- simulate_contemporaneous_deposits(
#'   n_deposits = 3,
#'   shared_window = c(2700, 2400),
#'   n_dates_per_deposit = 10,
#'   measurement_error = 30
#' )
#' @export
simulate_contemporaneous_deposits <- function(n_deposits,
                                               shared_window,
                                               n_dates_per_deposit,
                                               deposit_offsets = NULL,
                                               measurement_error = 30,
                                               calCurve = "intcal20",
                                               seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Validate inputs
  if (length(shared_window) != 2) {
    stop("shared_window must be a vector of length 2: c(start, end)")
  }

  if (shared_window[1] <= shared_window[2]) {
    stop("shared_window[1] (start) must be > shared_window[2] (end) in cal BP")
  }

  # Expand n_dates_per_deposit if single value
  if (length(n_dates_per_deposit) == 1) {
    n_dates_per_deposit <- rep(n_dates_per_deposit, n_deposits)
  }

  if (length(n_dates_per_deposit) != n_deposits) {
    stop("n_dates_per_deposit must be length 1 or equal to n_deposits")
  }

  # Set deposit offsets (for partial overlap)
  if (is.null(deposit_offsets)) {
    deposit_offsets <- rep(0, n_deposits)
  }

  if (length(deposit_offsets) != n_deposits) {
    stop("deposit_offsets must have length equal to n_deposits")
  }

  # Calculate deposit-specific windows
  deposit_windows <- matrix(0, n_deposits, 2)
  for (k in 1:n_deposits) {
    deposit_windows[k, ] <- shared_window + c(deposit_offsets[k], deposit_offsets[k])
  }

  # Load calibration curve
  cal_curve <- get_calibration_curve(calCurve)

  # Generate true calendar dates
  n_total <- sum(n_dates_per_deposit)
  true_cal_bp <- numeric(n_total)
  deposit_id <- integer(n_total)
  lab_codes <- character(n_total)

  idx <- 1
  for (k in 1:n_deposits) {
    n_k <- n_dates_per_deposit[k]
    indices <- idx:(idx + n_k - 1)

    # Sample uniformly from deposit window
    true_cal_bp[indices] <- runif(
      n_k,
      min = deposit_windows[k, 2],
      max = deposit_windows[k, 1]
    )

    deposit_id[indices] <- k
    lab_codes[indices] <- paste0("SIM-", LETTERS[k], "-", 1:n_k)

    idx <- idx + n_k
  }

  # Convert calendar dates to radiocarbon ages
  sim_data <- calendar_to_c14(
    true_cal_bp = true_cal_bp,
    measurement_error = measurement_error,
    cal_curve = cal_curve
  )

  # Create output data frame
  simulated_dates <- data.frame(
    lab_code = lab_codes,
    age = sim_data$c14_age,
    error = sim_data$c14_error,
    deposit = factor(deposit_id),
    true_cal_bp = true_cal_bp,
    stringsAsFactors = FALSE
  )

  # Return results
  result <- list(
    dates = simulated_dates,
    ground_truth = list(
      shared_window = shared_window,
      deposit_windows = deposit_windows,
      true_cal_bp = true_cal_bp,
      deposit_id = deposit_id,
      overlap_type = "contemporaneous"
    ),
    parameters = list(
      n_deposits = n_deposits,
      n_dates_per_deposit = n_dates_per_deposit,
      measurement_error = measurement_error,
      calCurve = calCurve,
      seed = seed
    )
  )

  class(result) <- c("simulated_c14", "list")
  return(result)
}

#' Simulate sequential deposits
#'
#' Generates radiocarbon dates from non-overlapping deposits with
#' defined temporal sequence.
#'
#' @param n_deposits Number of deposits to simulate
#' @param deposit_duration Duration of each deposit in calendar years (scalar or vector)
#' @param gap_between_deposits Gap between deposits in calendar years (scalar or vector)
#' @param earliest_start Starting date of first deposit in cal BP
#' @param n_dates_per_deposit Integer or vector of sample sizes per deposit
#' @param measurement_error Standard deviation of measurement error in 14C years (default 30)
#' @param calCurve Calibration curve to use (default "intcal20")
#' @param seed Random seed for reproducibility
#' @return List with simulated data and ground truth
#' @examples
#' sim <- simulate_sequential_deposits(
#'   n_deposits = 3,
#'   deposit_duration = 200,
#'   gap_between_deposits = 100,
#'   earliest_start = 3000,
#'   n_dates_per_deposit = 10
#' )
#' @export
simulate_sequential_deposits <- function(n_deposits,
                                          deposit_duration,
                                          gap_between_deposits,
                                          earliest_start,
                                          n_dates_per_deposit,
                                          measurement_error = 30,
                                          calCurve = "intcal20",
                                          seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Expand vectors if needed
  if (length(deposit_duration) == 1) {
    deposit_duration <- rep(deposit_duration, n_deposits)
  }

  if (length(gap_between_deposits) == 1) {
    gap_between_deposits <- rep(gap_between_deposits, n_deposits - 1)
  } else if (length(gap_between_deposits) != n_deposits - 1) {
    stop("gap_between_deposits must be length 1 or n_deposits - 1")
  }

  if (length(n_dates_per_deposit) == 1) {
    n_dates_per_deposit <- rep(n_dates_per_deposit, n_deposits)
  }

  # Calculate deposit windows
  deposit_windows <- matrix(0, n_deposits, 2)
  current_start <- earliest_start

  for (k in 1:n_deposits) {
    deposit_windows[k, 1] <- current_start  # start
    deposit_windows[k, 2] <- current_start - deposit_duration[k]  # end

    if (k < n_deposits) {
      current_start <- deposit_windows[k, 2] - gap_between_deposits[k]
    }
  }

  # Load calibration curve
  cal_curve <- get_calibration_curve(calCurve)

  # Generate true calendar dates
  n_total <- sum(n_dates_per_deposit)
  true_cal_bp <- numeric(n_total)
  deposit_id <- integer(n_total)
  lab_codes <- character(n_total)

  idx <- 1
  for (k in 1:n_deposits) {
    n_k <- n_dates_per_deposit[k]
    indices <- idx:(idx + n_k - 1)

    # Sample uniformly from deposit window
    true_cal_bp[indices] <- runif(
      n_k,
      min = deposit_windows[k, 2],
      max = deposit_windows[k, 1]
    )

    deposit_id[indices] <- k
    lab_codes[indices] <- paste0("SIM-", LETTERS[k], "-", 1:n_k)

    idx <- idx + n_k
  }

  # Convert to radiocarbon ages
  sim_data <- calendar_to_c14(
    true_cal_bp = true_cal_bp,
    measurement_error = measurement_error,
    cal_curve = cal_curve
  )

  # Create output data frame
  simulated_dates <- data.frame(
    lab_code = lab_codes,
    age = sim_data$c14_age,
    error = sim_data$c14_error,
    deposit = factor(deposit_id),
    true_cal_bp = true_cal_bp,
    stringsAsFactors = FALSE
  )

  # Return results
  result <- list(
    dates = simulated_dates,
    ground_truth = list(
      deposit_windows = deposit_windows,
      true_cal_bp = true_cal_bp,
      deposit_id = deposit_id,
      overlap_type = "sequential"
    ),
    parameters = list(
      n_deposits = n_deposits,
      deposit_duration = deposit_duration,
      gap_between_deposits = gap_between_deposits,
      n_dates_per_deposit = n_dates_per_deposit,
      measurement_error = measurement_error,
      calCurve = calCurve,
      seed = seed
    )
  )

  class(result) <- c("simulated_c14", "list")
  return(result)
}

#' Simulate partial overlap between deposits
#'
#' Generates radiocarbon dates from deposits with specified percentage
#' of temporal overlap.
#'
#' @param n_deposits Number of deposits to simulate
#' @param overlap_percentage Percentage of overlap (0-100)
#' @param deposit_duration Duration of each deposit in calendar years
#' @param earliest_start Starting date of first deposit in cal BP
#' @param n_dates_per_deposit Integer or vector of sample sizes per deposit
#' @param measurement_error Standard deviation of measurement error in 14C years (default 30)
#' @param calCurve Calibration curve to use (default "intcal20")
#' @param seed Random seed for reproducibility
#' @return List with simulated data and ground truth
#' @examples
#' sim <- simulate_partial_overlap(
#'   n_deposits = 2,
#'   overlap_percentage = 50,
#'   deposit_duration = 300,
#'   earliest_start = 2800,
#'   n_dates_per_deposit = 15
#' )
#' @export
simulate_partial_overlap <- function(n_deposits,
                                      overlap_percentage,
                                      deposit_duration,
                                      earliest_start,
                                      n_dates_per_deposit,
                                      measurement_error = 30,
                                      calCurve = "intcal20",
                                      seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Validate overlap percentage
  if (overlap_percentage < 0 || overlap_percentage > 100) {
    stop("overlap_percentage must be between 0 and 100")
  }

  # Expand vectors
  if (length(deposit_duration) == 1) {
    deposit_duration <- rep(deposit_duration, n_deposits)
  }

  if (length(n_dates_per_deposit) == 1) {
    n_dates_per_deposit <- rep(n_dates_per_deposit, n_deposits)
  }

  # Calculate deposit windows with specified overlap
  deposit_windows <- matrix(0, n_deposits, 2)
  deposit_windows[1, 1] <- earliest_start
  deposit_windows[1, 2] <- earliest_start - deposit_duration[1]

  for (k in 2:n_deposits) {
    # Calculate offset based on overlap percentage
    overlap_duration <- deposit_duration[k - 1] * (overlap_percentage / 100)
    offset <- deposit_duration[k - 1] - overlap_duration

    deposit_windows[k, 1] <- deposit_windows[k - 1, 1] - offset
    deposit_windows[k, 2] <- deposit_windows[k, 1] - deposit_duration[k]
  }

  # Load calibration curve
  cal_curve <- get_calibration_curve(calCurve)

  # Generate true calendar dates
  n_total <- sum(n_dates_per_deposit)
  true_cal_bp <- numeric(n_total)
  deposit_id <- integer(n_total)
  lab_codes <- character(n_total)

  idx <- 1
  for (k in 1:n_deposits) {
    n_k <- n_dates_per_deposit[k]
    indices <- idx:(idx + n_k - 1)

    # Sample uniformly from deposit window
    true_cal_bp[indices] <- runif(
      n_k,
      min = deposit_windows[k, 2],
      max = deposit_windows[k, 1]
    )

    deposit_id[indices] <- k
    lab_codes[indices] <- paste0("SIM-", LETTERS[k], "-", 1:n_k)

    idx <- idx + n_k
  }

  # Convert to radiocarbon ages
  sim_data <- calendar_to_c14(
    true_cal_bp = true_cal_bp,
    measurement_error = measurement_error,
    cal_curve = cal_curve
  )

  # Create output data frame
  simulated_dates <- data.frame(
    lab_code = lab_codes,
    age = sim_data$c14_age,
    error = sim_data$c14_error,
    deposit = factor(deposit_id),
    true_cal_bp = true_cal_bp,
    stringsAsFactors = FALSE
  )

  # Return results
  result <- list(
    dates = simulated_dates,
    ground_truth = list(
      deposit_windows = deposit_windows,
      true_cal_bp = true_cal_bp,
      deposit_id = deposit_id,
      overlap_percentage = overlap_percentage,
      overlap_type = "partial"
    ),
    parameters = list(
      n_deposits = n_deposits,
      deposit_duration = deposit_duration,
      overlap_percentage = overlap_percentage,
      n_dates_per_deposit = n_dates_per_deposit,
      measurement_error = measurement_error,
      calCurve = calCurve,
      seed = seed
    )
  )

  class(result) <- c("simulated_c14", "list")
  return(result)
}

#' Convert calendar dates to radiocarbon ages
#'
#' Models the radiocarbon dating process: applies calibration curve
#' and adds measurement error.
#'
#' @param true_cal_bp Numeric vector of true calendar dates in cal BP
#' @param measurement_error Standard deviation of measurement error (scalar or vector)
#' @param cal_curve Data frame with calibration curve (from get_calibration_curve)
#' @return List with simulated radiocarbon ages and errors
#' @export
calendar_to_c14 <- function(true_cal_bp, measurement_error, cal_curve) {

  n <- length(true_cal_bp)

  # Expand measurement error if scalar
  if (length(measurement_error) == 1) {
    measurement_error <- rep(measurement_error, n)
  }

  # Interpolate calibration curve at true calendar dates
  true_c14_age <- approx(
    x = cal_curve$cal_bp,
    y = cal_curve$c14_age,
    xout = true_cal_bp,
    rule = 2  # Use nearest value for out-of-range
  )$y

  cal_curve_error <- approx(
    x = cal_curve$cal_bp,
    y = cal_curve$error,
    xout = true_cal_bp,
    rule = 2
  )$y

  # Combine calibration curve uncertainty and measurement error
  # Total error is quadrature sum
  total_error <- sqrt(cal_curve_error^2 + measurement_error^2)

  # Add measurement error to get observed radiocarbon age
  observed_c14_age <- rnorm(n, mean = true_c14_age, sd = total_error)

  return(list(
    c14_age = observed_c14_age,
    c14_error = measurement_error,  # Report lab measurement error
    true_c14_age = true_c14_age,
    cal_curve_error = cal_curve_error
  ))
}

#' Calculate true overlap between deposits
#'
#' Computes the actual temporal overlap from ground truth data.
#'
#' @param ground_truth Ground truth list from simulation functions
#' @return Matrix of pairwise overlap proportions
#' @export
calculate_true_overlap <- function(ground_truth) {

  windows <- ground_truth$deposit_windows
  n_deposits <- nrow(windows)

  overlap_matrix <- matrix(0, n_deposits, n_deposits)

  for (i in 1:n_deposits) {
    for (j in 1:n_deposits) {
      if (i == j) {
        overlap_matrix[i, j] <- 1.0
      } else {
        # Calculate overlap
        overlap_start <- min(windows[i, 1], windows[j, 1])
        overlap_end <- max(windows[i, 2], windows[j, 2])

        # Actual overlapping interval
        actual_start <- min(windows[i, 1], windows[j, 1])
        actual_end <- max(windows[i, 2], windows[j, 2])

        # Check if there is overlap
        if (windows[i, 2] < windows[j, 1] || windows[j, 2] < windows[i, 1]) {
          # No overlap
          overlap_matrix[i, j] <- 0.0
        } else {
          # There is overlap
          overlap_start <- min(windows[i, 1], windows[j, 1])
          overlap_end <- max(windows[i, 2], windows[j, 2])

          # Total span
          total_span <- overlap_start - overlap_end

          # Overlapping span
          overlap_span <- min(windows[i, 1], windows[j, 1]) - max(windows[i, 2], windows[j, 2])

          overlap_matrix[i, j] <- overlap_span / total_span
        }
      }
    }
  }

  return(overlap_matrix)
}

#' Print summary of simulated data
#'
#' @param sim_result Output from simulation functions
#' @export
print.simulated_c14 <- function(sim_result) {
  cat("=== Simulated Radiocarbon Dataset ===\n\n")
  cat("Overlap type:", sim_result$ground_truth$overlap_type, "\n")
  cat("Number of deposits:", sim_result$parameters$n_deposits, "\n")
  cat("Total dates:", nrow(sim_result$dates), "\n")
  cat("Measurement error:", sim_result$parameters$measurement_error, "14C years\n\n")

  cat("Deposit windows (cal BP):\n")
  windows <- sim_result$ground_truth$deposit_windows
  for (k in 1:nrow(windows)) {
    cat(sprintf("  Deposit %d: %d - %d\n", k, windows[k, 1], windows[k, 2]))
  }
  cat("\n")

  cat("Dates per deposit:\n")
  print(table(sim_result$dates$deposit))
  cat("\n")
}

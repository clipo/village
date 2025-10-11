#' Contemporaneity Metrics
#'
#' Functions to compute temporal overlap and contemporaneity metrics
#' from Stan model posterior samples.
#'
#' @author Project Team
#' @date 2025-10-11

library(posterior)
library(dplyr)

#' Calculate pairwise overlap probability
#'
#' Computes probability that two deposits overlap in time based on
#' posterior samples of occupation boundaries.
#'
#' @param fit Fitted Stan model object
#' @param deposit_i Index of first deposit
#' @param deposit_j Index of second deposit
#' @return Numeric, probability of overlap (0-1)
#' @export
calculate_overlap_probability <- function(fit, deposit_i, deposit_j) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  # Get samples for specific deposits
  start_i <- theta_start[, deposit_i]
  end_i <- theta_end[, deposit_i]
  start_j <- theta_start[, deposit_j]
  end_j <- theta_end[, deposit_j]

  # Calculate overlap for each posterior sample
  # Deposits overlap if: end_i < start_j AND end_j < start_i
  # Equivalently: max(end_i, end_j) < min(start_i, start_j)
  overlaps <- (pmax(end_i, end_j) < pmin(start_i, start_j))

  # Probability is proportion of samples with overlap
  prob_overlap <- mean(overlaps)

  return(prob_overlap)
}

#' Calculate overlap duration
#'
#' Computes expected duration of overlap between two deposits.
#'
#' @param fit Fitted Stan model object
#' @param deposit_i Index of first deposit
#' @param deposit_j Index of second deposit
#' @param probs Quantiles to compute (default c(0.025, 0.5, 0.975))
#' @return List with mean, median, and quantiles of overlap duration
#' @export
calculate_overlap_duration <- function(fit, deposit_i, deposit_j,
                                        probs = c(0.025, 0.5, 0.975)) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  # Get samples for specific deposits
  start_i <- theta_start[, deposit_i]
  end_i <- theta_end[, deposit_i]
  start_j <- theta_start[, deposit_j]
  end_j <- theta_end[, deposit_j]

  # Calculate overlap duration for each posterior sample
  overlap_start <- pmin(start_i, start_j)
  overlap_end <- pmax(end_i, end_j)

  # Duration is positive only if there is overlap
  overlap_duration <- ifelse(overlap_end < overlap_start,
                              overlap_start - overlap_end,
                              0)

  # Compute summary statistics
  result <- list(
    mean = mean(overlap_duration),
    median = median(overlap_duration),
    sd = sd(overlap_duration),
    quantiles = quantile(overlap_duration, probs = probs),
    prob_overlap = mean(overlap_duration > 0)
  )

  return(result)
}

#' Create overlap probability matrix
#'
#' Computes pairwise overlap probabilities for all deposits.
#'
#' @param fit Fitted Stan model object
#' @param n_deposits Number of deposits
#' @return Matrix of pairwise overlap probabilities
#' @export
create_overlap_matrix <- function(fit, n_deposits) {

  overlap_matrix <- matrix(1, n_deposits, n_deposits)

  for (i in 1:(n_deposits - 1)) {
    for (j in (i + 1):n_deposits) {
      prob <- calculate_overlap_probability(fit, i, j)
      overlap_matrix[i, j] <- prob
      overlap_matrix[j, i] <- prob
    }
  }

  return(overlap_matrix)
}

#' Calculate occupation span metrics
#'
#' Computes duration and temporal range for each deposit.
#'
#' @param fit Fitted Stan model object
#' @param n_deposits Number of deposits
#' @param probs Quantiles to compute (default c(0.025, 0.16, 0.5, 0.84, 0.975))
#' @return Data frame with occupation metrics for each deposit
#' @export
calculate_occupation_spans <- function(fit, n_deposits,
                                        probs = c(0.025, 0.16, 0.5, 0.84, 0.975)) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  results <- data.frame(
    deposit = 1:n_deposits,
    duration_mean = NA_real_,
    duration_sd = NA_real_,
    start_mean = NA_real_,
    start_median = NA_real_,
    end_mean = NA_real_,
    end_median = NA_real_
  )

  # Add quantile columns
  for (p in probs) {
    results[[paste0("start_q", round(p * 100))]] <- NA_real_
    results[[paste0("end_q", round(p * 100))]] <- NA_real_
    results[[paste0("duration_q", round(p * 100))]] <- NA_real_
  }

  for (k in 1:n_deposits) {
    start_k <- theta_start[, k]
    end_k <- theta_end[, k]
    duration_k <- start_k - end_k

    results$duration_mean[k] <- mean(duration_k)
    results$duration_sd[k] <- sd(duration_k)
    results$start_mean[k] <- mean(start_k)
    results$start_median[k] <- median(start_k)
    results$end_mean[k] <- mean(end_k)
    results$end_median[k] <- median(end_k)

    # Quantiles
    start_q <- quantile(start_k, probs = probs)
    end_q <- quantile(end_k, probs = probs)
    duration_q <- quantile(duration_k, probs = probs)

    for (i in seq_along(probs)) {
      p_col <- paste0("q", round(probs[i] * 100))
      results[[paste0("start_", p_col)]][k] <- start_q[i]
      results[[paste0("end_", p_col)]][k] <- end_q[i]
      results[[paste0("duration_", p_col)]][k] <- duration_q[i]
    }
  }

  return(results)
}

#' Calculate multi-deposit contemporaneity
#'
#' Computes probability that all deposits overlap within a common window.
#'
#' @param fit Fitted Stan model object
#' @param n_deposits Number of deposits
#' @return List with contemporaneity probability and expected common window
#' @export
calculate_multi_deposit_contemporaneity <- function(fit, n_deposits) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  n_samples <- nrow(theta_start)

  # For each sample, find the common window (if any)
  common_window_exists <- logical(n_samples)
  common_start <- numeric(n_samples)
  common_end <- numeric(n_samples)
  common_duration <- numeric(n_samples)

  for (s in 1:n_samples) {
    # Common window is intersection of all deposit windows
    # Start is minimum of all starts
    # End is maximum of all ends
    window_start <- min(theta_start[s, ])
    window_end <- max(theta_end[s, ])

    # Check if all deposits overlap this window
    all_overlap <- all(theta_end[s, ] < window_start & theta_end[s, ] > window_end)

    common_window_exists[s] <- (window_end < window_start)
    common_start[s] <- window_start
    common_end[s] <- window_end
    common_duration[s] <- ifelse(common_window_exists[s],
                                   window_start - window_end,
                                   0)
  }

  result <- list(
    prob_all_overlap = mean(common_window_exists),
    common_duration_mean = mean(common_duration[common_window_exists]),
    common_duration_sd = sd(common_duration[common_window_exists]),
    common_window_start = mean(common_start[common_window_exists]),
    common_window_end = mean(common_end[common_window_exists])
  )

  return(result)
}

#' Calculate temporal ordering probability
#'
#' Computes probability of specific temporal ordering between deposits.
#'
#' @param fit Fitted Stan model object
#' @param ordering Integer vector specifying deposit order (earliest to latest)
#' @return Probability that deposits follow specified order
#' @export
calculate_ordering_probability <- function(fit, ordering) {

  n_deposits <- length(ordering)

  # Extract boundary samples
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  n_samples <- nrow(theta_end)

  # Check ordering for each sample
  correct_order <- logical(n_samples)

  for (s in 1:n_samples) {
    # Check if end dates follow specified order
    # (earlier deposits have larger end dates in BP)
    order_correct <- TRUE
    for (i in 1:(n_deposits - 1)) {
      dep_current <- ordering[i]
      dep_next <- ordering[i + 1]

      if (theta_end[s, dep_current] <= theta_end[s, dep_next]) {
        order_correct <- FALSE
        break
      }
    }
    correct_order[s] <- order_correct
  }

  prob_order <- mean(correct_order)

  return(prob_order)
}

#' Calculate total time span
#'
#' Computes the total time span encompassing all deposits.
#'
#' @param fit Fitted Stan model object
#' @param n_deposits Number of deposits
#' @param probs Quantiles to compute (default c(0.025, 0.5, 0.975))
#' @return List with time span statistics
#' @export
calculate_total_time_span <- function(fit, n_deposits,
                                       probs = c(0.025, 0.5, 0.975)) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  # Total span from earliest start to latest end
  earliest_start <- apply(theta_start, 1, max)
  latest_end <- apply(theta_end, 1, min)
  total_span <- earliest_start - latest_end

  result <- list(
    mean = mean(total_span),
    median = median(total_span),
    sd = sd(total_span),
    quantiles = quantile(total_span, probs = probs),
    earliest_start_mean = mean(earliest_start),
    latest_end_mean = mean(latest_end)
  )

  return(result)
}

#' Compute all contemporaneity metrics
#'
#' Convenience function that computes all relevant metrics.
#'
#' @param fit Fitted Stan model object
#' @param n_deposits Number of deposits
#' @return List with all metrics
#' @export
compute_all_metrics <- function(fit, n_deposits) {

  message("Computing contemporaneity metrics...")

  metrics <- list(
    overlap_matrix = create_overlap_matrix(fit, n_deposits),
    occupation_spans = calculate_occupation_spans(fit, n_deposits),
    multi_deposit = calculate_multi_deposit_contemporaneity(fit, n_deposits),
    total_span = calculate_total_time_span(fit, n_deposits)
  )

  message("Metrics computed successfully")

  return(metrics)
}

#' Print contemporaneity metrics summary
#'
#' @param metrics Output from compute_all_metrics()
#' @export
print_metrics_summary <- function(metrics) {

  cat("=== Contemporaneity Metrics Summary ===\n\n")

  cat("Pairwise Overlap Probabilities:\n")
  print(round(metrics$overlap_matrix, 3))
  cat("\n")

  cat("Occupation Spans:\n")
  print(metrics$occupation_spans[, c("deposit", "duration_mean", "start_median", "end_median")])
  cat("\n")

  cat("Multi-Deposit Contemporaneity:\n")
  cat("  Probability all deposits overlap: ",
      round(metrics$multi_deposit$prob_all_overlap, 3), "\n")
  cat("  Common window duration (mean): ",
      round(metrics$multi_deposit$common_duration_mean, 1), " years\n")
  cat("\n")

  cat("Total Time Span:\n")
  cat("  Mean: ", round(metrics$total_span$mean, 1), " years\n")
  cat("  95% CI: [", round(metrics$total_span$quantiles[1], 1), ", ",
      round(metrics$total_span$quantiles[3], 1), "]\n")
  cat("\n")
}

#' Hypothesis Tests for Temporal Relationships
#'
#' Functions to test specific hypotheses about site contemporaneity
#' and temporal ordering.
#'
#' @author Project Team
#' @date 2025-10-11

library(posterior)
library(dplyr)

#' Test directional precedence
#'
#' Tests whether site i definitively precedes site j (no overlap).
#'
#' @param fit Fitted Stan model object
#' @param site_i Index of potentially earlier site
#' @param site_j Index of potentially later site
#' @return List with probability of precedence and gap duration statistics
#' @export
test_precedence <- function(fit, site_i, site_j) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  # Site i precedes site j if end_i > start_j (BP scale, older/larger is earlier)
  start_i <- theta_start[, site_i]
  end_i <- theta_end[, site_i]
  start_j <- theta_start[, site_j]
  end_j <- theta_end[, site_j]

  # Check if i completely precedes j (i ends before j starts)
  precedes <- end_i > start_j

  # Calculate gap duration when precedence holds
  gap_duration <- ifelse(precedes, end_i - start_j, 0)

  result <- list(
    prob_precedes = mean(precedes),
    gap_mean = mean(gap_duration[precedes]),
    gap_median = median(gap_duration[precedes]),
    gap_sd = sd(gap_duration[precedes]),
    gap_95ci = quantile(gap_duration[precedes], c(0.025, 0.975))
  )

  return(result)
}

#' Test if overlap exceeds threshold
#'
#' Tests whether overlap duration exceeds a specified threshold.
#'
#' @param fit Fitted Stan model object
#' @param site_i Index of first site
#' @param site_j Index of second site
#' @param threshold Minimum overlap duration in years
#' @return List with probability and overlap statistics
#' @export
test_overlap_threshold <- function(fit, site_i, site_j, threshold = 100) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  start_i <- theta_start[, site_i]
  end_i <- theta_end[, site_i]
  start_j <- theta_start[, site_j]
  end_j <- theta_end[, site_j]

  # Calculate overlap duration
  overlap_start <- pmin(start_i, start_j)
  overlap_end <- pmax(end_i, end_j)
  overlap_duration <- ifelse(overlap_end < overlap_start,
                              overlap_start - overlap_end,
                              0)

  # Test if overlap exceeds threshold
  exceeds_threshold <- overlap_duration >= threshold

  result <- list(
    threshold = threshold,
    prob_exceeds = mean(exceeds_threshold),
    overlap_mean = mean(overlap_duration),
    overlap_median = median(overlap_duration),
    overlap_sd = sd(overlap_duration),
    overlap_95ci = quantile(overlap_duration, c(0.025, 0.975)),
    prob_any_overlap = mean(overlap_duration > 0)
  )

  return(result)
}

#' Test complete overlap (full contemporaneity)
#'
#' Tests whether one site's occupation is completely contained within another's.
#'
#' @param fit Fitted Stan model object
#' @param site_contained Index of potentially contained site
#' @param site_container Index of potentially containing site
#' @return List with probability and containment statistics
#' @export
test_complete_overlap <- function(fit, site_contained, site_container) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  start_contained <- theta_start[, site_contained]
  end_contained <- theta_end[, site_contained]
  start_container <- theta_start[, site_container]
  end_container <- theta_end[, site_container]

  # Complete overlap: contained site fully within container site
  # start_container >= start_contained (container starts earlier or same)
  # end_container <= end_contained (container ends later or same)
  complete_overlap <- (start_container >= start_contained) &
                      (end_container <= end_contained)

  result <- list(
    prob_complete_overlap = mean(complete_overlap),
    contained_duration_mean = mean(start_contained - end_contained),
    container_duration_mean = mean(start_container - end_container),
    overlap_fraction_mean = mean(
      (pmin(start_contained, start_container) - pmax(end_contained, end_container)) /
      (start_contained - end_contained)
    )
  )

  return(result)
}

#' Compute pairwise overlap durations
#'
#' Calculates overlap duration in years for all site pairs.
#'
#' @param fit Fitted Stan model object
#' @param site_names Character vector of site names
#' @return Data frame with pairwise overlap durations
#' @export
compute_pairwise_overlap_durations <- function(fit, site_names) {

  n_sites <- length(site_names)
  results <- data.frame()

  for (i in 1:(n_sites - 1)) {
    for (j in (i + 1):n_sites) {

      # Extract boundary samples
      theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
      theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

      start_i <- theta_start[, i]
      end_i <- theta_end[, i]
      start_j <- theta_start[, j]
      end_j <- theta_end[, j]

      # Calculate overlap
      overlap_start <- pmin(start_i, start_j)
      overlap_end <- pmax(end_i, end_j)
      overlap_duration <- ifelse(overlap_end < overlap_start,
                                  overlap_start - overlap_end,
                                  0)

      row <- data.frame(
        site_1 = site_names[i],
        site_2 = site_names[j],
        overlap_mean = mean(overlap_duration),
        overlap_median = median(overlap_duration),
        overlap_sd = sd(overlap_duration),
        overlap_95ci_lower = quantile(overlap_duration, 0.025),
        overlap_95ci_upper = quantile(overlap_duration, 0.975),
        prob_any_overlap = mean(overlap_duration > 0),
        prob_overlap_gt_50yr = mean(overlap_duration >= 50),
        prob_overlap_gt_100yr = mean(overlap_duration >= 100),
        prob_overlap_gt_200yr = mean(overlap_duration >= 200)
      )

      results <- rbind(results, row)
    }
  }

  return(results)
}

#' Test all pairwise precedence relationships
#'
#' Tests directional precedence for all site pairs.
#'
#' @param fit Fitted Stan model object
#' @param site_names Character vector of site names
#' @return Data frame with precedence probabilities
#' @export
test_all_precedence <- function(fit, site_names) {

  n_sites <- length(site_names)
  results <- data.frame()

  for (i in 1:n_sites) {
    for (j in 1:n_sites) {
      if (i == j) next

      prec <- test_precedence(fit, i, j)

      row <- data.frame(
        earlier_site = site_names[i],
        later_site = site_names[j],
        prob_precedence = prec$prob_precedes,
        gap_mean = ifelse(is.na(prec$gap_mean), 0, prec$gap_mean),
        gap_median = ifelse(is.na(prec$gap_median), 0, prec$gap_median),
        gap_95ci_lower = ifelse(is.na(prec$gap_95ci[1]), 0, prec$gap_95ci[1]),
        gap_95ci_upper = ifelse(is.na(prec$gap_95ci[2]), 0, prec$gap_95ci[2])
      )

      results <- rbind(results, row)
    }
  }

  return(results)
}

#' Test multiple ordering hypotheses
#'
#' Tests several possible temporal orderings and compares their probabilities.
#'
#' @param fit Fitted Stan model object
#' @param site_names Character vector of site names
#' @param orderings List of integer vectors specifying different orderings
#' @return Data frame with ordering probabilities
#' @export
test_multiple_orderings <- function(fit, site_names, orderings) {

  results <- data.frame()

  for (i in seq_along(orderings)) {
    ordering <- orderings[[i]]
    prob <- calculate_ordering_probability(fit, ordering)

    ordering_str <- paste(site_names[ordering], collapse = " -> ")

    row <- data.frame(
      hypothesis = paste0("H", i),
      ordering = ordering_str,
      probability = prob
    )

    results <- rbind(results, row)
  }

  # Add relative support
  results$relative_support <- results$probability / sum(results$probability)

  return(results)
}

#' Comprehensive hypothesis testing
#'
#' Runs all hypothesis tests and generates summary report.
#'
#' @param fit Fitted Stan model object
#' @param site_names Character vector of site names
#' @param output_dir Directory for output files
#' @export
run_comprehensive_tests <- function(fit, site_names, output_dir = "output/hypothesis_tests") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Running comprehensive hypothesis tests...\n")

  # 1. Pairwise overlap durations
  cat("  Computing overlap durations...\n")
  overlap_durations <- compute_pairwise_overlap_durations(fit, site_names)
  write.csv(overlap_durations,
            file.path(output_dir, "pairwise_overlap_durations.csv"),
            row.names = FALSE)

  # 2. Precedence relationships
  cat("  Testing precedence relationships...\n")
  precedence <- test_all_precedence(fit, site_names)
  write.csv(precedence,
            file.path(output_dir, "precedence_tests.csv"),
            row.names = FALSE)

  # 3. Complete overlap tests
  cat("  Testing complete overlap...\n")
  n_sites <- length(site_names)
  complete_overlap_results <- data.frame()

  for (i in 1:n_sites) {
    for (j in 1:n_sites) {
      if (i == j) next

      co <- test_complete_overlap(fit, i, j)

      row <- data.frame(
        contained_site = site_names[i],
        container_site = site_names[j],
        prob_complete_overlap = co$prob_complete_overlap,
        overlap_fraction = co$overlap_fraction_mean
      )

      complete_overlap_results <- rbind(complete_overlap_results, row)
    }
  }

  write.csv(complete_overlap_results,
            file.path(output_dir, "complete_overlap_tests.csv"),
            row.names = FALSE)

  cat("  Tests complete. Results saved to:", output_dir, "\n")

  return(list(
    overlap_durations = overlap_durations,
    precedence = precedence,
    complete_overlap = complete_overlap_results
  ))
}

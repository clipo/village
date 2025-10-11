#!/usr/bin/env Rscript
#' Run Comprehensive Hypothesis Tests
#'
#' This script runs detailed hypothesis tests on the fitted model,
#' including overlap durations, precedence relationships, and
#' ordering hypotheses.

library(cmdstanr)
library(dplyr)

# Source required functions
source("R/setup.R")
source("R/stan_utilities.R")
source("R/contemporaneity_metrics.R")
source("R/hypothesis_tests.R")

cat("=======================================================\n")
cat("Running Comprehensive Hypothesis Tests\n")
cat("=======================================================\n\n")

# Load the fitted model
cat("Loading fitted model...\n")
fit <- readRDS("output/sample_analysis/fit_partial_overlap.rds")

# Define site names (in order they appear in the model)
site_names <- c("Broome_Tech", "CCE_Site", "Chenango_Point", "Hill_Creek", "Thomas_Luckey")

cat("Sites analyzed:\n")
for (i in seq_along(site_names)) {
  cat(sprintf("  %d: %s\n", i, site_names[i]))
}
cat("\n")

# Run comprehensive tests
cat("=======================================================\n")
cat("HYPOTHESIS TESTS\n")
cat("=======================================================\n\n")

results <- run_comprehensive_tests(fit, site_names, output_dir = "output/hypothesis_tests")

# Print summaries
cat("\n=======================================================\n")
cat("OVERLAP DURATIONS (years)\n")
cat("=======================================================\n\n")

overlap_summary <- results$overlap_durations %>%
  select(site_1, site_2, overlap_mean, overlap_95ci_lower, overlap_95ci_upper, prob_any_overlap) %>%
  mutate(
    overlap_mean = round(overlap_mean, 0),
    overlap_95ci_lower = round(overlap_95ci_lower, 0),
    overlap_95ci_upper = round(overlap_95ci_upper, 0),
    prob_any_overlap = round(prob_any_overlap, 3)
  )

print(overlap_summary)

cat("\n=======================================================\n")
cat("PRECEDENCE RELATIONSHIPS\n")
cat("=======================================================\n")
cat("(Probability that earlier_site completely precedes later_site)\n\n")

precedence_summary <- results$precedence %>%
  filter(prob_precedence > 0.5) %>%
  arrange(desc(prob_precedence)) %>%
  select(earlier_site, later_site, prob_precedence, gap_mean, gap_95ci_lower, gap_95ci_upper) %>%
  mutate(
    prob_precedence = round(prob_precedence, 3),
    gap_mean = round(gap_mean, 0),
    gap_95ci_lower = round(gap_95ci_lower, 0),
    gap_95ci_upper = round(gap_95ci_upper, 0)
  )

if (nrow(precedence_summary) > 0) {
  print(precedence_summary)
} else {
  cat("No strong precedence relationships found (all probabilities < 0.5)\n")
}

cat("\n=======================================================\n")
cat("COMPLETE OVERLAP\n")
cat("=======================================================\n")
cat("(Probability that contained_site is fully within container_site)\n\n")

complete_overlap_summary <- results$complete_overlap %>%
  filter(prob_complete_overlap > 0.1) %>%
  arrange(desc(prob_complete_overlap)) %>%
  mutate(
    prob_complete_overlap = round(prob_complete_overlap, 3),
    overlap_fraction = round(overlap_fraction, 3)
  )

if (nrow(complete_overlap_summary) > 0) {
  print(complete_overlap_summary)
} else {
  cat("No strong complete overlap relationships found (all probabilities < 0.1)\n")
}

cat("\n=======================================================\n")
cat("SPECIFIC HYPOTHESES\n")
cat("=======================================================\n\n")

# Test specific ordering hypotheses
cat("Testing temporal ordering hypotheses:\n\n")

# Define hypotheses to test
# Indices: 1=Broome_Tech, 2=CCE_Site, 3=Chenango_Point, 4=Hill_Creek, 5=Thomas_Luckey

hypotheses <- list(
  c(4, 5, 1, 3, 2),  # H1: Hill_Creek -> Thomas_Luckey -> Broome_Tech -> Chenango_Point -> CCE_Site
  c(4, 1, 3, 5, 2),  # H2: Hill_Creek -> Broome_Tech -> Chenango_Point -> Thomas_Luckey -> CCE_Site
  c(4, 5, 3, 1, 2),  # H3: Hill_Creek -> Thomas_Luckey -> Chenango_Point -> Broome_Tech -> CCE_Site
  c(5, 4, 1, 3, 2),  # H4: Thomas_Luckey -> Hill_Creek -> Broome_Tech -> Chenango_Point -> CCE_Site
  c(1, 3, 5, 4, 2)   # H5: Broome_Tech -> Chenango_Point -> Thomas_Luckey -> Hill_Creek -> CCE_Site
)

ordering_results <- test_multiple_orderings(fit, site_names, hypotheses)
print(ordering_results)

cat("\n=======================================================\n")
cat("OVERLAP THRESHOLDS\n")
cat("=======================================================\n")
cat("Probability that overlap exceeds specific durations:\n\n")

# Test overlap thresholds for key site pairs
key_pairs <- list(
  c("Broome_Tech", "Thomas_Luckey", 1, 5),
  c("Broome_Tech", "Chenango_Point", 1, 3),
  c("Thomas_Luckey", "Chenango_Point", 5, 3),
  c("Hill_Creek", "CCE_Site", 4, 2)
)

threshold_results <- data.frame()

for (pair in key_pairs) {
  site_1_name <- pair[1]
  site_2_name <- pair[2]
  site_1_idx <- as.numeric(pair[3])
  site_2_idx <- as.numeric(pair[4])

  # Test 100-year threshold
  test_100 <- test_overlap_threshold(fit, site_1_idx, site_2_idx, threshold = 100)
  # Test 200-year threshold
  test_200 <- test_overlap_threshold(fit, site_1_idx, site_2_idx, threshold = 200)

  row <- data.frame(
    site_pair = paste(site_1_name, "-", site_2_name),
    prob_overlap_gt_100yr = round(test_100$prob_exceeds, 3),
    prob_overlap_gt_200yr = round(test_200$prob_exceeds, 3),
    mean_overlap = round(test_100$overlap_mean, 0)
  )

  threshold_results <- rbind(threshold_results, row)
}

print(threshold_results)

# Save ordering results
write.csv(ordering_results,
          "output/hypothesis_tests/ordering_hypotheses.csv",
          row.names = FALSE)

write.csv(threshold_results,
          "output/hypothesis_tests/overlap_thresholds.csv",
          row.names = FALSE)

cat("\n=======================================================\n")
cat("All results saved to: output/hypothesis_tests/\n")
cat("=======================================================\n")

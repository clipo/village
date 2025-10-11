#' Example Workflow: Complete Analysis Pipeline
#'
#' This script demonstrates a complete analysis workflow from data input
#' through model fitting, comparison, and interpretation.
#'
#' @author Project Team
#' @date 2025-10-11

# ==============================================================================
# SETUP
# ==============================================================================

# Load all required packages and functions
library(rcarbon)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(dplyr)
library(ggplot2)

# Source project functions
source("R/setup.R")
source("R/data_validation.R")
source("R/calibration.R")
source("R/simulation.R")
source("R/stan_utilities.R")
source("R/contemporaneity_metrics.R")
source("R/plotting.R")

# Verify setup
cat("=== Verifying Setup ===\n")
verify_stan_installation()

# Set random seed for reproducibility
set.seed(42)

# ==============================================================================
# EXAMPLE 1: SIMULATED DATA WITH KNOWN STRUCTURE
# ==============================================================================

cat("\n=== EXAMPLE 1: Simulated Partially Overlapping Deposits ===\n\n")

# Generate simulated data with 60% overlap
sim_data <- simulate_partial_overlap(
  n_deposits = 3,
  overlap_percentage = 60,
  deposit_duration = 250,
  earliest_start = 2800,
  n_dates_per_deposit = 10,
  measurement_error = 30,
  seed = 123
)

# Print simulation details
print(sim_data)

# Extract radiocarbon data
c14_data <- sim_data$dates[, c("lab_code", "age", "error", "deposit")]

# Validate data
cat("\n--- Data Validation ---\n")
validate_c14_data(c14_data)
print_c14_summary(c14_data)

# Calibrate
cat("\n--- Calibration ---\n")
cal_result <- calibrate_dates(
  ages = c14_data$age,
  errors = c14_data$error,
  lab_codes = c14_data$lab_code,
  calCurve = "intcal20",
  verbose = TRUE
)

# Plot calibrated distributions
cat("\n--- Plotting Calibrated Dates ---\n")
p1 <- plot_calibrated_dates(cal_result, c14_data)
print(p1)
ggsave("output/example1_calibrated_dates.png", p1, width = 10, height = 8)

# Prepare Stan data
cat("\n--- Preparing Stan Data ---\n")
stan_data <- prepare_stan_data(c14_data, cal_result)

# Compile models
cat("\n--- Compiling Stan Models ---\n")
model_contemp <- compile_stan_model("stan/contemporaneous_model.stan")
model_seq <- compile_stan_model("stan/sequential_model.stan")
model_partial <- compile_stan_model("stan/partial_overlap_model.stan")

# Fit all three models
cat("\n--- Fitting Models ---\n")

cat("\n1. Contemporaneous Model...\n")
fit_contemp <- fit_stan_model(
  model_contemp,
  stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 42
)

cat("\n2. Sequential Model...\n")
fit_seq <- fit_stan_model(
  model_seq,
  stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 42
)

cat("\n3. Partial Overlap Model...\n")
fit_partial <- fit_stan_model(
  model_partial,
  stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 42
)

# Compare models
cat("\n--- Model Comparison ---\n")
comparison <- compare_models_loo(
  contemporaneous = fit_contemp,
  sequential = fit_seq,
  partial = fit_partial,
  cores = 2
)

# Plot comparison
p2 <- plot_model_comparison(comparison)
print(p2)
ggsave("output/example1_model_comparison.png", p2, width = 8, height = 6)

# Identify best model
best_model <- names(comparison$weights)[which.max(comparison$weights)]
cat("\nBest model:", best_model, "\n")
cat("Model weight:", round(max(comparison$weights), 3), "\n")

# Analyze results using best model (typically partial overlap)
cat("\n--- Analyzing Results ---\n")

# Calculate all metrics
metrics <- compute_all_metrics(fit_partial, n_deposits = 3)
print_metrics_summary(metrics)

# Create diagnostic plots
cat("\n--- Creating Plots ---\n")
p3 <- plot_occupation_boundaries(fit_partial, n_deposits = 3)
print(p3)
ggsave("output/example1_boundaries.png", p3, width = 10, height = 8)

p4 <- plot_temporal_overlap(fit_partial, n_deposits = 3)
print(p4)
ggsave("output/example1_temporal_overlap.png", p4, width = 10, height = 6)

p5 <- plot_overlap_matrix(metrics$overlap_matrix)
print(p5)
ggsave("output/example1_overlap_matrix.png", p5, width = 8, height = 6)

p6 <- plot_posterior_predictive(fit_partial, c14_data)
print(p6)
ggsave("output/example1_ppc.png", p6, width = 10, height = 6)

# Compare to ground truth
cat("\n--- Validation Against Ground Truth ---\n")
true_windows <- sim_data$ground_truth$deposit_windows
est_spans <- calculate_occupation_spans(fit_partial, n_deposits = 3)

cat("\nTrue vs. Estimated Boundaries (cal BP):\n")
for (k in 1:3) {
  cat(sprintf("Deposit %d:\n", k))
  cat(sprintf("  True:      Start=%d, End=%d\n",
              round(true_windows[k, 1]), round(true_windows[k, 2])))
  cat(sprintf("  Estimated: Start=%d (95%% CI: %d-%d), End=%d (95%% CI: %d-%d)\n",
              round(est_spans$start_median[k]),
              round(est_spans$start_q2.5[k]),
              round(est_spans$start_q97.5[k]),
              round(est_spans$end_median[k]),
              round(est_spans$end_q2.5[k]),
              round(est_spans$end_q97.5[k])))
}

# ==============================================================================
# EXAMPLE 2: USER-PROVIDED DATA (Template)
# ==============================================================================

cat("\n\n=== EXAMPLE 2: Template for User Data ===\n\n")

# This is a template for analyzing your own data
# Replace with your actual data

# Create data frame with your radiocarbon dates
# my_data <- data.frame(
#   lab_code = c("OxA-12345", "OxA-12346", "Beta-54321", "Beta-54322"),
#   age = c(2450, 2480, 2520, 2550),
#   error = c(30, 35, 40, 30),
#   deposit = c("Deposit_A", "Deposit_A", "Deposit_B", "Deposit_B"),
#   material = c("charcoal", "charcoal", "bone", "bone"),  # optional
#   stringsAsFactors = FALSE
# )

# Then follow the same workflow:
# 1. validate_c14_data(my_data)
# 2. cal_result <- calibrate_dates(my_data$age, my_data$error, my_data$lab_code)
# 3. stan_data <- prepare_stan_data(my_data, cal_result)
# 4. Compile and fit models
# 5. Compare models
# 6. Analyze and interpret results

cat("To analyze your own data:\n")
cat("1. Create a data frame with columns: lab_code, age, error, deposit\n")
cat("2. Follow the workflow demonstrated in Example 1\n")
cat("3. Interpret results in archaeological context\n")

# ==============================================================================
# EXAMPLE 3: QUICK VALIDATION STUDY
# ==============================================================================

cat("\n\n=== EXAMPLE 3: Quick Validation Study ===\n\n")

cat("Running mini validation study (10 trials per condition)...\n")
cat("(For full validation, increase n_trials to 100)\n\n")

# Run small validation study
validation_results <- run_validation_study(
  n_trials = 10,  # Increase to 100 for full validation
  scenarios = c("contemporaneous", "partial"),
  sample_sizes = c(5, 10),
  n_deposits_range = c(2),
  measurement_errors = c(30),
  parallel = FALSE  # Set to TRUE with n_cores = 4 for speed
)

# Summarize
summary_df <- summarize_validation_results(validation_results)
print_validation_summary(summary_df)

# Plot
plots <- plot_validation_study(validation_results)
print(plots$accuracy)
ggsave("output/example3_validation_accuracy.png", plots$accuracy,
       width = 8, height = 6)

print(plots$coverage)
ggsave("output/example3_validation_coverage.png", plots$coverage,
       width = 8, height = 6)

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n\n=== WORKFLOW COMPLETE ===\n\n")
cat("All outputs saved to output/ directory\n")
cat("Generated files:\n")
cat("  - example1_calibrated_dates.png\n")
cat("  - example1_model_comparison.png\n")
cat("  - example1_boundaries.png\n")
cat("  - example1_temporal_overlap.png\n")
cat("  - example1_overlap_matrix.png\n")
cat("  - example1_ppc.png\n")
cat("  - example3_validation_accuracy.png\n")
cat("  - example3_validation_coverage.png\n")

cat("\nNext steps:\n")
cat("1. Review plots and model diagnostics\n")
cat("2. Interpret contemporaneity probabilities in archaeological context\n")
cat("3. Consider sensitivity analyses with different priors\n")
cat("4. For publication, run full validation study\n")
cat("5. Render full report: quarto render analysis.qmd\n")

cat("\nSession Info:\n")
print(sessionInfo())

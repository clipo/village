#' Analyze Your Radiocarbon Data
#'
#' Complete workflow for analyzing the radiocarbon dates in data/ directory.
#' This script loads your CSV file, groups by site, and analyzes contemporaneity.
#'
#' @author Project Team
#' @date 2025-10-11

# ==============================================================================
# SETUP
# ==============================================================================

cat("==========================================================\n")
cat("  Radiocarbon Dating Contemporaneity Analysis\n")
cat("  Your Data Analysis\n")
cat("==========================================================\n\n")

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
source("R/data_loading.R")
source("R/data_validation.R")
source("R/calibration.R")
source("R/stan_utilities.R")
source("R/contemporaneity_metrics.R")
source("R/plotting.R")

# Set random seed for reproducibility
set.seed(42)

# Create output directory if needed
if (!dir.exists("output/your_data")) {
  dir.create("output/your_data", recursive = TRUE)
}

# ==============================================================================
# LOAD AND EXPLORE DATA
# ==============================================================================

cat("\n=== Loading Data ===\n\n")

# Load your data file
# The function will automatically map columns:
#   site -> deposit
#   lab_no -> lab_code
#   c14_age -> age
#   c14_error -> error
data <- load_radiocarbon_data(
  "data/radiocarbon_dates.csv",
  deposit_col = "site"  # Use "site" column as deposit identifier
)

# Print detailed summary
print_data_summary(data)

# Check for outliers
cat("\n=== Checking for Outliers ===\n")
outliers <- check_outliers(data, sd_threshold = 3)

if (any(outliers$is_outlier, na.rm = TRUE)) {
  cat("\nPotential outliers detected:\n")
  print(outliers[outliers$is_outlier, c("lab_code", "deposit", "age", "z_score")])
  cat("\nReview these dates before proceeding.\n")
}

# ==============================================================================
# DATA FILTERING OPTIONS
# ==============================================================================

cat("\n=== Data Filtering ===\n\n")

# Option 1: Analyze all sites with enough dates
cat("Option 1: Analyzing sites with at least 3 dates each\n\n")

# Count dates per site
site_counts <- table(data$deposit)
cat("Sites with sample sizes:\n")
print(sort(site_counts, decreasing = TRUE))

# Filter to sites with minimum sample size
min_dates <- 3
sufficient_sites <- names(site_counts)[site_counts >= min_dates]

cat("\nSites with >= ", min_dates, " dates: ", length(sufficient_sites), "\n")

if (length(sufficient_sites) < 2) {
  stop("Need at least 2 sites with ", min_dates, "+ dates each for analysis")
}

# Option 2: Select specific sites for focused analysis
# Uncomment and modify to analyze specific sites:
# selected_sites <- c("Site_A", "Site_B", "Site_C")
# filtered_data <- data[data$deposit %in% selected_sites, ]

# For this example, let's analyze the top N sites by sample size
n_sites_to_analyze <- min(5, length(sufficient_sites))  # Analyze up to 5 sites

top_sites <- names(sort(site_counts[sufficient_sites], decreasing = TRUE))[1:n_sites_to_analyze]

cat("\nAnalyzing top ", n_sites_to_analyze, " sites:\n")
for (site in top_sites) {
  cat("  ", site, ": ", site_counts[site], " dates\n")
}

filtered_data <- data[data$deposit %in% top_sites, ]
filtered_data$deposit <- factor(filtered_data$deposit)  # Drop unused levels

cat("\nFinal dataset: ", nrow(filtered_data), " dates from ",
    length(unique(filtered_data$deposit)), " sites\n")

# ==============================================================================
# CALIBRATION
# ==============================================================================

cat("\n=== Calibrating Radiocarbon Dates ===\n\n")

cal_result <- calibrate_dates(
  ages = filtered_data$age,
  errors = filtered_data$error,
  lab_codes = filtered_data$lab_code,
  calCurve = "intcal20",
  verbose = TRUE
)

print(cal_result)

# Plot calibrated distributions
cat("\n--- Creating Calibration Plots ---\n")
p1 <- plot_calibrated_dates(cal_result, filtered_data)
print(p1)
ggsave("output/your_data/01_calibrated_dates.png", p1, width = 12, height = 10, dpi = 300)

# ==============================================================================
# PREPARE STAN DATA
# ==============================================================================

cat("\n=== Preparing Data for Stan ===\n\n")

stan_data <- prepare_stan_data(filtered_data, cal_result)

cat("Stan data structure:\n")
cat("  Number of dates:", stan_data$N, "\n")
cat("  Number of sites:", stan_data$n_deposits, "\n")
cat("  Mixture components:", stan_data$K, "\n")

# ==============================================================================
# COMPILE STAN MODELS
# ==============================================================================

cat("\n=== Compiling Stan Models ===\n\n")

model_contemp <- compile_stan_model("stan/contemporaneous_model.stan")
model_seq <- compile_stan_model("stan/sequential_model.stan")
model_partial <- compile_stan_model("stan/partial_overlap_model.stan")

# ==============================================================================
# FIT MODELS
# ==============================================================================

cat("\n=== Fitting Bayesian Models ===\n\n")

# Generate initial values
init_contemp <- generate_inits(stan_data, model_type = "contemporaneous")
init_partial <- generate_inits(stan_data, model_type = "partial")

cat("Fitting 3 models (this may take several minutes)...\n\n")

# Fit contemporaneous model
cat("1. Contemporaneous Model (single shared window)...\n")
fit_contemp <- fit_stan_model(
  model_contemp,
  stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 42,
  init = init_contemp,
  verbose = TRUE
)

# Fit sequential model
cat("\n2. Sequential Model (ordered, non-overlapping)...\n")
fit_seq <- fit_stan_model(
  model_seq,
  stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 42,
  init = init_partial,
  verbose = TRUE
)

# Fit partial overlap model
cat("\n3. Partial Overlap Model (flexible)...\n")
fit_partial <- fit_stan_model(
  model_partial,
  stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 42,
  init = init_partial,
  verbose = TRUE
)

# ==============================================================================
# MODEL COMPARISON
# ==============================================================================

cat("\n=== Model Comparison ===\n\n")

comparison <- compare_models_loo(
  contemporaneous = fit_contemp,
  sequential = fit_seq,
  partial = fit_partial,
  cores = 2
)

# Plot comparison
p2 <- plot_model_comparison(comparison)
print(p2)
ggsave("output/your_data/02_model_comparison.png", p2, width = 8, height = 6, dpi = 300)

# Determine best model
best_model_name <- names(comparison$weights)[which.max(comparison$weights)]
best_weight <- max(comparison$weights)

cat("\n=== Model Selection Results ===\n")
cat("Best model:", best_model_name, "\n")
cat("Model weight:", round(best_weight, 3), "\n\n")

if (best_weight > 0.9) {
  cat("Interpretation: Strong evidence for", best_model_name, "model\n")
} else if (best_weight > 0.7) {
  cat("Interpretation: Moderate evidence for", best_model_name, "model\n")
} else {
  cat("Interpretation: Weak evidence; consider model averaging\n")
}

# Select best fit for detailed analysis
if (best_model_name == "contemporaneous") {
  best_fit <- fit_contemp
} else if (best_model_name == "sequential") {
  best_fit <- fit_seq
} else {
  best_fit <- fit_partial
}

# ==============================================================================
# ANALYZE RESULTS
# ==============================================================================

cat("\n=== Analyzing Results ===\n\n")

n_deposits <- length(unique(filtered_data$deposit))
deposit_names <- levels(factor(filtered_data$deposit))

# Calculate all metrics
metrics <- compute_all_metrics(best_fit, n_deposits)
print_metrics_summary(metrics)

# Save metrics to file
sink("output/your_data/metrics_summary.txt")
cat("=== CONTEMPORANEITY ANALYSIS RESULTS ===\n")
cat("Dataset: ", basename(attr(data, "original_file")), "\n")
cat("Analysis date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Best model: ", best_model_name, " (weight = ", round(best_weight, 3), ")\n\n")
print_metrics_summary(metrics)
sink()

# ==============================================================================
# CREATE VISUALIZATIONS
# ==============================================================================

cat("\n=== Creating Visualizations ===\n\n")

# Occupation boundaries
p3 <- plot_occupation_boundaries(best_fit, n_deposits, deposit_names)
print(p3)
ggsave("output/your_data/03_occupation_boundaries.png", p3, width = 10, height = 8, dpi = 300)

# Temporal overlap
p4 <- plot_temporal_overlap(best_fit, n_deposits, deposit_names)
print(p4)
ggsave("output/your_data/04_temporal_overlap.png", p4, width = 10, height = 6, dpi = 300)

# Overlap matrix
p5 <- plot_overlap_matrix(metrics$overlap_matrix, deposit_names)
print(p5)
ggsave("output/your_data/05_overlap_matrix.png", p5, width = 8, height = 6, dpi = 300)

# Posterior predictive check
p6 <- plot_posterior_predictive(best_fit, filtered_data)
print(p6)
ggsave("output/your_data/06_posterior_predictive.png", p6, width = 10, height = 6, dpi = 300)

# Convergence diagnostics
cat("\n--- Creating Diagnostic Plots ---\n")
conv_plots <- plot_convergence_diagnostics(best_fit)
ggsave("output/your_data/07_trace_plots.png", conv_plots$trace, width = 10, height = 6, dpi = 300)
ggsave("output/your_data/08_rank_plots.png", conv_plots$rank, width = 10, height = 6, dpi = 300)

# ==============================================================================
# EXPORT RESULTS
# ==============================================================================

cat("\n=== Exporting Results ===\n\n")

# Export occupation span estimates
spans <- metrics$occupation_spans
spans$site <- deposit_names
write.csv(spans, "output/your_data/occupation_spans.csv", row.names = FALSE)

# Export overlap matrix
overlap_df <- as.data.frame(metrics$overlap_matrix)
colnames(overlap_df) <- deposit_names
overlap_df$site <- deposit_names
write.csv(overlap_df, "output/your_data/overlap_matrix.csv", row.names = FALSE)

# Export calibrated summaries
cal_export <- cal_result$summary
cal_export$site <- filtered_data$deposit
cal_export$material <- filtered_data$material
write.csv(cal_export, "output/your_data/calibrated_dates.csv", row.names = FALSE)

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("==========================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("==========================================================\n\n")

cat("Results saved to output/your_data/\n\n")

cat("Files created:\n")
cat("  01_calibrated_dates.png       - Calibrated radiocarbon distributions\n")
cat("  02_model_comparison.png       - Model comparison via LOO\n")
cat("  03_occupation_boundaries.png  - Posterior distributions of boundaries\n")
cat("  04_temporal_overlap.png       - Temporal overlap with uncertainty\n")
cat("  05_overlap_matrix.png         - Pairwise overlap probabilities\n")
cat("  06_posterior_predictive.png   - Model fit check\n")
cat("  07_trace_plots.png           - MCMC convergence traces\n")
cat("  08_rank_plots.png            - MCMC rank plots\n")
cat("  metrics_summary.txt          - Text summary of results\n")
cat("  occupation_spans.csv         - Occupation estimates (CSV)\n")
cat("  overlap_matrix.csv           - Overlap probabilities (CSV)\n")
cat("  calibrated_dates.csv         - Calibration results (CSV)\n\n")

cat("Key Findings:\n")
cat("  Analyzed ", nrow(filtered_data), " dates from ", n_deposits, " sites\n")
cat("  Best model: ", best_model_name, "\n")
cat("  Model weight: ", round(best_weight, 3), "\n")

if (n_deposits == 2) {
  cat("  Overlap probability: ", round(metrics$overlap_matrix[1, 2], 3), "\n")
}

cat("\nNext steps:\n")
cat("  1. Review all plots in output/your_data/\n")
cat("  2. Check convergence diagnostics (trace and rank plots)\n")
cat("  3. Examine overlap matrix for contemporaneity patterns\n")
cat("  4. Consider archaeological context when interpreting results\n")
cat("  5. For publication, render full report: quarto render analysis.qmd\n\n")

cat("Session Info:\n")
print(sessionInfo())

cat("\n==========================================================\n")

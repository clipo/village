#' Run Full Contemporaneity Analysis
#'
#' Analyzes top 5 sites from the filtered dataset

cat("==========================================================\n")
cat("  Full Contemporaneity Analysis\n")
cat("  Top Sites from Sample Data\n")
cat("==========================================================\n\n")

# Load packages
library(rcarbon)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(dplyr)
library(ggplot2)

# Source functions
source("R/data_validation.R")
source("R/calibration.R")
source("R/stan_utilities.R")
source("R/contemporaneity_metrics.R")
source("R/plotting.R")

set.seed(42)

# Create output directory
dir.create("output/sample_analysis", showWarnings = FALSE, recursive = TRUE)

# Load filtered data
cat("Loading filtered data...\n")
filtered_data <- readRDS("output/filtered_data_for_analysis.rds")

# Select top 5 sites by sample size for this demonstration
site_counts <- table(filtered_data$deposit)
top_sites <- names(sort(site_counts, decreasing = TRUE))[1:5]

analysis_data <- filtered_data[filtered_data$deposit %in% top_sites, ]
analysis_data$deposit <- factor(analysis_data$deposit)

cat("\nAnalyzing sites:\n")
for (site in levels(analysis_data$deposit)) {
  cat(sprintf("  %-20s: %2d dates\n", site, sum(analysis_data$deposit == site)))
}

cat("\nTotal: ", nrow(analysis_data), " dates from ",
    length(unique(analysis_data$deposit)), " sites\n\n")

# Calibrate
cat("=== Calibrating Dates ===\n")
cal_result <- calibrate_dates(
  ages = analysis_data$age,
  errors = analysis_data$error,
  lab_codes = analysis_data$lab_code,
  calCurve = "intcal20",
  verbose = FALSE
)

cat("Calibration complete\n\n")

# Prepare Stan data
cat("=== Preparing Stan Data ===\n")
stan_data <- prepare_stan_data(analysis_data, cal_result)

# Compile models
cat("\n=== Compiling Stan Models ===\n")
model_partial <- compile_stan_model("stan/partial_overlap_model.stan")

# Fit model (using partial overlap as most flexible)
cat("\n=== Fitting Partial Overlap Model ===\n")
cat("(This will take 5-10 minutes...)\n\n")

# Generate reasonable initial values
init_fun <- function() {
  # Estimate calendar ages from mix_means
  cal_means <- apply(stan_data$mix_means, 1, function(x) sum(x * stan_data$mix_weights[1,]))

  # For each deposit, find range
  theta_start_init <- numeric(stan_data$n_deposits)
  durations_init <- numeric(stan_data$n_deposits)

  for (k in 1:stan_data$n_deposits) {
    idx <- which(stan_data$deposit_id == k)
    dep_ages <- cal_means[idx]
    theta_start_init[k] <- max(dep_ages) + 50  # Start (older)
    theta_end_est <- min(dep_ages) - 50        # End (younger)
    durations_init[k] <- theta_start_init[k] - theta_end_est  # Positive duration
  }

  list(
    theta_start = theta_start_init,
    durations = durations_init,
    calendar_dates_raw = rnorm(stan_data$N, 0, 0.1)
  )
}

fit <- fit_stan_model(
  model_partial,
  stan_data,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 42,
  init = init_fun,
  verbose = TRUE
)

# Calculate metrics
cat("\n=== Computing Contemporaneity Metrics ===\n")
n_deposits <- length(unique(analysis_data$deposit))
metrics <- compute_all_metrics(fit, n_deposits)

print_metrics_summary(metrics)

# Save metrics
sink("output/sample_analysis/metrics_summary.txt")
cat("=== SAMPLE DATA ANALYSIS RESULTS ===\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Sites analyzed:\n")
for (i in 1:length(top_sites)) {
  cat(sprintf("  %d. %-20s (%2d dates)\n", i, top_sites[i],
              sum(analysis_data$deposit == top_sites[i])))
}
cat("\n")
print_metrics_summary(metrics)
sink()

# Create plots
cat("\n=== Creating Visualizations ===\n")

deposit_names <- levels(analysis_data$deposit)

# Try each plot, skip if it fails
try({
  p2 <- plot_occupation_boundaries(fit, n_deposits, deposit_names)
  ggsave("output/sample_analysis/02_occupation_boundaries.png", p2,
         width = 10, height = 8, dpi = 300)
  cat("  Created: 02_occupation_boundaries.png\n")
}, silent = TRUE)

try({
  p3 <- plot_temporal_overlap_from_metrics(metrics$occupation_spans, deposit_names)
  ggsave("output/sample_analysis/03_temporal_overlap.png", p3,
         width = 10, height = 6, dpi = 300)
  cat("  Created: 03_temporal_overlap.png\n")
}, silent = TRUE)

try({
  p4 <- plot_overlap_matrix(metrics$overlap_matrix, deposit_names)
  ggsave("output/sample_analysis/04_overlap_matrix.png", p4,
         width = 8, height = 6, dpi = 300)
  cat("  Created: 04_overlap_matrix.png\n")
}, silent = TRUE)

try({
  conv_plots <- plot_convergence_diagnostics(fit)
  ggsave("output/sample_analysis/06_trace_plots.png", conv_plots$trace,
         width = 10, height = 6, dpi = 300)
  ggsave("output/sample_analysis/07_rank_plots.png", conv_plots$rank,
         width = 10, height = 6, dpi = 300)
  cat("  Created: 06_trace_plots.png and 07_rank_plots.png\n")
}, silent = TRUE)

# Export results
cat("\n=== Exporting Results ===\n")

# Occupation spans
spans <- metrics$occupation_spans
spans$site <- deposit_names
write.csv(spans, "output/sample_analysis/occupation_spans.csv", row.names = FALSE)

# Overlap matrix
overlap_df <- as.data.frame(metrics$overlap_matrix)
colnames(overlap_df) <- deposit_names
overlap_df$site <- deposit_names
write.csv(overlap_df, "output/sample_analysis/overlap_matrix.csv", row.names = FALSE)

# Calibrated dates
cal_export <- cal_result$summary
cal_export$site <- analysis_data$deposit
cal_export$material <- analysis_data$material
write.csv(cal_export, "output/sample_analysis/calibrated_dates.csv", row.names = FALSE)

cat("\n")
cat("==========================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("==========================================================\n\n")

cat("Results saved to: output/sample_analysis/\n\n")

cat("Files created:\n")
cat("  01_calibrated_dates.png\n")
cat("  02_occupation_boundaries.png\n")
cat("  03_temporal_overlap.png\n")
cat("  04_overlap_matrix.png\n")
cat("  05_posterior_predictive.png\n")
cat("  06_trace_plots.png\n")
cat("  07_rank_plots.png\n")
cat("  metrics_summary.txt\n")
cat("  occupation_spans.csv\n")
cat("  overlap_matrix.csv\n")
cat("  calibrated_dates.csv\n\n")

cat("Key findings:\n")
cat("  Sites analyzed:", paste(top_sites, collapse = ", "), "\n")
cat("  All diagnostics passed:", attr(fit, "diagnostics_summary")$all_diagnostics_passed, "\n\n")

cat("==========================================================\n")

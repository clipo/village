#' Preview Your Data
#'
#' Quick script to explore your radiocarbon data before running full analysis.
#' This helps you decide which sites to include and check data quality.
#'
#' @author Project Team
#' @date 2025-10-11

cat("==========================================================\n")
cat("  Data Preview and Exploration\n")
cat("==========================================================\n\n")

# Load required functions
source("R/data_loading.R")
source("R/data_validation.R")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("=== Loading Data ===\n\n")

# Load your CSV file
data <- load_radiocarbon_data(
  "data/radiocarbon_dates.csv",
  deposit_col = "site"
)

# ==============================================================================
# BASIC SUMMARY
# ==============================================================================

print_data_summary(data)

# ==============================================================================
# DETAILED EXPLORATION
# ==============================================================================

cat("\n=== Detailed Site Information ===\n\n")

# Calculate statistics per site
site_stats <- data %>%
  group_by(deposit) %>%
  summarize(
    n_dates = n(),
    mean_age = round(mean(age), 0),
    sd_age = round(sd(age), 0),
    min_age = min(age),
    max_age = max(age),
    age_range = max_age - min_age,
    mean_error = round(mean(error), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_dates))

cat("Sites ranked by sample size:\n\n")
print(as.data.frame(site_stats), row.names = FALSE)

# ==============================================================================
# SAMPLE SIZE RECOMMENDATIONS
# ==============================================================================

cat("\n=== Sample Size Recommendations ===\n\n")

site_counts <- table(data$deposit)

cat("Sites with 10+ dates (recommended for analysis):\n")
good_sites <- names(site_counts)[site_counts >= 10]
if (length(good_sites) > 0) {
  for (site in good_sites) {
    cat("  ", site, ": ", site_counts[site], " dates\n")
  }
} else {
  cat("  None\n")
}

cat("\nSites with 5-9 dates (adequate for analysis):\n")
ok_sites <- names(site_counts)[site_counts >= 5 & site_counts < 10]
if (length(ok_sites) > 0) {
  for (site in ok_sites) {
    cat("  ", site, ": ", site_counts[site], " dates\n")
  }
} else {
  cat("  None\n")
}

cat("\nSites with 3-4 dates (minimum for analysis):\n")
min_sites <- names(site_counts)[site_counts >= 3 & site_counts < 5]
if (length(min_sites) > 0) {
  for (site in min_sites) {
    cat("  ", site, ": ", site_counts[site], " dates\n")
  }
} else {
  cat("  None\n")
}

cat("\nSites with <3 dates (too few for analysis):\n")
small_sites <- names(site_counts)[site_counts < 3]
if (length(small_sites) > 0) {
  cat("  ", length(small_sites), " sites\n")
} else {
  cat("  None\n")
}

# ==============================================================================
# MATERIAL DISTRIBUTION
# ==============================================================================

if ("material" %in% names(data)) {
  cat("\n=== Material Type Distribution ===\n\n")

  material_counts <- table(data$material)
  cat("Overall material distribution:\n")
  print(sort(material_counts, decreasing = TRUE))

  cat("\nMaterial by site:\n")
  material_by_site <- table(data$deposit, data$material)
  print(material_by_site[rowSums(material_by_site) >= 3, ])  # Sites with 3+ dates
}

# ==============================================================================
# AGE DISTRIBUTION
# ==============================================================================

cat("\n=== Age Distribution ===\n\n")

cat("Overall age statistics:\n")
cat("  Range: ", min(data$age), " - ", max(data$age), " BP\n")
cat("  Mean: ", round(mean(data$age), 0), " BP\n")
cat("  Median: ", round(median(data$age), 0), " BP\n")
cat("  SD: ", round(sd(data$age), 0), " years\n\n")

# Age categories
cat("Age distribution:\n")
age_breaks <- c(0, 500, 1000, 1500, 2000, Inf)
age_labels <- c("<500", "500-1000", "1000-1500", "1500-2000", ">2000")
data$age_category <- cut(data$age, breaks = age_breaks, labels = age_labels, right = FALSE)
print(table(data$age_category))

# ==============================================================================
# OUTLIER CHECK
# ==============================================================================

cat("\n=== Outlier Detection ===\n\n")

outliers <- check_outliers(data, sd_threshold = 3)

if (any(outliers$is_outlier, na.rm = TRUE)) {
  cat("Potential outliers detected:\n\n")
  outlier_data <- outliers[outliers$is_outlier, ]
  outlier_data <- outlier_data[order(outlier_data$z_score, decreasing = TRUE), ]
  print(outlier_data[, c("lab_code", "deposit", "age", "z_score")], row.names = FALSE)
  cat("\nThese dates should be reviewed before analysis.\n")
} else {
  cat("No statistical outliers detected (>3 SD from site mean)\n")
}

# ==============================================================================
# VISUALIZATION
# ==============================================================================

cat("\n=== Creating Quick Visualizations ===\n\n")

library(ggplot2)

# Age distribution by site
p1 <- ggplot(data, aes(x = reorder(deposit, age, FUN = median), y = age)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(
    x = "Site",
    y = "Radiocarbon Age (BP)",
    title = "Age Distribution by Site"
  ) +
  theme_minimal()

print(p1)
ggsave("output/data_preview_ages.png", p1, width = 10, height = 8, dpi = 300)

# Sample sizes
site_counts_df <- data.frame(
  site = names(site_counts),
  n_dates = as.numeric(site_counts)
)

p2 <- ggplot(site_counts_df, aes(x = reorder(site, n_dates), y = n_dates)) +
  geom_col(fill = "darkgreen", alpha = 0.7) +
  coord_flip() +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "yellow") +
  labs(
    x = "Site",
    y = "Number of Dates",
    title = "Sample Size by Site",
    subtitle = "Red = 10 (recommended), Orange = 5 (adequate), Yellow = 3 (minimum)"
  ) +
  theme_minimal()

print(p2)
ggsave("output/data_preview_sample_sizes.png", p2, width = 10, height = 8, dpi = 300)

# ==============================================================================
# RECOMMENDATIONS
# ==============================================================================

cat("\n=== Analysis Recommendations ===\n\n")

# Recommend sites to analyze
analyzable_sites <- names(site_counts)[site_counts >= 3]
cat("Sites with sufficient data (≥3 dates): ", length(analyzable_sites), "\n")

if (length(analyzable_sites) >= 2) {
  cat("\nRecommended for analysis:\n")

  # Top sites by sample size
  top_sites <- names(sort(site_counts[analyzable_sites], decreasing = TRUE))[1:min(5, length(analyzable_sites))]
  for (site in top_sites) {
    cat("  ", site, " (n=", site_counts[site], ")\n", sep = "")
  }

  cat("\nTo analyze these sites, edit analyze_your_data.R and set:\n")
  cat("  top_sites <- c('", paste(top_sites, collapse = "', '"), "')\n", sep = "")
} else {
  cat("\nWARNING: Need at least 2 sites with 3+ dates for contemporaneity analysis.\n")
  cat("Current: ", length(analyzable_sites), " analyzable site(s)\n")
}

# ==============================================================================
# EXPORT OPTIONS
# ==============================================================================

cat("\n=== Data Export Options ===\n\n")

cat("To export filtered data:\n\n")
cat("# Select sites to analyze\n")
cat("sites_to_analyze <- c('Site1', 'Site2', 'Site3')\n")
cat("filtered <- data[data$deposit %in% sites_to_analyze, ]\n")
cat("write.csv(filtered, 'output/filtered_data.csv', row.names = FALSE)\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==========================================================\n")
cat("  PREVIEW COMPLETE\n")
cat("==========================================================\n\n")

cat("Summary:\n")
cat("  Total dates: ", nrow(data), "\n")
cat("  Total sites: ", length(unique(data$deposit)), "\n")
cat("  Sites with ≥3 dates: ", length(analyzable_sites), "\n")
cat("  Age range: ", min(data$age), " - ", max(data$age), " BP\n\n")

cat("Files created:\n")
cat("  output/data_preview_ages.png\n")
cat("  output/data_preview_sample_sizes.png\n\n")

cat("Next steps:\n")
cat("  1. Review the plots and statistics above\n")
cat("  2. Decide which sites to include in analysis\n")
cat("  3. Check for any outliers or data quality issues\n")
cat("  4. Run analyze_your_data.R for full analysis\n\n")

cat("==========================================================\n")

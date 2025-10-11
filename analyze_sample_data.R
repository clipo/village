#' Analyze Sample Radiocarbon Data
#'
#' Analysis of the 157 radiocarbon dates, grouping by base site name
#' (extracting site from feature-specific identifiers)
#'

cat("==========================================================\n")
cat("  Sample Data Analysis\n")
cat("  157 Radiocarbon Dates\n")
cat("==========================================================\n\n")

# Load functions
library(dplyr)
library(stringr)
source("R/data_loading.R")
source("R/data_validation.R")

# Load data
cat("Loading data...\n")
data <- load_radiocarbon_data("data/radiocarbon_dates.csv", deposit_col = "site")

# Extract base site name (remove feature numbers and suffixes)
cat("\nExtracting base site names...\n")
data$base_site <- data$deposit

# Remove common feature patterns
data$base_site <- gsub("_F[0-9]+.*$", "", data$base_site)  # Remove _F123...
data$base_site <- gsub("_[0-9]+$", "", data$base_site)      # Remove trailing _1, _2
data$base_site <- gsub("_PAF.*$", "", data$base_site)       # Remove _PAF variants
data$base_site <- gsub("_House[0-9]+.*$", "", data$base_site)  # Remove _House1
data$base_site <- gsub("_Pit[0-9]+.*$", "", data$base_site)    # Remove _Pit4
data$base_site <- gsub("_post_mold.*$", "", data$base_site)    # Remove _post_mold
data$base_site <- gsub("_Trench_[0-9]+.*$", "", data$base_site) # Remove _Trench_44

# Manual fixes for complex names
data$base_site <- gsub("Blennerhassett_F23_[NS]", "Blennerhassett", data$base_site)
data$base_site <- gsub("Chenango_Point_South", "Chenango_Point", data$base_site)
data$base_site <- gsub("Chenango_Shores_East", "Chenango_Shores", data$base_site)
data$base_site <- gsub("JW_Wadsworth2", "JW_Wadsworth", data$base_site)
data$base_site <- gsub("Deposit_Airport1", "Deposit_Airport", data$base_site)
data$base_site <- gsub("Broome_Tech_PAF", "Broome_Tech", data$base_site)

# Show grouping results
cat("\nSite groupings:\n")
site_counts <- sort(table(data$base_site), decreasing = TRUE)
print(site_counts)

cat("\n=== Sites with Sufficient Data (>=3 dates) ===\n")
sufficient_sites <- names(site_counts)[site_counts >= 3]
cat("Number of sites: ", length(sufficient_sites), "\n\n")

for (site in sufficient_sites) {
  cat(sprintf("  %-30s: %2d dates\n", site, site_counts[site]))
}

# Use base_site for analysis
data$deposit <- factor(data$base_site)

# Filter to sites with at least 3 dates
filtered_data <- data[data$base_site %in% sufficient_sites, ]
filtered_data$deposit <- factor(filtered_data$deposit)

cat("\n=== Final Dataset for Analysis ===\n")
cat("Total dates:", nrow(filtered_data), "\n")
cat("Number of sites:", length(unique(filtered_data$deposit)), "\n")
cat("Sites:", paste(levels(filtered_data$deposit), collapse = ", "), "\n")

# Validate
cat("\nValidating data...\n")
validate_c14_data(filtered_data)

# Summary by site
cat("\n=== Summary by Site ===\n\n")
site_summary <- filtered_data %>%
  group_by(deposit) %>%
  summarize(
    n = n(),
    mean_age = round(mean(age)),
    sd_age = round(sd(age)),
    min_age = min(age),
    max_age = max(age),
    range = max_age - min_age,
    mean_error = round(mean(error), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

print(as.data.frame(site_summary), row.names = FALSE)

# Check for outliers
cat("\n=== Checking for Outliers ===\n")
outliers <- check_outliers(filtered_data, sd_threshold = 3)

if (any(outliers$is_outlier, na.rm = TRUE)) {
  cat("\nPotential outliers:\n")
  print(outliers[outliers$is_outlier, c("lab_code", "deposit", "age", "z_score")],
        row.names = FALSE)
}

# Material distribution
cat("\n=== Material Distribution ===\n")
print(table(filtered_data$material))

cat("\n=== Data Ready for Analysis ===\n")
cat("Proceed? The full analysis will take 10-15 minutes.\n")
cat("Run: source('run_full_analysis.R')\n")

# Save filtered data for analysis
saveRDS(filtered_data, "output/filtered_data_for_analysis.rds")
cat("\nFiltered data saved to: output/filtered_data_for_analysis.rds\n")

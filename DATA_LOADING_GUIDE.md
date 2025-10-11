# Data Loading Guide

Complete guide for loading and preparing your radiocarbon data.

## Quick Start

### Option 1: Preview Your Data First (Recommended)

```r
source("preview_data.R")
```

This will:
- Load your data with automatic column mapping
- Show summary statistics for each site
- Identify which sites have sufficient sample sizes
- Create visualizations
- Provide recommendations for analysis

### Option 2: Run Full Analysis

```r
source("analyze_your_data.R")
```

This will:
- Load and validate your data
- Select top sites by sample size
- Calibrate all dates
- Fit and compare models
- Calculate contemporaneity metrics
- Create all visualizations
- Export results

## Data File Formats

### Supported Formats

- **CSV** (`.csv`) - Comma-separated values
- **Excel** (`.xlsx`, `.xls`) - Microsoft Excel files

### Required Columns

Your data file must contain these columns (names can vary):

| Standard Name | Alternative Names | Description |
|---------------|------------------|-------------|
| `lab_code` | lab_no, labcode, sample_id | Laboratory identifier |
| `age` | c14_age, 14c_age, bp, date | Radiocarbon age (years BP) |
| `error` | c14_error, 14c_error, sd, uncertainty | Measurement error (1σ) |
| `deposit` | site, context, location, stratum | Deposit/site identifier |

### Optional Columns

| Column | Alternative Names | Description |
|--------|------------------|-------------|
| `material` | sample_type, type | Material dated (e.g., charcoal, bone) |
| `context` | archaeological_context, feature | Archaeological context |
| `reference` | ref, citation, source | Publication reference |

## Column Name Mapping

The system automatically detects column names, but you can also specify custom mappings.

### Automatic Detection

```r
# Automatically detects columns
data <- load_radiocarbon_data("data/your_file.csv")
```

The function will report what it found:
```
Column mapping detected/applied:
  lab_code     <- lab_no
  age          <- c14_age
  error        <- c14_error
  deposit      <- site
  material     <- material
  reference    <- reference
```

### Manual Column Mapping

If automatic detection fails or you want to specify:

```r
data <- load_radiocarbon_data(
  "data/your_file.csv",
  column_map = list(
    lab_code = "lab_no",
    age = "c14_age",
    error = "c14_error",
    deposit = "site"
  )
)
```

### Using Different Deposit Column

If you want to use a specific column for grouping (e.g., "site" instead of "deposit"):

```r
data <- load_radiocarbon_data(
  "data/your_file.csv",
  deposit_col = "site"  # Use "site" column as deposit identifier
)
```

## Loading Excel Files

```r
# Load first sheet (default)
data <- load_radiocarbon_data("data/dates.xlsx")

# Load specific sheet by name
data <- load_radiocarbon_data("data/dates.xlsx", sheet = "Dates")

# Load specific sheet by number
data <- load_radiocarbon_data("data/dates.xlsx", sheet = 2)
```

**Note**: Requires `readxl` package:
```r
install.packages("readxl")
```

## Complete Workflow

### 1. Load and Validate

```r
source("R/data_loading.R")
source("R/data_validation.R")

# Load data
data <- load_radiocarbon_data("data/radiocarbon_dates.csv")

# Validate
validate_c14_data(data)

# Print summary
print_data_summary(data)
```

### 2. Filter by Sample Size

```r
# Remove sites with < 3 dates
data_filtered <- load_and_prepare_data(
  "data/radiocarbon_dates.csv",
  min_dates_per_deposit = 3
)
```

### 3. Select Specific Sites

```r
# Load all data
data <- load_radiocarbon_data("data/radiocarbon_dates.csv")

# Select specific sites
sites_to_analyze <- c("Site_A", "Site_B", "Site_C")
data_filtered <- data[data$deposit %in% sites_to_analyze, ]

# Validate filtered data
validate_c14_data(data_filtered)
```

### 4. Export Standardized Data

```r
# Export to CSV
export_standardized_data(data, "output/standardized_dates.csv")
```

## Your Current Data

### File Structure

Your file `data/radiocarbon_dates.csv` has:
- **157 radiocarbon dates**
- **Columns**: site, material, lab_no, c14_age, c14_error, reference
- **Mapping**:
  - `site` → `deposit` (for analysis)
  - `lab_no` → `lab_code`
  - `c14_age` → `age`
  - `c14_error` → `error`
  - `material` → `material` (preserved)
  - `reference` → `reference` (preserved)

### Quick Analysis Commands

```r
# 1. Preview data
source("preview_data.R")

# 2. Analyze top 5 sites
source("analyze_your_data.R")

# 3. Or select specific sites
source("R/data_loading.R")
data <- load_radiocarbon_data("data/radiocarbon_dates.csv", deposit_col = "site")

# Choose sites
selected_sites <- c("Site1", "Site2", "Site3")
filtered_data <- data[data$deposit %in% selected_sites, ]

# Continue with analysis...
```

## Data Quality Checks

### Check for Issues

```r
source("R/data_validation.R")

# Load data
data <- load_radiocarbon_data("data/radiocarbon_dates.csv")

# Check for outliers
outliers <- check_outliers(data, sd_threshold = 3)
print(outliers[outliers$is_outlier, ])

# Summary by site
summarize_c14_data(data)
```

### Common Issues

1. **Duplicate lab codes**: Each date needs a unique identifier
2. **Missing values**: All required fields must have values
3. **Unrealistic ages**: <50 BP or >50000 BP
4. **Very small errors**: <10 years may cause issues
5. **Very large errors**: >500 years will have wide ranges

## Sample Size Guidelines

Based on validation studies:

| Sample Size | Suitability | Notes |
|-------------|-------------|-------|
| 1-2 dates | Insufficient | Cannot analyze contemporaneity |
| 3-4 dates | Minimum | Limited power, wide intervals |
| 5-9 dates | Adequate | Good for most analyses |
| 10+ dates | Recommended | Best precision and power |

**For contemporaneity analysis**: Need at least **2 deposits** with **≥3 dates each**

## Examples

### Example 1: Load CSV with Auto-Detection

```r
data <- load_radiocarbon_data("data/dates.csv")
print_data_summary(data)
```

### Example 2: Load Excel with Custom Mapping

```r
data <- load_radiocarbon_data(
  "data/dates.xlsx",
  sheet = "Sheet1",
  column_map = list(
    lab_code = "LabID",
    age = "Age_BP",
    error = "Error",
    deposit = "Context"
  )
)
```

### Example 3: Filter and Analyze

```r
# Load with filtering
data <- load_and_prepare_data(
  "data/dates.csv",
  min_dates_per_deposit = 5,  # At least 5 dates per site
  deposit_col = "site"
)

# Proceed with analysis
source("R/calibration.R")
cal_result <- calibrate_dates(data$age, data$error, data$lab_code)
```

### Example 4: Group by Custom Column

```r
# If your data has "region" and you want to compare regions
data <- load_radiocarbon_data("data/dates.csv")

# Use region as deposit
data$deposit <- data$region

# Validate and continue
validate_c14_data(data)
```

## Troubleshooting

### Column Not Found Error

```
Error: Mapped column 'lab_no' (for lab_code) not found in file.
Available columns: site, sample_id, age, error
```

**Solution**: Specify correct mapping:
```r
data <- load_radiocarbon_data(
  "data/dates.csv",
  column_map = list(lab_code = "sample_id")
)
```

### Missing Required Column

```
Error: Required columns not found or mapped: age
Available columns: site, c14_bp, error
```

**Solution**: Map the column:
```r
data <- load_radiocarbon_data(
  "data/dates.csv",
  column_map = list(age = "c14_bp")
)
```

### Excel Package Not Installed

```
Error: Package 'readxl' is required for Excel files
```

**Solution**:
```r
install.packages("readxl")
```

### Too Few Dates per Deposit

```
Warning: 10 deposits have fewer than 3 dates
```

**Solution**: Filter or select different deposits:
```r
# Option 1: Increase sample size filter
data <- load_and_prepare_data("data/dates.csv", min_dates_per_deposit = 3)

# Option 2: Manually select sites with enough dates
site_counts <- table(data$deposit)
good_sites <- names(site_counts)[site_counts >= 3]
filtered <- data[data$deposit %in% good_sites, ]
```

## Next Steps

After loading your data:

1. **Preview**: Run `preview_data.R` to explore
2. **Validate**: Check for outliers and data quality issues
3. **Select**: Choose which deposits to analyze
4. **Analyze**: Run `analyze_your_data.R` or custom workflow
5. **Interpret**: Review results in archaeological context

## Functions Reference

| Function | Purpose |
|----------|---------|
| `load_radiocarbon_data()` | Load and map columns |
| `load_and_prepare_data()` | Load, validate, and filter |
| `print_data_summary()` | Print formatted summary |
| `export_standardized_data()` | Save to CSV |
| `auto_detect_columns()` | Detect column names |
| `validate_c14_data()` | Validate data structure |
| `check_outliers()` | Identify statistical outliers |

## Additional Resources

- **README.md** - Complete project documentation
- **QUICKSTART.md** - Getting started guide
- **example_workflow.R** - Example with simulated data
- **analyze_your_data.R** - Template for your data

---

**Need Help?** Check README.md troubleshooting section or run `verify_installation.R`

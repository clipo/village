#' Data Loading Functions
#'
#' Functions to load radiocarbon data from CSV and Excel files
#' with flexible column name mapping.
#'
#' @author Project Team
#' @date 2025-10-11

library(dplyr)

#' Load radiocarbon data from file
#'
#' Loads radiocarbon dates from CSV or Excel files with flexible
#' column name mapping. Automatically detects file type and handles
#' various column naming conventions.
#'
#' @param file_path Path to CSV or Excel file
#' @param sheet Sheet name or number (for Excel files, default 1)
#' @param column_map Named list mapping standard names to file column names.
#'   Standard names are: lab_code, age, error, deposit, material, context, reference
#' @param auto_map Logical, attempt automatic column mapping (default TRUE)
#' @param deposit_col Column to use for deposit grouping. Can be "site" or other column name.
#'   If NULL, will look for "deposit", "site", or "context" (default NULL)
#' @return Data frame with standardized column names
#' @examples
#' # Load with automatic mapping
#' data <- load_radiocarbon_data("data/radiocarbon_dates.csv")
#'
#' # Load with custom column mapping
#' data <- load_radiocarbon_data(
#'   "data/my_dates.csv",
#'   column_map = list(
#'     lab_code = "lab_no",
#'     age = "c14_age",
#'     error = "c14_error",
#'     deposit = "site"
#'   )
#' )
#' @export
load_radiocarbon_data <- function(file_path,
                                   sheet = 1,
                                   column_map = NULL,
                                   auto_map = TRUE,
                                   deposit_col = NULL) {

  # Check file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Determine file type and load
  file_ext <- tolower(tools::file_ext(file_path))

  if (file_ext == "csv") {
    raw_data <- read.csv(file_path, stringsAsFactors = FALSE)
  } else if (file_ext %in% c("xlsx", "xls")) {
    # Load readxl if needed
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required for Excel files. Install with: install.packages('readxl')")
    }
    raw_data <- readxl::read_excel(file_path, sheet = sheet)
    raw_data <- as.data.frame(raw_data)
  } else {
    stop("Unsupported file type: ", file_ext, ". Use .csv, .xlsx, or .xls")
  }

  # If column_map not provided, create from auto-detection
  if (is.null(column_map)) {
    column_map <- auto_detect_columns(names(raw_data), deposit_col)
  } else {
    # Fill in auto-detected values for missing mappings
    auto_cols <- auto_detect_columns(names(raw_data), deposit_col)
    for (col_name in names(auto_cols)) {
      if (!col_name %in% names(column_map)) {
        column_map[[col_name]] <- auto_cols[[col_name]]
      }
    }
  }

  # Report column mapping
  message("Column mapping detected/applied:")
  for (std_name in names(column_map)) {
    message(sprintf("  %-12s <- %s", std_name, column_map[[std_name]]))
  }

  # Check required columns are mapped
  required <- c("lab_code", "age", "error", "deposit")
  missing <- setdiff(required, names(column_map))
  if (length(missing) > 0) {
    stop("Required columns not found or mapped: ", paste(missing, collapse = ", "),
         "\nAvailable columns: ", paste(names(raw_data), collapse = ", "))
  }

  # Check mapped columns exist in data
  for (std_name in names(column_map)) {
    file_col <- column_map[[std_name]]
    if (!file_col %in% names(raw_data)) {
      stop("Mapped column '", file_col, "' (for ", std_name, ") not found in file.",
           "\nAvailable columns: ", paste(names(raw_data), collapse = ", "))
    }
  }

  # Create standardized data frame
  std_data <- data.frame(
    lab_code = as.character(raw_data[[column_map$lab_code]]),
    age = as.numeric(raw_data[[column_map$age]]),
    error = as.numeric(raw_data[[column_map$error]]),
    deposit = as.character(raw_data[[column_map$deposit]]),
    stringsAsFactors = FALSE
  )

  # Add optional columns if mapped
  if ("material" %in% names(column_map)) {
    std_data$material <- as.character(raw_data[[column_map$material]])
  }

  if ("context" %in% names(column_map)) {
    std_data$context <- as.character(raw_data[[column_map$context]])
  }

  if ("reference" %in% names(column_map)) {
    std_data$reference <- as.character(raw_data[[column_map$reference]])
  }

  # Store original data and mapping as attributes
  attr(std_data, "original_file") <- file_path
  attr(std_data, "column_map") <- column_map
  attr(std_data, "load_date") <- Sys.time()

  message("\nLoaded ", nrow(std_data), " radiocarbon dates from ",
          length(unique(std_data$deposit)), " deposits")

  return(std_data)
}

#' Auto-detect column names
#'
#' Attempts to automatically map file column names to standard names
#' based on common naming conventions.
#'
#' @param file_columns Character vector of column names in file
#' @param deposit_col Preferred column for deposit (optional)
#' @return Named list with column mappings
#' @keywords internal
auto_detect_columns <- function(file_columns, deposit_col = NULL) {

  col_lower <- tolower(file_columns)
  mapping <- list()

  # Lab code patterns
  lab_patterns <- c("lab_code", "labcode", "lab_no", "labno", "lab_id",
                    "sample_id", "sample_no", "id")
  mapping$lab_code <- find_column_match(col_lower, lab_patterns, file_columns)

  # Age patterns
  age_patterns <- c("age", "c14_age", "c14age", "14c_age", "bp", "rcybp", "date")
  mapping$age <- find_column_match(col_lower, age_patterns, file_columns)

  # Error patterns
  error_patterns <- c("error", "c14_error", "c14error", "14c_error", "sd", "std", "uncertainty")
  mapping$error <- find_column_match(col_lower, error_patterns, file_columns)

  # Deposit patterns (flexible based on deposit_col)
  if (!is.null(deposit_col)) {
    # User specified which column to use
    if (deposit_col %in% file_columns) {
      mapping$deposit <- deposit_col
    } else {
      warning("Specified deposit_col '", deposit_col, "' not found in data. Using auto-detection.")
      deposit_patterns <- c("deposit", "site", "context", "location", "stratum", "layer")
      mapping$deposit <- find_column_match(col_lower, deposit_patterns, file_columns)
    }
  } else {
    # Auto-detect (prefer deposit, then site, then context)
    deposit_patterns <- c("deposit", "site", "context", "location", "stratum", "layer")
    mapping$deposit <- find_column_match(col_lower, deposit_patterns, file_columns)
  }

  # Optional columns
  material_patterns <- c("material", "sample_type", "type")
  material_match <- find_column_match(col_lower, material_patterns, file_columns, required = FALSE)
  if (!is.null(material_match)) {
    mapping$material <- material_match
  }

  context_patterns <- c("context", "archaeological_context", "feature", "stratum")
  context_match <- find_column_match(col_lower, context_patterns, file_columns, required = FALSE)
  if (!is.null(context_match)) {
    mapping$context <- context_match
  }

  reference_patterns <- c("reference", "ref", "citation", "source", "publication")
  reference_match <- find_column_match(col_lower, reference_patterns, file_columns, required = FALSE)
  if (!is.null(reference_match)) {
    mapping$reference <- reference_match
  }

  return(mapping)
}

#' Find matching column name
#'
#' Searches for column name matching any of the patterns.
#'
#' @param col_lower Lowercase column names
#' @param patterns Patterns to search for
#' @param file_columns Original column names
#' @param required Logical, whether column is required (default TRUE)
#' @return Column name or NULL
#' @keywords internal
find_column_match <- function(col_lower, patterns, file_columns, required = TRUE) {

  for (pattern in patterns) {
    # Exact match
    idx <- which(col_lower == pattern)
    if (length(idx) > 0) {
      return(file_columns[idx[1]])
    }

    # Partial match
    idx <- which(grepl(pattern, col_lower))
    if (length(idx) > 0) {
      return(file_columns[idx[1]])
    }
  }

  if (required) {
    stop("Could not find column matching patterns: ", paste(patterns, collapse = ", "),
         "\nAvailable columns: ", paste(file_columns, collapse = ", "))
  }

  return(NULL)
}

#' Load and prepare data for analysis
#'
#' Complete workflow: load file, validate, and optionally group deposits.
#'
#' @param file_path Path to data file
#' @param group_by Optional column to group/aggregate deposits by (e.g., "site")
#' @param min_dates_per_deposit Minimum dates required per deposit (default 3)
#' @param sheet Excel sheet (default 1)
#' @param column_map Custom column mapping (optional)
#' @return Validated and optionally filtered data frame
#' @export
load_and_prepare_data <- function(file_path,
                                   group_by = NULL,
                                   min_dates_per_deposit = 3,
                                   sheet = 1,
                                   column_map = NULL) {

  # Load data
  data <- load_radiocarbon_data(file_path, sheet = sheet, column_map = column_map)

  # Optional: Group deposits
  if (!is.null(group_by)) {
    if (!group_by %in% names(data)) {
      stop("Grouping column '", group_by, "' not found in data")
    }
    message("\nGrouping deposits by: ", group_by)
    data$deposit <- data[[group_by]]
  }

  # Validate
  message("\nValidating data...")
  source("R/data_validation.R")
  validate_c14_data(data)

  # Filter by minimum sample size
  deposit_counts <- table(data$deposit)
  small_deposits <- names(deposit_counts)[deposit_counts < min_dates_per_deposit]

  if (length(small_deposits) > 0) {
    message("\nWarning: ", length(small_deposits), " deposits have fewer than ",
            min_dates_per_deposit, " dates:")
    for (dep in small_deposits) {
      message("  ", dep, ": ", deposit_counts[dep], " dates")
    }

    message("\nRemoving deposits with < ", min_dates_per_deposit, " dates...")
    data <- data[!data$deposit %in% small_deposits, ]
    data$deposit <- factor(data$deposit)  # Drop unused levels

    message("Retained ", nrow(data), " dates from ",
            length(unique(data$deposit)), " deposits")
  }

  return(data)
}

#' Print data loading summary
#'
#' @param data Loaded data frame with attributes
#' @export
print_data_summary <- function(data) {

  cat("=== Radiocarbon Data Summary ===\n\n")

  # File info
  if (!is.null(attr(data, "original_file"))) {
    cat("Source file:", attr(data, "original_file"), "\n")
    cat("Loaded:", format(attr(data, "load_date"), "%Y-%m-%d %H:%M:%S"), "\n\n")
  }

  # Column mapping
  if (!is.null(attr(data, "column_map"))) {
    cat("Column mapping:\n")
    for (std_name in names(attr(data, "column_map"))) {
      cat(sprintf("  %-12s <- %s\n", std_name, attr(data, "column_map")[[std_name]]))
    }
    cat("\n")
  }

  # Data summary
  cat("Total dates:", nrow(data), "\n")
  cat("Number of deposits:", length(unique(data$deposit)), "\n\n")

  cat("Dates per deposit:\n")
  print(sort(table(data$deposit), decreasing = TRUE))
  cat("\n")

  cat("Age range:", min(data$age), "-", max(data$age), "BP\n")
  cat("Mean error:", round(mean(data$error), 1), "±",
      round(sd(data$error), 1), "years\n\n")

  if ("material" %in% names(data)) {
    cat("Material types:\n")
    print(sort(table(data$material), decreasing = TRUE))
    cat("\n")
  }

  if ("reference" %in% names(data)) {
    cat("References:\n")
    print(sort(table(data$reference), decreasing = TRUE))
    cat("\n")
  }
}

#' Export standardized data
#'
#' Save standardized data to CSV for future use.
#'
#' @param data Data frame with radiocarbon dates
#' @param output_path Path for output CSV file
#' @export
export_standardized_data <- function(data, output_path) {

  # Remove attributes before saving
  output_data <- data
  attributes(output_data) <- attributes(output_data)[c("names", "row.names", "class")]

  write.csv(output_data, output_path, row.names = FALSE)

  message("Exported ", nrow(output_data), " dates to: ", output_path)
}

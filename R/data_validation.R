#' Data Validation Functions for Radiocarbon Dating Analysis
#'
#' Functions to validate radiocarbon date data structures and catch common errors
#' before analysis. All functions return informative error messages.
#'
#' @author Project Team
#' @date 2025-10-11

#' Validate radiocarbon date data frame
#'
#' Checks that a data frame of radiocarbon dates meets all requirements
#' for downstream analysis, including proper structure, valid ranges,
#' and no duplicate laboratory codes.
#'
#' @param data Data frame with radiocarbon dates
#' @param require_columns Character vector of required column names
#' @return Invisible TRUE if valid, otherwise throws informative error
#' @examples
#' data <- data.frame(
#'   lab_code = c("OxA-1234", "OxA-1235"),
#'   age = c(2450, 2380),
#'   error = c(30, 25),
#'   deposit = c("A", "A")
#' )
#' validate_c14_data(data)
#' @export
validate_c14_data <- function(data,
                               require_columns = c("lab_code", "age", "error", "deposit")) {

  # Check that data is a data frame
  if (!is.data.frame(data)) {
    stop("Input must be a data frame. Received: ", class(data)[1])
  }

  # Check for required columns
  missing_cols <- setdiff(require_columns, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nAvailable columns: ", paste(names(data), collapse = ", "))
  }

  # Check for empty data
  if (nrow(data) == 0) {
    stop("Data frame contains no rows")
  }

  # Validate lab codes
  if (any(is.na(data$lab_code)) || any(data$lab_code == "")) {
    stop("Laboratory codes cannot be missing or empty")
  }

  # Check for duplicate lab codes
  duplicated_codes <- data$lab_code[duplicated(data$lab_code)]
  if (length(duplicated_codes) > 0) {
    stop("Duplicate laboratory codes found: ",
         paste(unique(duplicated_codes), collapse = ", "))
  }

  # Validate ages
  if (any(is.na(data$age))) {
    stop("Radiocarbon ages contain NA values at rows: ",
         paste(which(is.na(data$age)), collapse = ", "))
  }

  if (any(data$age < 0)) {
    stop("Radiocarbon ages must be non-negative. Negative values found at rows: ",
         paste(which(data$age < 0), collapse = ", "))
  }

  # Check reasonable date range (modern to limit of radiocarbon dating)
  if (any(data$age < 50)) {
    warning("Very young dates (<50 BP) found at rows: ",
            paste(which(data$age < 50), collapse = ", "),
            "\nThese may be modern samples or require special calibration.")
  }

  if (any(data$age > 50000)) {
    stop("Radiocarbon ages exceed 50,000 BP at rows: ",
         paste(which(data$age > 50000), collapse = ", "),
         "\nThese exceed the practical limit of radiocarbon dating.")
  }

  # Validate errors
  if (any(is.na(data$error))) {
    stop("Radiocarbon errors contain NA values at rows: ",
         paste(which(is.na(data$error)), collapse = ", "))
  }

  if (any(data$error <= 0)) {
    stop("Radiocarbon errors must be positive. Non-positive values found at rows: ",
         paste(which(data$error <= 0), collapse = ", "))
  }

  # Check for very small errors that may cause numerical issues
  if (any(data$error < 10)) {
    warning("Very small errors (<10 years) found at rows: ",
            paste(which(data$error < 10), collapse = ", "),
            "\nThese may cause numerical issues during calibration.")
  }

  # Check for unreasonably large errors
  if (any(data$error > 500)) {
    warning("Very large errors (>500 years) found at rows: ",
            paste(which(data$error > 500), collapse = ", "),
            "\nThese dates will have very wide calibrated ranges.")
  }

  # Validate deposit identifiers
  if (any(is.na(data$deposit))) {
    stop("Deposit identifiers contain NA values at rows: ",
         paste(which(is.na(data$deposit)), collapse = ", "))
  }

  # Check that each deposit has at least one date
  deposit_counts <- table(data$deposit)
  if (any(deposit_counts == 0)) {
    warning("Some deposits have no dates assigned")
  }

  # Check optional material type if present
  if ("material" %in% names(data)) {
    if (any(is.na(data$material) | data$material == "")) {
      warning("Material type contains missing values at rows: ",
              paste(which(is.na(data$material) | data$material == ""), collapse = ", "))
    }
  }

  message("Data validation passed: ",
          nrow(data), " dates from ",
          length(unique(data$deposit)), " deposits")

  return(invisible(TRUE))
}

#' Create standardized radiocarbon data structure
#'
#' Creates a properly formatted data frame for radiocarbon dates with
#' all required fields and validation.
#'
#' @param lab_code Character vector of laboratory codes
#' @param age Numeric vector of radiocarbon ages in years BP
#' @param error Numeric vector of standard errors in years
#' @param deposit Vector of deposit identifiers (character or numeric)
#' @param material Character vector of material types (optional)
#' @param context Character vector of archaeological context (optional)
#' @param stratigraphy Character vector of stratigraphic information (optional)
#' @return Data frame with standardized structure
#' @examples
#' c14_data <- create_c14_dataset(
#'   lab_code = c("OxA-1234", "OxA-1235", "Beta-5678"),
#'   age = c(2450, 2380, 2420),
#'   error = c(30, 25, 35),
#'   deposit = c("A", "A", "B"),
#'   material = c("charcoal", "charcoal", "bone")
#' )
#' @export
create_c14_dataset <- function(lab_code,
                                age,
                                error,
                                deposit,
                                material = NULL,
                                context = NULL,
                                stratigraphy = NULL) {

  # Check that vectors have same length
  n <- length(lab_code)
  if (length(age) != n || length(error) != n || length(deposit) != n) {
    stop("lab_code, age, error, and deposit must have the same length")
  }

  # Create base data frame
  data <- data.frame(
    lab_code = as.character(lab_code),
    age = as.numeric(age),
    error = as.numeric(error),
    deposit = deposit,
    stringsAsFactors = FALSE
  )

  # Add optional columns if provided
  if (!is.null(material)) {
    if (length(material) != n) {
      stop("material must have same length as lab_code")
    }
    data$material <- as.character(material)
  }

  if (!is.null(context)) {
    if (length(context) != n) {
      stop("context must have same length as lab_code")
    }
    data$context <- as.character(context)
  }

  if (!is.null(stratigraphy)) {
    if (length(stratigraphy) != n) {
      stop("stratigraphy must have same length as lab_code")
    }
    data$stratigraphy <- as.character(stratigraphy)
  }

  # Validate the created dataset
  validate_c14_data(data)

  return(data)
}

#' Check for potential outliers in radiocarbon dataset
#'
#' Identifies dates that are statistical outliers within their deposit
#' based on calibrated date ranges.
#'
#' @param data Data frame with radiocarbon dates
#' @param sd_threshold Number of standard deviations for outlier detection (default 3)
#' @return Data frame with outlier flags and statistics
#' @export
check_outliers <- function(data, sd_threshold = 3) {

  validate_c14_data(data)

  results <- data.frame(
    lab_code = data$lab_code,
    deposit = data$deposit,
    age = data$age,
    is_outlier = FALSE,
    z_score = NA_real_,
    stringsAsFactors = FALSE
  )

  # Check within each deposit
  for (dep in unique(data$deposit)) {
    dep_idx <- which(data$deposit == dep)

    if (length(dep_idx) > 2) {  # Need at least 3 dates to identify outliers
      dep_ages <- data$age[dep_idx]
      mean_age <- mean(dep_ages)
      sd_age <- sd(dep_ages)

      z_scores <- abs((dep_ages - mean_age) / sd_age)
      outliers <- z_scores > sd_threshold

      results$z_score[dep_idx] <- z_scores
      results$is_outlier[dep_idx] <- outliers
    }
  }

  n_outliers <- sum(results$is_outlier, na.rm = TRUE)
  if (n_outliers > 0) {
    warning("Detected ", n_outliers, " potential outliers:\n",
            paste(results$lab_code[results$is_outlier], collapse = ", "),
            "\nConsider reviewing these dates before analysis.")
  } else {
    message("No statistical outliers detected")
  }

  return(results)
}

#' Summarize radiocarbon dataset
#'
#' Provides summary statistics for a radiocarbon dating dataset
#'
#' @param data Data frame with radiocarbon dates
#' @return List with summary statistics
#' @export
summarize_c14_data <- function(data) {

  validate_c14_data(data)

  summary_list <- list(
    n_dates = nrow(data),
    n_deposits = length(unique(data$deposit)),
    deposits = unique(data$deposit),
    dates_per_deposit = table(data$deposit),
    age_range = range(data$age),
    mean_age = mean(data$age),
    median_age = median(data$age),
    mean_error = mean(data$error),
    median_error = median(data$error),
    error_range = range(data$error)
  )

  # Add material summary if available
  if ("material" %in% names(data)) {
    summary_list$materials <- table(data$material)
  }

  return(summary_list)
}

#' Print summary of radiocarbon dataset
#'
#' @param data Data frame with radiocarbon dates
#' @export
print_c14_summary <- function(data) {
  summ <- summarize_c14_data(data)

  cat("=== Radiocarbon Dataset Summary ===\n\n")
  cat("Total dates:", summ$n_dates, "\n")
  cat("Number of deposits:", summ$n_deposits, "\n")
  cat("Deposits:", paste(summ$deposits, collapse = ", "), "\n\n")

  cat("Dates per deposit:\n")
  print(summ$dates_per_deposit)
  cat("\n")

  cat("Age range:", summ$age_range[1], "-", summ$age_range[2], "BP\n")
  cat("Mean age:", round(summ$mean_age, 1), "BP\n")
  cat("Median age:", round(summ$median_age, 1), "BP\n\n")

  cat("Error statistics:\n")
  cat("  Mean error:", round(summ$mean_error, 1), "years\n")
  cat("  Median error:", round(summ$median_error, 1), "years\n")
  cat("  Error range:", summ$error_range[1], "-", summ$error_range[2], "years\n")

  if (!is.null(summ$materials)) {
    cat("\nMaterial types:\n")
    print(summ$materials)
  }

  cat("\n")
}

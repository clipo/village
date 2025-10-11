#' Unit Tests for Data Validation Functions
#'
#' @author Project Team
#' @date 2025-10-11

library(testthat)

source("../../R/data_validation.R")

test_that("validate_c14_data accepts valid data", {
  data <- data.frame(
    lab_code = c("OxA-1234", "OxA-1235"),
    age = c(2450, 2380),
    error = c(30, 25),
    deposit = c("A", "A")
  )

  expect_silent(validate_c14_data(data))
})

test_that("validate_c14_data rejects missing columns", {
  data <- data.frame(
    lab_code = c("OxA-1234", "OxA-1235"),
    age = c(2450, 2380)
  )

  expect_error(validate_c14_data(data), "Missing required columns")
})

test_that("validate_c14_data rejects negative ages", {
  data <- data.frame(
    lab_code = c("OxA-1234", "OxA-1235"),
    age = c(2450, -100),
    error = c(30, 25),
    deposit = c("A", "A")
  )

  expect_error(validate_c14_data(data), "must be non-negative")
})

test_that("validate_c14_data rejects ages > 50000 BP", {
  data <- data.frame(
    lab_code = c("OxA-1234", "OxA-1235"),
    age = c(2450, 55000),
    error = c(30, 25),
    deposit = c("A", "A")
  )

  expect_error(validate_c14_data(data), "exceed 50,000 BP")
})

test_that("validate_c14_data rejects non-positive errors", {
  data <- data.frame(
    lab_code = c("OxA-1234", "OxA-1235"),
    age = c(2450, 2380),
    error = c(30, 0),
    deposit = c("A", "A")
  )

  expect_error(validate_c14_data(data), "must be positive")
})

test_that("validate_c14_data warns about small errors", {
  data <- data.frame(
    lab_code = c("OxA-1234", "OxA-1235"),
    age = c(2450, 2380),
    error = c(30, 5),
    deposit = c("A", "A")
  )

  expect_warning(validate_c14_data(data), "Very small errors")
})

test_that("validate_c14_data rejects duplicate lab codes", {
  data <- data.frame(
    lab_code = c("OxA-1234", "OxA-1234"),
    age = c(2450, 2380),
    error = c(30, 25),
    deposit = c("A", "B")
  )

  expect_error(validate_c14_data(data), "Duplicate laboratory codes")
})

test_that("create_c14_dataset creates valid structure", {
  data <- create_c14_dataset(
    lab_code = c("OxA-1234", "OxA-1235"),
    age = c(2450, 2380),
    error = c(30, 25),
    deposit = c("A", "A")
  )

  expect_true(is.data.frame(data))
  expect_equal(nrow(data), 2)
  expect_equal(names(data), c("lab_code", "age", "error", "deposit"))
})

test_that("create_c14_dataset rejects mismatched lengths", {
  expect_error(
    create_c14_dataset(
      lab_code = c("OxA-1234", "OxA-1235"),
      age = c(2450),
      error = c(30, 25),
      deposit = c("A", "A")
    ),
    "same length"
  )
})

test_that("check_outliers identifies outliers correctly", {
  data <- data.frame(
    lab_code = c("OxA-1", "OxA-2", "OxA-3", "OxA-4"),
    age = c(2450, 2460, 2455, 3500),  # Last one is outlier
    error = c(30, 25, 28, 30),
    deposit = c("A", "A", "A", "A")
  )

  result <- check_outliers(data, sd_threshold = 3)

  expect_true(result$is_outlier[4])
  expect_false(result$is_outlier[1])
})

test_that("summarize_c14_data returns correct structure", {
  data <- data.frame(
    lab_code = c("OxA-1", "OxA-2", "OxA-3"),
    age = c(2450, 2460, 2455),
    error = c(30, 25, 28),
    deposit = c("A", "A", "B")
  )

  summ <- summarize_c14_data(data)

  expect_equal(summ$n_dates, 3)
  expect_equal(summ$n_deposits, 2)
  expect_true(is.numeric(summ$mean_age))
})

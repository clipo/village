#' Unit Tests for Calibration Functions
#'
#' @author Project Team
#' @date 2025-10-11

library(testthat)

source("../../R/calibration.R")
source("../../R/data_validation.R")

test_that("calibrate_dates returns correct structure", {
  ages <- c(2450, 2380)
  errors <- c(30, 25)

  result <- calibrate_dates(ages, errors)

  expect_true(inherits(result, "calibrated_c14"))
  expect_true(is.list(result$metadata))
  expect_true(is.data.frame(result$summary))
  expect_equal(nrow(result$summary), 2)
})

test_that("calibrate_dates rejects mismatched lengths", {
  ages <- c(2450, 2380)
  errors <- c(30)

  expect_error(calibrate_dates(ages, errors), "same length")
})

test_that("calibrate_dates rejects NA values", {
  ages <- c(2450, NA)
  errors <- c(30, 25)

  expect_error(calibrate_dates(ages, errors), "cannot contain NA")
})

test_that("calibrate_dates rejects non-positive errors", {
  ages <- c(2450, 2380)
  errors <- c(30, -5)

  expect_error(calibrate_dates(ages, errors), "must be positive")
})

test_that("extract_calibration_summary produces valid output", {
  ages <- c(2450, 2380)
  errors <- c(30, 25)

  result <- calibrate_dates(ages, errors)
  summary <- result$summary

  expect_true(all(!is.na(summary$median_cal_bp)))
  expect_true(all(!is.na(summary$hdi95_lower)))
  expect_true(all(!is.na(summary$hdi95_upper)))
  expect_true(all(summary$hdi95_lower > summary$hdi95_upper))  # BP convention
})

test_that("calculate_hdi produces valid intervals", {
  x <- 1:100
  density <- dnorm(x, mean = 50, sd = 10)

  hdi <- calculate_hdi(x, density, prob = 0.95)

  expect_equal(length(hdi), 2)
  expect_true(hdi[1] > hdi[2])  # BP convention: larger value first
})

test_that("get_calibration_curve loads IntCal20", {
  curve <- get_calibration_curve("intcal20")

  expect_true(is.data.frame(curve))
  expect_true("cal_bp" %in% names(curve))
  expect_true("c14_age" %in% names(curve))
  expect_true("error" %in% names(curve))
  expect_true(nrow(curve) > 0)
})

test_that("get_calibration_curve rejects invalid curve name", {
  expect_error(get_calibration_curve("invalid_curve"), "Unsupported calibration curve")
})

test_that("fit_gaussian_mixture_simple produces valid parameters", {
  x <- rnorm(100, mean = 50, sd = 10)
  weights <- rep(1/100, 100)

  mixture <- fit_gaussian_mixture_simple(x, weights, n_components = 2)

  expect_equal(length(mixture$weights), 2)
  expect_equal(length(mixture$means), 2)
  expect_equal(length(mixture$sds), 2)
  expect_true(abs(sum(mixture$weights) - 1.0) < 0.01)  # Weights sum to ~1
})

test_that("approximate_calibration_mixture works", {
  ages <- c(2450)
  errors <- c(30)
  result <- calibrate_dates(ages, errors)

  cal_date <- list(grids = list(result$calDates$grids[[1]]))
  mixture <- approximate_calibration_mixture(cal_date, n_components = 3)

  expect_equal(mixture$n_components, 3)
  expect_equal(length(mixture$weights), 3)
  expect_equal(length(mixture$means), 3)
  expect_equal(length(mixture$sds), 3)
})

test_that("prepare_stan_calibration_data produces mixture format", {
  ages <- c(2450, 2380)
  errors <- c(30, 25)
  result <- calibrate_dates(ages, errors)

  stan_data <- prepare_stan_calibration_data(result, method = "mixture", n_components = 3)

  expect_true("N" %in% names(stan_data))
  expect_true("K" %in% names(stan_data))
  expect_true("mix_weights" %in% names(stan_data))
  expect_true("mix_means" %in% names(stan_data))
  expect_true("mix_sds" %in% names(stan_data))

  expect_equal(stan_data$N, 2)
  expect_equal(stan_data$K, 3)
})

test_that("prepare_stan_calibration_data produces full curve format", {
  ages <- c(2450, 2380)
  errors <- c(30, 25)
  result <- calibrate_dates(ages, errors)

  stan_data <- prepare_stan_calibration_data(result, method = "full")

  expect_true("N_cal_curve" %in% names(stan_data))
  expect_true("cal_curve_bp" %in% names(stan_data))
  expect_true("cal_curve_c14" %in% names(stan_data))
  expect_true("cal_curve_error" %in% names(stan_data))
})

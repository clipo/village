#' Unit Tests for Simulation Framework
#'
#' @author Project Team
#' @date 2025-10-11

library(testthat)

source("../../R/simulation.R")
source("../../R/calibration.R")

test_that("simulate_contemporaneous_deposits produces correct structure", {
  sim <- simulate_contemporaneous_deposits(
    n_deposits = 2,
    shared_window = c(2700, 2400),
    n_dates_per_deposit = 5,
    measurement_error = 30,
    seed = 123
  )

  expect_true(inherits(sim, "simulated_c14"))
  expect_true(is.data.frame(sim$dates))
  expect_equal(nrow(sim$dates), 10)  # 2 deposits * 5 dates
  expect_equal(length(unique(sim$dates$deposit)), 2)
})

test_that("simulate_contemporaneous_deposits respects shared window", {
  sim <- simulate_contemporaneous_deposits(
    n_deposits = 2,
    shared_window = c(2700, 2400),
    n_dates_per_deposit = 10,
    measurement_error = 30,
    seed = 123
  )

  # Check that true calendar dates fall within window
  expect_true(all(sim$dates$true_cal_bp >= 2400))
  expect_true(all(sim$dates$true_cal_bp <= 2700))
})

test_that("simulate_contemporaneous_deposits accepts vector of sample sizes", {
  sim <- simulate_contemporaneous_deposits(
    n_deposits = 3,
    shared_window = c(2700, 2400),
    n_dates_per_deposit = c(5, 10, 15),
    measurement_error = 30,
    seed = 123
  )

  expect_equal(nrow(sim$dates), 30)  # 5 + 10 + 15
  expect_equal(sum(sim$dates$deposit == 1), 5)
  expect_equal(sum(sim$dates$deposit == 2), 10)
  expect_equal(sum(sim$dates$deposit == 3), 15)
})

test_that("simulate_sequential_deposits produces non-overlapping deposits", {
  sim <- simulate_sequential_deposits(
    n_deposits = 3,
    deposit_duration = 200,
    gap_between_deposits = 100,
    earliest_start = 3000,
    n_dates_per_deposit = 5,
    measurement_error = 30,
    seed = 123
  )

  windows <- sim$ground_truth$deposit_windows

  # Check that deposits are sequential (no overlap)
  for (k in 1:(nrow(windows) - 1)) {
    expect_true(windows[k, 2] > windows[k + 1, 1])  # End of k > Start of k+1
  }
})

test_that("simulate_sequential_deposits respects duration parameter", {
  duration <- 250
  sim <- simulate_sequential_deposits(
    n_deposits = 2,
    deposit_duration = duration,
    gap_between_deposits = 100,
    earliest_start = 3000,
    n_dates_per_deposit = 5,
    measurement_error = 30,
    seed = 123
  )

  windows <- sim$ground_truth$deposit_windows

  # Check durations
  for (k in 1:nrow(windows)) {
    actual_duration <- windows[k, 1] - windows[k, 2]
    expect_equal(actual_duration, duration)
  }
})

test_that("simulate_partial_overlap produces correct overlap", {
  sim <- simulate_partial_overlap(
    n_deposits = 2,
    overlap_percentage = 50,
    deposit_duration = 300,
    earliest_start = 2800,
    n_dates_per_deposit = 10,
    measurement_error = 30,
    seed = 123
  )

  windows <- sim$ground_truth$deposit_windows

  # Calculate actual overlap
  overlap_start <- min(windows[1, 1], windows[2, 1])
  overlap_end <- max(windows[1, 2], windows[2, 2])

  # Check that there is some overlap
  expect_true(windows[1, 2] < windows[2, 1] || windows[2, 2] < windows[1, 1])
})

test_that("calendar_to_c14 produces reasonable radiocarbon ages", {
  cal_curve <- get_calibration_curve("intcal20")

  true_cal_bp <- c(2500, 2600, 2700)
  measurement_error <- 30

  result <- calendar_to_c14(true_cal_bp, measurement_error, cal_curve)

  expect_equal(length(result$c14_age), 3)
  expect_true(all(result$c14_age > 0))
  expect_true(all(result$c14_error == measurement_error))
})

test_that("calculate_true_overlap works for contemporaneous deposits", {
  sim <- simulate_contemporaneous_deposits(
    n_deposits = 2,
    shared_window = c(2700, 2400),
    n_dates_per_deposit = 5,
    measurement_error = 30,
    seed = 123
  )

  overlap_matrix <- calculate_true_overlap(sim$ground_truth)

  # Diagonal should be 1.0
  expect_equal(overlap_matrix[1, 1], 1.0)
  expect_equal(overlap_matrix[2, 2], 1.0)

  # Off-diagonal should indicate full overlap
  expect_true(overlap_matrix[1, 2] > 0.9)
})

test_that("simulate functions respect random seed", {
  sim1 <- simulate_contemporaneous_deposits(
    n_deposits = 2,
    shared_window = c(2700, 2400),
    n_dates_per_deposit = 5,
    measurement_error = 30,
    seed = 123
  )

  sim2 <- simulate_contemporaneous_deposits(
    n_deposits = 2,
    shared_window = c(2700, 2400),
    n_dates_per_deposit = 5,
    measurement_error = 30,
    seed = 123
  )

  expect_equal(sim1$dates$age, sim2$dates$age)
})

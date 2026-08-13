source("../../R/susq_calibration.R")

test_that("the calendar grid matches the global constraint", {
  expect_equal(CAL_GRID$t0, 0); expect_equal(CAL_GRID$dt, 1)
  expect_equal(CAL_GRID$G, 2001L); expect_equal(length(CAL_GRID$t), 2001L)
})

test_that("IntCal20 is interpolated onto the grid without gaps", {
  cc <- intcal20_on_grid()
  expect_equal(nrow(cc), 2001L)
  expect_false(any(is.na(cc$c14))); expect_false(any(is.na(cc$err)))
  expect_true(all(cc$err > 0))
  expect_equal(cc$calBP, CAL_GRID$t)
})

test_that("cumulative_integral integrates a piecewise-linear function exactly", {
  # L(t) = t on a unit grid: integral from 0 to k is k^2 / 2, exactly
  # reproduced by the trapezoid rule on a linear function.
  L <- as.numeric(0:10)
  got <- cumulative_integral(L, dt = 1)
  expect_equal(got, (0:10)^2 / 2, tolerance = 1e-12)
  expect_equal(got[1], 0)
})

test_that("cumulative_integral is monotone non-decreasing for a likelihood", {
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  Phi <- cumulative_integral(L, CAL_GRID$dt)
  expect_true(all(diff(Phi) >= -1e-15))
})

test_that("calibrated density agrees with rcarbon", {
  # rcarbon is an independent implementation. Compare the normalised density
  # over the grid; agreement to 1e-3 in total variation is ample.
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  ours <- L / sum(L)

  rc <- rcarbon::calibrate(700, 30, calCurves = "intcal20",
                           normalised = TRUE, verbose = FALSE)
  g <- rc$grids[[1]]
  theirs <- rep(0, CAL_GRID$G)
  keep <- g$calBP >= 0 & g$calBP <= 2000
  theirs[match(g$calBP[keep], CAL_GRID$t)] <- g$PrDens[keep]
  theirs <- theirs / sum(theirs)

  expect_lt(0.5 * sum(abs(ours - theirs)), 1e-3)
})

test_that("exponential convolution shifts mass older and preserves total mass", {
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  Lk <- convolve_inbuilt(L, kappa = 60, dt = CAL_GRID$dt)
  expect_equal(convolve_inbuilt(L, 0, CAL_GRID$dt), L)
  # Convolution preserves total mass except for what the shift pushes past the
  # young end of the grid, which cannot be represented. That loss is one-sided
  # and, for the exponential kernel's far tail, around 3e-5 in relative terms.
  expect_lte(sum(Lk), sum(L))
  expect_equal(sum(Lk), sum(L), tolerance = 1e-4)
  # Inbuilt age means the sample grew before deposition, so the event-date
  # likelihood moves toward younger calendar dates, i.e. smaller cal BP.
  centroid <- function(v) sum(v * CAL_GRID$t) / sum(v)
  expect_lt(centroid(Lk), centroid(L))
  expect_gt(centroid(L) - centroid(Lk), 30)
  expect_lt(centroid(L) - centroid(Lk), 90)
})

test_that("convolution works at every kappa the model uses", {
  # The largest scale, 240 years, has an exponential tail longer than the grid.
  # An unclamped kernel asks for a shift past the end of the vector.
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  centroid <- function(v) sum(v * CAL_GRID$t) / sum(v)
  prev <- centroid(L)
  for (kappa in c(0, 15, 30, 60, 120, 240)) {
    Lk <- expect_no_error(convolve_inbuilt(L, kappa, CAL_GRID$dt))
    expect_length(Lk, CAL_GRID$G)
    expect_true(all(is.finite(Lk)))
    expect_true(all(Lk >= 0))
    # Larger inbuilt age moves the event date monotonically younger.
    if (kappa > 0) expect_lt(centroid(Lk), prev)
    prev <- centroid(Lk)
  }
})

test_that("pair_likelihood sharpens a wiggle-matched pair", {
  cc <- intcal20_on_grid()
  L1 <- calib_likelihood(870, 25, cc)
  L2 <- calib_likelihood(820, 25, cc)
  Lp <- pair_likelihood(L1, L2, offset = 17, dt = CAL_GRID$dt)
  sd_of <- function(v) { p <- v / sum(v); m <- sum(p * CAL_GRID$t)
                         sqrt(sum(p * (CAL_GRID$t - m)^2)) }
  expect_lt(sd_of(Lp), sd_of(L1))
  expect_true(all(Lp >= 0))
})

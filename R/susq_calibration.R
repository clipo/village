#' Calibration machinery. Builds, for each determination, the likelihood of the
#' observed radiocarbon age as a function of calendar date, on a fixed 1-year
#' grid, together with the exact integral of its piecewise-linear interpolant.
#' See spec section 3.3.
#'
#' Calendar time is cal BP throughout: larger values are older.

CAL_GRID <- list(t0 = 0, dt = 1, G = 2001L, t = 0:2000)

#' IntCal20 interpolated onto the analysis grid.
intcal20_on_grid <- function() {
  cc <- IntCal::ccurve(cc = 1)
  colnames(cc) <- c("calBP", "c14", "err")
  cc <- as.data.frame(cc)
  cc <- cc[cc$calBP <= max(CAL_GRID$t) + 50, ]
  data.frame(
    calBP = CAL_GRID$t,
    c14   = stats::approx(cc$calBP, cc$c14, xout = CAL_GRID$t, rule = 2)$y,
    err   = stats::approx(cc$calBP, cc$err, xout = CAL_GRID$t, rule = 2)$y
  )
}

#' Unnormalised p(y | t): the probability of the observed radiocarbon age given
#' a true calendar date, combining measurement and calibration-curve error.
#' `inflate` multiplies the measurement error, which is how the outlier variant
#' of a determination is built.
calib_likelihood <- function(c14_age, c14_error, curve, inflate = 1) {
  s <- sqrt((c14_error * inflate)^2 + curve$err^2)
  stats::dnorm(c14_age, mean = curve$c14, sd = s)
}

#' Convolve a likelihood with an exponential inbuilt-age distribution of mean
#' `kappa`. A sample with inbuilt age d grew d years before deposition, so the
#' event-date likelihood at t is the growth-date likelihood at t + d averaged
#' over d. Mass therefore moves toward smaller cal BP, that is toward younger
#' event dates. `kappa = 0` returns the input unchanged.
convolve_inbuilt <- function(L, kappa, dt) {
  if (kappa <= 0) return(L)
  G <- length(L)
  # Truncate the kernel where it no longer contributes, and never beyond the
  # grid: at kappa = 240 the 1-1e-8 quantile is about 4400 years, more than
  # twice the grid's span, and a shift of that size has nothing left to read.
  dmax <- ceiling(stats::qexp(1 - 1e-8, rate = 1 / kappa) / dt) * dt
  dmax <- min(dmax, (G - 1L) * dt)
  d <- seq(0, dmax, by = dt)
  w <- stats::dexp(d, rate = 1 / kappa)
  w <- w / sum(w)
  out <- numeric(G)
  for (k in seq_along(d)) {
    shift <- k - 1L
    if (shift == 0L) {
      out <- out + w[k] * L
    } else {
      idx <- seq_len(G - shift)
      out[idx] <- out[idx] + w[k] * L[idx + shift]
    }
  }
  out
}

#' Likelihood for a wiggle-matched pair. `offset` is how many years younger the
#' paired member is than the reference member, known from ring counting. The
#' returned curve is indexed by the calendar date of the reference member, so
#' the paired member's likelihood is read `offset` years toward the present.
pair_likelihood <- function(L_ref, L_other, offset, dt) {
  shift <- as.integer(round(offset / dt))
  G <- length(L_ref)
  shifted <- numeric(G)
  idx <- seq.int(shift + 1L, G)
  shifted[idx] <- L_other[idx - shift]
  L_ref * shifted
}

#' Exact integral of the piecewise-linear interpolant of L, from the grid start
#' to each knot. Element 1 is 0. This is the trapezoid rule, which is exact for
#' a piecewise-linear integrand, so the resulting Phi is piecewise quadratic and
#' continuously differentiable, which is what the Stan interpolant relies on.
cumulative_integral <- function(L, dt) {
  c(0, cumsum((utils::head(L, -1) + utils::tail(L, -1)) / 2 * dt))
}

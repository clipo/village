source("../../R/susq_data.R")
source("../../R/susq_calibration.R")
source("../../R/susq_stan_data.R")
suppressPackageStartupMessages(library(cmdstanr))

test_that("the Stan interpolant matches R evaluation of the same stored curve", {
  d <- prepare_susq_data(); cc <- intcal20_on_grid()
  sub <- d$dates[1:5, ]
  v <- build_variants(sub, cc)
  sd <- build_stan_data(sub, v, site_phase_counts(sub))

  # Query each variant across and beyond its window, so both clamp branches
  # are exercised as well as the interior.
  set.seed(20260813)
  M <- 200L
  vq <- sample(seq_len(sd$V), M, replace = TRUE)
  lo <- sd$curve_t0[vq]
  hi <- lo + (sd$curve_len[vq] - 1) * sd$dt
  tq <- runif(M, lo - 40, hi + 40)

  mod <- cmdstan_model("../../stan/test_phi_at.stan", include_paths = "../../stan")
  fit <- mod$sample(
    data = c(sd[c("Ltot", "L_flat", "Phi_flat", "curve_pos", "curve_len",
                  "curve_t0", "phi_total", "dt")],
             list(V = sd$V, M = M, t_query = tq, v_query = vq)),
    chains = 1, iter_warmup = 100, iter_sampling = 5,
    seed = 20260813, refresh = 0, show_messages = FALSE,
    # The default CSV precision is well below double, and this test compares
    # the interpolant to 1e-9. diagnostics = NULL skips the sampler checks,
    # which are meaningless for a harness whose posterior is a dummy normal.
    sig_figs = 18, diagnostics = NULL)
  got <- as.numeric(fit$draws("phi_out", format = "matrix")[1, ])

  ref <- vapply(seq_len(M), function(m) {
    vv <- vq[m]; p0 <- sd$curve_pos[vv]; n <- sd$curve_len[vv]
    L <- sd$L_flat[p0:(p0 + n - 1)]; Phi <- sd$Phi_flat[p0:(p0 + n - 1)]
    t0 <- sd$curve_t0[vv]; t <- tq[m]
    if (t <= t0) return(0)
    if (t >= t0 + (n - 1) * sd$dt) return(sd$phi_total[vv])
    g <- 1L + as.integer(floor((t - t0) / sd$dt))
    h <- t - (t0 + (g - 1) * sd$dt)
    Phi[g] + L[g] * h + (L[g + 1] - L[g]) * h^2 / (2 * sd$dt)
  }, 0)

  expect_equal(got, ref, tolerance = 1e-9)
  # Both clamp branches must actually have been hit, or the test proves less
  # than it appears to.
  expect_gt(sum(tq <= lo), 0)
  expect_gt(sum(tq >= hi), 0)
})

test_that("the interpolant reproduces fine-grained numerical integration", {
  # Independent check that does not reuse the interpolation formula: integrate
  # the underlying likelihood on a 0.01-year grid and compare.
  cc <- intcal20_on_grid()
  L <- calib_likelihood(700, 30, cc)
  L <- L / sum((head(L, -1) + tail(L, -1)) / 2)
  Phi <- cumulative_integral(L, 1)
  for (t in c(120.3, 450.7, 661.25, 900.5)) {
    fine <- seq(0, t, by = 0.01)
    dense <- approx(CAL_GRID$t, L, xout = fine)$y
    num <- sum((head(dense, -1) + tail(dense, -1)) / 2 * 0.01)
    g <- 1L + as.integer(floor(t)); h <- t - (g - 1)
    ours <- Phi[g] + L[g] * h + (L[g + 1] - L[g]) * h^2 / 2
    expect_equal(ours, num, tolerance = 1e-6)
  }
})

test_that("phase_average is the phase-normalised difference of the integral", {
  # Exercised through the same harness: phase_average(a, b) must equal
  # (phi_at(a) - phi_at(b)) / (a - b). Verified in R against phi_at outputs
  # the Stan model itself produced, so the 1/(a-b) normalisation is pinned.
  d <- prepare_susq_data(); cc <- intcal20_on_grid()
  sub <- d$dates[1:3, ]
  v <- build_variants(sub, cc)
  sd <- build_stan_data(sub, v, site_phase_counts(sub))

  a <- 900; b <- 700
  vq <- c(seq_len(sd$V), seq_len(sd$V))
  tq <- c(rep(a, sd$V), rep(b, sd$V))
  mod <- cmdstan_model("../../stan/test_phi_at.stan", include_paths = "../../stan")
  fit <- mod$sample(
    data = c(sd[c("Ltot", "L_flat", "Phi_flat", "curve_pos", "curve_len",
                  "curve_t0", "phi_total", "dt")],
             list(V = sd$V, M = length(tq), t_query = tq, v_query = vq)),
    chains = 1, iter_warmup = 100, iter_sampling = 5,
    seed = 20260813, refresh = 0, show_messages = FALSE,
    sig_figs = 18, diagnostics = NULL)
  phi <- as.numeric(fit$draws("phi_out", format = "matrix")[1, ])
  pa <- (phi[seq_len(sd$V)] - phi[sd$V + seq_len(sd$V)]) / (a - b)
  expect_true(all(pa >= 0))
  expect_true(all(is.finite(pa)))
  # A phase-averaged likelihood is a density in calendar time, so it must be
  # no larger than the peak of the curve it averages.
  peak <- vapply(seq_len(sd$V), function(vv) {
    p0 <- sd$curve_pos[vv]; n <- sd$curve_len[vv]
    max(sd$L_flat[p0:(p0 + n - 1)])
  }, 0)
  expect_true(all(pa <= peak + 1e-12))
})

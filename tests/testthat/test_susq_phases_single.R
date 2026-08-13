source("../../R/susq_data.R")
source("../../R/susq_calibration.R")
source("../../R/susq_stan_data.R")
source("../../R/susq_fit.R")
suppressPackageStartupMessages(library(posterior))

MODEL_PHASES <- "../../stan/susq_phases.stan"

#' Build a determination frame with the columns the pipeline expects.
sim_frame <- function(site, c14_age, c14_error) {
  n <- length(c14_age)
  data.frame(
    det_id = seq_len(n), site = site, context = as.character(seq_len(n)),
    lab_no = paste0(site, "-", seq_len(n)),
    c14_age = c14_age, c14_error = c14_error,
    material_class = "short_lived", class_id = 1L, reference = "sim",
    pair_lab_no = NA_character_, pair_c14_age = NA_real_,
    pair_c14_error = NA_real_, pair_offset = NA_real_, preprocess = "none",
    stringsAsFactors = FALSE)
}

test_that("the model recovers a known single occupation phase", {
  # Twelve short-lived determinations spread uniformly across a 60-year window
  # centred on 650 cal BP, with realistic errors. The posterior for the
  # boundaries must cover the truth.
  set.seed(20260813)
  cc <- intcal20_on_grid()
  truth <- list(a = 680, b = 620)
  t_true <- runif(12, truth$b, truth$a)
  mu <- approx(cc$calBP, cc$c14, xout = t_true)$y
  ce <- approx(cc$calBP, cc$err, xout = t_true)$y
  err <- rep(30, 12)
  obs <- round(rnorm(12, mu, sqrt(err^2 + ce^2)))

  sim <- sim_frame("SIM", obs, err)
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 12L, J = 1L))
  fit <- fit_susq(sd, dates = sim, model_file = MODEL_PHASES,
                  chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                  refresh = 0)

  # With a single phase the hierarchical duration prior is unidentified by
  # construction: only mu_d + sigma_d * log_dur_raw is determined, not its
  # three parts. That is expected here and disappears once 74 phases share the
  # hierarchy, so check convergence on the identified quantities.
  ident <- summarise_draws(
    fit$draws(c("phase_start_calBP", "phase_end_calBP", "phase_duration")),
    "rhat", "ess_bulk")
  expect_lt(max(ident$rhat), 1.01)
  expect_gt(min(ident$ess_bulk), 400)
  expect_equal(check_diagnostics(fit, strict = FALSE)$n_divergent, 0)

  s <- summarise_draws(fit$draws(c("phase_start_calBP", "phase_end_calBP")),
                       ~quantile(.x, c(0.025, 0.5, 0.975)))
  a_ci <- as.numeric(s[1, 2:4]); b_ci <- as.numeric(s[2, 2:4])
  expect_lt(a_ci[1], truth$a); expect_gt(a_ci[3], truth$a)
  expect_lt(b_ci[1], truth$b); expect_gt(b_ci[3], truth$b)
  # The phase must be oriented correctly in cal BP throughout.
  st <- fit$draws("phase_start_calBP", format = "matrix")
  en <- fit$draws("phase_end_calBP", format = "matrix")
  expect_true(all(st > en))
})

test_that("log_lik matches an independent R recomputation of the same mixture", {
  # The model block and generated quantities share variant_logprob(). If they
  # ever diverged, log_lik would stop describing the fitted model and every
  # LOO comparison built on it would be silently wrong. Recompute log_lik in R
  # from the posterior draws and the same stored curves, and compare.
  set.seed(11)
  cc <- intcal20_on_grid()
  # Include wood so the inbuilt-age simplex and a second outlier rate are
  # genuinely exercised, not just the degenerate short-lived branch.
  sim <- sim_frame("SIM", round(rnorm(10, 660, 40)), rep(30, 10))
  sim$material_class[6:10] <- "wood"; sim$class_id[6:10] <- 2L
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 10L, J = 1L))
  fit <- fit_susq(sd, dates = sim, model_file = MODEL_PHASES,
                  chains = 2, iter_warmup = 500, iter_sampling = 500,
                  refresh = 0, sig_figs = 18)

  a <- as.numeric(fit$draws("phase_start_calBP", format = "matrix")[, 1])
  b <- as.numeric(fit$draws("phase_end_calBP", format = "matrix")[, 1])
  rho <- fit$draws("rho", format = "matrix")
  q <- fit$draws("q", format = "matrix")
  ll_stan <- fit$draws("log_lik", format = "matrix")

  phi_r <- function(t, vv) {
    p0 <- sd$curve_pos[vv]; n <- sd$curve_len[vv]; t0 <- sd$curve_t0[vv]
    if (t <= t0) return(0)
    if (t >= t0 + (n - 1) * sd$dt) return(sd$phi_total[vv])
    L <- sd$L_flat[p0:(p0 + n - 1)]; Phi <- sd$Phi_flat[p0:(p0 + n - 1)]
    g <- 1L + as.integer(floor((t - t0) / sd$dt))
    h <- t - (t0 + (g - 1) * sd$dt)
    Phi[g] + L[g] * h + (L[g + 1] - L[g]) * h^2 / (2 * sd$dt)
  }

  draws_to_check <- sample(seq_along(a), 25)
  for (d in draws_to_check) {
    for (i in seq_len(sd$N)) {
      cls <- sd$class_id[i]
      vs <- sd$var_start[i] + seq_len(sd$n_var[i]) - 1L
      lp <- vapply(vs, function(vv) {
        lw <- if (sd$is_outlier[vv] == 1L) log(rho[d, cls])
              else if (cls == 1L) log1p(-rho[d, cls])
              else log1p(-rho[d, cls]) + log(q[d, sprintf("q[%d,%d]", cls - 1L,
                                                          sd$kappa_idx[vv])])
        A <- (phi_r(a[d], vv) - phi_r(b[d], vv)) / (a[d] - b[d])
        lw + log(A + 1e-300)
      }, 0)
      # as.numeric strips the draws_matrix dimnames that subsetting carries
      # along, which expect_equal would otherwise compare as attributes.
      expect_equal(as.numeric(ll_stan[d, i]), matrixStats::logSumExp(lp),
                   tolerance = 1e-8)
    }
  }
})

test_that("the prior predictive runs and produces plausible durations", {
  cc <- intcal20_on_grid()
  sim <- sim_frame("SIM", 650, 30)
  v <- build_variants(sim, cc)
  sd <- build_stan_data(sim, v, data.frame(site = "SIM", n = 1L, J = 1L))
  sd$prior_only <- 1L
  fit <- fit_susq(sd, dates = sim, model_file = MODEL_PHASES,
                  chains = 2, iter_warmup = 500, iter_sampling = 500,
                  refresh = 0)
  dur <- as.numeric(fit$draws("phase_duration", format = "matrix"))
  expect_gt(median(dur), 15); expect_lt(median(dur), 250)
  expect_lt(quantile(dur, 0.99), 5000)
})

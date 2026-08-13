#' Compile, initialise and fit the occupation models, and check the sampler.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })

SUSQ_SEED <- 20260813L

#' Rough calendar dates for initialisation only. This inverts the calibration
#' curve, which is not a calibration: the curve is non-monotonic, so a
#' radiocarbon age can map to several calendar dates and ties = mean averages
#' them. Good enough to start a chain near the data, and used for nothing else.
rough_calendar <- function(c14_age, curve) {
  stats::approx(curve$c14, curve$calBP, xout = c14_age, rule = 2, ties = mean)$y
}

#' Per-site centre and spread of the rough calendar dates.
site_windows <- function(stan_data, dates) {
  cc <- intcal20_on_grid()
  approx_cal <- rough_calendar(dates$c14_age, cc)
  idx <- as.character(seq_len(stan_data$S))
  centre <- as.numeric(tapply(approx_cal, stan_data$site, stats::median)[idx])
  span <- as.numeric(tapply(approx_cal, stan_data$site,
                            function(x) diff(range(x)))[idx])
  centre[is.na(centre)] <- stats::median(approx_cal)
  span[is.na(span)] <- 0
  # A determination lying wholly outside a phase makes the phase average
  # exactly zero at both boundaries, so it contributes a constant with no
  # gradient and nothing pulls the boundary back toward it. Start each phase
  # wide enough to contain its own site's spread and let the data tighten it.
  list(centre = centre, width = pmax(1.5 * span + 100, 150))
}

#' Inits for the phase-mixture model (stan/susq_phases.stan).
#' Spacing simplex placing J ordered midpoints across the calendar range,
#' clustered around a site's own data. Element k of the returned simplex is the
#' step from midpoint k-1 to midpoint k, as a fraction of the range.
spacing_simplex <- function(J, centre, spread, cal_min, cal_max) {
  rng <- cal_max - cal_min
  targets <- if (J == 1) centre else
    seq(centre - spread / 2, centre + spread / 2, length.out = J)
  targets <- pmin(pmax(targets, cal_min + 20), cal_max - 20)
  targets <- sort(targets)
  steps <- c(targets[1] - cal_min, diff(targets), cal_max - targets[J])
  steps <- pmax(steps, rng * 1e-3)
  steps / sum(steps)
}

init_phases <- function(stan_data, dates, chains) {
  w <- site_windows(stan_data, dates)
  sigma_d0 <- 0.5
  mu_d0 <- stan_data$mu_d_prior_mean
  # dur = exp(mu_d + sigma_d * log_dur_raw), so invert for the target width.
  raw_for <- function(target) (log(target) - mu_d0) / sigma_d0
  site_of_phase <- rep(seq_len(stan_data$S), times = stan_data$n_phase)
  cmin <- stan_data$cal_min; cmax <- stan_data$cal_max

  lapply(seq_len(chains), function(k) {
    jitter <- stats::rnorm(stan_data$S, 0, 25)
    ctr <- pmin(pmax(w$centre + jitter, cmin + 60), cmax - 60)
    out <- list(
      log_dur_raw = raw_for(w$width[site_of_phase]) +
                    stats::rnorm(stan_data$P, 0, 0.2),
      mu_d = mu_d0, sigma_d = sigma_d0, rho = rep(0.1, 3))
    if (stan_data$n_site_J1 > 0)
      out$mid1 <- ctr[stan_data$sites_J1]
    if (stan_data$n_site_J2 > 0)
      out$w2 <- lapply(stan_data$sites_J2, function(s)
        spacing_simplex(2, ctr[s], max(w$width[s], 200), cmin, cmax))
    if (stan_data$n_site_J3 > 0)
      out$w3 <- lapply(stan_data$sites_J3, function(s)
        spacing_simplex(3, ctr[s], max(w$width[s], 300), cmin, cmax))
    if (stan_data$n_site_J4 > 0)
      out$w4 <- lapply(stan_data$sites_J4, function(s)
        spacing_simplex(4, ctr[s], max(w$width[s], 400), cmin, cmax))
    out
  })
}

#' Inits for the sequential model (stan/susq_sequential.stan), whose single
#' parameter vector must be strictly increasing in cal BP. Sites arrive already
#' ordered youngest first, so the boundaries are laid out in that order and
#' then forced strictly increasing, which the ordered declaration requires.
init_sequential <- function(stan_data, dates, chains) {
  w <- site_windows(stan_data, dates)
  S <- stan_data$S

  lapply(seq_len(chains), function(k) {
    half <- pmax(w$width, 60) / 2
    b <- w$centre - half
    a <- w$centre + half
    bound <- as.numeric(rbind(b, a))            # b_1, a_1, b_2, a_2, ...
    bound <- bound + stats::rnorm(length(bound), 0, 5)
    bound <- sort(bound)
    # Enforce a strict gap so the ordered constraint is satisfied at t = 0.
    bound <- cummax(bound + seq_along(bound) * 1e-3)
    bound <- pmin(pmax(bound, stan_data$cal_min + 1), stan_data$cal_max - 1)
    bound <- cummax(bound + seq_along(bound) * 1e-6)
    list(bound = bound, mu_d = stan_data$mu_d_prior_mean, sigma_d = 0.5,
         rho = rep(0.1, 3))
  })
}

fit_susq <- function(stan_data, dates,
                     model_file = "stan/susq_phases.stan",
                     chains = 4, iter_warmup = 2000, iter_sampling = 2000,
                     adapt_delta = 0.95, max_treedepth = 12,
                     seed = SUSQ_SEED, refresh = 200, sig_figs = NULL) {
  if (is.null(stan_data$prior_only)) stan_data$prior_only <- 0L
  payload <- stan_data[setdiff(names(stan_data), "site_names")]

  mod <- cmdstan_model(model_file, include_paths = dirname(model_file))

  # Pick inits from the model's own parameter block rather than from the file
  # name, so a renamed or added model cannot silently receive wrong inits.
  pars <- mod$variables()$parameters
  init <- if ("bound" %in% names(pars)) {
    init_sequential(stan_data, dates, chains)
  } else {
    init_phases(stan_data, dates, chains)
  }

  mod$sample(
    data = payload,
    chains = chains,
    parallel_chains = min(chains, max(1L, parallel::detectCores() - 1L)),
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    seed = seed, refresh = refresh, sig_figs = sig_figs, init = init
  )
}

#' Sampler diagnostics against the thresholds in spec section 6.3.
check_diagnostics <- function(fit, strict = TRUE) {
  pars <- intersect(c("mid", "dur", "mu_d", "sigma_d", "rho"),
                    fit$metadata()$stan_variables)
  s <- posterior::summarise_draws(fit$draws(pars), "rhat", "ess_bulk", "ess_tail")
  dg <- fit$diagnostic_summary(quiet = TRUE)

  out <- data.frame(
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(s$ess_tail, na.rm = TRUE),
    n_divergent = sum(dg$num_divergent),
    n_max_treedepth = sum(dg$num_max_treedepth),
    ebfmi_min = min(dg$ebfmi)
  )
  if (strict) {
    if (out$max_rhat >= 1.01) stop("R-hat ", round(out$max_rhat, 4), " >= 1.01")
    if (out$min_ess_bulk < 400) stop("bulk ESS ", round(out$min_ess_bulk), " < 400")
    if (out$min_ess_tail < 400) stop("tail ESS ", round(out$min_ess_tail), " < 400")
    if (out$n_divergent > 0) stop(out$n_divergent, " divergent transitions")
    if (out$ebfmi_min < 0.3) stop("E-BFMI ", round(out$ebfmi_min, 3), " < 0.3")
  }
  out
}

#' Compile, initialise and fit the occupation models, and check the sampler.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })

SUSQ_SEED <- 20260813L

#' Initialise each site's oldest phase midpoint near the calibrated dates that
#' belong to it. Without this, a phase can start far from all data, where every
#' phase-average is zero and the gradient carries no information.
init_from_data <- function(stan_data, dates, chains) {
  cc <- intcal20_on_grid()
  # Inverting the calibration curve is only a rough starting point, not a
  # calibration. The curve is non-monotonic, so a radiocarbon age can map to
  # several calendar dates; ties = mean averages them and keeps the warning
  # out of otherwise pristine test output.
  approx_cal <- stats::approx(cc$c14, cc$calBP, xout = dates$c14_age,
                              rule = 2, ties = mean)$y
  by_site <- tapply(approx_cal, stan_data$site, stats::median)
  centre <- as.numeric(by_site[as.character(seq_len(stan_data$S))])
  centre[is.na(centre)] <- stats::median(approx_cal)

  # A determination lying wholly outside a phase contributes a constant with
  # zero gradient, because the phase average is then exactly zero at both
  # boundaries. Nothing pulls the boundary back toward that determination, so
  # a phase initialised narrower than its own dates can strand the chain.
  # Start each phase wide enough to contain the site's calibrated spread and
  # let the likelihood tighten it.
  spread <- tapply(approx_cal, stan_data$site, function(x) diff(range(x)))
  span <- as.numeric(spread[as.character(seq_len(stan_data$S))])
  span[is.na(span)] <- 0
  init_dur_site <- pmax(1.5 * span + 100, 150)

  sigma_d0 <- 0.5
  mu_d0 <- stan_data$mu_d_prior_mean
  # dur = exp(mu_d + sigma_d * log_dur_raw), so invert for the target width.
  raw_for <- function(target) (log(target) - mu_d0) / sigma_d0

  # Every phase of a site starts at that site's width.
  site_of_phase <- rep(seq_len(stan_data$S), times = stan_data$n_phase)
  init_dur <- init_dur_site[site_of_phase]

  lapply(seq_len(chains), function(k) {
    list(
      mid_first   = pmin(pmax(centre + stats::rnorm(stan_data$S, 0, 25),
                              stan_data$cal_min + 50), stan_data$cal_max - 50),
      log_gap     = stats::rnorm(max(stan_data$P - stan_data$S, 0), log(150), 0.2),
      log_dur_raw = raw_for(init_dur) + stats::rnorm(stan_data$P, 0, 0.2),
      mu_d        = mu_d0,
      sigma_d     = sigma_d0,
      rho         = rep(0.1, 3)
    )
  })
}

fit_susq <- function(stan_data, dates,
                     model_file = "stan/susq_phases.stan",
                     chains = 4, iter_warmup = 2000, iter_sampling = 2000,
                     adapt_delta = 0.95, max_treedepth = 12,
                     seed = SUSQ_SEED, refresh = 200, sig_figs = NULL) {
  if (is.null(stan_data$prior_only)) stan_data$prior_only <- 0L
  payload <- stan_data[setdiff(names(stan_data), "site_names")]

  stan_dir <- dirname(model_file)
  mod <- cmdstan_model(model_file, include_paths = stan_dir)
  mod$sample(
    data = payload,
    chains = chains,
    parallel_chains = min(chains, max(1L, parallel::detectCores() - 1L)),
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    seed = seed, refresh = refresh, sig_figs = sig_figs,
    init = init_from_data(stan_data, dates, chains)
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

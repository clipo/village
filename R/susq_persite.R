#' Fit sites independently and recombine.
#'
#' Why: the joint model couples every site through the shared outlier rates
#' `rho` and inbuilt-age simplexes `q`. Bisection showed the cost of that
#' coupling directly, with one phase per site and no duration hierarchy, so
#' that nothing else could be responsible:
#'
#'   Otsiningo_Market alone  (3 short-lived)      0.0% divergent, R-hat 1.00
#'   Lords_Wells alone       (4 maize)            0.0% divergent, R-hat 1.00
#'   Street alone            (8 wood)             0.0% divergent, R-hat 1.00
#'   Port Dickinson alone    (4 wood)             0.0% divergent, R-hat 1.00
#'   Roundtop alone          (6 short + 5 wood)   0.0% divergent, R-hat 1.00
#'   4 sites                                     26.7% divergent, R-hat 1.50
#'   10 sites                                    87.8% divergent, R-hat 1.73
#'
#' Every site samples perfectly alone, including wood-only and mixed-material
#' sites, and divergences scale with the number of sites. Every phase position
#' depends on rho and q, and rho and q depend on every site, which ties the
#' whole assemblage into one strongly coupled multimodal posterior.
#'
#' What this costs: the inbuilt-age and outlier parameters are no longer pooled
#' across sites. A site carrying both short-lived and wood determinations still
#' estimates its own offset from its own contrast; a site carrying only one
#' material class returns the prior for that offset, and is reported as such.
#' Since the design already recorded that the wood offset rests on only
#' Roundtop and Bates, the pooling was carrying little weight to begin with.
#' Short-lived determinations are unaffected either way, being fixed at
#' kappa = 0.
#'
#' What this preserves: within each site, the full model is unchanged. Multiple
#' occupation phases, marginalised phase membership, marginalised calendar
#' dates, inbuilt age, outliers and the uniform-phase normalisation all still
#' apply exactly as specified.

#' Fit every site on its own. Returns one entry per site holding the fit, the
#' Stan data, and the sampler diagnostics.
fit_all_sites <- function(dates, curve, phase_counts,
                          model_file = "stan/susq_phases.stan",
                          mu_d_prior = log(60), chains = 4,
                          iter_warmup = 1000, iter_sampling = 1000,
                          verbose = TRUE) {
  sites <- sort(unique(dates$site))
  out <- vector("list", length(sites))
  names(out) <- sites

  for (s in sites) {
    sub <- dates[dates$site == s, , drop = FALSE]
    sub$det_id <- seq_len(nrow(sub))
    v <- build_variants(sub, curve)
    pc <- phase_counts[phase_counts$site == s, ]
    pc$n <- nrow(sub)
    sd <- build_stan_data(sub, v, pc, mu_d_prior = mu_d_prior,
                          hier_duration = 0L)
    fit <- fit_susq(sd, dates = sub, model_file = model_file, chains = chains,
                    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                    refresh = 0)
    dg <- check_diagnostics(fit, strict = FALSE)
    if (verbose) {
      cat(sprintf("%-34s n=%2d J=%d  div=%3d  rhat=%.3f  ESS=%5.0f\n",
                  s, nrow(sub), pc$J, dg$n_divergent, dg$max_rhat,
                  dg$min_ess_bulk))
    }
    out[[s]] <- list(site = s, fit = fit, stan_data = sd,
                     diagnostics = dg, n = nrow(sub), J = pc$J)
  }
  out
}

#' Site-by-year occupancy across all sites, drawn from the independent fits.
#'
#' The fits are independent, so draw d of one site carries no relationship to
#' draw d of another. Combining them index-wise is therefore valid: the joint
#' posterior over sites factorises, which is exactly the assumption that makes
#' fitting them separately correct in the first place.
occupancy_persite <- function(site_fits, report_dt = 10, n_draws = NULL) {
  t_grid <- seq(0, 2000, by = report_dt)
  sites <- names(site_fits)
  D <- min(vapply(site_fits, function(x)
    nrow(x$fit$draws("phase_start_calBP", format = "matrix")), 0L))
  if (!is.null(n_draws)) D <- min(D, n_draws)

  O <- array(FALSE, dim = c(D, length(sites), length(t_grid)))
  for (si in seq_along(sites)) {
    f <- site_fits[[si]]$fit
    st <- f$draws("phase_start_calBP", format = "matrix")[seq_len(D), , drop = FALSE]
    en <- f$draws("phase_end_calBP", format = "matrix")[seq_len(D), , drop = FALSE]
    ac <- f$draws("active", format = "matrix")[seq_len(D), , drop = FALSE]
    for (p in seq_len(ncol(st))) {
      cover <- outer(st[, p], t_grid, ">=") & outer(en[, p], t_grid, "<=")
      O[, si, ] <- O[, si, ] | (cover & (ac[, p] == 1))
    }
  }
  list(t = t_grid, O = O, sites = sites, dt = report_dt)
}

#' Derived quantities computed in R from the posterior draws.
#' Keeping these out of Stan keeps the fit output small and lets the same draws
#' answer questions posed after the fact, such as restricting the count to
#' well-dated sites. See spec section 4.
#'
#' Every function here takes an occupancy object: a list with `t` (reporting
#' grid in cal BP), `O` (draws x sites x grid logical array), `sites`, and
#' `dt`. Both the joint fit and the per-site fits produce that same shape.

suppressPackageStartupMessages(library(posterior))

#' Site-by-calendar-year occupancy for every draw, from raw draw matrices.
#' A site is occupied at t when at least one of its active phases covers t.
occupancy_from_arrays <- function(start, end, active, stan_data, report_dt = 10) {
  t_grid <- seq(0, 2000, by = report_dt)
  D <- nrow(start); S <- stan_data$S; G <- length(t_grid)
  O <- array(FALSE, dim = c(D, S, G))
  for (s in seq_len(S)) {
    p_idx <- stan_data$phase_start[s] + seq_len(stan_data$n_phase[s]) - 1L
    for (p in p_idx) {
      cover <- outer(start[, p], t_grid, ">=") & outer(end[, p], t_grid, "<=")
      O[, s, ] <- O[, s, ] | (cover & matrix(active[, p] == 1, D, G))
    }
  }
  list(t = t_grid, O = O, sites = stan_data$site_names, dt = report_dt)
}

occupancy_draws <- function(fit, stan_data, report_dt = 10) {
  occupancy_from_arrays(
    fit$draws("phase_start_calBP", format = "matrix"),
    fit$draws("phase_end_calBP", format = "matrix"),
    fit$draws("active", format = "matrix"),
    stan_data, report_dt)
}

#' Posterior number of simultaneously occupied sites at each calendar year.
count_curve <- function(occ, sites = NULL) {
  keep <- if (is.null(sites)) seq_along(occ$sites) else which(occ$sites %in% sites)
  counts <- apply(occ$O[, keep, , drop = FALSE], c(1, 3), sum)
  q <- function(p) apply(counts, 2, stats::quantile, probs = p, names = FALSE)
  data.frame(
    calBP = occ$t, AD = 1950 - occ$t,
    mean = colMeans(counts),
    median = q(0.5),
    lo50 = q(0.25), hi50 = q(0.75),
    lo95 = q(0.025), hi95 = q(0.975)
  )
}

#' Years of joint occupation for every pair of sites.
#' Reported as a duration with a credible interval, and as the probability of
#' at least 25 years of overlap. A plain indicator of any overlap saturates at
#' 1 for any pair of fuzzy boundaries and carries no information.
pairwise_overlap <- function(occ, site_names = NULL) {
  if (is.null(site_names)) site_names <- occ$sites
  S <- length(site_names); rows <- list(); k <- 0L
  for (i in seq_len(S - 1)) for (j in (i + 1):S) {
    years <- rowSums(occ$O[, i, ] & occ$O[, j, ]) * occ$dt
    k <- k + 1L
    rows[[k]] <- data.frame(
      site_a = site_names[i], site_b = site_names[j],
      median_years = stats::median(years),
      lo95 = unname(stats::quantile(years, 0.025)),
      hi95 = unname(stats::quantile(years, 0.975)),
      p_ge_25 = mean(years >= 25),
      p_any = mean(years > 0),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

#' Per-phase summary from independently fitted sites, with boundary and
#' duration quantiles conditional on the phase being active. An inactive phase
#' has no archaeological meaning, so summarising it unconditionally would mix
#' real episodes with empty ones.
phase_summary_persite <- function(site_fits) {
  rows <- list(); k <- 0L
  for (nm in names(site_fits)) {
    f <- site_fits[[nm]]$fit
    st <- f$draws("phase_start_calBP", format = "matrix")
    en <- f$draws("phase_end_calBP", format = "matrix")
    du <- f$draws("phase_duration", format = "matrix")
    ac <- f$draws("active", format = "matrix")
    for (p in seq_len(ncol(st))) {
      sel <- ac[, p] == 1
      k <- k + 1L
      qs <- function(v) if (sum(sel) < 20) rep(NA_real_, 3) else
        unname(stats::quantile(v[sel], c(0.025, 0.5, 0.975)))
      a <- qs(st[, p]); b <- qs(en[, p]); dd <- qs(du[, p])
      rows[[k]] <- data.frame(
        site = nm, phase = p, n = site_fits[[nm]]$n,
        p_active = mean(sel),
        start_lo = a[1], start_med = a[2], start_hi = a[3],
        end_lo = b[1], end_med = b[2], end_hi = b[3],
        dur_lo = dd[1], dur_med = dd[2], dur_hi = dd[3],
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

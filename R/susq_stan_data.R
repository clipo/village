#' Build the ragged likelihood-variant structure and assemble Stan data.
#' See spec sections 3.3 to 3.6.
#'
#' Each determination contributes several likelihood variants: one per inbuilt
#' age scale, plus one outlier variant. Each variant is stored only over the
#' stretch of calendar grid carrying meaningful mass, which keeps the payload
#' handed to CmdStan around 5 MB rather than 25 MB.

KAPPAS <- c(0, 15, 30, 60, 120, 240)
OUTLIER_INFLATE <- 5
WINDOW_REL_TOL <- 1e-10
WINDOW_PAD <- 25L

#' Trim a curve to the stretch carrying meaningful mass, and return the window
#' together with the exact cumulative integral inside it.
window_curve <- function(L, dt, t_grid) {
  thresh <- max(L) * WINDOW_REL_TOL
  keep <- which(L > thresh)
  lo <- max(1L, min(keep) - WINDOW_PAD)
  hi <- min(length(L), max(keep) + WINDOW_PAD)
  Lw <- L[lo:hi]
  Phi <- cumulative_integral(Lw, dt)
  list(L = Lw, Phi = Phi, t0 = t_grid[lo], len = length(Lw),
       phi_total = Phi[length(Phi)])
}

build_variants <- function(dates, curve) {
  curves <- list(); index <- list(); vid <- 0L

  for (i in seq_len(nrow(dates))) {
    row <- dates[i, ]
    base <- calib_likelihood(row$c14_age, row$c14_error, curve)
    out  <- calib_likelihood(row$c14_age, row$c14_error, curve,
                             inflate = OUTLIER_INFLATE)

    # A wiggle-matched pair contributes the product of its members'
    # likelihoods with the known ring offset applied, in both forms.
    if (!is.na(row$pair_offset)) {
      p_base <- calib_likelihood(row$pair_c14_age, row$pair_c14_error, curve)
      p_out  <- calib_likelihood(row$pair_c14_age, row$pair_c14_error, curve,
                                 inflate = OUTLIER_INFLATE)
      base <- pair_likelihood(base, p_base, row$pair_offset, CAL_GRID$dt)
      out  <- pair_likelihood(out,  p_out,  row$pair_offset, CAL_GRID$dt)
    }

    # One normaliser for every variant of this determination, so the relative
    # scale between variants survives. Normalising each variant separately
    # would destroy the evidence distinguishing an outlier from a normal
    # determination. See spec section 3.4.
    norm <- sum((utils::head(base, -1) + utils::tail(base, -1)) / 2 * CAL_GRID$dt)
    stopifnot(norm > 0)
    base <- base / norm
    out  <- out  / norm

    n_kappa <- if (row$material_class == "short_lived") 1L else length(KAPPAS)
    for (k in seq_len(n_kappa)) {
      Lk <- convolve_inbuilt(base, KAPPAS[k], CAL_GRID$dt)
      vid <- vid + 1L
      curves[[vid]] <- window_curve(Lk, CAL_GRID$dt, CAL_GRID$t)
      index[[vid]] <- data.frame(variant_id = vid, det_id = row$det_id,
                                 kappa_idx = k, is_outlier = 0L)
    }
    vid <- vid + 1L
    curves[[vid]] <- window_curve(out, CAL_GRID$dt, CAL_GRID$t)
    index[[vid]] <- data.frame(variant_id = vid, det_id = row$det_id,
                               kappa_idx = 0L, is_outlier = 1L)
  }
  list(curves = curves, index = do.call(rbind, index))
}

build_stan_data <- function(dates, variants, phase_counts,
                            mu_d_prior = log(60),
                            model = c("multi", "indep", "single", "sequential")) {
  model <- match.arg(model)

  sites <- sort(unique(dates$site))
  site_id <- match(dates$site, sites)
  J <- phase_counts$J[match(sites, phase_counts$site)]

  if (model == "indep") J <- rep(1L, length(sites))

  if (model == "sequential") {
    # One phase per site, and the sites are ordered youngest first so that the
    # model's single ordered boundary vector runs in the same direction. The
    # ordering is chosen from the data, which is the ordering most favourable
    # to the serial hypothesis. See spec section 5.
    J <- rep(1L, length(sites))
    cc <- intcal20_on_grid()
    approx_cal <- stats::approx(cc$c14, cc$calBP, xout = dates$c14_age, rule = 2)$y
    med <- tapply(approx_cal, dates$site, stats::median)
    sites <- names(sort(med))
    site_id <- match(dates$site, sites)
  }

  if (model == "single") {
    sites <- "ALL"; site_id <- rep(1L, nrow(dates)); J <- 1L
  }

  S <- length(J); P <- sum(J)
  phase_start <- as.integer(cumsum(c(1L, utils::head(J, -1))))

  idx <- variants$index[order(variants$index$det_id, variants$index$variant_id), ]
  n_var <- as.integer(table(factor(idx$det_id, levels = dates$det_id)))
  var_start <- as.integer(cumsum(c(1L, utils::head(n_var, -1))))

  cvs <- variants$curves[idx$variant_id]
  curve_len <- vapply(cvs, function(x) as.integer(x$len), 1L)
  curve_pos <- as.integer(cumsum(c(1L, utils::head(curve_len, -1))))

  list(
    N = nrow(dates), S = S, P = P,
    site = site_id,
    phase_start = phase_start, n_phase = as.integer(J),
    class_id = as.integer(dates$class_id),

    V = nrow(idx),
    var_start = var_start, n_var = n_var,
    kappa_idx = as.integer(idx$kappa_idx),
    is_outlier = as.integer(idx$is_outlier),

    Ltot = sum(curve_len),
    L_flat = unlist(lapply(cvs, `[[`, "L"), use.names = FALSE),
    Phi_flat = unlist(lapply(cvs, `[[`, "Phi"), use.names = FALSE),
    curve_pos = curve_pos, curve_len = curve_len,
    curve_t0 = vapply(cvs, function(x) as.numeric(x$t0), 0),
    phi_total = vapply(cvs, function(x) as.numeric(x$phi_total), 0),
    dt = CAL_GRID$dt,

    n_kappa = length(KAPPAS),
    cal_min = min(CAL_GRID$t), cal_max = max(CAL_GRID$t),

    # as.integer() matters even when these are empty: a zero-length vector of
    # the wrong type serialises to JSON in a form CmdStan rejects for an int
    # array, and every configuration except "multi" leaves both of them empty.
    sites_J2 = as.integer(which(J == 2L)), n_site_J2 = sum(J == 2L),
    sites_J4 = as.integer(which(J == 4L)), n_site_J4 = sum(J == 4L),

    mu_d_prior_mean = mu_d_prior,
    site_names = sites
  )
}

#' Choose how many occupation phases each site gets.
#'
#' The model infers which phase a determination belongs to, but the number of
#' candidate phases is fixed before fitting. Leaving spare components that no
#' determination occupies makes the posterior geometry intractable: on
#' simulated data with two known episodes, four candidate phases produced 78
#' to 95 per cent divergent transitions at every Dirichlet concentration and
#' step size tried, while two candidate phases produced none. So the candidate
#' count is set from the data, and verified afterwards by refitting the
#' well-dated sites with one extra phase and comparing by LOO.
#'
#' This is a data-dependent model-selection step and is reported as such.

#' Posterior mean and standard deviation of each determination's calibrated
#' distribution, computed from the same likelihood the model uses.
calibrated_moments <- function(dates, curve) {
  t <- curve$calBP
  out <- vapply(seq_len(nrow(dates)), function(i) {
    L <- calib_likelihood(dates$c14_age[i], dates$c14_error[i], curve)
    p <- L / sum(L)
    m <- sum(p * t)
    c(mean = m, sd = sqrt(sum(p * (t - m)^2)))
  }, c(mean = 0, sd = 0))
  data.frame(det_id = dates$det_id, site = dates$site,
             cal_mean = out["mean", ], cal_sd = out["sd", ])
}

#' Number of occupation episodes a site's dates actually separate into.
#'
#' Complete-linkage clustering on the calibrated means, cut at `gap_years`.
#' The cut height is deliberately larger than typical calibration uncertainty,
#' around 40 to 60 calendar years for this assemblage, so only genuinely
#' separated groups split rather than every wiggle in the curve.
#'
#' The binding requirement is that no candidate phase starts out empty, not
#' that each holds many determinations, so the cluster count is capped only by
#' `max_J` and by the number of determinations. A phase resting on one or two
#' dates has a prior-driven duration, which is already true of the single-date
#' sites, and every result is reported with and without those sites.
estimate_phase_counts <- function(dates, curve, gap_years = 150,
                                  max_J = 4L) {
  cm <- calibrated_moments(dates, curve)
  sites <- sort(unique(dates$site))

  rows <- lapply(sites, function(s) {
    x <- cm$cal_mean[cm$site == s]
    n <- length(x)
    cap <- min(max_J, n)
    if (n == 1L) {
      k <- 1L
    } else {
      h <- stats::hclust(stats::dist(x), method = "complete")
      k <- length(unique(stats::cutree(h, h = gap_years)))
      k <- min(k, cap)
    }
    data.frame(site = s, n = n, J = as.integer(k),
               cap = as.integer(cap),
               cal_span = round(diff(range(x))), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[order(out$site), ]
}

#' Figures for the contemporaneity analysis. See spec section 7.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

#' Grey bands marking the IntCal20 reversals, in AD.
reversal_bands <- function() {
  d <- data.frame(xmin = 1950 - CURVE_REVERSALS$start_calBP,
                  xmax = 1950 - CURVE_REVERSALS$end_calBP)
  geom_rect(data = d, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey91")
}

#' Figure 1. Posterior number of simultaneously occupied sites through time.
plot_count_through_time <- function(cnt, cnt_subset = NULL,
                                    cnt_nopalimpsest = NULL,
                                    xlim = c(600, 1800)) {
  lines <- cbind(cnt, series = "all 34 deposits")
  if (!is.null(cnt_subset))
    lines <- rbind(lines, cbind(cnt_subset, series = "deposits with 3+ dates"))
  if (!is.null(cnt_nopalimpsest))
    lines <- rbind(lines, cbind(cnt_nopalimpsest,
                                series = "excluding multi-component deposits"))
  lines$series <- factor(lines$series, levels = c(
    "all 34 deposits", "deposits with 3+ dates",
    "excluding multi-component deposits"))

  ggplot(cnt, aes(AD)) +
    reversal_bands() +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = OKABE_ITO[6], alpha = 0.18) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = OKABE_ITO[6], alpha = 0.35) +
    geom_line(data = lines, aes(y = median, colour = series, linetype = series),
              linewidth = 0.75) +
    scale_colour_manual(values = c(OKABE_ITO[6], OKABE_ITO[2], OKABE_ITO[7]),
                        name = NULL) +
    scale_linetype_manual(values = c("solid", "22", "44"), name = NULL) +
    annotate("text", x = xlim[1] + 30, y = Inf, hjust = 0, vjust = 1.6,
             size = 2.6, colour = "grey35", family = susq_font(),
             label = "grey bands: IntCal20 reversals, where calibration is multimodal") +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    scale_x_continuous(breaks = seq(600, 1800, 200)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    guides(colour = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
    labs(x = "Year AD", y = "Simultaneously occupied deposits",
         caption = "Ribbons show 50% and 95% credible intervals for all 34 deposits.") +
    theme_susq() +
    theme(plot.caption = element_text(size = 7, colour = "grey35", hjust = 0))
}

#' Figure 2. Occupation intervals, ordered by median start date.
plot_occupation_intervals <- function(ps) {
  d <- ps[!is.na(ps$start_med), ]
  d$ad_start <- 1950 - d$start_med; d$ad_end <- 1950 - d$end_med
  d$ad_lo <- 1950 - d$start_lo;     d$ad_hi <- 1950 - d$end_hi
  d$kind <- ifelse(d$site %in% PALIMPSEST_SITES, "multi-component (not one occupation)",
            ifelse(d$n <= 2, "duration prior-driven (1-2 dates)", "single occupation"))
  d <- d[order(d$start_med), ]
  d$label <- factor(sprintf("%s (n=%d)", d$site, d$n),
                    levels = sprintf("%s (n=%d)", d$site, d$n))

  ggplot(d, aes(y = label, colour = kind)) +
    reversal_bands() +
    geom_segment(aes(x = ad_lo, xend = ad_hi, yend = label),
                 linewidth = 0.9, alpha = 0.35) +
    geom_segment(aes(x = ad_start, xend = ad_end, yend = label),
                 linewidth = 2.6) +
    scale_colour_manual(values = c(
      "single occupation" = OKABE_ITO[6],
      "duration prior-driven (1-2 dates)" = OKABE_ITO[3],
      "multi-component (not one occupation)" = OKABE_ITO[7]), name = NULL) +
    guides(colour = guide_legend(nrow = 3)) +
    labs(x = "Year AD", y = NULL) +
    theme_susq() + theme(axis.text.y = element_text(size = 7))
}

#' Figure 3. Probability that a pair of deposits shared at least 25 years.
#' This is the figure that answers which deposits are contemporaneous.
plot_contemporaneity_matrix <- function(po, order_sites) {
  d <- rbind(po, transform(po, site_a = po$site_b, site_b = po$site_a))
  d$site_a <- factor(d$site_a, levels = order_sites)
  d$site_b <- factor(d$site_b, levels = rev(order_sites))
  d$band <- cut(d$p_ge_25, c(-0.01, 0.05, 0.5, 0.95, 1.01),
                labels = c("not contemporaneous (P<0.05)", "unlikely (0.05-0.5)",
                           "likely (0.5-0.95)", "contemporaneous (P>0.95)"))
  ggplot(d, aes(site_a, site_b, fill = band)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = c(OKABE_ITO[6], "#BBD9EE", "#F6D6A8", OKABE_ITO[7]),
                      name = NULL, drop = FALSE) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = NULL) + coord_equal() +
    theme_susq() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5.5),
          panel.grid = element_blank(), axis.line = element_blank())
}

#' Figure 4. Median years of overlap for every pair.
plot_overlap_years <- function(po, order_sites) {
  d <- rbind(po, transform(po, site_a = po$site_b, site_b = po$site_a))
  d$site_a <- factor(d$site_a, levels = order_sites)
  d$site_b <- factor(d$site_b, levels = rev(order_sites))
  d$show <- ifelse(d$median_years == 0, NA, d$median_years)
  ggplot(d, aes(site_a, site_b, fill = show)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_viridis_c(name = "Median overlap (years)", na.value = "grey93",
                         option = "viridis") +
    labs(x = NULL, y = NULL) + coord_equal() +
    theme_susq() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5.5),
          panel.grid = element_blank(), axis.line = element_blank())
}

#' Figure 5. Duration posteriors against the prior, so it is visible which
#' deposits are constrained by their own dates and which return the prior.
plot_duration_posteriors <- function(site_fits, ps, mu_d_prior = log(60),
                                     sigma_d = 0.6) {
  long <- do.call(rbind, lapply(names(site_fits), function(nm) {
    du <- as.numeric(site_fits[[nm]]$fit$draws("phase_duration", format = "matrix")[, 1])
    data.frame(site = nm, n = site_fits[[nm]]$n, duration = du)
  }))
  ord <- tapply(long$duration, long$site, median)
  long$label <- factor(sprintf("%s (n=%d)", long$site, long$n),
                       levels = sprintf("%s (n=%d)",
                                        names(sort(ord)),
                                        vapply(names(sort(ord)),
                                               function(s) long$n[long$site == s][1], 0)))
  prior_q <- qlnorm(c(0.025, 0.5, 0.975), mu_d_prior, sigma_d)
  ggplot(long, aes(duration, label)) +
    geom_vline(xintercept = prior_q, linetype = c("22", "solid", "22"),
               colour = "grey45", linewidth = 0.35) +
    ggridges::geom_density_ridges(scale = 2.2, fill = OKABE_ITO[4],
                                  colour = NA, alpha = 0.8,
                                  rel_min_height = 0.02) +
    scale_x_continuous(trans = "log10", breaks = c(10, 30, 100, 300, 1000),
                       limits = c(5, 2000)) +
    labs(x = "Occupation duration (years, log scale)", y = NULL) +
    theme_susq() + theme(axis.text.y = element_text(size = 7))
}

#' Figure 6. The underlying evidence: calibrated determinations by deposit,
#' coloured by material class.
plot_calibrated_dates <- function(dates, curve, order_sites) {
  rows <- lapply(seq_len(nrow(dates)), function(i) {
    L <- calib_likelihood(dates$c14_age[i], dates$c14_error[i], curve)
    L <- L / max(L)
    keep <- L > 0.02
    data.frame(row = i, site = dates$site[i],
               material_class = dates$material_class[i],
               AD = 1950 - curve$calBP[keep], dens = L[keep])
  })
  df <- do.call(rbind, rows)
  df$site <- factor(df$site, levels = order_sites)
  ggplot(df, aes(AD, factor(row), height = dens, fill = material_class)) +
    reversal_bands() +
    ggridges::geom_ridgeline(scale = 2.2, colour = NA, alpha = 0.85) +
    facet_grid(rows = vars(site), scales = "free_y", space = "free_y",
               switch = "y") +
    scale_fill_manual(values = c(short_lived = OKABE_ITO[4],
                                 wood = OKABE_ITO[7],
                                 indeterminate = OKABE_ITO[3]), name = NULL) +
    coord_cartesian(xlim = c(200, 1800)) +
    labs(x = "Year AD", y = NULL) + theme_susq() +
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 5, hjust = 1),
          panel.spacing.y = unit(0.05, "lines"))
}

#' Figure 7. Sampler diagnostics per deposit.
plot_diagnostics <- function(dg) {
  d <- dg
  d$site <- factor(d$site, levels = d$site[order(d$max_rhat)])
  a <- ggplot(d, aes(max_rhat, site)) +
    geom_vline(xintercept = 1.01, linetype = "22", colour = OKABE_ITO[7]) +
    geom_point(aes(colour = n_divergent > 0), size = 1.4) +
    scale_colour_manual(values = c(`FALSE` = OKABE_ITO[6], `TRUE` = OKABE_ITO[7]),
                        labels = c("none", "some"), name = "divergences") +
    labs(x = "R-hat", y = NULL) + theme_susq() +
    theme(axis.text.y = element_text(size = 6))
  b <- ggplot(d, aes(pmax(n_divergent, 0.5), site)) +
    geom_vline(xintercept = 1, linetype = "22", colour = "grey60") +
    geom_point(colour = OKABE_ITO[7], size = 1.4) +
    scale_x_continuous(trans = "log10") +
    labs(x = "Divergent transitions of 4000", y = NULL) + theme_susq() +
    theme(axis.text.y = element_blank())
  patchwork::wrap_plots(a, b, nrow = 1, widths = c(1.5, 1)) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

#' Figure 8. Time slices: probability each deposit was occupied in each
#' interval. This is the non-spatial form of the time-slice map, and it is
#' built from the same per-slice occupancy the map uses, so the two always
#' agree.
slice_occupancy <- function(occ, slice_years = 50, ad_range = c(900, 1750)) {
  edges <- seq(ad_range[1], ad_range[2], by = slice_years)
  ad <- 1950 - occ$t
  rows <- list(); k <- 0L
  for (e in seq_len(length(edges) - 1L)) {
    inside <- ad >= edges[e] & ad < edges[e + 1L]
    if (!any(inside)) next
    for (si in seq_along(occ$sites)) {
      # Occupied in the slice if occupied at any reporting year inside it.
      p <- mean(apply(occ$O[, si, inside, drop = FALSE], 1, any))
      k <- k + 1L
      rows[[k]] <- data.frame(site = occ$sites[si], slice_start = edges[e],
                              slice_end = edges[e + 1L], p_occupied = p,
                              stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

plot_time_slices <- function(sl, order_sites) {
  sl$site <- factor(sl$site, levels = rev(order_sites))
  sl$lab <- sprintf("%d-%d", sl$slice_start, sl$slice_end)
  sl$lab <- factor(sl$lab, levels = unique(sl$lab[order(sl$slice_start)]))
  ggplot(sl, aes(lab, site, fill = p_occupied)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_viridis_c(name = "P(occupied)", limits = c(0, 1),
                         option = "mako", direction = -1) +
    labs(x = "Year AD", y = NULL) + theme_susq() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          panel.grid = element_blank(), axis.line = element_blank())
}

#' Figure 9. Time-slice map. Requires a coordinates table with columns
#' `site`, `easting`, `northing` (or `lon`, `lat`), which the radiocarbon
#' spreadsheet does not contain. Point area shows the probability the deposit
#' was occupied during the slice, so a reader can see which deposits are
#' plausibly contemporaneous and where they sit relative to one another.
plot_time_slice_map <- function(sl, coords, slice_years = 50) {
  xy <- if (all(c("easting", "northing") %in% names(coords)))
    data.frame(site = coords$site, x = coords$easting, y = coords$northing)
  else data.frame(site = coords$site, x = coords$lon, y = coords$lat)

  d <- merge(sl, xy, by = "site")
  missing <- setdiff(unique(sl$site), xy$site)
  if (length(missing))
    warning("no coordinates for: ", paste(missing, collapse = ", "))
  d$lab <- factor(sprintf("AD %d-%d", d$slice_start, d$slice_end),
                  levels = unique(sprintf("AD %d-%d", d$slice_start, d$slice_end)[
                    order(d$slice_start)]))

  ggplot(d, aes(x, y)) +
    geom_point(data = xy, inherit.aes = FALSE, aes(x, y),
               colour = "grey85", size = 0.9) +
    geom_point(aes(size = p_occupied, alpha = p_occupied),
               colour = OKABE_ITO[6]) +
    facet_wrap(~lab) +
    scale_size_continuous(range = c(0.4, 4), limits = c(0, 1),
                          name = "P(occupied)") +
    scale_alpha_continuous(range = c(0.15, 0.9), limits = c(0, 1),
                           guide = "none") +
    coord_equal() +
    labs(x = NULL, y = NULL) + theme_susq() +
    theme(axis.text = element_blank(), panel.grid = element_blank())
}

#' Plotting Functions for Radiocarbon Dating Analysis
#'
#' Visualization functions for calibrated dates, posterior distributions,
#' and model diagnostics.
#'
#' @author Project Team
#' @date 2025-10-11

library(ggplot2)
library(bayesplot)
library(dplyr)
library(tidyr)

#' Plot calibrated radiocarbon dates
#'
#' Creates density plots of calibrated radiocarbon dates by deposit.
#'
#' @param cal_result Calibrated dates from calibrate_dates()
#' @param c14_data Original radiocarbon data with deposit assignments
#' @param fill_by Character, color fill by "deposit" or NULL (default "deposit")
#' @return ggplot object
#' @export
plot_calibrated_dates <- function(cal_result, c14_data, fill_by = "deposit") {

  # Extract calibrated densities
  n_dates <- length(cal_result$calDates$grids)

  plot_data <- data.frame()

  for (i in 1:n_dates) {
    grid <- cal_result$calDates$grids[[i]]

    temp_df <- data.frame(
      cal_bp = grid$calBP,
      density = grid$PrDens,
      lab_code = cal_result$summary$lab_code[i],
      deposit = c14_data$deposit[i]
    )

    plot_data <- rbind(plot_data, temp_df)
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = cal_bp, y = density))

  if (!is.null(fill_by) && fill_by == "deposit") {
    p <- p + geom_area(aes(fill = deposit), alpha = 0.5, position = "identity")
  } else {
    p <- p + geom_area(alpha = 0.5)
  }

  p <- p +
    facet_wrap(~lab_code, ncol = 2, scales = "free_y") +
    scale_x_reverse() +  # BP convention
    labs(
      x = "Calendar age (cal BP)",
      y = "Probability density",
      title = "Calibrated Radiocarbon Dates"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8)
    )

  return(p)
}

#' Plot occupation boundaries
#'
#' Visualizes posterior distributions of occupation boundaries for each deposit.
#'
#' @param fit Fitted Stan model object
#' @param n_deposits Number of deposits
#' @param deposit_names Character vector of deposit names (optional)
#' @return ggplot object
#' @export
plot_occupation_boundaries <- function(fit, n_deposits, deposit_names = NULL) {

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  if (is.null(deposit_names)) {
    deposit_names <- paste("Deposit", 1:n_deposits)
  }

  # Prepare data for plotting
  plot_data <- data.frame()

  for (k in 1:n_deposits) {
    temp_df <- data.frame(
      value = c(theta_start[, k], theta_end[, k]),
      boundary = rep(c("Start", "End"), each = nrow(theta_start)),
      deposit = deposit_names[k]
    )
    plot_data <- rbind(plot_data, temp_df)
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = value, fill = boundary)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~deposit, ncol = 1, scales = "free_y") +
    scale_x_reverse() +  # BP convention
    labs(
      x = "Calendar age (cal BP)",
      y = "Posterior density",
      title = "Occupation Boundaries",
      fill = "Boundary"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  return(p)
}

#' Plot temporal overlap from metrics
#'
#' Helper function that takes pre-computed occupation spans
#'
#' @param occupation_spans Data frame with occupation span statistics
#' @param deposit_names Character vector of deposit names
#' @return ggplot object
#' @export
plot_temporal_overlap_from_metrics <- function(occupation_spans, deposit_names) {
  # Simple duration-based plot from metrics
  plot_data <- data.frame(
    deposit = deposit_names,
    start_median = occupation_spans$start_median,
    end_median = occupation_spans$end_median,
    duration = occupation_spans$duration_mean
  )

  p <- ggplot(plot_data, aes(y = deposit)) +
    geom_segment(aes(x = end_median, xend = start_median, yend = deposit),
                 linewidth = 4, color = "steelblue", alpha = 0.7) +
    geom_point(aes(x = start_median), size = 3, shape = 21, fill = "darkblue") +
    geom_point(aes(x = end_median), size = 3, shape = 21, fill = "darkred") +
    scale_x_reverse() +
    labs(
      x = "Calendar Years BP",
      y = "Site",
      title = "Estimated Occupation Spans"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 12)
    )

  return(p)
}

#' Plot temporal overlap
#'
#' Creates a visual representation of deposit occupation ranges with uncertainty.
#'
#' @param fit Fitted Stan model object
#' @param n_deposits Number of deposits
#' @param deposit_names Character vector of deposit names (optional)
#' @param probs Credible interval probabilities (default c(0.68, 0.95))
#' @return ggplot object
#' @export
plot_temporal_overlap <- function(fit, n_deposits, deposit_names = NULL,
                                   probs = c(0.68, 0.95)) {

  if (is.null(deposit_names)) {
    deposit_names <- paste("Deposit", 1:n_deposits)
  }

  # Extract boundary summaries
  occupation_spans <- calculate_occupation_spans(fit, n_deposits, probs = c(0.025, 0.16, 0.5, 0.84, 0.975))

  # Prepare data
  plot_data <- data.frame(
    deposit = deposit_names,
    start_median = occupation_spans$start_median,
    start_q16 = occupation_spans$start_q16,
    start_q84 = occupation_spans$start_q84,
    start_q2.5 = occupation_spans$start_q2.5,
    start_q97.5 = occupation_spans$start_q97.5,
    end_median = occupation_spans$end_median,
    end_q16 = occupation_spans$end_q16,
    end_q84 = occupation_spans$end_q84,
    end_q2.5 = occupation_spans$end_q2.5,
    end_q97.5 = occupation_spans$end_q97.5
  )

  # Create plot
  p <- ggplot(plot_data, aes(y = deposit)) +
    # 95% CI
    geom_segment(aes(x = end_q2.5, xend = start_q97.5, yend = deposit),
                 linewidth = 1, alpha = 0.3) +
    # 68% CI
    geom_segment(aes(x = end_q16, xend = start_q84, yend = deposit),
                 linewidth = 3, alpha = 0.6) +
    # Median
    geom_point(aes(x = start_median), size = 3, shape = 21, fill = "blue") +
    geom_point(aes(x = end_median), size = 3, shape = 21, fill = "red") +
    scale_x_reverse() +  # BP convention
    labs(
      x = "Calendar age (cal BP)",
      y = "Deposit",
      title = "Occupation Ranges with Uncertainty",
      subtitle = "Dark band: 68% CI, Light band: 95% CI\nBlue: Start, Red: End"
    ) +
    theme_minimal()

  return(p)
}

#' Plot overlap probability matrix
#'
#' Heatmap of pairwise overlap probabilities.
#'
#' @param overlap_matrix Matrix from create_overlap_matrix()
#' @param deposit_names Character vector of deposit names (optional)
#' @return ggplot object
#' @export
plot_overlap_matrix <- function(overlap_matrix, deposit_names = NULL) {

  n_deposits <- nrow(overlap_matrix)

  if (is.null(deposit_names)) {
    deposit_names <- paste("Deposit", 1:n_deposits)
  }

  # Convert to long format
  plot_data <- expand.grid(
    deposit_i = 1:n_deposits,
    deposit_j = 1:n_deposits
  )
  plot_data$probability <- as.vector(overlap_matrix)
  plot_data$deposit_i_name <- deposit_names[plot_data$deposit_i]
  plot_data$deposit_j_name <- deposit_names[plot_data$deposit_j]

  # Create heatmap
  p <- ggplot(plot_data, aes(x = deposit_i_name, y = deposit_j_name, fill = probability)) +
    geom_tile() +
    geom_text(aes(label = round(probability, 2)), color = "white", size = 4) +
    scale_fill_gradient(low = "darkred", high = "darkgreen", limits = c(0, 1)) +
    labs(
      x = "Deposit",
      y = "Deposit",
      title = "Pairwise Overlap Probabilities",
      fill = "P(Overlap)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_equal()

  return(p)
}

#' Plot posterior predictive check
#'
#' Compares observed data to posterior predictive distribution.
#'
#' @param fit Fitted Stan model object
#' @param c14_data Original radiocarbon data
#' @param n_samples Number of posterior samples to use (default 50)
#' @return ggplot object
#' @export
plot_posterior_predictive <- function(fit, c14_data, n_samples = 50) {

  # Extract posterior predictive samples
  c14_rep <- tryCatch({
    extract_posterior(fit, pars = "c14_age_rep", format = "matrix")
  }, error = function(e) {
    # If extraction fails, create a simple placeholder plot
    message("Could not extract posterior predictive samples")
    return(NULL)
  })

  if (is.null(c14_rep) || nrow(c14_rep) == 0) {
    # Create placeholder plot
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "Posterior Predictive Check\n(Data extraction in progress)",
               size = 8) +
      theme_void()
    return(p)
  }

  # Sample subset
  if (nrow(c14_rep) > n_samples) {
    sample_idx <- sample(1:nrow(c14_rep), n_samples)
    c14_rep <- c14_rep[sample_idx, ]
  }

  # Convert matrix to long format
  c14_rep_long <- as.data.frame(c14_rep)
  c14_rep_long$sample_id <- 1:nrow(c14_rep_long)
  c14_rep_long <- tidyr::pivot_longer(c14_rep_long,
                                       cols = -sample_id,
                                       names_to = "observation",
                                       values_to = "age")

  # Handle both c14_age and age column names
  age_col <- if("c14_age" %in% names(c14_data)) c14_data$c14_age else c14_data$age

  # Create density plot
  p <- ggplot() +
    geom_density(data = c14_rep_long, aes(x = age, group = sample_id),
                 alpha = 0.05, color = "blue", linewidth = 0.3) +
    geom_density(data = data.frame(age = age_col), aes(x = age),
                 linewidth = 1.2, color = "black") +
    labs(
      x = "Radiocarbon age (BP)",
      y = "Density",
      title = "Posterior Predictive Check",
      subtitle = "Black: Observed data, Blue: Posterior predictions"
    ) +
    theme_minimal()

  return(p)
}

#' Plot contemporaneity visualization
#'
#' Creates a timeline showing occupation spans with uncertainty bands,
#' colored by contemporaneity groups and showing overlap relationships.
#'
#' @param fit Fitted Stan model object
#' @param site_names Character vector of site names
#' @param overlap_matrix Matrix of pairwise overlap probabilities
#' @param overlap_threshold Threshold for defining contemporaneity (default 0.95)
#' @return ggplot object
#' @export
plot_contemporaneity_timeline <- function(fit, site_names, overlap_matrix,
                                          overlap_threshold = 0.95) {

  n_sites <- length(site_names)

  # Extract boundary samples
  theta_start <- extract_posterior(fit, pars = "theta_start", format = "matrix")
  theta_end <- extract_posterior(fit, pars = "theta_end", format = "matrix")

  # Compute summaries for each site
  timeline_data <- data.frame()

  for (i in 1:n_sites) {
    start_samples <- theta_start[, i]
    end_samples <- theta_end[, i]

    row <- data.frame(
      site = site_names[i],
      site_num = i,
      start_median = median(start_samples),
      start_q16 = quantile(start_samples, 0.16),
      start_q84 = quantile(start_samples, 0.84),
      start_q025 = quantile(start_samples, 0.025),
      start_q975 = quantile(start_samples, 0.975),
      end_median = median(end_samples),
      end_q16 = quantile(end_samples, 0.16),
      end_q84 = quantile(end_samples, 0.84),
      end_q025 = quantile(end_samples, 0.025),
      end_q975 = quantile(end_samples, 0.975)
    )

    timeline_data <- rbind(timeline_data, row)
  }

  # Identify contemporaneity groups using hierarchical clustering
  # Convert overlap matrix to distance matrix
  dist_matrix <- as.dist(1 - overlap_matrix)
  hc <- hclust(dist_matrix, method = "complete")

  # Cut tree to identify groups (sites with >threshold overlap)
  height_cutoff <- 1 - overlap_threshold
  groups <- cutree(hc, h = height_cutoff)
  timeline_data$group <- as.factor(groups)

  # Order sites by median start date (earliest at top)
  timeline_data <- timeline_data[order(-timeline_data$start_median), ]
  timeline_data$site_ordered <- factor(timeline_data$site,
                                       levels = timeline_data$site)

  # Create the plot
  p <- ggplot(timeline_data, aes(y = site_ordered, color = group, fill = group)) +
    # 95% CI bands
    geom_segment(aes(x = start_q975, xend = end_q025,
                     y = site_ordered, yend = site_ordered),
                 linewidth = 6, alpha = 0.2) +
    # 68% CI bands (1 sigma)
    geom_segment(aes(x = start_q84, xend = end_q16,
                     y = site_ordered, yend = site_ordered),
                 linewidth = 8, alpha = 0.4) +
    # Median span
    geom_segment(aes(x = start_median, xend = end_median,
                     y = site_ordered, yend = site_ordered),
                 linewidth = 2, alpha = 0.9) +
    # Median points for start/end
    geom_point(aes(x = start_median), size = 3, shape = 21,
               fill = "white", stroke = 1.5) +
    geom_point(aes(x = end_median), size = 3, shape = 21,
               fill = "white", stroke = 1.5) +
    scale_x_reverse(breaks = seq(0, 1200, 100)) +
    scale_color_brewer(palette = "Set1", name = "Contemporaneity\nGroup") +
    scale_fill_brewer(palette = "Set1", name = "Contemporaneity\nGroup") +
    labs(
      x = "Calendar Years BP",
      y = "",
      title = "Site Occupation Timelines with Contemporaneity Groups",
      subtitle = "Lines show median spans; bands show 68% (dark) and 95% (light) credible intervals"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      axis.text.y = element_text(size = 11, face = "bold")
    )

  return(p)
}

#' Plot convergence diagnostics
#'
#' Creates trace plots and rank plots for key parameters.
#'
#' @param fit Fitted Stan model object
#' @param pars Parameters to plot (default c("theta_start", "theta_end"))
#' @return List of ggplot objects
#' @export
plot_convergence_diagnostics <- function(fit, pars = c("theta_start", "theta_end")) {

  draws <- extract_posterior(fit, pars = pars, format = "draws")

  plots <- list()

  # Trace plots
  plots$trace <- mcmc_trace(draws, pars = pars) +
    labs(title = "Trace Plots") +
    theme_minimal()

  # Rank plots
  plots$rank <- mcmc_rank_overlay(draws, pars = pars) +
    labs(title = "Rank Plots") +
    theme_minimal()

  return(plots)
}

#' Plot model comparison
#'
#' Visualizes model comparison results from LOO.
#'
#' @param comparison_result Output from compare_models_loo()
#' @return ggplot object
#' @export
plot_model_comparison <- function(comparison_result) {

  comp_df <- as.data.frame(comparison_result$comparison)
  comp_df$model <- rownames(comp_df)

  weights_df <- data.frame(
    model = names(comparison_result$weights),
    weight = as.numeric(comparison_result$weights)
  )

  # Create barplot of model weights
  p <- ggplot(weights_df, aes(x = reorder(model, weight), y = weight)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = round(weight, 3)), vjust = -0.5) +
    ylim(0, 1) +
    labs(
      x = "Model",
      y = "LOO Weight",
      title = "Model Comparison via LOO",
      subtitle = "Higher weight indicates better model fit"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}

#' Plot validation results
#'
#' Visualizes simulation-based validation results.
#'
#' @param validation_results Data frame with validation results
#' @return ggplot object
#' @export
plot_validation_results <- function(validation_results) {

  # Assuming validation_results has columns:
  # scenario, sample_size, correct_model_selected, overlap_bias, etc.

  p1 <- ggplot(validation_results, aes(x = sample_size, y = correct_model_selected, color = scenario)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
    ylim(0, 1) +
    labs(
      x = "Sample size per deposit",
      y = "Proportion correct",
      title = "Model Selection Accuracy",
      subtitle = "Dashed line: 70% threshold"
    ) +
    theme_minimal()

  return(p1)
}

#' Create comprehensive diagnostic plot
#'
#' Combines multiple diagnostic plots into a single figure.
#'
#' @param fit Fitted Stan model object
#' @param c14_data Original radiocarbon data
#' @param n_deposits Number of deposits
#' @return List of ggplot objects
#' @export
create_diagnostic_plots <- function(fit, c14_data, n_deposits) {

  plots <- list()

  plots$boundaries <- plot_occupation_boundaries(fit, n_deposits)
  plots$overlap <- plot_temporal_overlap(fit, n_deposits)
  plots$ppc <- plot_posterior_predictive(fit, c14_data)

  message("Diagnostic plots created successfully")

  return(plots)
}

#' Save all plots to files
#'
#' Convenience function to save plots to output directory.
#'
#' @param plots List of ggplot objects
#' @param output_dir Directory to save plots (default "output/plots")
#' @param width Plot width in inches (default 10)
#' @param height Plot height in inches (default 8)
#' @export
save_plots <- function(plots, output_dir = "output/plots", width = 10, height = 8) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  for (plot_name in names(plots)) {
    filename <- file.path(output_dir, paste0(plot_name, ".png"))
    ggsave(filename, plots[[plot_name]], width = width, height = height, dpi = 300)
    message("Saved: ", filename)
  }
}

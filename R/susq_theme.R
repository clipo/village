#' Shared figure style. No titles; descriptive information belongs in captions.

suppressPackageStartupMessages({ library(ggplot2); library(systemfonts) })

OKABE_ITO <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

susq_font <- function() {
  if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "sans"
}

theme_susq <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = susq_font()) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      axis.line = element_line(linewidth = 0.35, colour = "grey20"),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey20"),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      legend.position = "bottom", legend.key.height = unit(0.7, "lines"),
      legend.title = element_text(size = base_size - 1),
      strip.text = element_text(face = "plain", hjust = 0)
    )
}

susq_save <- function(plot, file, width = 7, height = 4.5) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(file, plot, device = ragg::agg_png, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
  message("wrote ", file)
  invisible(file)
}

#' Intervals where IntCal20 reverses, measured directly from the curve over the
#' study range. Calibration is multimodal inside these and dating resolution
#' collapses, so they are marked on every figure with a calendar axis.
CURVE_REVERSALS <- data.frame(
  start_calBP = c(411, 886, 613, 1233, 1463, 1041),
  end_calBP   = c(344, 837, 577, 1199, 1429, 1014)
)

#' Sites whose single fitted interval is wide enough that it is describing a
#' palimpsest rather than one occupation. Flagged on every figure where their
#' width would otherwise be read as a village lifespan.
PALIMPSEST_SITES <- c("CCE_Site", "Chenango_Point", "Roundtop", "Broome_Tech")

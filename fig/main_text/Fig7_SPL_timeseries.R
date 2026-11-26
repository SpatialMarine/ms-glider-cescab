
# ============================================================
# Ocean gliders for PAM
# ============================================================
# Greta Jankauskaite | @gretajan97
# ============================================================
# Fig7_SPL_timeseries.R
# ============================================================

# Calculate daily SPL exceedance levels for the cleaned 63 Hz, 125 Hz, and 8 kHz data. 
# Save time-series panels.

source("setup.R")

# ------------------------------------------------------------
# 1. Load cleaned data and define output directory
# ------------------------------------------------------------
input_file <- file.path(
  output_dir, "noise", "timeseries",
  "CESCAB_SPL_63_125_8k.rds"
)

data <- readRDS(input_file)

outdir <- file.path(output_dir, "noise", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Calculate SPL exceedance levels
# ------------------------------------------------------------

# Ln is the SPL exceeded by n% of estimates within each time bin.
calculate_el <- function(data, band, period, min_samples) {
  if (!band %in% c("Hz63", "Hz125", "Hz8k")) {
    stop("Unknown band: ", band)
  }
  
  if (!period %in% c("hour", "day")) {
    stop("Unknown period: ", period)
  }
  
  group_time <- lubridate::floor_date(data$Time, unit = period)
  
  band_data <- data.frame(
    Time = group_time,
    SPL = data[[band]]
  )
  
  band_data |>
    dplyr::filter(
      !is.na(Time),
      is.finite(SPL)
    ) |>
    dplyr::group_by(Time) |>
    dplyr::filter(dplyr::n() >= min_samples) |>
    dplyr::summarise(
      n = dplyr::n(),
      L99 = as.numeric(stats::quantile(SPL, 0.01, type = 5)),
      L90 = as.numeric(stats::quantile(SPL, 0.10, type = 5)),
      L75 = as.numeric(stats::quantile(SPL, 0.25, type = 5)),
      L50 = as.numeric(stats::quantile(SPL, 0.50, type = 5)),
      L25 = as.numeric(stats::quantile(SPL, 0.75, type = 5)),
      L10 = as.numeric(stats::quantile(SPL, 0.90, type = 5)),
      L1 = as.numeric(stats::quantile(SPL, 0.99, type = 5)),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# 3. Define time-series plot
# ------------------------------------------------------------
plot_el <- function(stats_data, title) {
  plot_theme <- if (exists("theme_article", mode = "function")) {
    theme_article()
  } else {
    ggplot2::theme_bw()
  }
  
  ggplot2::ggplot(stats_data, ggplot2::aes(x = Time)) +
    # Pink shading: middle 80% of SPL estimates.
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = L90, ymax = L10),
      fill = "#E16462FF",
      alpha = 0.6
    ) +
    # Blue shading: interquartile range.
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = L75, ymax = L25),
      fill = "lightskyblue1",
      alpha = 0.9
    ) +
    # Dashed lines: upper and lower exceedance levels.
    ggplot2::geom_line(
      ggplot2::aes(y = L1),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = L99),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    # Dark blue line: median.
    ggplot2::geom_line(
      ggplot2::aes(y = L50),
      colour = "royalblue4",
      linewidth = 0.8
    ) +
    ggplot2::scale_x_datetime(
      date_breaks = "3 days",
      date_labels = "%m/%d",
      timezone = "UTC"
    ) +
    ggplot2::labs(
      title = title,
      x = "Date (UTC)",
      y = "SPL (dB re 1 µPa)"
    ) +
    plot_theme +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 15)
    )
}

# ------------------------------------------------------------
# 4. Define manuscript panels
# ------------------------------------------------------------
panels <- list(
  # Fig4b_hourly_63Hz = list(
  #   band = "Hz63", period = "hour", min_samples = 200L,
  #   title = "Hourly 63-Hz 1/3-octave band SPL"
  # ),
  # Fig4c_hourly_125Hz = list(
  #   band = "Hz125", period = "hour", min_samples = 200L,
  #   title = "Hourly 125-Hz 1/3-octave band SPL"
  # ),
  Fig4d_daily_63Hz = list(
    band = "Hz63", period = "day", min_samples = 800L,
    title = "Daily 63-Hz 1/3-octave band SPL"
  ),
  Fig4e_daily_125Hz = list(
    band = "Hz125", period = "day", min_samples = 800L,
    title = "Daily 125-Hz 1/3-octave band SPL"
  ),
  FigS4a_daily_8kHz = list(
    band = "Hz8k", period = "day", min_samples = 800L,
    title = "Daily 8-kHz 1/3-octave band SPL"
  )
)

# ------------------------------------------------------------
# 5. Save panels
# ------------------------------------------------------------
for (panel_name in names(panels)) {
  config <- panels[[panel_name]]
  values <- calculate_el(data, config$band, config$period, config$min_samples)
  if (!nrow(values)) stop("No valid time bins for: ", panel_name)
  
  plot <- plot_el(values, config$title)
  
  ggplot2::ggsave(
    filename = file.path(outdir, paste0(panel_name, ".png")),
    plot = plot, width = 8, height = 3, dpi = 300
  )
  
  message("Saved ", panel_name, " (", nrow(values), " time bins)")
}

# ------------------------------------------------------------
# 6. Save common legend
# ------------------------------------------------------------
grDevices::png(
  file.path(outdir, "SPL_legend.png"),
  width = 2,
  height = 1.5,
  units = "in",
  res = 300
)

grid::grid.newpage()

grid::grid.rect(
  gp = grid::gpar(fill = "white", col = "grey65")
)

# Median and dashed exceedance-level lines.
grid::grid.segments(
  x0 = 0.10,
  y0 = c(0.82, 0.65, 0.52),
  x1 = 0.32,
  y1 = c(0.82, 0.65, 0.52),
  gp = grid::gpar(
    col = c("royalblue4", "black", "black"),
    lty = c(1, 2, 2),
    lwd = c(1.8, 1.5, 1.5)
  )
)

# Shaded ranges with the same transparency as the panels.
grid::grid.rect(
  x = 0.21,
  y = 0.35,
  width = 0.22,
  height = 0.08,
  gp = grid::gpar(
    fill = grDevices::adjustcolor("#E16462FF", alpha.f = 0.6),
    col = NA
  )
)

# Pink background beneath the blue interquartile ribbon.
grid::grid.rect(
  x = 0.21,
  y = 0.21,
  width = 0.22,
  height = 0.08,
  gp = grid::gpar(
    fill = grDevices::adjustcolor("#E16462FF", alpha.f = 0.6),
    col = NA
  )
)

grid::grid.rect(
  x = 0.21,
  y = 0.21,
  width = 0.22,
  height = 0.08,
  gp = grid::gpar(
    fill = grDevices::adjustcolor("lightskyblue1", alpha.f = 0.9),
    col = NA
  )
)

grid::grid.text(
  label = c(
    "L50 (median)",
    "L1",
    "L99",
    "L90–L10",
    "L75–L25"
  ),
  x = 0.40,
  y = c(0.82, 0.65, 0.52, 0.35, 0.21),
  just = "left",
  gp = grid::gpar(cex = 0.9)
)

grDevices::dev.off()

# ------------------------------------------------------------
# 7. Calculate daily L50 ranges for Results
# ------------------------------------------------------------
daily_63 <- calculate_el(
  data, "Hz63", "day", min_samples = 800L
)

daily_125 <- calculate_el(
  data, "Hz125", "day", min_samples = 800L
)

daily_L50_ranges <- data.frame(
  Band = c("63 Hz", "125 Hz"),
  Minimum_L50 = c(
    min(daily_63$L50),
    min(daily_125$L50)
  ),
  Maximum_L50 = c(
    max(daily_63$L50),
    max(daily_125$L50)
  )
)

print(daily_L50_ranges, digits = 6)

# ------------------------------------------------------------
# 8. Calculate hourly L50 ranges for Results
# ------------------------------------------------------------
# hourly_63 <- calculate_el(
#   data, "Hz63", "hour", min_samples = 200L
# )
#
# hourly_125 <- calculate_el(
#   data, "Hz125", "hour", min_samples = 200L
# )
#
# hourly_L50_ranges <- data.frame(
#   Band = c("63 Hz", "125 Hz"),
#   Minimum_L50 = c(
#     min(hourly_63$L50),
#     min(hourly_125$L50)
#   ),
#   Maximum_L50 = c(
#     max(hourly_63$L50),
#     max(hourly_125$L50)
#   )
# )
#
# print(hourly_L50_ranges, digits = 6)
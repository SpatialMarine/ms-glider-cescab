# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# 02_plot_CESCAB_noise_timeseries_FINAL.R

# Calculate hourly and daily SPL exceedance levels from the cleaned
# 63 Hz, 125 Hz, and 8 kHz data. Save the corresponding manuscript
# time-series panels and the values plotted in each panel.
# ============================================================

source("setup.R")

# ------------------------------------------------------------
# 1. Load packages and cleaned data
# ------------------------------------------------------------

input_file <- file.path(output_dir, "noise", "timeseries", "CESCAB_SPL_63_125_8k.rds")
data <- readRDS(input_file)

outdir <- file.path(output_dir, "noise", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Calculate hourly and daily SPL exceedance levels (ELs)
# ------------------------------------------------------------
calculate_el <- function(data, band, period, min_samples) {
  if (!band %in% c("Hz63", "Hz125", "Hz8k")) stop("Unknown band: ", band)
  if (!period %in% c("hour", "day")) stop("Unknown period: ", period)

  group_time <- lubridate::floor_date(data$Time, unit = period)
  band_data <- data.frame(Time = group_time, SPL = data[[band]])

  # Exclude missing/non-finite values, require enough samples, then compute
  # percentiles for every retained bin. L50 is the median (50 EL).
  band_data |>
    dplyr::filter(is.finite(SPL)) |>
    dplyr::group_by(Time) |>
    dplyr::filter(dplyr::n() >= min_samples) |>
    dplyr::summarise(
      n = dplyr::n(),
      EL99 = as.numeric(stats::quantile(SPL, 0.01, type = 5)),
      EL90 = as.numeric(stats::quantile(SPL, 0.10, type = 5)),
      EL75 = as.numeric(stats::quantile(SPL, 0.25, type = 5)),
      L50  = as.numeric(stats::quantile(SPL, 0.50, type = 5)),
      EL25 = as.numeric(stats::quantile(SPL, 0.75, type = 5)),
      EL10 = as.numeric(stats::quantile(SPL, 0.90, type = 5)),
      EL1  = as.numeric(stats::quantile(SPL, 0.99, type = 5)),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# 3. Plot the same layers and colours as the submitted EL plots
# ------------------------------------------------------------
plot_el <- function(stats_data, title) {
  # The original plotting code uses theme_article() from the user's setup.
  # Use it when loaded; theme_bw() makes this script runnable on its own.
  plot_theme <- if (exists("theme_article", mode = "function")) {
    theme_article()
  } else {
    ggplot2::theme_bw()
  }

  ggplot2::ggplot(stats_data, ggplot2::aes(x = Time)) +
    # Pink ribbon: 10–90 exceedance levels (10th–90th percentiles).
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = EL90, ymax = EL10),
      fill = "#E16462FF", alpha = 0.6
    ) +
    # Light blue ribbon: 25–75 exceedance levels (interquartile range).
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = EL75, ymax = EL25),
      fill = "lightskyblue1", alpha = 0.9
    ) +
    # Dashed lines: the 1 EL (loud) and 99 EL (quiet).
    ggplot2::geom_line(ggplot2::aes(y = EL1),
                       linetype = "dashed", linewidth = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = EL99),
                       linetype = "dashed", linewidth = 0.6) +
    # Dark blue line: the median (50 EL).
    ggplot2::geom_line(ggplot2::aes(y = L50),
                       colour = "royalblue4", linewidth = 0.8) +
    ggplot2::scale_x_datetime(date_breaks = "3 days", date_labels = "%m/%d",
                               timezone = "UTC") +
    ggplot2::labs(title = title, x = "Date (UTC)", y = "SPL (dB re 1 µPa)") +
    plot_theme +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 15)
    )
}

# ------------------------------------------------------------
# 4. Define only the panels used in the manuscript
# ------------------------------------------------------------
panels <- list(
  Fig4b_hourly_63Hz = list(
    band = "Hz63", period = "hour", min_samples = 200L,
    title = "Hourly 63-Hz 1/3-octave band SPL"
  ),
  Fig4c_hourly_125Hz = list(
    band = "Hz125", period = "hour", min_samples = 200L,
    title = "Hourly 125-Hz 1/3-octave band SPL"
  ),
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
# 5. Save each plot alongside its exact plotted values
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
  utils::write.csv(values, file.path(outdir, paste0(panel_name, ".csv")),
                   row.names = FALSE)
  message("Saved ", panel_name, " (", nrow(values), " time bins)")
}


# ------------------------------------------------------------
# 5. Save a common legend for the plots
# ------------------------------------------------------------
png(file.path(outdir, "SPL_legend.png"),
    width = 2, height = 1.5, units = "in", res = 300)

grid::grid.newpage()
grid::grid.rect(gp = grid::gpar(fill = "white", col = "grey65"))
grid::grid.segments(.12, c(.82, .65, .52), .34, c(.82, .65, .52),
                    gp = grid::gpar(col = c("royalblue4", "black", "black"),
                                    lty = c(1, 2, 2), lwd = 1.5))
grid::grid.rect(x = .23, y = c(.35, .21), width = .22, height = .08,
                gp = grid::gpar(fill = c("#E16462FF", "lightskyblue1"), col = NA))
grid::grid.text(c("median", "1 EL", "99 EL", "10–90 EL", "25–75 EL"),
                x = .43, y = c(.82, .65, .52, .35, .21),
                just = "left", gp = grid::gpar(cex = .9))


dev.off()

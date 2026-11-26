
# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# Fig6a_track_SPLs.R
#
# Match per-file median SPLs in the 63 Hz and 125 Hz bands
# to the nearest glider positions within five minutes. Inspect
# SPL distributions and save maps showing noise levels along
# the glider track over bathymetry.
# ============================================================

source("setup.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
nc_file <- file.path(
  input_dir,
  "glider_data",
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)

noise_file <- file.path(
  output_dir, "noise", "timeseries",
  "CESCAB_SPL_63_125_8k.rds"
)

outdir <- file.path(output_dir, "noise", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Import track and calculate per-recording median SPLs
# ------------------------------------------------------------
track_nc <- ncdf4::nc_open(nc_file)

track <- tibble::tibble(
  lon = as.numeric(ncdf4::ncvar_get(track_nc, "longitude")),
  lat = as.numeric(ncdf4::ncvar_get(track_nc, "latitude")),
  time = as.numeric(ncdf4::ncvar_get(track_nc, "time"))
)

ncdf4::nc_close(track_nc)

track <- track %>%
  filter(
    is.finite(lon),
    is.finite(lat),
    is.finite(time),
    time > 0
  ) %>%
  mutate(
    time = as.POSIXct(
      time,
      origin = "1970-01-01",
      tz = "UTC"
    )
  ) %>%
  arrange(time)

# Load cleaned SPL samples.
noise <- readRDS(noise_file)

# Calculate one median per recording, separately for each band.
band_63Hz <- noise %>%
  filter(!is.na(FileStart), is.finite(Hz63)) %>%
  group_by(FileStart) %>%
  summarise(
    Hz63_median = median(Hz63),
    .groups = "drop"
  ) %>%
  mutate(ValueTime = FileStart)

band_125Hz <- noise %>%
  filter(!is.na(FileStart), is.finite(Hz125)) %>%
  group_by(FileStart) %>%
  summarise(
    Hz125_median = median(Hz125),
    .groups = "drop"
  ) %>%
  mutate(ValueTime = FileStart)

# ------------------------------------------------------------
# 3. Match SPLs to the nearest glider positions
# ------------------------------------------------------------
track <- data.table::as.data.table(track)

track <- track[
  is.finite(lon) & is.finite(lat) & !is.na(time)
]

data.table::setorder(track, time)

# Preserve the actual track timestamp separately from the join key.
track[, track_time := time]
data.table::setkey(track, time)

tol_sec <- 5 * 60

match_noise_to_track <- function(band_data, spl_column) {
  
  band <- data.table::as.data.table(band_data)
  
  band <- band[
    !is.na(ValueTime) & is.finite(get(spl_column))
  ]
  
  band[, recording_time := ValueTime]
  
  matched <- track[
    band,
    on = .(time = ValueTime),
    roll = "nearest"
  ]
  
  matched[, time_diff_sec := abs(
    as.numeric(difftime(
      track_time,
      recording_time,
      units = "secs"
    ))
  )]
  
  matched <- matched[
    !is.na(time_diff_sec) & time_diff_sec <= tol_sec
  ]
  
  result <- matched[, .(
    lon,
    lat,
    closest_data_time = track_time,
    MedianTime = recording_time,
    time_diff_sec,
    SPL_median = get(spl_column)
  )]
  
  data.table::setnames(result, "SPL_median", spl_column)
  
  message(
    spl_column, ": ",
    nrow(result), " of ", nrow(band),
    " valid recordings matched within five minutes."
  )
  
  result
}

noise_63 <- match_noise_to_track(
  band_63Hz,
  "Hz63_median"
)

noise_125 <- match_noise_to_track(
  band_125Hz,
  "Hz125_median"
)

# ------------------------------------------------------------
# 4. Inspect SPL distributions
# ------------------------------------------------------------
# Values above the 99th percentile are omitted only from the
# trimmed histogram. All matched values remain available for maps.

inspect_spl <- function(values, band_label) {
  
  values <- round(values[is.finite(values)], 1)
  
  if (length(values) == 0L) {
    stop("No valid SPL values available for ", band_label, ".")
  }
  
  hist(
    values,
    breaks = "FD",
    col = "grey90",
    border = "grey40",
    main = paste(band_label, "per-recording medians"),
    xlab = "SPL (dB re 1 µPa)"
  )
  
  q99 <- unname(quantile(values, 0.99))
  trimmed <- values[values <= q99]
  
  n_trimmed <- sum(values > q99)
  pct_trimmed <- 100 * n_trimmed / length(values)
  
  hist(
    trimmed,
    breaks = "FD",
    col = "grey90",
    border = "grey40",
    main = paste(
      band_label,
      "per-recording medians (≤ 99th percentile)"
    ),
    xlab = "SPL (dB re 1 µPa)"
  )
  
  abline(v = q99, col = "red", lty = 2, lwd = 2)
  abline(v = median(trimmed), lty = 2, lwd = 2)
  abline(v = mean(trimmed), lty = 3, lwd = 2)
  
  legend(
    "topleft",
    legend = c(
      sprintf(
        "q99 = %.2f dB (trimmed %d, %.1f%%)",
        q99, n_trimmed, pct_trimmed
      ),
      "Median (dashed)",
      "Mean (dotted)"
    ),
    bty = "n"
  )
  
  invisible(range(trimmed))
}

db_limits_63 <- inspect_spl(
  noise_63$Hz63_median,
  "63 Hz"
)

db_limits_125 <- inspect_spl(
  noise_125$Hz125_median,
  "125 Hz"
)

# ------------------------------------------------------------
# 5. Import spatial layers
# ------------------------------------------------------------

# Bathymetry
bathy <- terra::rast(
  file.path(
    carto_dir,
    "EMODnet-bathymetry_2024",
    "bathymetry_dtm_2024_bb7b_a70e_b1f5.nc"
  )
)

bathy_crop <- terra::crop(bathy, study_area)

bathy_df <- as.data.frame(
  bathy_crop,
  xy = TRUE,
  na.rm = TRUE
)

names(bathy_df) <- c("x", "y", "value")

bathy_df <- bathy_df %>%
  filter(value < 1)

# Bathymetry contours
c_file <- paste0(
  carto_dir,
  "/EMODnet_Bathymetry_2022_contours/Contours_2022.shp"
)

contours <- st_read(c_file)

contour_1000 <- subset(
  contours,
  Elevation == "1000.00"
)

cont_int_1000 <- st_intersection(
  contour_1000,
  study_area
)

# Spain boundaries
spain <- sf::st_read(
  file.path(
    carto_dir,
    "spain_boundaries_GISData_MAPOG",
    "spain_Spain_Country_Boundary.shp"
  ),
  quiet = TRUE
) %>%
  sf::st_transform(4326)

# Cabrera Archipelago National Park
cabrera_MPA_pr <- sf::st_read(
  file.path(
    carto_dir,
    "MPA",
    "Limite_Archipielago_Cabrera.shp"
  ),
  quiet = TRUE
) %>%
  sf::st_transform(4326)

# ------------------------------------------------------------
# 6. Create SPL maps
# ------------------------------------------------------------
# bathy_df, cabrera_pr and spain are supplied by setup.R.
# Values outside the colour limits are shown at the nearest
# scale endpoint; recordings are not removed from the maps.

make_spl_map <- function(
    noise_data,
    spl_column,
    band_label,
    colour_limits,
    colour_breaks
) {
  
  ggplot() +
    geom_raster(
      data = bathy_df,
      aes(x = x, y = y, fill = value)
    ) +
    scale_fill_gradient(
      low = "royalblue4",
      high = "lightskyblue1",
      na.value = "white",
      breaks = c(-2800, 0),
      labels = c("-2800", "0"),
      name = "Depth (m)",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    geom_sf(
      data = cabrera_MPA_pr,
      color = NA,
      fill = "darkblue",
      alpha = 0.2
    ) +
    geom_sf(
      data = spain,
      fill = "lightgrey",
      colour = "black",
      linewidth = 0.1
    ) +
    geom_sf(
      data = cont_int_1000,
      colour = "grey1",
      linewidth = 0.07
    ) +
    geom_path(
      data = track,
      aes(x = lon, y = lat),
      colour = "black",
      linewidth = 0.8,
      lineend = "round"
    ) +
    geom_point(
      data = noise_data,
      aes(
        x = lon,
        y = lat,
        colour = .data[[spl_column]]
      ),
      size = 2.2,
      stroke = 0
    ) +
    scale_color_viridis_c(
      option = "plasma",
      name = paste0(
        "Median SPL (", band_label, ")\n",
        "(dB re 1 µPa)"
      ),
      limits = colour_limits,
      breaks = colour_breaks,
      oob = scales::squish,
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5
      )
    ) +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(
        size = 15,
        margin = margin(t = 7)
      ),
      axis.title.y = element_text(
        size = 15,
        margin = margin(r = 7)
      ),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.4
      ),
      axis.ticks.length = unit(0.08, "cm"),
      legend.direction = "horizontal",
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      )
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        title.vjust = 0.5,
        direction = "horizontal"
      )
    ) +
    coord_sf(
      xlim = c(2.12, 3.47),
      ylim = c(38.22, 39.35)
    )
}

spl_63 <- noise_63$Hz63_median
spl_63 <- spl_63[is.finite(spl_63)]

cap_99 <- unname(quantile(spl_63, probs = 0.99))

c(
  percentile_99_dB = cap_99,
  n_above_cap = sum(spl_63 > cap_99),
  percent_above_cap = 100 * mean(spl_63 > cap_99)
) # 87.741534 99th percentile

map_noise_63 <- make_spl_map(
  noise_data = noise_63,
  spl_column = "Hz63_median",
  band_label = "63 Hz",
  colour_limits = c(75, 88),
  colour_breaks = c(75, 82, 88)
)



spl_125 <- noise_125$Hz125_median
spl_125 <- spl_125[is.finite(spl_125)]

cap_99 <- unname(quantile(spl_125, probs = 0.99))

c(
  percentile_99_dB = cap_99,
  n_above_cap = sum(spl_125 > cap_99),
  percent_above_cap = 100 * mean(spl_125 > cap_99)
) # 88.340625

map_noise_125 <- make_spl_map(
  noise_data = noise_125,
  spl_column = "Hz125_median",
  band_label = "125 Hz",
  colour_limits = c(75, 88),
  colour_breaks = c(75, 82, 88)
)

# print(map_noise_63)
# print(map_noise_125)

# ------------------------------------------------------------
# 7. Save the maps
# ------------------------------------------------------------
ggsave(
  filename = file.path(outdir, "track_63Hz_median_SPL.png"),
  plot = map_noise_63,
  width = 9,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(outdir, "track_125Hz_median_SPL.png"),
  plot = map_noise_125,
  width = 9,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ------------------------------------------------------------
# 8. Get stats of per-file median SPL ranges per band
# ------------------------------------------------------------
tibble::tibble(
  band = c("63 Hz", "125 Hz"),
  lowest_dB = c(min(noise_63$Hz63_median), min(noise_125$Hz125_median)),
  highest_dB = c(max(noise_63$Hz63_median), max(noise_125$Hz125_median))
)

# band   lowest_dB highest_dB
#  63 Hz       74.8       98.3
# 125 Hz      76.6      103. 

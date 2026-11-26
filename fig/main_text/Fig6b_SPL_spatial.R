# ============================================================
# Ocean gliders for PAM
# ============================================================
# Greta Jankauskaite | @gretajan97
# ============================================================
# Fig6b_SPL_spatial.R
#
# Calculate per-file median SPLs in the 63 Hz and 125 Hz bands.
# Match recordings to glider positions and CTD depths within
# five minutes and extract bathymetry at the matched positions.
# Save square violin plots across bathymetric, glider-depth,
# and MPA-status intervals.
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

bathy_file <- file.path(
  carto_dir,
  "EMODnet-bathymetry_2024",
  "bathymetry_dtm_2024_bb7b_a70e_b1f5.nc"
)

outdir <- file.path(output_dir, "noise", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Import glider positions and CTD depth
# ------------------------------------------------------------

nc <- ncdf4::nc_open(nc_file)

track <- tibble::tibble(
  lon = as.numeric(ncdf4::ncvar_get(nc, "longitude")),
  lat = as.numeric(ncdf4::ncvar_get(nc, "latitude")),
  time = as.numeric(ncdf4::ncvar_get(nc, "time")),
  depth = as.numeric(ncdf4::ncvar_get(nc, "depth_ctd"))
)

ncdf4::nc_close(nc)

track <- track %>%
  filter(is.finite(time), time > 0) %>%
  mutate(
    time = as.POSIXct(
      time,
      origin = "1970-01-01",
      tz = "UTC"
    )
  ) %>%
  arrange(time)

# Prepare positions and depths separately because valid depth
# measurements are not available at every position timestamp.
positions <- track %>%
  filter(is.finite(lon), is.finite(lat)) %>%
  select(time, lon, lat) %>%
  distinct(time, .keep_all = TRUE)

depths <- track %>%
  filter(is.finite(depth), depth >= 0) %>%
  select(time, depth) %>%
  distinct(time, .keep_all = TRUE)


# ------------------------------------------------------------
# 3. Calculate per-file median SPLs
# ------------------------------------------------------------

noise <- readRDS(noise_file)

file_spl <- noise %>%
  select(FileStart, Hz63, Hz125) %>%
  pivot_longer(
    cols = c(Hz63, Hz125),
    names_to = "band",
    values_to = "SPL"
  ) %>%
  filter(!is.na(FileStart), is.finite(SPL)) %>%
  group_by(FileStart, band) %>%
  summarise(
    SPL_median = median(SPL),
    .groups = "drop"
  ) %>%
  arrange(FileStart, band)


# ------------------------------------------------------------
# 4. Match recordings to positions and depths
# ------------------------------------------------------------

tol_sec <- 5 * 60

match_covariates <- function(recordings, covariates) {
  
  covariates <- data.table::as.data.table(covariates)
  recordings <- data.table::as.data.table(recordings)
  
  covariates[, matched_time := time]
  recordings[, recording_time := FileStart]
  
  data.table::setkey(covariates, time)
  
  matched <- covariates[
    recordings,
    on = .(time = FileStart),
    roll = "nearest"
  ]
  
  matched[, time_diff_sec := abs(
    as.numeric(difftime(
      matched_time,
      recording_time,
      units = "secs"
    ))
  )]
  
  matched <- matched[
    is.finite(time_diff_sec) & time_diff_sec <= tol_sec
  ]
  
  as.data.frame(matched)
}

position_matched <- match_covariates(
  file_spl,
  positions
) %>%
  transmute(
    FileStart = recording_time,
    band,
    SPL_median,
    lon,
    lat
  )

depth_matched <- match_covariates(
  file_spl,
  depths
) %>%
  transmute(
    FileStart = recording_time,
    band,
    SPL_median,
    depth
  )


# ------------------------------------------------------------
# 5. Extract bathymetry at matched positions
# ------------------------------------------------------------

bathy <- terra::rast(bathy_file)

position_points <- terra::vect(
  position_matched,
  geom = c("lon", "lat"),
  crs = "EPSG:4326"
)

position_points <- terra::project(
  position_points,
  terra::crs(bathy)
)

bathy_values <- terra::extract(
  bathy,
  position_points,
  method = "bilinear"
)

bathy_matched <- position_matched %>%
  mutate(bathy = bathy_values[[2]]) %>%
  filter(is.finite(bathy), bathy < 0)


# ------------------------------------------------------------
# 6. Define glider-depth intervals
# ------------------------------------------------------------

# Retain the original depth intervals, including the upper endpoint.
depth_breaks <- c(seq(20, 920, by = 100), 940)

depth_labels <- paste0(
  head(depth_breaks, -1),
  "–",
  tail(depth_breaks, -1)
)

df_depth <- depth_matched %>%
  mutate(
    bin = cut(
      depth,
      breaks = depth_breaks,
      labels = depth_labels,
      include.lowest = TRUE,
      right = TRUE
    )
  ) %>%
  filter(!is.na(bin))


# ------------------------------------------------------------
# 7. Define bathymetric intervals
# ------------------------------------------------------------

# Group depths shallower than 100 m into the shallowest interval.
bathy_matched <- bathy_matched %>%
  mutate(bathy_for_bin = pmin(bathy, -100))

deep <- floor(
  min(bathy_matched$bathy_for_bin) / 100
) * 100

shallow <- -100
span <- shallow - deep

candidate_steps <- c(100, 200, 250, 300, 400, 500, 600)

valid_steps <- candidate_steps[
  ceiling(span / candidate_steps) <= 10
]

bin_step <- if (length(valid_steps)) {
  min(valid_steps)
} else {
  max(candidate_steps)
}

bathy_breaks <- seq(deep, shallow, by = bin_step)

if (tail(bathy_breaks, 1) < shallow) {
  bathy_breaks <- c(bathy_breaks, shallow)
}

bathy_labels <- paste0(
  tail(bathy_breaks, -1),
  " – ",
  head(bathy_breaks, -1)
)

df_bathy <- bathy_matched %>%
  mutate(
    bin = cut(
      bathy_for_bin,
      breaks = bathy_breaks,
      labels = bathy_labels,
      include.lowest = TRUE,
      right = TRUE
    ),
    bin = factor(
      bin,
      levels = rev(bathy_labels)
    )
  ) %>%
  filter(!is.na(bin))


# ------------------------------------------------------------
# 8. Define the violin-plot function
# ------------------------------------------------------------

make_violin <- function(plot_data, x_label, band_label) {
  
  ggplot(
    plot_data,
    aes(x = bin, y = SPL_median)
  ) +
    geom_violin(
      fill = "lightskyblue1",
      alpha = 0.9,
      colour = "black",
      linewidth = 0.4,
      scale = "width",
      trim = TRUE
    ) +
    geom_boxplot(
      width = 0.15,
      fill = "royalblue4",
      alpha = 0.7,
      colour = "black",
      linewidth = 0.4,
      outlier.shape = NA
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 2.2,
      colour = "black"
    ) +
    labs(
      x = x_label,
      y = paste0(
        band_label,
        " 1/3-octave band SPL (dB re 1 µPa)"
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 13
      ),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(
        size = 15,
        margin = margin(t = 16)
      ),
      axis.title.y = element_text(
        size = 15,
        margin = margin(r = 16)
      )
    )
}


# ------------------------------------------------------------
# 9. Assign MPA status to matched recording positions
# ------------------------------------------------------------

mpa <- sf::st_read(
  file.path(carto_dir, "MPA", "Limite_Archipielago_Cabrera.shp"),
  quiet = TRUE
) %>%
  sf::st_make_valid()

recording_points <- sf::st_as_sf(
  position_matched,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
) %>%
  sf::st_transform(sf::st_crs(mpa))

df_mpa <- position_matched %>%
  mutate(
    bin = factor(
      ifelse(
        lengths(sf::st_intersects(recording_points, mpa)) > 0,
        "Inside",
        "Outside"
      ),
      levels = c("Outside", "Inside")
    )
  )


# ------------------------------------------------------------
# 10. Export square plots for each frequency band
# ------------------------------------------------------------

bands <- c(
  Hz63 = "63 Hz",
  Hz125 = "125 Hz"
)

for (band_name in names(bands)) {
  
  band_label <- unname(bands[band_name])
  
  datasets <- list(
    bathy = filter(df_bathy, band == band_name),
    glider_depth = filter(df_depth, band == band_name),
    MPA = filter(df_mpa, band == band_name)
  )
  
  x_labels <- c(
    bathy = "Bathymetry (m)",
    glider_depth = "Glider depth (m)",
    MPA = "MPA status"
  )
  
  if (!nrow(datasets$glider_depth) || !nrow(datasets$bathy)) {
    stop("No valid plotting data for ", band_label, ".")
  }
  
  # Same SPL scale for all three plots within each frequency band.
  spl_range <- range(
    unlist(lapply(datasets, function(d) d$SPL_median)),
    na.rm = TRUE
  )
  
  y_limits <- c(
    floor(spl_range[1] / 5) * 5,
    ceiling(spl_range[2] / 5) * 5
  )
  
  for (v in names(datasets)) {
    
    p <- make_violin(
      datasets[[v]],
      x_label = unname(x_labels[v]),
      band_label = band_label
    ) +
      scale_y_continuous(
        breaks = seq(y_limits[1], y_limits[2], by = 5),
        expand = expansion(mult = c(0, 0.02))
      ) +
      coord_cartesian(ylim = y_limits) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 11
        ),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(
          size = 12,
          margin = margin(t = 10)
        ),
        axis.title.y = element_text(
          size = 12,
          margin = margin(r = 10)
        ),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    # Invisible second-line labels reserve similar space below this panel.
    if (v == "MPA") {
      p <- p +
        scale_x_discrete(
          labels = c(
            Outside = "Outside\n ",
            Inside = "Inside\n "
          )
        ) +
        theme(
          axis.text.x = element_text(
            angle = 0,
            hjust = 0.5,
            size = 11
          )
        )
    }
    
    print(p)
    
    ggsave(
      filename = file.path(
        outdir,
        paste0(
          gsub(" ", "", band_label),
          "-violin-", v, "-square.png"
        )
      ),
      plot = p,
      width = 6,
      height = 6,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }
  
  message(
    band_label, ": ",
    nrow(datasets$glider_depth),
    " recordings in glider-depth plots; ",
    nrow(datasets$bathy),
    " recordings in bathymetry plots."
  )
}


# ------------------------------------------------------------
# 11. Calculate and print summary statistics
# ------------------------------------------------------------

panel_medians <- dplyr::bind_rows(
  list(
    bathymetry = df_bathy,
    glider_depth = df_depth,
    MPA = df_mpa
  ),
  .id = "panel"
) |>
  dplyr::group_by(panel, band, bin) |>
  dplyr::summarise(
    n_recordings = dplyr::n(),
    median_SPL = median(SPL_median),
    .groups = "drop"
  )

print(panel_medians, n = Inf)

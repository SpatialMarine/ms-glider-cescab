
# ============================================================
# Ocean gliders for PAM
# ============================================================
# Greta Jankauskaite | @gretajan97
# ============================================================
# FigS2_vessel_density_glider_track.R
# ============================================================

# Extract monthly vessel-density values along the CESCAB glider track.
# August positions are matched to the August vessel-density raster,
# and September positions are matched to the September raster.

source("setup.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
nc_dir <- file.path(input_dir, "glider_data")
traffic_dir <- file.path(input_dir, "marine_traffic")

outdir <- file.path(output_dir, "noise", "fig")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

nc_file <- file.path(
  nc_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)

traffic_files <- list(
  "08" = file.path(
    traffic_dir,
    "GMTDS_2024_08-2024_09_All",
    "grid_float_All_2024_08_converted.tif"
  ),
  "09" = file.path(
    traffic_dir,
    "GMTDS_2024_08-2024_09_All",
    "grid_float_All_2024_09_converted.tif"
  )
)

# ------------------------------------------------------------
# 2. Import glider track
# ------------------------------------------------------------
track_nc <- nc_open(nc_file)

track <- tibble(
  longitude = as.numeric(ncvar_get(track_nc, "longitude")),
  latitude  = as.numeric(ncvar_get(track_nc, "latitude")),
  time      = as.POSIXct(
    ncvar_get(track_nc, "time"),
    origin = "1970-01-01",
    tz = "UTC"
  )
)

nc_close(track_nc)

track <- track %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    !is.na(time)
  ) %>%
  mutate(
    month = format(time, "%m")
  ) %>%
  arrange(time)

# ------------------------------------------------------------
# 3. Function to extract vessel density from one monthly raster
# ------------------------------------------------------------
extract_density <- function(track_df, raster_file) {
  
  traffic_raster <- rast(raster_file)
  
  # Global Marine Traffic rasters are already in WGS84 after conversion.
  # If not, uncomment the project line below.
  # traffic_raster <- terra::project(traffic_raster, "EPSG:4326")
  
  traffic_crop <- terra::crop(
    traffic_raster,
    study_area,
    snap = "out"
  )
  
  track_sf <- st_as_sf(
    track_df,
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )
  
  track_sf <- st_transform(track_sf, crs(traffic_crop))
  
  extracted <- terra::extract(
    traffic_crop,
    terra::vect(track_sf)
  )
  
  track_df %>%
    mutate(
      vessel_density = extracted[[2]]
    )
}

# ------------------------------------------------------------
# 4. Extract monthly vessel density along glider path
# ------------------------------------------------------------
track_density_df <- purrr::map_dfr(
  names(traffic_files),
  function(month_id) {
    
    track_month <- track %>%
      filter(month == month_id)
    
    if (nrow(track_month) == 0) {
      return(tibble())
    }
    
    extract_density(
      track_df = track_month,
      raster_file = traffic_files[[month_id]]
    )
  }
) %>%
  filter(!is.na(vessel_density)) %>%
  arrange(time)

max_vessel_density <- max(track_density_df$vessel_density, na.rm = TRUE)
print(max_vessel_density)
# Max vessel density along the glider track was 3.3 h km⁻² month⁻¹

# ------------------------------------------------------------
# 5. Plot vessel density along glider path through time
# ------------------------------------------------------------
traffic_plot <- ggplot(
  track_density_df,
  aes(x = time, y = vessel_density)
) +
  geom_line(
    linewidth = 0.7,
    colour = "royalblue4"
  ) +
  labs(
    x = "Date",
    y = expression("Vessel density (hr/km"^{2}*")")
  ) +
  scale_x_datetime(
    date_breaks = "3 days",
    date_labels = "%m/%d"
  ) +
  theme_minimal(base_size = 11) +
  theme_article() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 13.5)
  )

#traffic_plot

# ------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------
ggsave(
  filename = file.path(outdir, "CESCAB_vessel_density_along_glider_track.png"),
  plot = traffic_plot,
  width = 8,
  height = 3,
  dpi = 300
)



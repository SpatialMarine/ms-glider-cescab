
# ============================================================
# Ocean gliders for PAM
# ============================================================
# Greta Jankauskaite | @gretajan97
# ============================================================
# FigS3_vessel_density_map.R
# ============================================================

# Plot August 2024 vessel density around the CESCAB glider track.
# Calculate the study-area 99th-percentile threshold and display cells exceeding this threshold in grey. 
# Calculate maximum vessel densities before and after exclusion across the study area and Cabrera MPA.

source("setup.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
nc_dir <- file.path(input_dir, "glider_data")
traffic_dir <- file.path(input_dir, "marine_traffic")

outdir <- file.path(output_dir, "noise", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

nc_file <- file.path(
  nc_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)

traffic_file <- file.path(
  traffic_dir,
  "GMTDS_2024_08-2024_09_All",
  "grid_float_All_2024_08_converted.tif"
)

mpa_file <- file.path(
  carto_dir,
  "MPA",
  "Limite_Archipielago_Cabrera.shp"
)

spain_file <- file.path(
  carto_dir,
  "spain_boundaries_GISData_MAPOG",
  "spain_Spain_Country_Boundary.shp"
)

# ------------------------------------------------------------
# 2. Import glider track
# ------------------------------------------------------------
track_nc <- ncdf4::nc_open(nc_file)

track <- tibble::tibble(
  longitude = as.numeric(ncdf4::ncvar_get(track_nc, "longitude")),
  latitude = as.numeric(ncdf4::ncvar_get(track_nc, "latitude")),
  time = as.POSIXct(
    as.numeric(ncdf4::ncvar_get(track_nc, "time")),
    origin = "1970-01-01",
    tz = "UTC"
  )
)

ncdf4::nc_close(track_nc)

track <- track %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    !is.na(time)
  ) %>%
  arrange(time)

# ------------------------------------------------------------
# 3. Import vessel-density raster
# ------------------------------------------------------------
traffic_aug <- terra::rast(traffic_file)

traffic_crop <- terra::crop(
  traffic_aug,
  study_area,
  snap = "out"
)

# Treat source NA cells as zero recorded vessel density.
traffic_crop[is.na(traffic_crop)] <- 0

# ------------------------------------------------------------
# 4. Import land and protected-area boundaries
# ------------------------------------------------------------
spain <- sf::st_read(spain_file, quiet = TRUE)
cabrera_MPA <- sf::st_read(mpa_file, quiet = TRUE)

spain_plot <- sf::st_transform(spain, 4326)
cabrera_plot <- sf::st_transform(cabrera_MPA, 4326)

cabrera_vect <- cabrera_MPA %>%
  sf::st_transform(terra::crs(traffic_crop)) %>%
  terra::vect()

# ------------------------------------------------------------
# 5. Calculate the 99th-percentile threshold
# ------------------------------------------------------------
v <- terra::values(traffic_crop, mat = FALSE)
v <- v[is.finite(v)]

upper_threshold <- as.numeric(
  quantile(v, probs = 0.99)
)

exclusion_summary <- data.frame(
  percentile = 99,
  threshold = upper_threshold,
  cells_excluded = sum(v > upper_threshold),
  percent_excluded = 100 * mean(v > upper_threshold)
)

print(exclusion_summary)

# Exclude cells exceeding the threshold.
traffic_clean <- terra::ifel(
  traffic_crop > upper_threshold,
  NA,
  traffic_crop
)

# ------------------------------------------------------------
# 6. Calculate original and filtered maximum densities
# ------------------------------------------------------------
traffic_cabrera <- traffic_crop %>%
  terra::crop(cabrera_vect) %>%
  terra::mask(cabrera_vect)

traffic_cabrera_clean <- traffic_clean %>%
  terra::crop(cabrera_vect) %>%
  terra::mask(cabrera_vect)

maximum_density <- function(raster) {
  terra::global(raster, "max", na.rm = TRUE)[1, 1]
}

maximum_summary <- data.frame(
  area = c("Study area", "Cabrera MPA"),
  original_max = c(
    maximum_density(traffic_crop),
    maximum_density(traffic_cabrera)
  ),
  filtered_max = c(
    maximum_density(traffic_clean),
    maximum_density(traffic_cabrera_clean)
  )
)

print(maximum_summary)
# Study area: 22.2 h km⁻² month⁻¹
# Cabrera MPA 20.8 h km⁻² month⁻¹

# ------------------------------------------------------------
# 7. Prepare raster cells for plotting
# ------------------------------------------------------------
traffic_plot <- traffic_crop

if (!terra::same.crs(traffic_plot, "EPSG:4326")) {
  traffic_plot <- terra::project(
    traffic_plot,
    "EPSG:4326",
    method = "near"
  )
}

traffic_df <- as.data.frame(
  traffic_plot,
  xy = TRUE,
  na.rm = TRUE
)

names(traffic_df)[3] <- "vessel_density"

retained_df <- traffic_df %>%
  filter(vessel_density <= upper_threshold)

excluded_df <- traffic_df %>%
  filter(vessel_density > upper_threshold)

# ------------------------------------------------------------
# 8. Plot vessel density and glider track
# ------------------------------------------------------------
density_breaks <- seq(0, 20, by = 5)
density_labels <- as.character(density_breaks)

map_xlim <- range(track$longitude, na.rm = TRUE) + c(-0.05, 0.05)
map_ylim <- range(track$latitude, na.rm = TRUE) + c(-0.05, 0.05)

map_traffic <- ggplot() +
  geom_raster(
    data = retained_df,
    aes(x = x, y = y, fill = vessel_density)
  ) +
  scale_fill_viridis_c(
    name = expression(
      "Vessel density (h km"^{-2}*" month"^{-1}*")"
    ),
    option = "viridis",
    begin = 0.3,
    end = 1,
    limits = c(0, upper_threshold),
    breaks = density_breaks,
    labels = density_labels,
    na.value = "white"
  ) +
  geom_raster(
    data = excluded_df,
    aes(x = x, y = y),
    fill = "grey90",
    show.legend = FALSE
  ) +
  geom_sf(
    data = cabrera_plot,
    colour = "turquoise3",
    linewidth = 0.3,
    fill = NA
  ) +
  geom_sf(
    data = spain_plot,
    fill = "grey80",
    colour = "black",
    linewidth = 0.1
  ) +
  geom_path(
    data = track,
    aes(x = longitude, y = latitude),
    colour = "black",
    linewidth = 0.6,
    lineend = "round"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.2
  ) +
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering()
  ) +
  coord_sf(
    crs = sf::st_crs(4326),
    default_crs = sf::st_crs(4326),
    xlim = map_xlim,
    ylim = map_ylim,
    expand = FALSE
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
    axis.text = element_text(size = 14),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.4
    ),
    axis.ticks.length = grid::unit(0.08, "cm"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  ) +
  guides(
    fill = guide_colourbar(
      order = 1,
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(8, "cm"),
      barheight = grid::unit(0.4, "cm")
    )
  )

#print(map_traffic)

# ------------------------------------------------------------
# 9. Save map and summary tables
# ------------------------------------------------------------
ggsave(
  filename = file.path(
    outdir,
    "CESCAB_vessel_traffic_august2024.png"
  ),
  plot = map_traffic,
  width = 9,
  height = 9,
  dpi = 300
)


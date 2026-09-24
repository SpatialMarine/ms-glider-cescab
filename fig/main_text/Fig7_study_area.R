# ============================================================
# Ocean gliders for PAM
# ============================================================
# Greta Jankauskaite | @gretajan97
# ============================================================
# 01_glider_track_plot.R
# ============================================================

# Plot glider track within the study area.

source("setup.R")


# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------

indir <- file.path(output_dir, "cetaceans")

outdir <- file.path(output_dir, "cetaceans", "fig")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

nc_dir   <- file.path(input_dir, "glider_data")

# ------------------------------------------------------------
# 2. Import glider track
# ------------------------------------------------------------

nc_file <- file.path(
  nc_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)

track_nc <- nc_open(nc_file)

track <- tibble(
  longitude = as.numeric(ncvar_get(track_nc, "longitude")),
  latitude  = as.numeric(ncvar_get(track_nc, "latitude")),
  time      = as.numeric(ncvar_get(track_nc, "time"))
)

nc_close(track_nc)

track <- track %>%
  mutate(
    time = as.POSIXct(
      time,
      origin = "1970-01-01",
      tz = "UTC"
    )
  ) %>%
  arrange(time)


# ------------------------------------------------------------
# 3. Import layers for plotting
# ------------------------------------------------------------

# Bathymetry
bathy_file <- file.path(
  carto_dir,
  "EMODnet-bathymetry_2024",
  "bathymetry_dtm_2024_bb7b_a70e_b1f5.nc"
)

bathy <- rast(bathy_file)


# Bathymetry contours
cont_file <- file.path(
  carto_dir,
  "EMODnet_Bathymetry_2022_contours",
  "Contours_2022.shp"
)

contours <- st_read(cont_file)
contour_1000 <- subset(contours, Elevation == "1000.00")


# Spain boundaries
spain_file <- file.path(
  carto_dir,
  "spain_boundaries_GISData_MAPOG",
  "spain_Spain_Country_Boundary.shp"
)

spain <- st_read(spain_file)


# Cabrera Archipelago MPA
mpa_file <- file.path(
  carto_dir,
  "MPA",
  "Limite_Archipielago_Cabrera.shp"
)

cabrera_MPA <- st_read(mpa_file)
cabrera_MPA_pr <- st_transform(cabrera_MPA, crs(spain))


# World boundaries
world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)


# Western Mediterranean
westMed_file <- file.path(
  carto_dir,
  "eea-msfd",
  "western_mediterranean.shp"
)

westMed <- terra::vect(westMed_file)


# ------------------------------------------------------------
# 4. Prepare layers for glider track map
# ------------------------------------------------------------

# Crop bathymetry to CESCAB study area
bathy_crop <- crop(bathy, study_area)

# Convert bathymetry to dataframe
bathy_df <- as.data.frame(bathy_crop, xy = TRUE) %>%
  drop_na()

colnames(bathy_df) <- c("x", "y", "value")

# Keep depths below sea level
bathy_df <- bathy_df %>%
  filter(value < 1)


# Crop 1000 m contour to CESCAB study area
cont_int_1000 <- st_intersection(
  contour_1000,
  study_area
)


# ------------------------------------------------------------
# 5. Plot glider track
# ------------------------------------------------------------

glider_track_plot <- ggplot() +
  geom_raster(
    data = bathy_df,
    aes(x = x, y = y, fill = value)
  ) +
  scale_fill_gradient(
    low = "royalblue4",
    high = "lightskyblue1",
    na.value = "white",
    breaks = c(0, -2800),
    labels = c("0", "-2800"),
    name = "Depth (m)"
  ) +
  geom_sf(
    data = cabrera_MPA_pr,
    color = "black",
    fill = "darkblue",
    alpha = 0.2
  ) +
  geom_sf(
    data = spain,
    fill = "lightgrey",
    colour = "black",
    size = 0.1
  ) +
  geom_sf(
    data = cont_int_1000,
    colour = "grey1",
    linewidth = 0.07
  ) +
  geom_point(
    data = track,
    aes(x = longitude, y = latitude),
    color = "black",
    size = 1,
    shape = 21,
    fill = "black"
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 7, r = 0, b = 0, l = 0)
    ),
    axis.title.y = element_text(
      size = 13,
      margin = margin(t = 0, r = 7, b = 0, l = 0)
    ),
    legend.direction = "horizontal",
    legend.position = c(0.9, 0.06),
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
    ),
    size = guide_legend(
      title.position = "top"
    )
  ) +
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.2
  ) +
  ggspatial::annotation_north_arrow(
    location = "tl",
    which_north = "true",
    style = north_arrow_fancy_orienteering()
  ) +
  coord_sf(
    xlim = c(2.12, 3.47),
    ylim = c(38.22, 39.35)
  )

# Save
ggsave(
  filename = file.path(outdir, "study_area_track.png"),
  plot = glider_track_plot,
  width = 9,
  height = 11,
  bg = "transparent"
)


# ------------------------------------------------------------
# 6. Prepare layers for inset map
# ------------------------------------------------------------

# Crop bathymetry to Western Mediterranean
bathy_inset <- crop(bathy, westMed)

# Reduce resolution for inset map
bathy_inset <- aggregate(
  bathy_inset,
  fact = 10,
  fun = mean,
  na.rm = TRUE
)

# Convert bathymetry to dataframe
bathy_inset_df <- as.data.frame(
  bathy_inset,
  xy = TRUE,
  na.rm = FALSE
)

colnames(bathy_inset_df) <- c("x", "y", "value")

# Replace missing cells and keep depths below sea level
bathy_inset_df$value[is.na(bathy_inset_df$value)] <- 0

bathy_inset_df <- bathy_inset_df %>%
  filter(value < 1)

# Western Mediterranean extent
westMed_ext <- ext(westMed)


# ------------------------------------------------------------
# 7. Plot inset map
# ------------------------------------------------------------

inset_map <- ggplot() +
  geom_raster(
    data = bathy_inset_df,
    aes(x = x, y = y, fill = value)
  ) +
  scale_fill_gradient(
    low = "royalblue4",
    high = "lightskyblue1",
    na.value = "white"
  ) +
  geom_sf(
    data = world,
    fill = "lightgrey",
    colour = "grey50",
    linewidth = 0.1
  ) +
  geom_sf(
    data = study_area,
    fill = NA,
    colour = "black",
    linewidth = 0.5
  ) +
  coord_sf(
    xlim = c(westMed_ext[1], westMed_ext[2]),
    ylim = c(westMed_ext[3], westMed_ext[4]),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    legend.position = "none"
  )

# Save
ggsave(
  filename = file.path(outdir, "study_area_inset.png"),
  plot = inset_map,
  width = 10,
  height = 5,
  bg = "transparent"
)

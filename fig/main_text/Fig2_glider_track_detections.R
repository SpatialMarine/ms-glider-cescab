
# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# 01_glider_track_detection_plots.R
#
# Assign selection-level cetacean detections to the glider track and
# map their spatial distribution. annotations spanning multiple FLAC
# files are clipped to actual 30-s recording windows to exclude
# duty-cycle gaps before being matched to glider track timestamps.
# ============================================================

source("setup.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
indir    <- file.path(output_dir, "cetaceans", "final-annotation-table")

outdir <- file.path(output_dir, "cetaceans", "fig")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

flac_dir <- file.path(input_dir, "1_CESCAB-raw-flac")

nc_dir   <- file.path(input_dir, "glider_data")

# ------------------------------------------------------------
# 2. Import final annotation table
# ------------------------------------------------------------
annotations <- read.delim(
  file.path(indir, "2_CESCAB-annotation-table.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 3. Import glider track
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
# 4. Create complete screened FLAC file list
# ------------------------------------------------------------
flacs <- tibble(
  file = list.files(
    flac_dir,
    pattern = "\\.flac$",
    full.names = FALSE
  )
) %>%
  mutate(
    fileStartUTC = ymd_hms(
      str_remove(
        str_extract(file, "CESCAB_\\d{8}_\\d{6}"),
        "^CESCAB_"
      ),
      tz = "UTC"
    ),
    fileEndUTC = fileStartUTC + seconds(30)
  ) %>%
  filter(!is.na(fileStartUTC)) %>%
  distinct(file, .keep_all = TRUE) %>%
  arrange(fileStartUTC)

# ------------------------------------------------------------
# 5. Prepare selection-level detection intervals
# ------------------------------------------------------------
detection_intervals <- annotations %>%
  mutate(
    detStartUTC = parse_date_time(
      `Begin Date Time`,
      orders = c("ymd HMS", "ymd IMS"),
      tz = "UTC"
    ),
    detEndUTC = parse_date_time(
      paste(`End Date`, `End Clock Time`),
      orders = c("ymd HMS", "ymd IMS"),
      tz = "UTC"
    )
  ) %>%
  transmute(
    species = Species,
    signalType = `Signal Type`,
    detStartUTC,
    detEndUTC
  ) %>%
  filter(
    !is.na(species),
    !is.na(detStartUTC),
    !is.na(detEndUTC),
    detEndUTC > detStartUTC
  )

# ------------------------------------------------------------
# 6. Clip annotations to actual 30-s FLAC recording windows
# ------------------------------------------------------------
# This prevents duty-cycle gaps between consecutive recordings
# from being treated as periods containing detections.

detection_intervals_recorded <- flacs %>%
  left_join(
    detection_intervals,
    join_by(
      fileStartUTC < detEndUTC,
      fileEndUTC > detStartUTC
    )
  ) %>%
  filter(!is.na(species)) %>%
  mutate(
    detStartRecordedUTC = pmax(fileStartUTC, detStartUTC),
    detEndRecordedUTC   = pmin(fileEndUTC, detEndUTC)
  ) %>%
  filter(detEndRecordedUTC > detStartRecordedUTC) %>%
  select(
    file,
    species,
    signalType,
    detStartRecordedUTC,
    detEndRecordedUTC
  )

# ------------------------------------------------------------
# 7. Assign recorded detections to glider track
# ------------------------------------------------------------
# A glider-track point is marked only when its timestamp falls
# within the recorded portion of an annotated selection.

track_marked <- track %>%
  left_join(
    detection_intervals_recorded,
    join_by(
      time >= detStartRecordedUTC,
      time <= detEndRecordedUTC
    )
  ) %>%
  mutate(
    del_detection = species == "Delphinidae",
    phy_detection = species == "Physeter macrocephalus",
    gra_detection = species == "Grampus griseus"
  ) %>%
  group_by(longitude, latitude, time) %>%
  summarise(
    del_detection = any(del_detection, na.rm = TRUE),
    phy_detection = any(phy_detection, na.rm = TRUE),
    gra_detection = any(gra_detection, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    del_detection = replace_na(del_detection, FALSE),
    phy_detection = replace_na(phy_detection, FALSE),
    gra_detection = replace_na(gra_detection, FALSE)
  ) %>%
  filter(
    !is.na(longitude),
    !is.na(latitude)
  ) %>%
  arrange(time)

# ------------------------------------------------------------
# 8. Import spatial layers for plotting
# ------------------------------------------------------------

# Bathymetry
b_file <- paste0(
  carto_dir,
  "/EMODnet-bathymetry_2024/bathymetry_dtm_2024_bb7b_a70e_b1f5.nc"
)

bathy <- terra::rast(b_file)
bathy_crop <- terra::crop(bathy, study_area)

bathy_df <- as.data.frame(
  bathy_crop,
  xy = TRUE,
  na.rm = TRUE
)

colnames(bathy_df) <- c("x", "y", "value")

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

# World boundaries
w_file <- paste0(
  carto_dir,
  "/world_boundaries_shp/world.shp"
)

world <- st_read(w_file)

# Spain boundaries
s_file <- paste0(
  carto_dir,
  "/spain_boundaries_GISData_MAPOG/spain_Spain_Country_Boundary.shp"
)

spain <- st_read(s_file)

# Cabrera Archipelago National Park
mpa_file <- paste0(
  carto_dir,
  "/MPA/pn_cabrera.gpkg"
)

cabrera_MPA <- st_read(mpa_file)

cabrera_MPA_pr <- st_transform(
  cabrera_MPA,
  crs(spain)
)

# ------------------------------------------------------------
# 9. Create common map
# ------------------------------------------------------------
base_map <- ggplot() +
  geom_raster(
    data = bathy_df,
    aes(
      x = x,
      y = y,
      fill = value
    )
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
    color = NA,
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

# ------------------------------------------------------------
# 10. Plot sperm whale detections
# ------------------------------------------------------------
map_phy <- base_map +
  geom_point(
    data = track_marked %>%
      filter(!phy_detection),
    aes(
      x = longitude,
      y = latitude
    ),
    color = "black",
    fill = "black",
    size = 0.3,
    shape = 21
  ) +
  geom_point(
    data = track_marked %>%
      filter(phy_detection),
    aes(
      x = longitude,
      y = latitude
    ),
    color = "#FCA636FF",
    fill = "#FCA636FF",
    size = 1.2,
    shape = 21
  )

map_phy

ggsave(
  filename = file.path(outdir, "CESCAB_track_phy_det.png"),
  plot = map_phy,
  width = 9,
  height = 11,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Plot Delphinidae and Risso's dolphin detections
# ------------------------------------------------------------
map_del_gra <- base_map +
  geom_point(
    data = track_marked %>%
      filter(
        !del_detection,
        !gra_detection
      ),
    aes(
      x = longitude,
      y = latitude
    ),
    color = "black",
    fill = "black",
    size = 0.3,
    shape = 21
  ) +
  geom_point(
    data = track_marked %>%
      filter(del_detection),
    aes(
      x = longitude,
      y = latitude
    ),
    color = "#E16462FF",
    fill = "#E16462FF",
    size = 1.2,
    shape = 21
  ) +
  geom_point(
    data = track_marked %>%
      filter(gra_detection),
    aes(
      x = longitude,
      y = latitude
    ),
    color = "#A9C46C",
    fill = "#A9C46C",
    size = 1.2,
    shape = 21
  )

map_del_gra

ggsave(
  filename = file.path(outdir, "CESCAB_track_del_gra_det.png"),
  plot = map_del_gra,
  width = 9,
  height = 11,
  dpi = 300
)

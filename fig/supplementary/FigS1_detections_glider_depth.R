# ============================================================
# Ocean gliders for PAM
# ============================================================
# Greta Jankauskaite | @gretajan97
# ============================================================
# FigS1_detections_glider_depth.R
# ============================================================

# Plot the glider dive profile with SELECTION-LEVEL cetacean
# detections.
#
# Annotations spanning multiple files are clipped to the actual
# 30 s FLAC recording windows. Therefore, duty-cycle gaps are
# not shown as detections.
#
# Glider depth is interpolated at the exact beginning and end
# of each recorded selection interval so that short annotations
# are retained even if no NetCDF timestamp falls inside them.

source("setup.R")


# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------

flac_dir <- file.path(input_dir, "1_CESCAB-raw-flac")
nc_dir   <- file.path(input_dir, "glider_data")

indir <- file.path(output_dir, "cetaceans", "final-annotation-table")

outdir <- file.path(output_dir, "cetaceans", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Import final selection table
# ------------------------------------------------------------

annotations <- read.delim(
  file.path(indir, "2_CESCAB-annotation-table.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Import glider track and depth profile
# ------------------------------------------------------------´

nc_file <- file.path(
  nc_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)

track_nc <- nc_open(nc_file)

track <- tibble(
  time = as.numeric(ncvar_get(track_nc, "time")),
  depth_ctd = as.numeric(ncvar_get(track_nc, "depth_ctd"))
)

nc_close(track_nc)

track <- track %>%
  mutate(
    time = as.POSIXct(
      time,
      origin = "1970-01-01",
      tz = "UTC"
    ),
    depth = -abs(depth_ctd)
  ) %>%
  filter(
    !is.na(time),
    !is.na(depth),
    is.finite(depth)
  ) %>%
  distinct(time, .keep_all = TRUE) %>%
  arrange(time)

stopifnot(
  inherits(track$time, "POSIXct"),
  identical(attr(track$time, "tzone"), "UTC")
)


# ------------------------------------------------------------
# 4. Create complete screened FLAC file list
# ------------------------------------------------------------

flacs <- tibble(
  file = list.files(
    flac_dir,
    pattern = "\\.flac$",
    full.names = FALSE,
    ignore.case = TRUE
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

stopifnot(
  inherits(flacs$fileStartUTC, "POSIXct"),
  identical(attr(flacs$fileStartUTC, "tzone"), "UTC")
)


# ------------------------------------------------------------
# 5. Prepare selection-level detection intervals
# ------------------------------------------------------------

detection_intervals <- annotations %>%
  mutate(
    original_selection_id = row_number(),
    
    detStartUTC = parse_date_time(
      `Begin Date Time`,
      orders = c(
        "ymd HMS",
        "ymd HMSOS",
        "ymd IMS",
        "ymd IMSOS"
      ),
      tz = "UTC"
    ),
    
    detEndUTC = parse_date_time(
      paste(`End Date`, `End Clock Time`),
      orders = c(
        "ymd HMS",
        "ymd HMSOS",
        "ymd IMS",
        "ymd IMSOS"
      ),
      tz = "UTC"
    )
  ) %>%
  transmute(
    original_selection_id,
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
# 6. Clip annotations to actual FLAC windows
# ------------------------------------------------------------

detection_intervals_recorded <- flacs %>%
  left_join(
    detection_intervals,
    join_by(
      fileStartUTC < detEndUTC,
      fileEndUTC > detStartUTC
    ),
    relationship = "many-to-many"
  ) %>%
  filter(!is.na(species)) %>%
  mutate(
    detStartRecordedUTC = pmax(fileStartUTC, detStartUTC),
    detEndRecordedUTC   = pmin(fileEndUTC, detEndUTC)
  ) %>%
  filter(detEndRecordedUTC > detStartRecordedUTC) %>%
  arrange(original_selection_id, detStartRecordedUTC) %>%
  mutate(recorded_interval_id = row_number()) %>%
  select(
    recorded_interval_id,
    original_selection_id,
    file,
    species,
    signalType,
    detStartRecordedUTC,
    detEndRecordedUTC
  )


# ------------------------------------------------------------
# 7. Interpolate depth at selection boundaries
# ------------------------------------------------------------

track_time_numeric <- as.numeric(track$time)

interpolate_depth <- function(datetime) {
  approx(
    x = track_time_numeric,
    y = track$depth,
    xout = as.numeric(datetime),
    rule = 1
  )$y
}

selection_segments <- lapply(
  seq_len(nrow(detection_intervals_recorded)),
  function(i) {
    
    interval <- detection_intervals_recorded[i, ]
    
    start_time <- interval$detStartRecordedUTC
    end_time   <- interval$detEndRecordedUTC
    
    start_depth <- interpolate_depth(start_time)
    end_depth   <- interpolate_depth(end_time)
    
    if (is.na(start_depth) || is.na(end_depth)) return(NULL)
    
    internal_track <- track %>%
      filter(
        time > start_time,
        time < end_time
      ) %>%
      select(time, depth)
    
    bind_rows(
      tibble(
        time = start_time,
        depth = start_depth
      ),
      internal_track,
      tibble(
        time = end_time,
        depth = end_depth
      )
    ) %>%
      arrange(time) %>%
      mutate(
        recorded_interval_id = interval$recorded_interval_id,
        original_selection_id = interval$original_selection_id,
        species = interval$species,
        signalType = interval$signalType,
        file = interval$file
      )
  }
) %>%
  bind_rows()


# ------------------------------------------------------------
# 8. Define species colours
# ------------------------------------------------------------

species_colors <- c(
  "Delphinidae" = "#E16462FF",
  "Grampus griseus" = "#A9C46C",
  "Physeter macrocephalus" = "#FCA636FF"
)


# ------------------------------------------------------------
# 9. Plot glider depth profile with detections
# ------------------------------------------------------------

glider_depth_detections <- ggplot() +
  
  geom_line(
    data = track,
    aes(
      x = time,
      y = depth
    ),
    colour = "grey22",
    linewidth = 0.55,
    lineend = "round"
  ) +
  
  geom_line(
    data = selection_segments,
    aes(
      x = time,
      y = depth,
      colour = species,
      group = recorded_interval_id
    ),
    linewidth = 2.5,
    lineend = "round"
  ) +
  
  scale_colour_manual(
    values = species_colors,
    name = "Species",
    na.translate = FALSE
  ) +
  
  scale_x_datetime(
    timezone = "UTC",
    breaks = seq(
      as.POSIXct("2024-08-02 00:00:00", tz = "UTC"),
      max(track$time),
      by = "3 days"
    ),
    date_labels = "%m/%d",
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  
  coord_cartesian(
    xlim = c(
      as.POSIXct("2024-08-02 00:00:00", tz = "UTC"),
      max(track$time)
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(-1000, 0, by = 200),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    x = "Date (UTC)",
    y = "Glider depth (m)"
  ) +
  
  theme_bw(base_size = 16) +
  
  theme(
    panel.grid = element_blank(),
    
    axis.title.x = element_text(
      size = 16,
      margin = margin(t = 7)
    ),
    
    axis.title.y = element_text(
      size = 16,
      margin = margin(r = 7)
    ),
    
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 16
    ),
    
    axis.text.y = element_text(
      size = 16
    ),
    
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.4
    ),
    
    axis.ticks.length = unit(0.08, "cm"),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold", 
      size = 16
    ),
    
    legend.text = element_text(
      size = 16
    ),
    
    plot.margin = margin(
      20, 40, 20, 20
    )
  )

#glider_depth_detections


ggsave(
  filename = file.path(
    outdir,
    "detections-glider-depth.png"
  ),
  plot = glider_depth_detections,
  width = 40,
  height = 14,
  units = "cm",
  dpi = 300,
  bg = "white"
)


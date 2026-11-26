
# ============================================================
# Ocean gliders for PAM
# ============================================================
# Greta Jankauskaite | @gretajan97
# ============================================================
# FigS1_detections_glider_depth.R
# ============================================================

# Plot the glider dive profile with file-level cetacean presence.
# One point per species present in each 30-s recording.

source("setup.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------

indir <- file.path(output_dir, "cetaceans", "detection-stats")
nc_dir <- file.path(input_dir, "glider_data")

outdir <- file.path(output_dir, "cetaceans", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Import file-level presence/absence table
# ------------------------------------------------------------

files_presAbs <- read.csv(
  file.path(indir, "file_presence_absence.csv"),
  stringsAsFactors = FALSE
) %>%
  select(-any_of("X")) %>%
  mutate(
    fileStartUTC = ymd_hms(fileStartUTC, tz = "UTC"),
    fileEndUTC = ymd_hms(fileEndUTC, tz = "UTC"),
    fileMidUTC = fileStartUTC + (fileEndUTC - fileStartUTC) / 2
  ) %>%
  arrange(fileStartUTC)


# ------------------------------------------------------------
# 3. Import glider track and depth profile
# ------------------------------------------------------------

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
    is.finite(depth)
  ) %>%
  distinct(time, .keep_all = TRUE) %>%
  arrange(time)


# ------------------------------------------------------------
# 4. Match closest depth to each recording midpoint
# ------------------------------------------------------------

nearest_value <- function(time, value, target) {
  
  ok <- is.finite(as.numeric(time)) & is.finite(value)
  
  x <- as.numeric(time[ok])
  y <- value[ok]
  q <- as.numeric(target)
  
  ans <- rep(NA_real_, length(q))
  
  if (!length(x)) return(ans)
  
  valid <- which(is.finite(q))
  i <- findInterval(q[valid], x)
  
  before <- pmax(i, 1L)
  after <- pmin(i + 1L, length(x))
  
  nearest <- ifelse(
    abs(q[valid] - x[before]) <= abs(q[valid] - x[after]),
    before,
    after
  )
  
  ans[valid] <- y[nearest]
  
  ans
}

files_presAbs$depth <- nearest_value(
  track$time,
  track$depth,
  files_presAbs$fileMidUTC
)


# ------------------------------------------------------------
# 5. Prepare file-level detection points
# ------------------------------------------------------------

# Keep Delphinidae and Grampus separate, as in the original plot.
detection_points <- files_presAbs %>%
  select(
    file,
    fileMidUTC,
    depth,
    pres_Del,
    pres_Gg,
    pres_Pm
  ) %>%
  pivot_longer(
    cols = c(pres_Del, pres_Gg, pres_Pm),
    names_to = "species",
    values_to = "presence"
  ) %>%
  filter(
    presence == 1,
    !is.na(fileMidUTC),
    is.finite(depth)
  ) %>%
  mutate(
    species = factor(
      species,
      levels = c("pres_Del", "pres_Gg", "pres_Pm"),
      labels = c(
        "Delphinidae",
        "Grampus griseus",
        "Physeter macrocephalus"
      )
    )
  ) %>%
  distinct(file, species, .keep_all = TRUE)



# ------------------------------------------------------------
# 7. Prepare depth data for both plots
# ------------------------------------------------------------

# Combine unidentified delphinid and Risso's dolphin signals.
# Each recording counts once per signal group.
# Preserve missing values rather than treating them as absences.

depth_data <- files_presAbs %>%
  transmute(
    file,
    depth_m = abs(depth),
    sperm_clicks = as.logical(pres_Pm),
    delphinid_clicks = as.logical(del_clicks) |
      as.logical(Gg_clicks),
    delphinid_whistles = as.logical(del_whistles) |
      as.logical(Gg_whistles)
  ) %>%
  filter(
    is.finite(depth_m),
    between(depth_m, 0, 1000)
  ) %>%
  pivot_longer(
    cols = c(
      sperm_clicks,
      delphinid_clicks,
      delphinid_whistles
    ),
    names_to = "group",
    values_to = "presence"
  ) %>%
  mutate(presence = as.integer(presence)) %>%
  filter(!is.na(presence)) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "delphinid_clicks",
        "delphinid_whistles",
        "sperm_clicks"
      ),
      labels = c(
        "Delphinid echolocation clicks",
        "Delphinid whistles",
        "Sperm whale echolocation clicks"
      )
    ),
    presabs = factor(
      presence,
      levels = c(0, 1),
      labels = c("Absent", "Present")
    )
  ) %>%
  distinct(file, group, .keep_all = TRUE)

depth_colours <- c(
  "Sperm whale echolocation clicks" = "#FCA636FF",
  "Delphinid echolocation clicks" = "#7D2D3F",
  "Delphinid whistles" = "#E8B2A3"
)

depth_theme <- theme_bw(base_size = 13) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    panel.spacing = grid::unit(1, "lines"),
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey65"
    ),
    strip.text = element_text(size = 12),
    axis.text = element_text(colour = "gray20"),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )


# ------------------------------------------------------------
# 8. Density plots: presence versus absence
# ------------------------------------------------------------

# Plot absence curves in grey and presence curves in group colours.
# Density is normalised separately for each curve.
# Use a common smoothing bandwidth of 50 m.

depth_density <- ggplot() +
  geom_density(
    data = filter(depth_data, presence == 0),
    aes(
      x = depth_m,
      fill = "Absent",
      colour = "Absent"
    ),
    alpha = 0.4,
    linewidth = 0.6,
    bw = 50,
    trim = TRUE
  ) +
  geom_density(
    data = filter(depth_data, presence == 1),
    aes(
      x = depth_m,
      fill = group,
      colour = group
    ),
    alpha = 0.5,
    linewidth = 0.6,
    bw = 50,
    trim = TRUE
  ) +
  facet_wrap(
    ~ group,
    nrow = 1,
    scales = "free_y",
    labeller = label_wrap_gen(width = 24)
  ) +
  scale_fill_manual(
    values = c(
      "Absent" = "grey65",
      depth_colours
    ),
    breaks = "Absent"
  ) +
  scale_colour_manual(
    values = c(
      "Absent" = "grey65",
      depth_colours
    ),
    breaks = "Absent"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1000, by = 200),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(xlim = c(0, 1000)) +
  scale_y_continuous(
    limits = c(0, NA),
    breaks = scales::breaks_pretty(n = 4),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Glider depth (m)",
    y = "Density",
    caption = "Coloured curves: present; grey curves: absent."
  ) +
  depth_theme +
  theme(
    legend.position = "none",
    plot.caption = element_text(
      hjust = 0.5,
      size = 11,
      margin = margin(t = 10)
    )
  )

print(depth_density)


# ------------------------------------------------------------
# 9. Calculate detection percentages in 100 m depth bins
# ------------------------------------------------------------

depth_breaks <- seq(0, 1000, by = 100)

depth_labels <- paste0(
  head(depth_breaks, -1),
  "–",
  tail(depth_breaks, -1)
)

# Bins are [0,100), [100,200), etc.
# The final bin is 900–1000 m and includes exactly 1000 m.

depth_binned <- depth_data %>%
  mutate(
    bin_id = findInterval(
      depth_m,
      depth_breaks,
      rightmost.closed = TRUE
    ),
    depth_bin = factor(
      depth_labels[bin_id],
      levels = depth_labels
    )
  ) %>%
  group_by(group, depth_bin, .drop = FALSE) %>%
  summarise(
    n_recordings = n(),
    n_detections = sum(presence == 1),
    .groups = "drop"
  ) %>%
  mutate(
    detection_pct = if_else(
      n_recordings > 0,
      100 * n_detections / n_recordings,
      NA_real_
    )
  )


# ------------------------------------------------------------
# 10. Combined depth-percentage figure
# ------------------------------------------------------------

depth_binned <- depth_binned %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Delphinid echolocation clicks",
        "Delphinid whistles",
        "Sperm whale echolocation clicks"
      )
    )
  )

depth_percentage <- ggplot(
  depth_binned,
  aes(
    x = depth_bin,
    y = detection_pct,
    fill = group
  )
) +
  geom_col(
    width = 0.9,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ group,
    nrow = 1,
    scales = "fixed",
    axes = "all_y",
    axis.labels = "all_y"
  ) +
  scale_fill_manual(
    values = c(
      "Delphinid echolocation clicks" = "#7D2D3F",
      "Delphinid whistles" = "#E8B2A3",
      "Sperm whale echolocation clicks" = "#FCA636FF"
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    drop = FALSE,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Glider depth (m)",
    y = "Percentage of data with target signals"
  ) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 15),
    axis.text = element_text(
      size = 14,
      colour = "gray20"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "gray40"),
    panel.spacing = grid::unit(1.5, "lines"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  )

#print(depth_percentage)


# ------------------------------------------------------------
# 11. Save combined figure
# ------------------------------------------------------------

ggsave(
  file.path(outdir, "detections-glider-depth-percentage-100m.png"),
  plot = depth_percentage,
  width = 15,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

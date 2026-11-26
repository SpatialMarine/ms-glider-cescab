
# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# Fig4_diel_patterns_plots.R
# ============================================================

source("setup.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
indir <- file.path(output_dir, "cetaceans", "detection-stats")

outdir <- file.path(output_dir, "cetaceans", "fig")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------
# 2. Import file-level presence/absence table
# ------------------------------------------------------------
files_presAbs <- read.csv(
  file.path(indir, "file_presence_absence.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # The timestamps stored in the CSV represent UTC
    fileStartUTC = lubridate::ymd_hms(fileStartUTC, tz = "UTC"),
    fileEndUTC   = lubridate::ymd_hms(fileEndUTC,   tz = "UTC")
  )

# Import day/night boundaries derived from deployment-median sunrise and sunset times 
# calculated using "suncalc" package and rounded to the nearest hour.
diel_boundaries <- read.csv(
  file.path(
    output_dir, "cetaceans", "glm-clean", "scratch", "diel_boundaries.csv"
  )
)

day_start <- diel_boundaries$day_start
night_start <- diel_boundaries$night_start

# ------------------------------------------------------------
# 3. Convert recording midpoints from UTC to local time
# ------------------------------------------------------------
files_presAbs <- files_presAbs %>%
  mutate(
    fileMidUTC = fileStartUTC + (fileEndUTC - fileStartUTC) / 2,
    fileMidLocal = lubridate::with_tz(
      fileMidUTC,
      tzone = "Europe/Madrid"
    ),
    hourLocal = lubridate::hour(fileMidLocal)
  )

# ------------------------------------------------------------
# 4. Calculate hourly recording effort
# ------------------------------------------------------------
effort_hour <- files_presAbs %>%
  count(hourLocal, name = "total_files") %>%
  complete(hourLocal = 0:23, fill = list(total_files = 0)) %>%
  mutate(
    observation_minutes = total_files * 30 / 60
  )

effort_summary <- effort_hour %>%
  summarise(
    mean_obs = mean(observation_minutes),
    sd_obs   = sd(observation_minutes),
    min_obs  = min(observation_minutes),
    max_obs  = max(observation_minutes)
  )

effort_summary

# ------------------------------------------------------------
# 5. Prepare hourly detection percentages by species group
# ------------------------------------------------------------
detected_hour <- bind_rows(
  files_presAbs %>%
    group_by(hourLocal) %>%
    summarise(
      detected_files = sum(pres_Del | pres_Gg),
      .groups = "drop"
    ) %>%
    mutate(group = "Delphinidae"),
  
  files_presAbs %>%
    group_by(hourLocal) %>%
    summarise(
      detected_files = sum(pres_Pm),
      .groups = "drop"
    ) %>%
    mutate(group = "Physeter macrocephalus")
) %>%
  complete(
    group = c("Delphinidae", "Physeter macrocephalus"),
    hourLocal = 0:23,
    fill = list(detected_files = 0)
  )

plot_data <- detected_hour %>%
  left_join(effort_hour, by = "hourLocal") %>%
  mutate(
    percentage = ifelse(
      total_files > 0,
      100 * detected_files / total_files,
      NA_real_
    ),
    hour = factor(sprintf("%02d", hourLocal),
                  levels = sprintf("%02d", 0:23))
  )


# ------------------------------------------------------------
# 7. Plot sperm whale diel detections
# ------------------------------------------------------------
sperm_wh_p <- ggplot(
  filter(plot_data, group == "Physeter macrocephalus"),
  aes(x = hour, y = percentage)
) +
  annotate(
    "rect",
    xmin = 0, xmax = day_start + 1,
    ymin = -Inf, ymax = Inf,
    fill = "grey90", alpha = 0.7
  ) +
  annotate(
    "rect",
    xmin = night_start + 1, xmax = 25,
    ymin = -Inf, ymax = Inf,
    fill = "grey90", alpha = 0.7
  ) +
  geom_bar(stat = "identity", fill = "#FCA636FF", alpha = 0.9) +
  labs(
    x = "Hour of day",
    y = "Percentage of data with target signals"
  ) +
  scale_x_discrete(
    drop = FALSE,
    expand = c(0, 0),
    breaks = sprintf("%02d", seq(0, 23, by = 2))
  ) +
  scale_y_continuous(
    limits = c(0, 35),
    breaks = seq(0, 35, by = 5)
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14, color = "gray20"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "gray40"),
    legend.position = "none"
  )


ggsave(
  file.path(outdir, "sp_wh_diel.png"),
  sperm_wh_p,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Prepare hourly detection percentages by signal type
# ------------------------------------------------------------
detected_hour_signalType <- bind_rows(
  files_presAbs %>%
    group_by(hourLocal) %>%
    summarise(
      detected_files = sum(del_clicks | Gg_clicks),
      .groups = "drop"
    ) %>%
    mutate(signalType = "clicks"),
  
  files_presAbs %>%
    group_by(hourLocal) %>%
    summarise(
      detected_files = sum(del_whistles | Gg_whistles),
      .groups = "drop"
    ) %>%
    mutate(signalType = "whistles")
) %>%
  complete(
    hourLocal = 0:23,
    signalType = c("clicks", "whistles"),
    fill = list(detected_files = 0)
  )

plot_data_signalType <- detected_hour_signalType %>%
  left_join(effort_hour, by = "hourLocal") %>%
  mutate(
    percentage = ifelse(
      total_files > 0,
      100 * detected_files / total_files,
      NA_real_
    ),
    hour = factor(sprintf("%02d", hourLocal),
                  levels = sprintf("%02d", 0:23))
  )


# ------------------------------------------------------------
# 9. Plot clicks separately
# ------------------------------------------------------------
click_p <- ggplot(
  filter(plot_data_signalType, signalType == "clicks"),
  aes(x = hour, y = percentage)
) +
  annotate(
    "rect",
    xmin = 0, xmax = day_start + 1,
    ymin = -Inf, ymax = Inf,
    fill = "grey90", alpha = 0.7
  ) +
  annotate(
    "rect",
    xmin = night_start + 1, xmax = 25,
    ymin = -Inf, ymax = Inf,
    fill = "grey90", alpha = 0.7
  ) +
  geom_bar(stat = "identity", fill = "#7D2D3F", alpha = 0.9) +
  labs(
    x = "Hour of day",
    y = "Percentage of data with target signals"
  ) +
  scale_x_discrete(
    drop = FALSE,
    expand = c(0, 0),
    breaks = sprintf("%02d", seq(0, 23, by = 2))
  ) +
  scale_y_continuous(
    limits = c(0, 35),
    breaks = seq(0, 35, by = 5) # standardise for figure
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14, color = "gray20"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "gray40"),
    legend.position = "none"
  )


ggsave(
  file.path(outdir, "click_diel.png"),
  click_p,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Plot whistles separately
# ------------------------------------------------------------
whistle_p <- ggplot(
  filter(plot_data_signalType, signalType == "whistles"),
  aes(x = hour, y = percentage)
) +
  annotate(
    "rect",
    xmin = 0, xmax = day_start + 1,
    ymin = -Inf, ymax = Inf,
    fill = "grey90", alpha = 0.7
  ) +
  annotate(
    "rect",
    xmin = night_start + 1, xmax = 25,
    ymin = -Inf, ymax = Inf,
    fill = "grey90", alpha = 0.7
  ) +
  geom_bar(stat = "identity", fill = "#E8B2A3", alpha = 0.9) +
  labs(
    x = "Hour of day",
    y = "Percentage of data with target signals"
  ) +
  scale_x_discrete(
    drop = FALSE,
    expand = c(0, 0),
    breaks = sprintf("%02d", seq(0, 23, by = 2))
  ) +
  scale_y_continuous(
    limits = c(0, 35),
    breaks = seq(0, 35, by = 5)
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14, color = "gray20"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "gray40"),
    legend.position = "none"
  )


ggsave(
  file.path(outdir, "whistles_diel.png"),
  whistle_p,
  width = 7,
  height = 6,
  dpi = 300
)



# ============================================================
# Ocean gliders for PAM
# ============================================================
# Greta Jankauskaite | @gretajan97
# ============================================================
# Fig6_daily_detection_percentages.R
# ============================================================

# Calculate daily file-level detection percentages per species.
# Uses "file_presence_absence.csv" created in "analysis/cetaceans/01_detection_percentages.R".

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
files_presAbs <- read.csv(file.path(indir, "file_presence_absence.csv"))

# ------------------------------------------------------------
# 3. Add local time and date
# ------------------------------------------------------------
files_presAbs <- files_presAbs %>%
  mutate(
    fileStartLocal = with_tz(fileStartUTC, "Europe/Madrid"),
    hourLocal = hour(fileStartLocal),
    Date = as.Date(fileStartLocal)
  )

# ------------------------------------------------------------
# 4. Calculate daily effort
# ------------------------------------------------------------
effort_day <- files_presAbs %>%
  count(Date, name = "total_files")

# ------------------------------------------------------------
# 5. Calculate detected files per day by species
# ------------------------------------------------------------
detected_day <- bind_rows(
  files_presAbs %>%
    group_by(Date) %>%
    summarise(detected_files = sum(pres_Pm), .groups = "drop") %>%
    mutate(species = "Physeter macrocephalus"),
  
  files_presAbs %>%
    group_by(Date) %>%
    summarise(detected_files = sum(pres_Del), .groups = "drop") %>%
    mutate(species = "Delphinidae"),
  
  files_presAbs %>%
    group_by(Date) %>%
    summarise(detected_files = sum(pres_Gg), .groups = "drop") %>%
    mutate(species = "Grampus griseus")
)

# ------------------------------------------------------------
# 6. Create complete daily percentage table
# ------------------------------------------------------------
daily_counts <- expand.grid(
  species = c(
    "Physeter macrocephalus",
    "Delphinidae",
    "Grampus griseus"
  ),
  Date = seq(min(files_presAbs$Date), max(files_presAbs$Date), by = "day"),
  stringsAsFactors = FALSE
) %>%
  left_join(detected_day, by = c("species", "Date")) %>%
  left_join(effort_day, by = "Date") %>%
  mutate(
    detected_files = replace_na(detected_files, 0),
    percentage = ifelse(
      total_files > 0,
      100 * detected_files / total_files,
      NA_real_
    ),
    species = factor(
      species,
      levels = c(
        "Grampus griseus",
        "Delphinidae",
        "Physeter macrocephalus"
      )
    )
  )

# ------------------------------------------------------------
# 7. Prepare data for plotting
# ------------------------------------------------------------
daily_counts_plot <- daily_counts %>%
  filter(Date >= as.Date("2024-08-02")) %>%
  mutate(
    perc_bin = case_when(
      percentage < 10 ~ "<10%",
      percentage >= 10 & percentage < 30 ~ "10–30%",
      percentage >= 30 ~ "30–55%"
    ),
    perc_bin = factor(
      perc_bin,
      levels = c("<10%", "10–30%", "30–55%")
    )
  )

taxa_lines <- expand.grid(
  species = levels(daily_counts$species),
  Date = seq(
    as.Date("2024-08-02"),
    max(daily_counts$Date),
    by = "day"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    species = factor(species, levels = levels(daily_counts$species))
  )

max_daily_percentage <- max(daily_counts$percentage, na.rm = TRUE)
print(max_daily_percentage)

# ------------------------------------------------------------
# 8. Plot daily detection percentages
# ------------------------------------------------------------
daily_detection_plot <- ggplot() +
  geom_line(
    data = taxa_lines,
    aes(x = Date, y = species, group = species),
    color = "black",
    linewidth = 0.3
  ) +
  geom_point(
    data = daily_counts_plot %>% filter(percentage > 0),
    aes(x = Date, y = species, color = species, size = perc_bin),
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "Delphinidae" = "#E16462FF",
      "Physeter macrocephalus" = "#FCA636FF",
      "Grampus griseus" = "#A9C46C"
    ),
    name = "Species"
  ) +
  scale_size_manual(
    name = "Percentage (%)",
    values = c(
      "<10%" = 3,
      "10–30%" = 6,
      "30–55%" = 9
    ),
    drop = FALSE
  ) +
  scale_x_date(
    limits = c(as.Date("2024-08-02"), max(daily_counts$Date)),
    breaks = seq(
      as.Date("2024-08-02"),
      max(daily_counts$Date),
      by = "3 days"
    ),
    date_labels = "%m/%d",
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(x = "Date", y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", margin = margin(r = 20)),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 40, 20, 20)
  )

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------
ggsave(
  filename = file.path(outdir, "CESCAB_daily_detection_percentages_species.png"),
  plot = daily_detection_plot,
  width = 40,
  height = 7,
  units = "cm",
  dpi = 300
)


# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# Fig8_wind_8k_SPL_regression.R
#
# Calculate median 8-kHz third-octave SPLs per recording and match
# them to the nearest wind observations. Plot SPL against wind
# speed and fit linear regressions using log10-transformed wind
# speed, including a separate model for wind speeds ≥3 m/s.
# ============================================================

source("setup.R")

outdir <- paste0(output_dir, "/noise/fig")

# import data

# Import cleaned SPL samples.
noise <- readRDS(
  file.path(output_dir, "noise", "timeseries", "CESCAB_SPL_63_125_8k.rds")
)

# Calculate one median 8-kHz SPL per 30-second recording.
noise8k <- noise %>%
  filter(is.finite(Hz8k)) %>%
  group_by(FileStart) %>%
  summarise(
    Hz8k_median = median(Hz8k),
    .groups = "drop"
  ) %>%
  rename(ValueTime = FileStart)


# wind
wind <- read.csv(paste0(output_dir, "/wind-speed.csv"), colClasses = "character")

# fix missing time in dates
date_only <- nchar(wind$time) == 10
wind$time[date_only] <- paste0(wind$time[date_only], " 00:00:00")

# convert formats
wind$time <- as.POSIXct(wind$time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
wind$ws <- as.numeric(wind$ws)


# # --- JOIN DATA (nearest time) ---
# 
# data <- noise8k %>%
#   fuzzyjoin::difference_inner_join(
#     wind,
#     by = c("ValueTime" = "time"),
#     max_dist = as.difftime(1, units = "hours")
#   ) %>%
#   group_by(ValueTime) %>%
#   slice_min(abs(ValueTime - time)) %>%
#   ungroup()
# 
# # line data: keep your existing bins, but make line start at 2 m/s
# line_data <- bin_stats %>%
#   filter(ws_center >= 2) %>%
#   { 
#     add_row(
#       .,
#       ws_center = 2,
#       median_noise = .$median_noise[1],
#       .before = 1
#     )
#   }

# Match each recording to the nearest wind observation within one hour.
# If two observations are equally close, retain the earlier one.
data <- noise8k %>%
  fuzzyjoin::difference_inner_join(
    wind,
    by = c("ValueTime" = "time"),
    max_dist = as.difftime(1, units = "hours")
  ) %>%
  arrange(ValueTime, time) %>%
  group_by(ValueTime) %>%
  slice_min(
    abs(as.numeric(difftime(ValueTime, time, units = "secs"))),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

# Keep observations within the plotted wind-speed range.
df <- data %>%
  filter(ws >= 2, ws <= 10.1)

# Calculate median SPL within 1 m/s wind-speed bins.
bin_stats <- df %>%
  mutate(ws_center = floor(ws) + 0.5) %>%
  group_by(ws_center) %>%
  summarise(
    median_noise = median(Hz8k_median),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(ws_center)

# Values for the dark blue line.
line_data <- bin_stats


# plot
p <- ggplot(df %>% filter(ws >= 2), aes(x = ws, y = Hz8k_median)) +
  
  geom_point(shape = 16, size = 2, color = "lightskyblue1", alpha = 0.9) +
  
  geom_line(
    data = line_data,
    aes(x = ws_center, y = median_noise),
    color = "royalblue4",
    linewidth = 1.2
  ) +
  
  scale_x_log10(
    limits = c(2, 10.1),
    breaks = c(2, 5, 10),
    labels = c("2", "5", "10"),
    minor_breaks = 2:10
  ) +
  
  scale_y_continuous(
    breaks = seq(70, 110, 5),
    minor_breaks = seq(70, 110, 1)
  ) +
  
  labs(
    x = "Wind speed [m/s]",
    y = "8 kHz 1/3-octave band SPL [dB re 1 µPa]"
  ) +
  
  theme_article() +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey70", linetype = "dotted"),
    panel.grid.minor = element_line(linewidth = 0.4, colour = "grey65", linetype = "dotted"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15)
  )

#p

# save
#ggsave(filename = file.path(outdir, "noise-vs-wind_medians.png"), plot = p, width = 6, height = 6, dpi = 300)





model <- lm(
  Hz8k_median ~ log10(ws),
  data = data %>% filter(is.finite(ws), ws > 0)
)

summary(model)


# Data used for the regression.
reg_data <- data %>%
  filter(is.finite(ws), ws > 0, is.finite(Hz8k_median))

model <- lm(Hz8k_median ~ log10(ws), data = reg_data)

# Predictions from the fitted model.
prediction <- data.frame(
  ws = 10^seq(
    log10(min(reg_data$ws)),
    log10(max(reg_data$ws)),
    length.out = 200
  )
)
prediction$SPL <- predict(model, newdata = prediction)

# Equation and R² from the fitted model.
equation <- sprintf(
  "SPL = %.2f + %.2f log10(U)\nR² = %.3f",
  coef(model)[1],
  coef(model)[2],
  summary(model)$r.squared
)

p <- ggplot(reg_data, aes(x = ws, y = Hz8k_median)) +
  geom_point(
    shape = 16, size = 2,
    colour = "lightskyblue1", alpha = 0.9
  ) +
  geom_line(
    data = prediction,
    aes(x = ws, y = SPL),
    colour = "royalblue4",
    linewidth = 1.2
  ) +
  annotate(
    "text",
    x = min(reg_data$ws), y = Inf,
    label = equation,
    hjust = 0, vjust = 1.3,
    size = 4.5
  ) +
  scale_x_log10() +
  scale_y_continuous(
    breaks = seq(70, 110, 5),
    minor_breaks = seq(70, 110, 1)
  ) +
  labs(
    x = "Wind speed (m/s)",
    y = "8 kHz 1/3-octave band SPL (dB re 1 µPa)"
  ) +
  theme_article() +
  theme(
    panel.grid.major = element_line(
      linewidth = 0.3, colour = "grey70", linetype = "dotted"
    ),
    panel.grid.minor = element_line(
      linewidth = 0.4, colour = "grey65", linetype = "dotted"
    ),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15)
  )

print(p)

ggsave(
  filename = file.path(outdir, "noise-vs-wind_regression.png"),
  plot = p,
  width = 6, height = 6, dpi = 300
)



model_3 <- lm(
  Hz8k_median ~ log10(ws),
  data = data %>%
    filter(is.finite(ws), ws >= 3, is.finite(Hz8k_median))
)

summary(model_3)





















# ============================================================
# Ocean gliders for PAM
# ============================================================
#
# Greta Jankauskaite | @gretajan97
#
# ============================================================
# FigS7_sound_speed_profile.R
# ============================================================
#
# This script calculates a mission-averaged vertical sound-speed
# profile from the temperature, salinity, and pressure measurements
# collected by the CESCAB Slocum glider.
#
# Practical Salinity is first converted to Absolute Salinity using
# the TEOS-10 Gibbs SeaWater functions. Sound speed is then calculated
# for every valid CTD observation and averaged within 5 m depth bins.
#
# The resulting profile is plotted against depth and saved as a
# high-resolution PNG. The depth bin containing the minimum mean
# sound speed is also reported.
#

# ============================================================
# 1. Packages
# ============================================================

library(gsw)

# ============================================================
# 2. Input and output paths
# ============================================================

input_dir <- paste0(
  "C:/Users/greta/SML Dropbox/Dream Team/data/campaign/",
  "240801-CESCAB/glider/oceanographic&track"
)

nc_file <- file.path(
  input_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)

outdir <- file.path(output_dir, "noise", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile <- file.path(
  outdir,
  "sound_speed_profile.png"
)

if (!file.exists(nc_file)) {
  stop("NetCDF file not found: ", nc_file)
}


# ============================================================
# 3. Read CTD variables from the NetCDF file
# ============================================================

nc <- nc_open(nc_file)

longitude   <- ncvar_get(nc, "longitude")
latitude    <- ncvar_get(nc, "latitude")
depth_ctd   <- ncvar_get(nc, "depth_ctd")
pressure    <- ncvar_get(nc, "pressure")
salinity    <- ncvar_get(nc, "salinity")
temperature <- ncvar_get(nc, "temperature")

nc_close(nc)


# ============================================================
# 4. Create and clean the CTD dataset
# ============================================================

ctd <- data.frame(
  longitude   = as.numeric(longitude),
  latitude    = as.numeric(latitude),
  depth_ctd   = as.numeric(depth_ctd),
  pressure    = as.numeric(pressure),
  salinity    = as.numeric(salinity),
  temperature = as.numeric(temperature)
) %>%
  filter(
    is.finite(longitude),
    is.finite(latitude),
    is.finite(depth_ctd),
    is.finite(pressure),
    is.finite(salinity),
    is.finite(temperature),
    depth_ctd >= 0,
    pressure >= 0,
    salinity >= 0
  )

if (nrow(ctd) == 0) {
  stop("No valid CTD observations remained after filtering.")
}

message(
  "Valid CTD observations: ",
  format(nrow(ctd), big.mark = ",")
)

message(
  "Maximum CTD depth: ",
  round(max(ctd$depth_ctd), 1),
  " m"
)


# ============================================================
# 5. Calculate sound speed
# ============================================================

ctd <- ctd %>%
  mutate(
    # Convert Practical Salinity to Absolute Salinity (g/kg)
    absolute_salinity = gsw_SA_from_SP(
      SP        = salinity,
      p         = pressure,
      longitude = longitude,
      latitude  = latitude
    ),
    
    # Calculate sound speed from Absolute Salinity,
    # in-situ temperature, and sea pressure
    sound_speed = gsw_sound_speed_t_exact(
      SA = absolute_salinity,
      t  = temperature,
      p  = pressure
    )
  ) %>%
  filter(is.finite(sound_speed))


# ============================================================
# 6. Calculate the mean profile in 5 m depth bins
# ============================================================

bin_m <- 5

profile <- ctd %>%
  mutate(
    depth_bin = floor(depth_ctd / bin_m) * bin_m
  ) %>%
  group_by(depth_bin) %>%
  summarise(
    sound_speed = mean(sound_speed),
    n_observations = n(),
    .groups = "drop"
  ) %>%
  arrange(depth_bin)

if (nrow(profile) == 0) {
  stop("The binned sound-speed profile could not be calculated.")
}


# ============================================================
# 7. Identify the sound-speed minimum
# ============================================================

sound_speed_minimum <- profile %>%
  slice_min(
    order_by = sound_speed,
    n = 1,
    with_ties = FALSE
  )

message(
  "Minimum mean sound speed: ",
  round(sound_speed_minimum$sound_speed, 2),
  " m/s at the ",
  sound_speed_minimum$depth_bin,
  " m depth bin."
)

print(sound_speed_minimum)


# ============================================================
# 8. Plot the sound-speed profile
# ============================================================

maximum_plot_depth <- ceiling(
  max(profile$depth_bin) / 100
) * 100

p <- ggplot(
  profile,
  aes(
    x = sound_speed,
    y = depth_bin
  )
) +
  geom_path(linewidth = 0.9) +
  scale_y_reverse(
    limits = c(maximum_plot_depth, 0),
    breaks = seq(
      from = 0,
      to = maximum_plot_depth,
      by = 200
    ),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    position = "top",
    breaks = scales::breaks_pretty(n = 6)
  ) +
  labs(
    x = expression("Sound speed (m s"^{-1}*")"),
    y = "Depth (m)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3),
    axis.title.x = element_text(
      margin = margin(b = 10)
    ),
    axis.title.y = element_text(
      margin = margin(r = 10)
    )
  )

print(p)


# ============================================================
# 9. Save the figure
# ============================================================

ggsave(
  filename = outfile,
  plot = p,
  width = 6,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

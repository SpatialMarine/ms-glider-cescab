# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# FigS6_ERA5_wind_speed_timeseries.R
#
# Calculate 10 m wind speed from ERA5 u and v components and
# extract hourly values along the glider track using the nearest
# recorded positions and bilinear spatial interpolation.
# Save the wind-speed time series as a CSV and plot.
# ============================================================

source("setup.R")

# import track------------------------------------------------------------------------------
# set input dir
indir <- "C:/Users/greta/SML Dropbox/data/campaign/240801-CESCAB/glider/oceanographic&track/"

nc_dir   <- file.path(input_dir, "glider_data")

outdir <- file.path(output_dir, "noise", "fig")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

### import file
nc_file <- file.path(
  nc_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)

track <- nc_open(nc_file)

# get latitude and longitude 
lon <- ncvar_get(track, "longitude")
lat <- ncvar_get(track, "latitude")

nc_close(track)

# combine lat and lon into a data frame
coords <- data.frame(lon = as.numeric(lon), lat = as.numeric(lat))


# get time (Albert said use "time", not "time_ctd")
time <- as.data.frame(ncvar_get(track, "time")) 
names(time) <- c("time")
# remove the 0s, make them NAs
time <- time %>%
  mutate(
    .num = suppressWarnings(as.numeric(as.character(time))),
    time = ifelse(!is.na(.num) & .num == 0, NA, time)
  ) %>%
  dplyr::select(-.num)
# Convert the time to UTC
time$time <- as.POSIXct(time$time, origin = "1970-01-01", tz = "UTC")

# combine 
data <- cbind(coords, time)
# filter out rows with NA in either coords or time
track_df <- data %>% 
  filter(!is.na(lat) & !is.na(lon) & !is.na(time)) # remove all rows with NAs from lat, lon and time (time is sampled every 8-9 gps points or 6 minutes)

# import wind data

# read the GRIB file----------------------------------------------------------
input_dir <- paste0(input_dir, "/wind-speed")
file_path <- file.path(input_dir, "ae443b3896c3b5f7977a9fac1bbe9fc5.grib")
data_aug <- rast(file_path)

file_path <- file.path(input_dir, "cee3fe406ed2498efbdfaad22d6cea54.grib")
data_sep <- rast(file_path)

# Calculate u and v for august
nlay <- nlyr(data_aug)

# Odd layers: 1,3,5,... = u
u <- data_aug[[seq(1, nlay, by = 2)]]

# Even layers: 2,4,6,... = v
v <- data_aug[[seq(2, nlay, by = 2)]]

# check that timestamps align
time(u)[1:5]
time(v)[1:5]

# compute wind speed from u and v components 
ws_aug <- sqrt(u^2 + v^2)

# Carry over time (terra should keep it, but just in case)
time(ws_aug) <- time(u)

# Give nicer names if you like
names(ws_aug) <- paste0("WS10_", format(time(ws_aug), "%Y%m%d_%H"))

# print created wind speed layer
ws_aug


# September now
nlay <- nlyr(data_sep)

# Odd layers: 1,3,5,... = u
u <- data_sep[[seq(1, nlay, by = 2)]]

# Even layers: 2,4,6,... = v
v <- data_sep[[seq(2, nlay, by = 2)]]

# check that timestamps align
time(u)[1:5]
time(v)[1:5]


# compute wind speed from u and v components 
ws_sep <- sqrt(u^2 + v^2)

# Carry over time (terra should keep it, but just in case)
time(ws_sep) <- time(u)

# Give nicer names if you like
names(ws_sep) <- paste0("WS10_", format(time(ws_sep), "%Y%m%d_%H"))

# print created wind speed layer
ws_sep


# combine into one wind speed raster for the entire CESCAB campaign
ws_cescab <- c(ws_aug, ws_sep)
ws_cescab


# Move on to extracting wind speed values by taking closest timestamp coordinates and extracting-----------------------

# grab the ERA5 timestamps
t_ws <- time(ws_cescab)
head(t_ws)
tail(t_ws)

# make sure times ok - has to be posixct and UTC
str(track_df$time)
attr(track_df$time, "tzone") # all good

# set start and end times of CESCAB mission
mission_start <- as.POSIXct("2024-08-01 14:00:00", tz = "UTC") # "2024-08-01 14:04:27 UTC" (first wind speed at 3PM) 
# but we'll set it to 14:00:00 cause track is available 14:04:27 and first audio file 14:17, dont want to lose that
mission_end   <- max(track_df$time)  # "2024-09-02 08:52:42 UTC" (last wind speed at 8AM)

# Keep only ERA5 hours within mission
keep <- t_ws >= mission_start & t_ws <= mission_end

ws_mission <- ws_cescab[[keep]]  # subset of layers
t_mission  <- t_ws[keep]         # matching times

ws_mission
head(t_mission)
tail(t_mission)


# Make sure track_df is a data.table keyed by time
setDT(track_df)
setkey(track_df, time)

# For each ERA5 time in t_mission, find the nearest glider position
# use nearest because times in "t_mission" are taken from wind speed layer, exactly 1722560400 seconds since 1970 → 01:00:00.000
# Some glider timestamps are exact seconds, Some are slightly off (e.g. +0.848 s), R compares them at microsecond level, so %in% / match() give a lot of NA
track_hourly <- track_df[.(t_mission), roll = "nearest"] 


head(track_hourly)
nrow(track_hourly)    # should be equal to length(t_mission) = nlyr(ws_mission)


# Extract wind speed for the glider track-------------------------------------------

# points along the glider track, one per hour
pts <- vect(track_hourly[, .(lon, lat)], crs = crs(ws_mission))

n <- nlyr(ws_mission)   # should be 762

ws_vals <- sapply(seq_len(n), function(i) {
  # extract() returns a data.frame: ID + value
  terra::extract(ws_mission[[i]], pts[i, ], method = "bilinear")[1, 2]
})

ws_ts <- data.frame(
  time = t_mission,          # or track_hourly$time (same)
  lon  = track_hourly$lon,
  lat  = track_hourly$lat,
  ws   = ws_vals             # 10 m wind speed [m/s]
)
head(ws_ts)

# save csv
# Create the full file path
output_dir <- paste(main_dir, "output", sep="/")
outfile <- file.path(output_dir, "wind-speed.csv")
# Save as CSV
write.csv(ws_ts, outfile, row.names = FALSE)


# plot time wind speed time series
wind <- ggplot(ws_ts, aes(x = time, y = ws)) +
  geom_line(colour = "royalblue4",  linewidth = 0.8) +
  labs(
    x = "Time (UTC)",
    y = "10 m wind speed (m/s)",
    title = "ERA 5 Wind Speed"
  ) +
  theme_bw()
# save plot
ggsave(filename = file.path(output_dir, "/noise/fig/wind-speed-time-series-dark.png"), plot = wind, width = 8, height = 4, dpi = 300)






# Average hourly ERA5 wind speeds along the track within each UTC day.
wind_daily <- ws_ts %>%
  dplyr::filter(is.finite(ws)) %>%
  dplyr::group_by(date = as.Date(time, tz = "UTC")) %>%
  dplyr::summarise(
    mean_ws = mean(ws),
    n_hours = dplyr::n(),
    .groups = "drop"
  )

wind_daily_plot <- ggplot(wind_daily, aes(x = date, y = mean_ws)) +
  geom_line(colour = "royalblue4", linewidth = 0.8) +
  scale_x_date(date_breaks = "3 days", date_labels = "%m/%d") +
  labs(
    x = "Date (UTC)",
    y = "Mean 10 m wind speed (m/s)"
  ) +
  theme_bw()

#print(wind_daily_plot)

ggsave(
  filename = file.path(output_dir, "/noise/fig/wind-speed-time-series-daily.png"),
  plot = wind_daily_plot,
  width = 8,
  height = 4,
  dpi = 300
)

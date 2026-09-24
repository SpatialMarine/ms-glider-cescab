
# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# 02_diel_patterns_statistics.R

# Analyse diel variation in cetacean acoustic presence by assigning each
# 30-s FLAC file to daytime or nighttime based on local sunrise/sunset times,
# and compare file-level detection frequencies between periods.
# ============================================================

source("setup.R")
source("analysis/cetaceans/02_detection_statistics/fun/fun_diel_statistics.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
indir <- file.path(output_dir, "cetaceans")
outdir <- file.path(output_dir, "cetaceans")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------
# 2. Import file-level presence/absence table
# ------------------------------------------------------------
files_presAbs <- read.csv(
  file.path(indir, "files_presAbs.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    # CSV imports datetime columns as character.
    # These timestamps represent UTC, so parse them explicitly as UTC.
    fileStartUTC = lubridate::ymd_hms(fileStartUTC, tz = "UTC"),
    fileEndUTC   = lubridate::ymd_hms(fileEndUTC,   tz = "UTC")
  )

# ------------------------------------------------------------
# 3. Convert UTC timestamps to local Madrid time
# ------------------------------------------------------------
files_presAbs <- files_presAbs %>%
  mutate(
    fileStartLocal = lubridate::with_tz(
      fileStartUTC,
      tzone = "Europe/Madrid"
    ),
    hourLocal = lubridate::hour(fileStartLocal)
  )

# ------------------------------------------------------------
# 3.1. Check
# ------------------------------------------------------------
# Check datetime classes and time zones
class(files_presAbs$fileStartUTC)
attr(files_presAbs$fileStartUTC, "tzone")

class(files_presAbs$fileStartLocal)
attr(files_presAbs$fileStartLocal, "tzone")

head(
  files_presAbs %>%
    dplyr::select(file, fileStartUTC, fileStartLocal, hourLocal)
)

# ------------------------------------------------------------
# 4. Define study area centroid for sunrise/sunset calculations
# ------------------------------------------------------------
nc_file <- file.path(paste0(input_dir, "/glider_data/dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"))
track_nc <- nc_open(nc_file)

lon <- mean(ncvar_get(track_nc, "longitude"), na.rm = TRUE)
lat <- mean(ncvar_get(track_nc, "latitude"), na.rm = TRUE)

nc_close(track_nc)

# ------------------------------------------------------------
# 5. Calculate sunrise/sunset times for the campaign
# ------------------------------------------------------------
dayNight_times <- suncalc::getSunlightTimes(
  date = as.Date(files_presAbs$fileStartLocal),
  lat = lat,
  lon = lon,
  keep = c("sunrise", "sunset"),
  tz = "Europe/Madrid"
) %>%
  distinct(date, .keep_all = TRUE)

median_sunrise <- median(dayNight_times$sunrise)
median_sunset  <- median(dayNight_times$sunset)
median_sunrise
median_sunset
# Median sunrise is 07:06, so 07:00 is used as the sunrise hour. 
# Median sunset is 20:42, so 21:00 is used as the sunset hour.

day_hour <- 7
night_hour <- 21

# ------------------------------------------------------------
# 6. Add day/night period
# ------------------------------------------------------------
files_presAbs <- files_presAbs %>%
  mutate(
    period = ifelse(
      hourLocal >= day_hour & hourLocal < night_hour,
      "day",
      "night"
    )
  )

# ------------------------------------------------------------
# 7. Prepare binary detection variables
# ------------------------------------------------------------
files_presAbs <- files_presAbs %>%
  mutate(
    pres_Pm_num = as.integer(pres_Pm),
    pres_del_total = as.integer(pres_Del | pres_Gg),
    clicks_num = as.integer(del_clicks | Gg_clicks),
    whistles_num = as.integer(del_whistles | Gg_whistles)
  )

# ------------------------------------------------------------
# 8. Run diel statistics
# ------------------------------------------------------------
diel_stats <- bind_rows(
  run_diel_stats(files_presAbs, "pres_Pm_num", "Physeter macrocephalus"),
  run_diel_stats(files_presAbs, "pres_del_total", "Delphinidae + Grampus griseus"),
  run_diel_stats(files_presAbs, "clicks_num", "Clicks"),
  run_diel_stats(files_presAbs, "whistles_num", "Whistles")
)

diel_stats

# ------------------------------------------------------------------
# Summary of diel activity patterns
#
# - Physeter macrocephalus detections were significantly more frequent
#   during daytime (9.56%) than nighttime (5.55%)
#   (Wilcoxon rank-sum test, p < 0.001).
#
# - Delphinidae (including Grampus griseus due to the low number of
#   species-specific detections) were significantly more frequently
#   detected during nighttime (26.2%) than daytime (12.9%)
#   (Wilcoxon rank-sum test, p < 0.001).
#
# - Delphinid echolocation clicks were significantly more frequent
#   during nighttime (18.3%) than daytime (3.88%)
#   (Wilcoxon rank-sum test, p < 0.001).
#
# - Delphinid whistles were slightly more frequent during nighttime
#   (11.6%) than daytime (9.90%). 
#   (Wilcoxon rank-sum test, p = 0.032).
# ------------------------------------------------------------------

# ------------------------------------------------------------
# 9. Save diel statistics
# ------------------------------------------------------------
write.csv(diel_stats, file.path(outdir, "diel_patterns_statistics.csv"))


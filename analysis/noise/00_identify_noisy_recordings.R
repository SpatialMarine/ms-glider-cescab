
# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97
#
# ============================================================
# 00_identify_noisy_recordings.R
# ============================================================
#
# This script documents how recordings affected by glider self-noise or surfacing noise were identified for exclusion from the underwater-noise time series.
#
# In the raw acoustic dataset available from the associated Zenodo repository (10.5281/zenodo.20702261, folder "1_CESCAB-raw-flac"), the 432 affected FLAC files have already been flagged in their filenames. 
# Therefore, running this script is not required to reproduce the analyses presented in the manuscript. 
# It is provided solely for transparency and does not delete, move, rename, or modify any files.
#
# Two categories of noisy files were excluded:
#
#   1. Files overlapping battery or oil-pump movement, together with additional files containing manually identified glider noise.
#
#   2. Files overlapping shallow glider periods (<20 m depth), where surfacing and wave or wind action may affect the recordings.
#
# Some files contain both noise types. In total, 432 unique files contain at least one of these noise types and are excluded from the underwater-noise time series.

source("setup.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
nc_dir <- file.path(input_dir, "glider_data")

flac_dir <- file.path(input_dir, "1_CESCAB-raw-flac") # this folder in Zenodo repository already includes renamed files. 

outdir <- file.path(output_dir, "noise")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

L0_file <- file.path(
  nc_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L0_2024-08-01_data_dt.nc"
)

L1_file <- file.path(
  nc_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)


# ------------------------------------------------------------
# 3. Import file list
# ------------------------------------------------------------
flac_dt <- data.table(
  file = list.files(flac_dir, pattern = "\\.flac$", full.names = FALSE)
)

flac_dt[, filePath := file.path(flac_dir, file)]

flac_dt[, file_start := ymd_hms(
  str_remove(str_extract(file, "CESCAB_\\d{8}_\\d{6}"), "^CESCAB_"),
  tz = "UTC"
)]

flac_dt[, file_end := file_start + seconds(30)]

setorder(flac_dt, file_start)

# ============================================================
# Identify files with pump and battery self-noise
# These glider operations produce the most noise
# ============================================================

# ------------------------------------------------------------
# 4. Import L0 glider operations log
# ------------------------------------------------------------
L0 <- nc_open(L0_file) # use raw (L0) delayed-time product which contains glider operation logs

operations <- data.table(
  c_batt_pos             = ncvar_get(L0, "c_battpos"),
  m_is_battpos_moving   = ncvar_get(L0, "m_is_battpos_moving"),
  c_de_oil_vol          = ncvar_get(L0, "c_de_oil_vol"),
  m_is_de_pump_moving   = ncvar_get(L0, "m_is_de_pump_moving"),
  m_present_time        = as.POSIXct(
    ncvar_get(L0, "m_present_time"),
    origin = "1970-01-01",
    tz = "UTC"
  ),
  m_depth               = ncvar_get(L0, "m_depth")
)

nc_close(L0)

# ------------------------------------------------------------
# 5. Identify pump/battery movement timestamps
# ------------------------------------------------------------

# Only pitch-battery and oil-pump movement are used here because these
# were identified as the main contributors to glider self-noise.
motion <- operations[
  (m_is_battpos_moving == 1 | m_is_de_pump_moving == 1) &
    !is.na(m_present_time),
  .(
    batt = as.integer(m_is_battpos_moving == 1),
    pump = as.integer(m_is_de_pump_moving == 1)
  ),
  by = m_present_time
][
  ,
  .(
    batt = as.integer(any(batt == 1)),
    pump = as.integer(any(pump == 1))
  ),
  by = m_present_time
][
  order(m_present_time)
]

# ------------------------------------------------------------
# 6. Match movement timestamps to 30 s WAV files
# ------------------------------------------------------------
setkey(motion, m_present_time)

idx  <- motion[.(flac_dt$file_start), which = TRUE, roll = -Inf]
cand <- motion[idx]

pump_batt <- copy(flac_dt)[
  ,
  `:=`(
    m_present_time      = cand$m_present_time,
    m_is_battpos_moving = cand$batt,
    m_is_de_pump_moving = cand$pump
  )
][
  is.na(m_present_time) | m_present_time >= file_end,
  `:=`(
    m_present_time      = as.POSIXct(NA),
    m_is_battpos_moving = NA_integer_,
    m_is_de_pump_moving = NA_integer_
  )
][
  ,
  remove := !is.na(m_present_time)
]

# Also flag files outside the available glider metadata period.
first_timestamp <- min(operations$m_present_time, na.rm = TRUE)
last_timestamp  <- max(operations$m_present_time, na.rm = TRUE)

pump_batt[
  file_end <= first_timestamp | file_start >= last_timestamp,
  remove := TRUE
]


# Add manually identified glider-noise files
# These files showed unusually high sound levels in the time series generated by "analysis/noise/02_create_timeseries.m". 
# Their spectrograms were reviewed manually and confirmed to contain glider noise.
# They may correspond to glider operations not captured by the variables used above.

manual_dir <- file.path(input_dir, "noise/00_manual-glider-noise-det.csv")

if (file.exists(manual_dir)) {
  
  manual <- fread(manual_dir)
  
  manual[, file_start := dmy_hms(file_start, tz = "UTC")]
  
  pump_batt[
    file_start %in% manual[glider_noise == TRUE, file_start],
    remove := TRUE
  ]
}

# filter to files containing pumping and battery movement noise
pump_batt <- pump_batt[remove == TRUE]

# check
pump_batt

# fwrite(
#   pump_batt,
#   file.path(outdir, "glider-noise-files.csv")
# )

# ============================================================
# Surfacing / shallow-depth noise
# Files recorded under 20 m depth contain significant surfacing noise
# ============================================================

# ------------------------------------------------------------
# 7. Import L1 depth data
# ------------------------------------------------------------
L1 <- nc_open(L1_file) # use processed (L1) delayed-time product (most accurate)

depth_data <- data.table(
  time      = ncvar_get(L1, "time_ctd"),
  depth_ctd = ncvar_get(L1, "depth_ctd"), # Use CTD depth where available; fall back to depth (less precise) for missing depth_ctd at beginning of the mission
  depth     = ncvar_get(L1, "depth"),
  lon       = ncvar_get(L1, "longitude"),
  lat       = ncvar_get(L1, "latitude")
)

nc_close(L1)

depth_data[time == 0, time := NA_real_]

depth_data[, time := as.POSIXct(
  time,
  origin = "1970-01-01",
  tz = "UTC"
)]

depth_data <- depth_data[
  !is.na(time) &
    !(is.na(depth_ctd) & is.na(depth))
]

# ------------------------------------------------------------
# 8. Define shallow intervals
# ------------------------------------------------------------

threshold <- 20
min_overlap_seconds <- 1

depth_data[, depth_used := fifelse(!is.na(depth_ctd), depth_ctd, depth)]

setorder(depth_data, time)

depth_data[, shallow := depth_used < threshold]
depth_data[, run_id := rleid(shallow)]

shallow_intervals <- depth_data[
  shallow == TRUE,
  .(
    interval_start = min(time),
    interval_end   = max(time),
    run_min_depth  = min(depth_used, na.rm = TRUE)
  ),
  by = run_id
]

shallow_intervals[!is.finite(run_min_depth), run_min_depth := NA_real_]

# ------------------------------------------------------------
# 9. Match shallow intervals to FLAC files
# ------------------------------------------------------------

surfacing <- copy(flac_dt)

surfacing[, `:=`(
  file_start = file_start,
  file_end   = file_end
)]

shallow_intervals[, `:=`(
  int_start = interval_start,
  int_end   = interval_end
)]

setkey(surfacing, file_start, file_end)
setkey(shallow_intervals, int_start, int_end)

overlaps <- foverlaps(
  surfacing,
  shallow_intervals,
  by.x = c("file_start", "file_end"),
  by.y = c("int_start", "int_end"),
  nomatch = 0L
)

if (nrow(overlaps) > 0) {
  
  overlaps[, overlap_start := pmax(file_start, interval_start)]
  overlaps[, overlap_end   := pmin(file_end, interval_end)]
  
  overlaps[, overlap_sec := as.numeric(
    overlap_end - overlap_start,
    units = "secs"
  )]
  
  overlaps <- overlaps[overlap_sec > 0]
  
  surfacing <- overlaps[
    ,
    .(
      overlap_total = sum(overlap_sec, na.rm = TRUE),
      interval_start = interval_start[which.max(overlap_sec)],
      interval_end   = interval_end[which.max(overlap_sec)],
      run_id         = run_id[which.max(overlap_sec)],
      run_min_depth  = run_min_depth[which.max(overlap_sec)]
    ),
    by = .(filePath, file_start, file_end)
  ][
    overlap_total >= min_overlap_seconds
  ]
  
  surfacing[, remove_for_shallow := TRUE]
  
} else {
  
  surfacing <- data.table(
    filePath = character(),
    file_start = as.POSIXct(character(), tz = "UTC"),
    file_end = as.POSIXct(character(), tz = "UTC"),
    overlap_total = numeric(),
    interval_start = as.POSIXct(character(), tz = "UTC"),
    interval_end = as.POSIXct(character(), tz = "UTC"),
    run_id = integer(),
    run_min_depth = numeric(),
    remove_for_shallow = logical()
  )
}

# format dates for CSV
dt_cols <- c("file_start", "file_end", "interval_start", "interval_end")

surfacing[, (dt_cols) := lapply(
  .SD,
  function(x) format(x, "%Y-%m-%d %H:%M:%S")
), .SDcols = dt_cols]

# check
surfacing

# fwrite(
#   surfacing,
#   file.path(outdir, "surfacing-files.csv")
# )

# ------------------------------------------------------------
# 10. Summary
# ------------------------------------------------------------

# check if any files contain both glider operations and surfacing noise
common_files <- intersect(
  surfacing$filePath,
  pump_batt$filePath
) 

cat("Pump/battery files: ", nrow(pump_batt), "\n")
cat("Surfacing files:    ", nrow(surfacing), "\n")
cat("Common files:       ", length(common_files), "\n")
cat("Unique files:       ", length(unique(c(
  pump_batt$filePath,
  surfacing$filePath
))), "\n")

# Unique files:        432

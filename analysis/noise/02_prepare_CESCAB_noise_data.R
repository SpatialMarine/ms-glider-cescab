# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# 01_prepare_CESCAB_noise_data.R

# Read the MATLAB third-octave SPL time series, retain the 63 Hz, 125 Hz,
# and 8 kHz bands, and remove files marked as surfacing or glider noise
# by their FLAC filenames. Save cleaned SPL samples and an exclusion audit.
# ============================================================

source("setup.R")

# ------------------------------------------------------------
# 1. Package and path checks
# ------------------------------------------------------------
if (!requireNamespace("R.matlab", quietly = TRUE)) {
  stop("Install R.matlab first: install.packages('R.matlab')")
}

mat_file <- file.path(output_dir, "noise", "timeseries", "CESCAB_timeseries_to.mat")
flac_dir <- file.path(input_dir, "1_CESCAB-raw-flac")

outdir <- file.path(output_dir, "noise", "timeseries")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Read recording times and noise flags from FLAC filenames
# ------------------------------------------------------------
flac_names <- list.files(
  flac_dir,
  pattern = "\\.flac$",
  recursive = TRUE,
  ignore.case = TRUE
)

timestamp_text <- sub(
  "^CESCAB_([0-9]{8}_[0-9]{6}).*$",
  "\\1",
  basename(flac_names),
  ignore.case = TRUE
)

files <- data.frame(
  FileStart = as.POSIXct(
    timestamp_text,
    format = "%Y%m%d_%H%M%S",
    tz = "UTC"
  ),
  FileName = basename(flac_names),
  # Recordings containing surfacing or glider self-noise are flagged in their filenames.
  # A recording is retained only if its complete filename follows:
  # CESCAB_yyyymmdd_HHMMSS.flac
  #
  # Files containing an additional suffix are excluded, for example:
  #   CESCAB_20240801_121022_surface.flac
  #   CESCAB_20240801_124022_glider_noise.flac
  #   CESCAB_20240801_145522_surface_glider_noise.flac
  Excluded = grepl(
    "_(surface|glider_noise|surface_glider_noise)\\.flac$",
    flac_names,
    ignore.case = TRUE
  )
)

files <- files[order(files$FileStart), ]
row.names(files) <- NULL

message("FLAC files found: ", nrow(files)) # 6147
message("Marked for exclusion: ", sum(files$Excluded)) # 432
message("Expected to retain: ", sum(!files$Excluded)) # 5715


# ------------------------------------------------------------
# 3. Read the MATLAB third-octave time series
# ------------------------------------------------------------
# R.matlab::readMat() replaces underscores in MATLAB variable names with dots:
#   spectro_to -> spectro.to, freq_to -> freq.to, time_axis -> time.axis.
mat <- R.matlab::readMat(mat_file)

spl <- mat$spectro.to #  Sound Pressure Level (SPL) values: frequency bands × measurement times
freq <- as.numeric(mat$freq.to) # centre frequency (Hz) of each row in "spl"
time_axis <- as.numeric(mat$time.axis) # # MATLAB timestamp of each column in "spl"

# MATLAB datenums are days since MATLAB's origin; 719529 corresponds to
# 1970-01-01. The MATLAB producer already added each recording's start_date.
# Consequently no extra anchor-time adjustment is applied.
sample_time <- as.POSIXct((time_axis - 719529) * 86400,
                          origin = "1970-01-01", tz = "UTC")

message("Sample timerange: ", length(sample_time), " (",
        format(sample_time[1], "%Y-%m-%d %H:%M:%OS3", tz = "UTC"),
        " to ", format(tail(sample_time, 1), "%Y-%m-%d %H:%M:%OS3", tz = "UTC"), " UTC)")

# ------------------------------------------------------------
# 4. Select the three third-octave band centres
# ------------------------------------------------------------
target_hz <- c(Hz63 = 63, Hz125 = 125, Hz8k = 8000)
band_index <- vapply(target_hz, function(x) which.min(abs(freq - x)), integer(1))

# ------------------------------------------------------------
# 5. Match SPL samples to their 30-second FLAC recordings
# ------------------------------------------------------------
file_index <- findInterval(
  as.numeric(sample_time),
  as.numeric(files$FileStart)
)

in_file <- file_index > 0L
in_file[in_file] <- (
  as.numeric(sample_time[in_file]) <
    as.numeric(files$FileStart[file_index[in_file]]) + 30
)
stopifnot(all(in_file))

samples_per_file <- tabulate(file_index, nbins = nrow(files))
message(
  "SPL samples per FLAC file: ",
  min(samples_per_file), "–", max(samples_per_file)
)

#SPL samples per FLAC file: 58–59

# ------------------------------------------------------------
# 6. Keep only samples from unflagged FLAC files
# ------------------------------------------------------------
keep <- !files$Excluded[file_index]
clean_samples <- data.frame(
  Time = sample_time[keep],
  FileStart = files$FileStart[file_index[keep]],
  Hz63 = as.numeric(spl[band_index["Hz63"], keep]),
  Hz125 = as.numeric(spl[band_index["Hz125"], keep]),
  Hz8k = as.numeric(spl[band_index["Hz8k"], keep])
)
message("Clean SPL samples retained: ", nrow(clean_samples),
        " from ", sum(!files$Excluded), " recordings")

# ------------------------------------------------------------
# 7. Save the intermediate data and the exclusion audit
# ------------------------------------------------------------
saveRDS(
  clean_samples,
  file.path(outdir, "CESCAB_SPL_63_125_8k.rds"),
  compress = "xz"
)

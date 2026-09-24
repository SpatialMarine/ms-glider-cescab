
# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# 01_detection_percentages.R

# Calculate file-level cetacean acoustic presence and detection percentages 
# from the final annotation table. 
# Create and export file-level presence/absence and summary tables.
# ============================================================

source("setup.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
indir <- file.path(output_dir, "cetaceans/final-annotation-table")

flac_dir <- file.path(input_dir, "1_CESCAB-raw-flac")

outdir <- file.path(output_dir, "cetaceans", "detection-stats")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------
# 2. Import final detection table
# ------------------------------------------------------------
annotations <- read.delim(
  file.path(indir, "2_CESCAB-annotation-table.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 3. Create complete screened FLAC file list
# ------------------------------------------------------------
flacs <- tibble(
  file = list.files(flac_dir, pattern = "\\.flac$", full.names = FALSE)
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
  distinct(file, .keep_all = TRUE) %>%
  arrange(fileStartUTC)

total_files <- nrow(flacs)
total_hours <- total_files * 30 / 3600

# ------------------------------------------------------------
# 4. Prepare detection intervals
# ------------------------------------------------------------
annotations <- annotations %>%
  mutate(
    beginDateTime = parse_date_time(
      `Begin Date Time`,
      orders = c("ymd HMS", "ymd IMS"),
      tz = "UTC"
    ),
    endDateTime = parse_date_time(
      paste(`End Date`, `End Clock Time`),
      orders = c("ymd HMS", "ymd IMS"),
      tz = "UTC"
    )
  )

detection_intervals <- annotations %>%
  transmute(
    species = Species,
    signalType = `Signal Type`,
    detStartUTC = beginDateTime,
    detEndUTC = endDateTime
  ) %>%
  filter(
    !is.na(detStartUTC),
    !is.na(detEndUTC),
    detEndUTC >= detStartUTC
  )

# ------------------------------------------------------------
# 5. Create file-level presence/absence table
# ------------------------------------------------------------
files_presAbs <- flacs %>%
  left_join(
    detection_intervals,
    join_by(
      fileStartUTC < detEndUTC,
      fileEndUTC > detStartUTC
    )
  ) %>%
  group_by(file, fileStartUTC, fileEndUTC) %>%
  summarise(
    presence = any(!is.na(species)),
    pres_Pm = any(species == "Physeter macrocephalus", na.rm = TRUE),
    pres_Del = any(species == "Delphinidae", na.rm = TRUE),
    pres_Gg = any(species == "Grampus griseus", na.rm = TRUE),
    del_whistles = any(species == "Delphinidae" & signalType == "whistles", na.rm = TRUE),
    del_clicks = any(species == "Delphinidae" & signalType == "clicks", na.rm = TRUE),
    Gg_whistles = any(species == "Grampus griseus" & signalType == "whistles", na.rm = TRUE),
    Gg_clicks = any(species == "Grampus griseus" & signalType == "clicks", na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Summarise detection percentages
# ------------------------------------------------------------
detection_files <- tibble(
  Category = c(
    "Any cetacean presence",
    "Physeter macrocephalus",
    "Delphinidae",
    "Grampus griseus",
    "Delphinidae whistles",
    "Delphinidae clicks",
    "Grampus griseus whistles",
    "Grampus griseus clicks"
  ),
  Detected_files = c(
    nrow(filter(files_presAbs, presence)),
    nrow(filter(files_presAbs, pres_Pm)),
    nrow(filter(files_presAbs, pres_Del)),
    nrow(filter(files_presAbs, pres_Gg)),
    nrow(filter(files_presAbs, del_whistles)),
    nrow(filter(files_presAbs, del_clicks)),
    nrow(filter(files_presAbs, Gg_whistles)),
    nrow(filter(files_presAbs, Gg_clicks))
  )
) %>%
  mutate(
    Percentage = 100 * Detected_files / total_files
  )

detection_files

# summary                Detected_files Percentage
#   1 Any cetacean presence              1570     25.5  
# 2 Physeter macrocephalus              485      7.89 
# 3 Delphinidae                        1098     17.9  
# 4 Grampus griseus                      36      0.586
# 5 Delphinidae whistles                652     10.6  
# 6 Delphinidae clicks                  573      9.32 
# 7 Grampus griseus whistles              0      0    
# 8 Grampus griseus clicks               36      0.586

# ------------------------------------------------------------
# 7. Save outputs needed by later scripts
# ------------------------------------------------------------
write.csv(files_presAbs, file.path(outdir, "file_presence_absence.csv"))
write.csv(detection_files, file.path(outdir, "detection_percentage_summary.csv"), row.names = FALSE)

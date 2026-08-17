
# ==============================================================================
# Ocean gliders for PAM
# ==============================================================================

# Greta Jankauskaite | @gretajan97

# ==============================================================================
# 00_inter-annotator-comparison.R

# Inter-annotator agreement: 
# file-level presence/absence comparison between two full-dataset annotators
# ==============================================================================

source("setup.R")
source("analysis/cetaceans/01_inter-annotator/fun/fun_inter-annotator.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
indir <- file.path(input_dir, "CESCAB-annotation-tables")
flac_dir <- file.path(input_dir, "1_CESCAB-raw-flac")

# ------------------------------------------------------------
# 2. Import Raven selection tables
# ------------------------------------------------------------
annotations_a1 <- read_sel("annotator_1_selection_table.txt", "A1")
annotations_a2 <- read_sel("annotator_2_selection_table.txt", "A2")

# -------------------------------------------------------------------
# 3. Reclassify Risso's dolphin as Delphinidae for agreement analysis
# -------------------------------------------------------------------
annotations_a1 <- annotations_a1 %>%
  mutate(
    Species = if_else(
      Species == "Grampus griseus",
      "Delphinidae",
      Species
    )
  )

# ------------------------------------------------------------
# 4. Create complete screened FLAC file list
# ------------------------------------------------------------
flacs <- tibble(
  file = list.files(flac_dir, pattern = "\\.flac$", full.names = FALSE)
) %>%
  mutate(
    core_id = str_extract(file, "CESCAB_\\d{8}_\\d{6}"),
    timestamp_str = str_remove(core_id, "^CESCAB_"),
    fileStartUTC = ymd_hms(timestamp_str, tz = "UTC"),
    fileEndUTC = fileStartUTC + seconds(30),
    fileStartLocal = with_tz(fileStartUTC, "Europe/Madrid"),
    Date = as.Date(fileStartLocal)
  ) %>%
  distinct(file, .keep_all = TRUE) %>%
  arrange(fileStartUTC) %>%
  mutate(file_index = row_number())

file_lookup <- flacs %>%
  dplyr::select(file, file_index)

# ------------------------------------------------------------
# 5. Expand selections to all FLAC files they span
# ------------------------------------------------------------
det_a1 <- get_files(annotations_a1)
det_a2 <- get_files(annotations_a2)

# ------------------------------------------------------------
# 6. Create file-level presence/absence matrix
# ------------------------------------------------------------
pres <- bind_rows(det_a1, det_a2) %>%
  make_presence(c("A1", "A2"))

# ------------------------------------------------------------
# 7. Calculate Cohen's kappa
# ------------------------------------------------------------
cats <- c(
  "all",
  "sperm_whale",
  "delphinid",
  "delphinid_clicks",
  "delphinid_whistles"
)

kappa <- map_dfr(
  cats,
  ~ calc_kappa(pres, "A1", "A2", .x)
)

kappa

# ------------------------------------------------------------
# Inter-annotator agreement summary
# ------------------------------------------------------------
# Across all cetacean detections (delphinid echolocation clicks and whistles, and sperm whale echolocation clicks), 
# annotators agreed on 97.3% of files, with a Cohen's kappa of 0.909, indicating "almost perfect" agreement following McHugh (2012).

# Reference: McHugh, M. L. Interrater reliability: the kappa statistic. Biochemia Medica 22, 276–282 (2012).


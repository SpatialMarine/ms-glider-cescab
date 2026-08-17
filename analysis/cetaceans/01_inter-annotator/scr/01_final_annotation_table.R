
# ============================================================
# Ocean gliders for PAM
# ============================================================

# Greta Jankauskaite | @gretajan97

# ============================================================
# 01_final_annotation_table.R

# Create final CESCAB selection table combined across annotators 
# ============================================================

source("setup.R")
source("analysis/cetaceans/01_inter-annotator/fun/fun_inter-annotator.R")

# ------------------------------------------------------------
# 1. Define paths
# ------------------------------------------------------------
indir <- file.path(input_dir, "CESCAB-annotation-tables")
outdir <- file.path(output_dir, "cetaceans/final-annotation-table")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

flac_dir <- file.path(input_dir, "1_CESCAB-raw-flac")

# ------------------------------------------------------------
# 2. Import Raven selection tables
# ------------------------------------------------------------
annotations_a1 <- read_table("annotator_1_table.txt", "A1")
annotations_a2 <- read_table("annotator_2_table.txt", "A2")
annotations_a3 <- read_table("annotator_3_table.txt", "A3")

# ------------------------------------------------------------
# 3. Create complete screened FLAC file list
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
# 4. Expand annotations to all FLAC files they span
# ------------------------------------------------------------
det_a1 <- get_files(annotations_a1)
det_a2 <- get_files(annotations_a2)
det_a3 <- get_files(annotations_a3)

# ------------------------------------------------------------
# 5. Identify files already covered by A1
# ------------------------------------------------------------
covered_files <- det_a1 %>%
  distinct(file)

# ------------------------------------------------------------
# 6. Add A2 annotations only if they contain at least one new file not in A1 table
# ------------------------------------------------------------
add_a2_ids <- det_a2 %>%
  anti_join(covered_files, by = "file") %>%
  distinct(row_id)

add_a2 <- annotations_a2 %>%
  semi_join(add_a2_ids, by = "row_id")

# Update covered files after adding A2
covered_files <- bind_rows(
  covered_files,
  det_a2 %>% semi_join(add_a2_ids, by = "row_id") %>% distinct(file)
) %>%
  distinct(file)

# ------------------------------------------------------------
# 7. Add A3 selections only if they contain at least one file not covered by A1 or added A2
# ------------------------------------------------------------
add_a3_ids <- det_a3 %>%
  anti_join(covered_files, by = "file") %>%
  distinct(row_id)

add_a3 <- annotations_a3 %>%
  semi_join(add_a3_ids, by = "row_id")

# ------------------------------------------------------------
# 8. Combine selections and renumber chronologically
# ------------------------------------------------------------
final_annotations <- bind_rows(
  annotations_a1,
  add_a2,
  add_a3
) %>%
  arrange(`Begin Time (s)`) %>%
  mutate(Selection = row_number()) %>%
  dplyr::select(
    -row_id,
    -any_of(c("View", "Channel"))
  )
# 1626 annotations.

# ------------------------------------------------------------
# 9. Export final Zenodo table
# ------------------------------------------------------------
write.table(final_annotations, file = file.path(outdir, "2_CESCAB-annotation-table.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


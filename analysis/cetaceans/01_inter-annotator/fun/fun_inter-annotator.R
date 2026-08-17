
#-----------------------------------------------------------------------------------------
# fun_inter-annotator.R          Common functions for Raven annotation and agreement analyses
#-----------------------------------------------------------------------------------------
# read_sel                  Import Raven selection tables and add annotator ID
# get_files                 Expand selections to all FLAC files they span (Raven selection tables store only the start and end files, omitting any files in between)
# make_presence             Create file-level presence/absence matrix
# calc_kappa                Calculate file-level agreement statistics and Cohen's kappa 
#-----------------------------------------------------------------------------------------

# ------------------------------------------------------------
# Import Raven selection table
# ------------------------------------------------------------
read_table <- function(file, id) {
  
  read.delim(
    file.path(indir, file),
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    mutate(
      `Begin Date Time` = format(
        parse_date_time(
          `Begin Date Time`,
          orders = "ymd HMS",
          tz = "UTC"
        ),
        "%Y-%m-%d %H:%M:%OS4"
      ),
      
      endDateTime = parse_date_time(
        paste(`End Date`, `End Clock Time`),
        orders = "ymd HMS",
        tz = "UTC"
      ),
      
      `End Date` = format(
        endDateTime,
        "%Y-%m-%d"
      ),
      
      `End Clock Time` = format(
        endDateTime,
        "%H:%M:%OS4"
      ),
      
      annotator = id,
      row_id = row_number()
    ) %>%
    select(-endDateTime)
}

# ------------------------------------------------------------
# Expand selections to all FLAC files they span
# ------------------------------------------------------------
get_files <- function(sel) {
  sel %>%
    left_join(file_lookup, by = c("Begin File" = "file")) %>%
    rename(i1 = file_index) %>%
    left_join(file_lookup, by = c("End File" = "file")) %>%
    rename(i2 = file_index) %>%
    filter(!is.na(i1), !is.na(i2)) %>%
    mutate(
      file_index = map2(
        pmin(i1, i2),
        pmax(i1, i2),
        ~ seq(.x, .y)
      )
    ) %>%
    unnest(file_index) %>%
    left_join(file_lookup, by = "file_index") %>%
    select(annotator, row_id, file, Species, `Signal Type`)
}

# ------------------------------------------------------------
# Create file-level presence/absence matrix
# ------------------------------------------------------------
make_presence <- function(det, annotators) {
  
  cats <- c(
    "all",
    "sperm_whale",
    "delphinid",
    "delphinid_clicks",
    "delphinid_whistles"
  )
  
  det %>%
    mutate(
      all = 1L,
      sperm_whale = as.integer(Species == "Physeter macrocephalus"),
      delphinid = as.integer(Species == "Delphinidae"),
      delphinid_clicks = as.integer(
        Species == "Delphinidae" &
          `Signal Type` == "clicks"
      ),
      delphinid_whistles = as.integer(
        Species == "Delphinidae" &
          `Signal Type` == "whistles"
      )
    ) %>%
    group_by(annotator, file) %>%
    summarise(
      across(all_of(cats), ~ as.integer(any(.x == 1))),
      .groups = "drop"
    ) %>%
    right_join(
      expand_grid(
        annotator = annotators,
        file = flacs$file
      ),
      by = c("annotator", "file")
    ) %>%
    mutate(
      across(all_of(cats), ~ replace_na(.x, 0L))
    )
}

# ------------------------------------------------------------
# Calculate Cohen's kappa
# ------------------------------------------------------------
calc_kappa <- function(pres, a, b, cat) {
  
  x <- pres %>%
    filter(annotator %in% c(a, b)) %>%
    select(annotator, file, value = all_of(cat)) %>%
    pivot_wider(
      names_from = annotator,
      values_from = value,
      values_fill = 0
    ) %>%
    arrange(file)
  
  v1 <- x[[a]]
  v2 <- x[[b]]
  
  tibble(
    comparison = paste(a, "vs", b),
    category = cat,
    n_files = nrow(x),
    both_present = sum(v1 == 1 & v2 == 1),
    both_absent = sum(v1 == 0 & v2 == 0),
    disagreement = sum(v1 != v2),
    agreement_percent = round(mean(v1 == v2) * 100, 1),
    cohens_kappa = round(kappa2(data.frame(v1, v2))$value, 3),
    A1_present = sum(v1 == 1),
    A1_absent = sum(v1 == 0),
    A2_present = sum(v2 == 1),
    A2_absent = sum(v2 == 0)
  )
}


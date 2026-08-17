
#-----------------------------------------------------------------------------------------
# fun_diel_detection_statistics.R
# Common functions for diel activity analyses
#-----------------------------------------------------------------------------------------

# run_diel_stats
# Compare day and night detection proportions using a chi-squared
# test of independence and summarise detection percentages
#-----------------------------------------------------------------------------------------

run_diel_stats <- function(data, response_var, label) {
  
  # Summarise detections by period
  summary <- data %>%
    group_by(period) %>%
    summarise(
      n_files = n(),
      detected_files = sum(.data[[response_var]], na.rm = TRUE),
      detection_percent = 100 * mean(.data[[response_var]], na.rm = TRUE),
      .groups = "drop"
    )
  
  day <- summary %>%
    filter(period == "day")
  
  night <- summary %>%
    filter(period == "night")
  
  # Construct 2 x 2 contingency table:
  # rows = day/night
  # columns = detected/not detected
  contingency_table <- matrix(
    c(
      day$detected_files,
      day$n_files - day$detected_files,
      night$detected_files,
      night$n_files - night$detected_files
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      period = c("day", "night"),
      detection = c("detected", "not_detected")
    )
  )
  
  # Pearson's chi-squared test of independence
  test <- chisq.test(
    contingency_table,
    correct = FALSE
  )
  
  tibble(
    category = label,
    day_detected = day$detected_files,
    day_total = day$n_files,
    percent_day = day$detection_percent,
    night_detected = night$detected_files,
    night_total = night$n_files,
    percent_night = night$detection_percent,
    p_value = test$p.value,
    statistic = unname(test$statistic)
  )
}
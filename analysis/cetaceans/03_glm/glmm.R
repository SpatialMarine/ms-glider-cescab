
# ============================================================
# Ocean gliders for PAM
# Greta Jankauskaite | @gretajan97
# glmm.R
# File-level binomial GLMMs 
# ============================================================

source("setup.R")

library(diveMove)
library(glmmTMB)
library(suncalc)
library(patchwork)

# 1. Paths ----------------------------------------------------

indir <- file.path(output_dir, "cetaceans", "detection-stats")
nc_dir <- file.path(input_dir, "glider_data")
outdir <- file.path(output_dir, "cetaceans", "glm-clean")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Supporting files stay separate from manuscript-ready outputs.
scratch_dir <- file.path(outdir, "scratch")
dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)


mpa_file <- file.path(
  carto_dir, "MPA", "Limite_Archipielago_Cabrera.shp"
)

bathy_file <- file.path(
  carto_dir, "EMODnet-bathymetry_2024",
  "bathymetry_dtm_2024_bb7b_a70e_b1f5.nc"
)

nc_file <- file.path(
  nc_dir,
  "dep0008_sl3u1064_scb-sl3u1064_L1_2024-08-01_data_dt.nc"
)


# 2. Read the presence/absence table ---------------------------

files_presAbs <- read.csv(
  file.path(indir, "file_presence_absence.csv"),
  stringsAsFactors = FALSE
) %>%
  select(-any_of("X")) %>%
  mutate(
    fileStartUTC = ymd_hms(fileStartUTC, tz = "UTC"),
    fileEndUTC = ymd_hms(fileEndUTC, tz = "UTC")
  ) %>%
  arrange(fileStartUTC)

# Include Grampus in Delphinidae and calculate recording midpoints.
files_presAbs <- files_presAbs %>%
  mutate(
    del_whistles_all = del_whistles | Gg_whistles,
    del_clicks_all = del_clicks | Gg_clicks,
    fileMidUTC = fileStartUTC + (fileEndUTC - fileStartUTC) / 2,
    fileMidLocal = with_tz(fileMidUTC, "Europe/Madrid"),
    hourLocal = hour(fileMidLocal),
    dateLocal = factor(
      as.Date(fileMidLocal, tz = "Europe/Madrid")
    )
  )


# 3. Read the glider track -------------------------------------

track_nc <- nc_open(nc_file)

track_data <- tibble(
  longitude = as.numeric(ncvar_get(track_nc, "longitude")),
  latitude = as.numeric(ncvar_get(track_nc, "latitude")),
  time = as.numeric(ncvar_get(track_nc, "time")),
  depth = as.numeric(ncvar_get(track_nc, "depth_ctd"))
)

nc_close(track_nc)

track_data <- track_data %>%
  mutate(
    time = as.POSIXct(time, origin = "1970-01-01", tz = "UTC"),
    depth = abs(depth)
  ) %>%
  mutate(
    longitude = if_else(is.finite(longitude) & abs(longitude) <= 180,
                        longitude, NA_real_),
    latitude = if_else(is.finite(latitude) & abs(latitude) <= 90,
                       latitude, NA_real_)
  ) %>%
  filter(!is.na(time)) %>%
  distinct(time, .keep_all = TRUE) %>%
  arrange(time)


# Calculate sunrise and sunset at the glider track centroid.
positions <- track_data %>%
  filter(is.finite(longitude), is.finite(latitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)


centroid <- positions %>%
  st_transform(32631) %>%
  st_union() %>%
  st_centroid() %>%
  st_transform(4326) %>%
  st_coordinates()

sun_times <- getSunlightTimes(
  date = sort(unique(
    as.Date(files_presAbs$fileMidLocal, tz = "Europe/Madrid")
  )),
  lat = centroid[1, "Y"],
  lon = centroid[1, "X"],
  keep = c("sunrise", "sunset"),
  tz = "Europe/Madrid"
)

# Round deployment-median sunrise and sunset to the nearest hour.
local_hour <- function(x) {
  hour(x) + minute(x) / 60 + second(x) / 3600
}

day_start <- round(median(local_hour(sun_times$sunrise)))
night_start <- round(median(local_hour(sun_times$sunset)))

# save for figures
write.csv(
  data.frame(
    day_start = day_start,
    night_start = night_start
  ),
  file.path(scratch_dir, "diel_boundaries.csv"),
  row.names = FALSE
)

# assign each recording to day or night using its local midpoint hour.
files_presAbs <- files_presAbs %>%
  mutate(
    diel = factor(
      if_else(
        hourLocal >= day_start & hourLocal < night_start,
        "day",
        "night"
      ),
      levels = c("day", "night")
    )
  )


# 4. Match measurements to each recording midpoint -------------

# Match each variable separately to its nearest valid measurement.
# Return NA if no valid measurement is within 15 seconds.

max_distance_s <- 15

nearest_value <- function(time, value, target,
                          max_distance = max_distance_s) {
  
  ok <- is.finite(as.numeric(time)) & is.finite(value)
  
  x <- as.numeric(time[ok])
  y <- value[ok]
  q <- as.numeric(target)
  
  ans <- rep(NA_real_, length(q))
  
  if (!length(x)) return(ans)
  
  valid <- which(is.finite(q))
  i <- findInterval(q[valid], x)
  
  before <- pmax(i, 1L)
  after <- pmin(i + 1L, length(x))
  
  nearest <- ifelse(
    abs(q[valid] - x[before]) <= abs(q[valid] - x[after]),
    before,
    after
  )
  
  keep <- abs(q[valid] - x[nearest]) <= max_distance
  
  ans[valid[keep]] <- y[nearest[keep]]
  
  ans
}

model_data <- files_presAbs

for (v in c(
  "longitude", "latitude", "depth"
)) {
  
  model_data[[v]] <- nearest_value(
    track_data$time,
    track_data[[v]],
    model_data$fileMidUTC
  )
}

# Recover missing position/depth matches within 180 seconds for any file,
# only within valid temporal coverage. Never borrow measurements hours away.

gap_tolerance_s <- 180



for (v in c("longitude", "latitude", "depth")) {
  
  valid <- is.finite(track_data[[v]])
  source_time <- track_data$time[valid]
  if (!length(source_time)) next
  
  rows <- which(
    is.na(model_data[[v]]) &
      model_data$fileMidUTC >= min(source_time) &
      model_data$fileMidUTC <= max(source_time)
  )
  
  model_data[[v]][rows] <- nearest_value(
    source_time,
    track_data[[v]][valid],
    model_data$fileMidUTC[rows],
    max_distance = gap_tolerance_s
  )
}


# 5. MPA status and bathymetry at the glider position ------------

cabrera_MPA <- st_read(mpa_file, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(4326)

bathy <- rast(bathy_file)

valid_position <- which(
  is.finite(model_data$longitude) &
    is.finite(model_data$latitude) &
    abs(model_data$longitude) <= 180 &
    abs(model_data$latitude) <= 90
)

model_data$MPA <- NA_character_
model_data$bathy <- NA_real_

if (length(valid_position)) {
  
  points <- st_as_sf(
    model_data[valid_position, ],
    coords = c("longitude", "latitude"),
    crs = 4326
  )
  
  model_data$MPA[valid_position] <- ifelse(
    lengths(st_intersects(points, cabrera_MPA)) > 0,
    "inside",
    "outside"
  )
  
  points_bathy <- project(vect(points), crs(bathy))
  
  elevation <- terra::extract(
    bathy, points_bathy, ID = FALSE
  )[[1]]
  
  # Convert negative seabed elevation to positive water depth.
  model_data$bathy[valid_position] <- ifelse(
    elevation <= 0, -elevation, NA_real_
  )
}

model_data$MPA <- factor(
  model_data$MPA,
  levels = c("outside", "inside")
)


# 6. Identify dives using diveMove ------------------------------

# Identify dives deeper than 20 m, including descent and ascent.
# CTD depths are assumed to be calibrated; apply zero offset.
# Process separate track segments to avoid identifying dives across gaps.

max_gap_s <- 120
dive_threshold_m <- 20

dive_track <- track_data %>%
  filter(is.finite(depth)) %>%
  mutate(
    track_segment = cumsum(
      c(TRUE, diff(as.numeric(time)) > max_gap_s)
    )
  )


dive_track$dive_id <- NA_integer_
dive_track$partial_dive <- FALSE

offset_id <- 0L

for (seg in unique(dive_track$track_segment)) {
  
  rows <- which(dive_track$track_segment == seg)
  d <- dive_track[rows, ]
  
  # Require both surface and diving measurements.
  if (nrow(d) < 4 ||
      !any(d$depth <= dive_threshold_m) ||
      !any(d$depth > dive_threshold_m)) next
  
  # Set shallow readings to zero for dive identification only.
  tdr <- diveMove::createTDR(
    time = d$time,
    depth = ifelse(d$depth <= dive_threshold_m, 0, d$depth),
    dtime = median(diff(as.numeric(d$time))),
    file = paste0("CESCAB_segment_", seg)
  )
  
  calibrated <- diveMove::calibrateDepth(
    tdr,
    dive.thr = dive_threshold_m,
    zoc.method = "offset",
    offset = 0,
    wet.thr = 0,
    interp.wet = FALSE,
    dive.model = "smooth.spline"
  )
  
  ids <- as.integer(as.character(
    diveMove::getDAct(calibrated, "dive.id")
  ))
  
  ids[ids <= 0L] <- NA_integer_
  
  partial <- rep(FALSE, length(ids))
  
  # Flag dives cut off at either end of the segment.
  if (!is.na(ids[1])) {
    partial[which(ids == ids[1])] <- TRUE
  }
  
  if (!is.na(ids[length(ids)])) {
    partial[which(ids == ids[length(ids)])] <- TRUE
  }
  
  if (any(!is.na(ids))) {
    ids[!is.na(ids)] <- ids[!is.na(ids)] + offset_id
    offset_id <- max(ids, na.rm = TRUE)
  }
  
  dive_track$dive_id[rows] <- ids
  dive_track$partial_dive[rows] <- partial
}

# Total dives identified.
n_distinct(dive_track$dive_id, na.rm = TRUE)


# Assign shallow intervals to the surrounding dive cycles.
# Between dives, split at the last shallowest measurement.
# Assign shallow segment edges to their adjacent identified dive.

for (seg in unique(dive_track$track_segment)) {
  
  rows <- which(dive_track$track_segment == seg)
  d <- dive_track[rows, ]
  
  missing_id <- is.na(d$dive_id)
  
  runs <- split(
    seq_len(nrow(d)),
    cumsum(c(TRUE, diff(missing_id) != 0))
  )
  
  for (idx in runs) {
    
    if (!all(missing_id[idx])) next
    if (!all(d$depth[idx] <= dive_threshold_m)) next
    
    previous_row <- min(idx) - 1L
    next_row <- max(idx) + 1L
    
    previous_id <- if (previous_row >= 1L) {
      d$dive_id[previous_row]
    } else {
      NA_integer_
    }
    
    next_id <- if (next_row <= nrow(d)) {
      d$dive_id[next_row]
    } else {
      NA_integer_
    }
    
    if (is.na(previous_id) && is.na(next_id)) next
    
    # Leading shallow interval: assign to the next dive.
    if (is.na(previous_id)) {
      d$dive_id[idx] <- next_id
      d$partial_dive[idx] <- d$partial_dive[next_row]
      next
    }
    
    # Trailing shallow interval: assign to the previous dive.
    if (is.na(next_id)) {
      d$dive_id[idx] <- previous_id
      d$partial_dive[idx] <- d$partial_dive[previous_row]
      next
    }
    
    # Between dives: split at the last shallowest measurement.
    shallowest <- min(d$depth[idx])
    turn <- max(idx[d$depth[idx] == shallowest])
    
    d$dive_id[idx[idx <= turn]] <- previous_id
    d$dive_id[idx[idx > turn]] <- next_id
    
    d$partial_dive[idx[idx <= turn]] <-
      d$partial_dive[previous_row]
    
    d$partial_dive[idx[idx > turn]] <-
      d$partial_dive[next_row]
  }
  
  dive_track$dive_id[rows] <- d$dive_id
  dive_track$partial_dive[rows] <- d$partial_dive
}


# Assign dive information from the nearest valid depth measurement.
# Use the same 180-second maximum and temporal coverage restriction.

x <- as.numeric(dive_track$time)
q <- as.numeric(model_data$fileMidUTC)

valid <- which(is.finite(q))
i <- findInterval(q[valid], x)

before <- pmax(i, 1L)
after <- pmin(i + 1L, length(x))

nearest <- ifelse(
  abs(q[valid] - x[before]) <= abs(q[valid] - x[after]),
  before,
  after
)

tolerance <- rep(gap_tolerance_s, length(valid))

keep <- abs(q[valid] - x[nearest]) <= tolerance &
  q[valid] >= min(x) &
  q[valid] <= max(x) &
  is.finite(model_data$depth[valid])

file_rows <- valid[keep]
track_rows <- nearest[keep]

model_data$dive_id <- NA_integer_
model_data$partial_dive <- NA

model_data$dive_id[file_rows] <- dive_track$dive_id[track_rows]

model_data$partial_dive[file_rows] <-
  dive_track$partial_dive[track_rows]

# Preserve the full dive sequence.
model_data$dive_id <- factor(
  model_data$dive_id,
  levels = seq_len(offset_id)
)

# Group the AR(1) correlation structure by day of year.
model_data$julian_day <- factor(
  yday(with_tz(model_data$fileMidUTC, "Europe/Madrid"))
)

model_data$AR_group <- model_data$julian_day


# Identify recordings with all required model variables.
eligible_no_env <- complete.cases(model_data[, c(
  "depth", "bathy", "MPA", "diel", "dive_id", "AR_group"
)])

# Number of recordings retained.
sum(eligible_no_env)

# Inspect excluded recordings.
excluded_recordings <- model_data %>%
  filter(!eligible_no_env) %>%
  select(
    file, fileMidUTC, longitude, latitude,
    depth, bathy, MPA, diel, dive_id, AR_group
  )

print(
  tibble::as_tibble(excluded_recordings),
  n = Inf,
  width = Inf
)

# Show the remaining excluded files.
model_data %>%
  filter(!eligible_no_env) %>%
  select(
    file, fileMidUTC, longitude, latitude, depth,
    bathy, MPA, diel, dive_id, AR_group, partial_dive
  ) %>%
  as_tibble() %>%
  print(n = Inf, width = Inf)


# ============================================================
# Final modelling and reporting workflow
# ============================================================
# Same covariates in all models; compare only their functional forms.
# Model-based Wald intervals; no resampling or bootstrap.
# The script does not assume which response needs quadratic depth.

# ------------------------------------------------------------
# 7. Model settings and analysis data
# ------------------------------------------------------------

confidence_level <- 0.95
# Intervals and tests are conditional on the selected model form.

responses <- c("pres_Pm", "del_clicks_all", "del_whistles_all")
response_labels <- c(pres_Pm = "Sperm whale",
                     del_clicks_all = "Delphinid echolocation clicks",
                     del_whistles_all = "Delphinid whistles")
colours <- c("Sperm whale" = "#FCA636FF",
             "Delphinid echolocation clicks" = "#7D2D3F", "Delphinid whistles" = "#E8B2A3")

as_binary <- function(x, name) {
  if (is.logical(x)) return(as.integer(x))
  if (is.character(x)) {
    x[x == "TRUE"] <- "1"
    x[x == "FALSE"] <- "0"
  }
  x <- suppressWarnings(as.numeric(x))
  x
}
for (response in responses) {
  model_data[[response]] <- as_binary(model_data[[response]], response)
}
common_eligible <- eligible_no_env & complete.cases(model_data[, responses])
core_data <- model_data[common_eligible, ] %>% arrange(fileMidUTC)
core_data$AR_group <- droplevels(core_data$AR_group)
# Retain the full dive_id factor sequence defined above: missing levels
# must not be dropped, because AR(1) distances follow that sequence.
core_data$date_local <- as.Date(core_data$fileMidUTC, tz = "Europe/Madrid")

scaling <- data.frame(variable = c("depth", "bathy"),
                      mean = c(mean(core_data$depth), mean(core_data$bathy)),
                      sd = c(sd(core_data$depth), sd(core_data$bathy)))
for (v in scaling$variable) {
  ss <- scaling[scaling$variable == v, ]
  core_data[[paste0(v, "_z")]] <- (core_data[[v]] - ss$mean) / ss$sd
}
# Squared terms are computed from centred/scaled depth inside formulas.
# Full quadratic = z + I(z^2). The linear term is never removed.

write.csv(core_data, file.path(scratch_dir, "core_model_data.csv"), row.names = FALSE)
write.csv(scaling, file.path(scratch_dir, "covariate_scaling.csv"), row.names = FALSE)
write.csv(model_data[!common_eligible, ], file.path(scratch_dir, "excluded_files.csv"),
          row.names = FALSE)

candidate_dir <- file.path(scratch_dir, "candidate_models")
dir.create(candidate_dir, showWarnings = FALSE)

# Protect cached fits against input or setting changes.

AIC_tolerance <- 2
cache_settings <- list(
  script_version = "final-depth-bathy-wald-v1",
  data = core_data[, c("file", "fileMidUTC", "depth_z", "bathy_z", "MPA", "diel",
                       "dive_id", "AR_group", "date_local", responses)],
  scaling = scaling,
  AIC_tolerance = AIC_tolerance,
  glmmTMB_version = as.character(packageVersion("glmmTMB")),
  TMB_version = as.character(packageVersion("TMB"))
)
settings_path <- file.path(scratch_dir, "cache_settings.rds")
if (file.exists(settings_path)) {
  if (!identical(readRDS(settings_path), cache_settings)) {
    stop("Inputs/settings/package changed. Choose a NEW outdir to avoid stale caches.")
  }
} else saveRDS(cache_settings, settings_path)

# ------------------------------------------------------------
# 8. Fit and compare four functional-form candidates per signal
# ------------------------------------------------------------
# LL: both linear; QL: quadratic glider depth only;
# LQ: quadratic bathymetry only; QQ: both quadratic.
# Every candidate contains MPA, diel and the SAME dive/day AR(1).
# Compare ML AIC on identical observations. This is shape comparison,
# not covariate inclusion/exclusion. No dredge or stepwise selection.
# Among candidates within 2 AIC units of the lowest AIC, choose the
# fewest quadratic terms. Break remaining ties by lowest AIC.
# This is a pragmatic parsimonious rule, not a significance test.

forms <- data.frame(
  form = c("LL", "QL", "LQ", "QQ"),
  depth_quadratic = c(FALSE, TRUE, FALSE, TRUE),
  bathy_quadratic = c(FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
fixed_formula <- function(response = NULL, depth_q = FALSE, bathy_q = FALSE) {
  rhs <- c("depth_z", if (depth_q) "I(depth_z^2)",
           "bathy_z", if (bathy_q) "I(bathy_z^2)", "MPA", "diel")
  as.formula(paste(if (is.null(response)) "" else response,
                   "~", paste(rhs, collapse = " + ")))
}
full_formula <- function(response, depth_q, bathy_q) {
  update(fixed_formula(response, depth_q, bathy_q),
         . ~ . + ar1(dive_id + 0 | AR_group))
}

ar_stats <- function(m) {
  block <- VarCorr(m)$cond$AR_group
  if (is.null(block)) stop("Missing AR_group covariance block.")
  sd_ar <- as.numeric(attr(block, "stddev")[1])
  rho <- as.numeric(attr(block, "correlation")[1, 2])
  c(SD = sd_ar, variance = sd_ar^2, rho = rho)
}
fit_checked <- function(d, response, depth_q, bathy_q) {
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    if (length(unique(d[[response]])) != 2L) stop("Response has only one class.")
    if (n_distinct(d$MPA) < 2L || n_distinct(d$diel) < 2L) {
      stop("Resample lacks a categorical predictor level.")
    }
    m <- glmmTMB(
      full_formula(response, depth_q, bathy_q), data = d,
      family = binomial(link = "logit"), REML = FALSE, na.action = na.fail,
      control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000))
    )
    if (!isTRUE(m$fit$convergence == 0) || !isTRUE(m$sdr$pdHess)) {
      stop("Convergence/Hessian checks failed.")
    }
    beta <- fixef(m)$cond
    V <- as.matrix(vcov(m)$cond)
    if (any(!is.finite(beta)) || any(!is.finite(V)) || any(diag(V) <= 0)) {
      stop("Invalid coefficient estimates/covariance.")
    }
    # Do not accept a silent rank adjustment (different fixed-effect set).
    expected <- colnames(model.matrix(fixed_formula(NULL, depth_q, bathy_q), d))
    if (!identical(names(beta), expected)) stop("Fixed-effect design is rank deficient/reordered.")
    ap <- ar_stats(m)
    if (any(!is.finite(ap)) || !is.finite(AIC(m))) stop("Invalid variance or AIC.")
    list(ok = TRUE, model = m, beta = beta, V = V, ar = ap,
         extreme = any(abs(beta) > 20) || any(sqrt(diag(V)) > 10) || ap[["SD"]] > 10)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(ok = FALSE, error = conditionMessage(e)))
  result$warnings <- paste(unique(warnings), collapse = " | ")
  result
}

candidate_models <- list()
comparison_rows <- list()
for (response in responses) {
  candidate_models[[response]] <- list()
  for (j in seq_len(nrow(forms))) {
    ff <- forms[j, ]
    path <- file.path(candidate_dir, paste0(response, "_", ff$form, ".rds"))
    result <- if (file.exists(path)) readRDS(path) else {
      rr <- fit_checked(core_data, response, ff$depth_quadratic, ff$bathy_quadratic)
      saveRDS(rr, path); rr
    }
    good <- isTRUE(result$ok)
    if (good) candidate_models[[response]][[ff$form]] <- result$model
    comparison_rows[[paste(response, ff$form)]] <- data.frame(
      response, form = ff$form, depth_quadratic = ff$depth_quadratic,
      bathy_quadratic = ff$bathy_quadratic,
      n_quadratic = as.integer(ff$depth_quadratic) + as.integer(ff$bathy_quadratic),
      converged = good, AIC = if (good) AIC(result$model) else NA_real_,
      warning = result$warnings,
      error = if (good) "" else result$error,
      stringsAsFactors = FALSE
    )
  }
}
comparison <- bind_rows(comparison_rows)
write.csv(comparison, file.path(scratch_dir, "functional_form_comparison.csv"), row.names = FALSE)
# Do not choose forms on an incomplete comparison without review.
if (any(!comparison$converged)) {
  stop("A functional-form candidate failed. Inspect functional_form_comparison.csv before selection.")
}
comparison <- comparison %>% group_by(response) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>% ungroup()
choices <- comparison %>% filter(delta_AIC <= AIC_tolerance + 1e-8) %>%
  arrange(response, n_quadratic, AIC) %>% group_by(response) %>% slice(1) %>% ungroup()
comparison$selected <- paste(comparison$response, comparison$form) %in%
  paste(choices$response, choices$form)
write.csv(comparison, file.path(scratch_dir, "functional_form_comparison.csv"), row.names = FALSE)
write.csv(choices, file.path(scratch_dir, "selected_functional_forms.csv"), row.names = FALSE)
print(comparison)

models <- setNames(vector("list", length(responses)), responses)
shape <- list()
for (response in responses) {
  ch <- choices[choices$response == response, ]
  models[[response]] <- candidate_models[[response]][[ch$form]]
  shape[[response]] <- list(depth_q = ch$depth_quadratic, bathy_q = ch$bathy_quadratic)
  saveRDS(models[[response]], file.path(scratch_dir, paste0(response, "_final_model.rds")))
}

# Fixed-effect design matrices for predictions.
fixed_X <- function(d, response) {
  sp <- shape[[response]]
  model.matrix(fixed_formula(NULL, sp$depth_q, sp$bathy_q), d)
}

# ------------------------------------------------------------
# 9. Prediction grids and fixed-effect model summaries
# ------------------------------------------------------------
# Continuous variables: full observed range, without percentile trimming. Other predictor at
# its mean; MPA outside and diel day. Plot categorical contrasts at
# both continuous means. All random effects are set to zero.
# These are reference-condition probabilities, not prevalence estimates.

predictor_labels <- c(depth = "Glider depth (m)", bathy = "Bathymetry (m)",
                      MPA = "MPA status", diel = "Time of day")
newdata_by_predictor <- list()
for (v in names(predictor_labels)) {
  values <- if (v %in% c("depth", "bathy")) {
    lower <- min(core_data[[v]])
    upper <- max(core_data[[v]])
    seq(lower, upper, length.out = 200)
  } else levels(core_data[[v]])
  nd <- core_data[rep(1L, length(values)), , drop = FALSE]
  nd$depth_z <- 0; nd$bathy_z <- 0
  nd$MPA <- factor(rep("outside", nrow(nd)), levels = levels(core_data$MPA))
  nd$diel <- factor(rep("day", nrow(nd)), levels = levels(core_data$diel))
  if (v %in% c("depth", "bathy")) {
    ss <- scaling[scaling$variable == v, ]
    nd[[paste0(v, "_z")]] <- (values - ss$mean) / ss$sd
  } else nd[[v]] <- factor(values, levels = levels(core_data[[v]]))
  newdata_by_predictor[[v]] <- nd
}
pred_meta <- bind_rows(lapply(names(predictor_labels), function(v) {
  nd <- newdata_by_predictor[[v]]
  data.frame(predictor = v,
             x = if (v %in% c("depth", "bathy")) {
               ss <- scaling[scaling$variable == v, ]
               nd[[paste0(v, "_z")]] * ss$sd + ss$mean
             } else seq_len(nrow(nd)),
             category = if (v %in% c("depth", "bathy")) NA_character_ else as.character(nd[[v]]))
}))
X_predictions <- setNames(lapply(responses, function(response) {
  do.call(rbind, lapply(newdata_by_predictor, fixed_X, response = response))
}), responses)

# Exact point estimates and model-based uncertainty (retained for audit).
# Table and figure use matching model-based Wald 95% intervals.
coefficient_rows <- list(); metric_rows <- list(); prediction_point <- list(); joint_rows <- list()
for (response in responses) {
  m <- models[[response]]; beta <- fixef(m)$cond
  V <- as.matrix(vcov(m)$cond); se <- sqrt(diag(V))
  X <- fixed_X(core_data, response)
  X <- X[, names(beta), drop = FALSE]
  eta <- as.numeric(X %*% beta)
  ap <- ar_stats(m)
  var_fixed <- var(eta)
  var_random <- ap[["variance"]]
  var_distribution <- pi^2 / 3
  total <- var_fixed + var_random + var_distribution
  # Latent-scale variance decomposition with theoretical Bernoulli-logit
  # residual variance pi^2/3. AR marginal variance is the mean diagonal
  # of Z Sigma Z'; for one homogeneous AR effect per file, it is sigma_AR^2.
  # ICC below is conditional on fixed covariates for SAME dive/day.
  # For different dives k levels apart within one day: ICC * rho^k.
  # Across different day groups the model assumes independent AR effects.
  metric_rows[[response]] <- data.frame(response,
                                        observations = nobs(m), presences = sum(core_data[[response]] == 1),
                                        absences = sum(core_data[[response]] == 0), dives = n_distinct(core_data$dive_id),
                                        days = n_distinct(core_data$date_local), AR1_SD = ap[["SD"]],
                                        AR1_variance = var_random, AR1_rho = ap[["rho"]],
                                        residual_variance_latent = var_distribution,
                                        ICC_same_dive_day_latent = var_random / (var_random + var_distribution),
                                        marginal_R2_latent = var_fixed / total,
                                        conditional_R2_latent = (var_fixed + var_random) / total,
                                        AIC = AIC(m), convergence = m$fit$convergence, pdHess = m$sdr$pdHess)
  coefficient_rows[[response]] <- data.frame(response, term = names(beta),
                                             estimate = as.numeric(beta), SE_log_odds = as.numeric(se),
                                             odds_ratio = exp(as.numeric(beta)), Wald_lower95 = exp(beta - qnorm(.975) * se),
                                             Wald_upper95 = exp(beta + qnorm(.975) * se), Wald_z = as.numeric(beta / se),
                                             Wald_p = as.numeric(2 * pnorm(abs(beta / se), lower.tail = FALSE)), row.names = NULL)
  XP <- X_predictions[[response]][, names(beta), drop = FALSE]
  pp <- pred_meta
  pp$response <- response
  pp$probability <- plogis(as.numeric(XP %*% beta))
  prediction_point[[response]] <- pp
  # Overall Wald tests for each covariate, including both terms together
  # when quadratic. No term is removed on the basis of these tests.
  term_groups <- list(depth = intersect(c("depth_z", "I(depth_z^2)"), names(beta)),
                      bathy = intersect(c("bathy_z", "I(bathy_z^2)"), names(beta)),
                      MPA = "MPAinside", diel = "dielnight")
  joint_rows[[response]] <- bind_rows(lapply(names(term_groups), function(v) {
    idx <- match(term_groups[[v]], names(beta)); b <- beta[idx]
    w <- as.numeric(crossprod(b, solve(V[idx, idx, drop = FALSE], b)))
    data.frame(response, predictor = v, Wald_chisq = w, df = length(idx),
               Wald_p = pchisq(w, length(idx), lower.tail = FALSE))
  }))
}
coefficients <- bind_rows(coefficient_rows)
model_metrics <- bind_rows(metric_rows)
joint_tests <- bind_rows(joint_rows)
# Audit outputs are shared across responses; no separate coefficient tables.
write.csv(coefficients, file.path(scratch_dir, "coefficient_audit_Wald.csv"), row.names = FALSE)
write.csv(model_metrics, file.path(scratch_dir, "model_metrics.csv"), row.names = FALSE)
write.csv(joint_tests, file.path(scratch_dir, "joint_covariate_tests.csv"), row.names = FALSE)

rm(candidate_models)
invisible(gc())

# ------------------------------------------------------------
# 10. Model-based confidence intervals for partial effects
# ------------------------------------------------------------
# Compute uncertainty on the logit scale, then transform to probability.
# Fixed-effect predictions set the fitted AR effects to zero.
# These are pointwise intervals, conditional on the selected functional form.

coefficient_results <- coefficients
prediction_rows <- lapply(responses, function(response) {
  m <- models[[response]]
  beta <- fixef(m)$cond
  V <- as.matrix(vcov(m)$cond)
  XP <- X_predictions[[response]][, names(beta), drop = FALSE]
  eta <- as.numeric(XP %*% beta)
  se <- sqrt(pmax(0, rowSums((XP %*% V) * XP)))
  point <- prediction_point[[response]]
  critical <- qnorm((1 + confidence_level) / 2)
  point$lower <- plogis(eta - critical * se)
  point$upper <- plogis(eta + critical * se)
  point
})
prediction_data <- bind_rows(prediction_rows)
prediction_data$signal <- factor(unname(response_labels[prediction_data$response]),
                                 levels = unname(response_labels))
write.csv(prediction_data, file.path(scratch_dir, "partial_effect_predictions.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 11. Predictor checks and recording counts
# ------------------------------------------------------------
# Fit convergence/Hessian were checked above. Save outcome counts by factor
# so sparse presence/absence cells can be inspected for each signal.
# These checks do not automatically accept or reject a model.

cell_counts <- bind_rows(lapply(responses, function(response) {
  bind_rows(lapply(c("MPA", "diel"), function(v) {
    core_data %>% group_by(level = .data[[v]]) %>%
      summarise(recordings = n(), presences = sum(.data[[response]] == 1),
                absences = sum(.data[[response]] == 0), .groups = "drop") %>%
      mutate(response = response, predictor = v, .before = 1)
  }))
}))
write.csv(cell_counts, file.path(scratch_dir, "presence_absence_by_factor.csv"), row.names = FALSE)
write.csv(data.frame(depth_bathymetry_correlation = cor(core_data$depth_z, core_data$bathy_z)),
          file.path(scratch_dir, "continuous_predictor_correlation.csv"), row.names = FALSE)
vif_rows <- bind_rows(lapply(responses, function(response) {
  X <- fixed_X(core_data, response)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  bind_rows(lapply(seq_len(ncol(X)), function(j) {
    r2 <- summary(lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared
    data.frame(response, term = colnames(X)[j], VIF = 1 / (1 - r2))
  }))
}))
# Squared and linear terms are related: their separate VIFs need joint interpretation.
write.csv(vif_rows, file.path(scratch_dir, "collinearity_VIF.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 12. Concise manuscript table: selected models only
# ------------------------------------------------------------
# S5 format: OR, 95% CI and p for each of the three selected models.
# Include only coefficients present in the selected model, not candidates.
# A full quadratic has BOTH z and z^2: these are two coefficients of one
# selected relationship. A dash means that coefficient is not in the model.
# No form labels, duplicate overall tests, rho, residual variance or AIC.
# Random-effect variance, ICC and R2 remain, as in the S5 example.
# Test statistics/df remain in audit files.

# Refresh labels here too, allowing sections 12-13 to run on loaded results.
response_labels <- c(pres_Pm = "Sperm whale",
                     del_clicks_all = "Delphinid echolocation clicks",
                     del_whistles_all = "Delphinid whistles")
colours <- c("Sperm whale" = "#FCA636FF",
             "Delphinid echolocation clicks" = "#7D2D3F",
             "Delphinid whistles" = "#C97878")  # original color was E8B2A3 but too light
prediction_data$signal <- factor(unname(response_labels[prediction_data$response]),
                                 levels = unname(response_labels))
predictor_labels[["diel"]] <- "Time of day"

term_order <- c("(Intercept)", "depth_z", "I(depth_z^2)",
                "bathy_z", "I(bathy_z^2)", "MPAinside", "dielnight")
term_labels <- c("(Intercept)" = "Intercept",
                 depth_z = "Glider depth", "I(depth_z^2)" = "Glider depth^2",
                 bathy_z = "Bathymetry", "I(bathy_z^2)" = "Bathymetry^2",
                 MPAinside = "MPA (inside)", dielnight = "Time of day (night)")
fmt <- function(x, digits = 2) {
  if (!is.finite(x)) return(as.character(x))
  if (x != 0 && (abs(x) < .01 || abs(x) >= 10000)) {
    formatC(x, format = "e", digits = 2)
  } else formatC(x, format = "f", digits = digits)
}
fmt_p <- function(p) {
  if (is.na(p)) "" else if (p < .001) "<0.001" else formatC(p, format = "f", digits = 3)
}
used_terms <- term_order[term_order %in% coefficient_results$term]
summary_labels <- c("Random effects", "AR(1) variance", "ICC",
                    "Marginal R2 / Conditional R2")
row_names <- c(unname(term_labels[used_terms]), summary_labels)
paper_table <- data.frame(Predictors = row_names, stringsAsFactors = FALSE)
for (response in responses) {
  lab <- response_labels[[response]]
  or_col <- ci_col <- p_col <- rep("", length(row_names))
  for (i in seq_along(used_terms)) {
    rr <- coefficient_results[coefficient_results$response == response &
                                coefficient_results$term == used_terms[i], ]
    if (!nrow(rr)) {
      or_col[i] <- ci_col[i] <- p_col[i] <- "-"
      next
    }
    stopifnot(nrow(rr) == 1L)
    or_col[i] <- fmt(rr$odds_ratio)
    ci_col[i] <- paste(fmt(rr$Wald_lower95), fmt(rr$Wald_upper95), sep = " - ")
    p_col[i] <- fmt_p(rr$Wald_p)
  }
  mm <- model_metrics[model_metrics$response == response, ]
  stopifnot(nrow(mm) == 1L)
  or_col[length(used_terms) + seq_along(summary_labels)] <- c(
    "", fmt(mm$AR1_variance), fmt(mm$ICC_same_dive_day_latent),
    paste(fmt(mm$marginal_R2_latent, 3), fmt(mm$conditional_R2_latent, 3), sep = " / "))
  paper_table[[paste(lab, "OR", sep = " | ")]] <- or_col
  paper_table[[paste(lab, "95% CI", sep = " | ")]] <- ci_col
  paper_table[[paste(lab, "p", sep = " | ")]] <- p_col
}
# ASCII punctuation in CSV avoids corrupted dash/superscript characters in Excel.
write.csv(paper_table, file.path(outdir, "paper_results_table.csv"),
          row.names = FALSE, na = "", fileEncoding = "UTF-8")

# Supporting numbers for writing the Results, NOT another manuscript table.
# The test statistic varies. Coefficient chi-square tests have 1 df;
# overall quadratic tests have 2 df. Do not mix these two tests in text.
results_for_text <- joint_tests %>% left_join(
  model_metrics %>% select(response, observations, presences, absences),
  by = "response")
write.csv(results_for_text, file.path(outdir, "statistics_for_results_text.csv"), row.names = FALSE)

escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
# Same compact formatting for the combined and individual model tables.
write_paper_table <- function(table, included_responses, filename) {
  stopifnot(ncol(table) == 1L + 3L * length(included_responses))
  # Keep every exported column to a single quantity. Summary values have
  # their own block and CSV; they never appear in the OR/CI/p columns.
  summary <- table[table$Predictors %in% setdiff(summary_labels, "Random effects"),
                   c(1L, 2L + 3L * (seq_along(included_responses) - 1L)), drop = FALSE]
  names(summary) <- c("Statistic", unname(response_labels[included_responses]))
  table <- table[!table$Predictors %in% summary_labels, , drop = FALSE]
  html <- c('<!DOCTYPE html><html><head><meta charset="utf-8"><style>',
            'body{font-family:Arial,sans-serif;margin:24px;font-size:12px}table{border-collapse:collapse;width:100%}',
            'th,td{padding:6px 8px;text-align:center}th:first-child,td:first-child{text-align:left}',
            'thead{border-top:1px solid #333;border-bottom:1px solid #333}p{line-height:1.5}',
            '.section{border-top:1px solid #333}.r2{border-top:1px solid #333;border-bottom:1px solid #333}',
            '</style></head><body><table><thead><tr><th rowspan="2">Predictor</th>',
            paste0('<th colspan="3">', escape_html(unname(response_labels[included_responses])),
                   '</th>', collapse = ""),
            '</tr><tr>', paste(rep('<th>OR</th><th>95% CI</th><th>p</th>',
                                   length(included_responses)), collapse = ""),
            '</tr></thead><tbody>')
  for (i in seq_len(nrow(table))) {
    values <- escape_html(unlist(table[i, ], use.names = FALSE))
    values[1] <- gsub("^2", "<sup>2</sup>", values[1], fixed = TRUE)
    if (table$Predictors[i] == "Marginal R2 / Conditional R2") {
      values[1] <- "Marginal R<sup>2</sup> / Conditional R<sup>2</sup>"
    }
    term <- names(term_labels)[match(table$Predictors[i], unname(term_labels))]
    if (!is.na(term)) for (j in seq_along(included_responses)) {
      rr <- coefficient_results[coefficient_results$response == included_responses[j] &
                                  coefficient_results$term == term, ]
      if (nrow(rr) && is.finite(rr$Wald_p) && rr$Wald_p < .05) {
        at <- 1L + 3L * j
        values[at] <- paste0("<b>", values[at], "</b>")
      }
    }
    cls <- if (table$Predictors[i] == "Random effects") ' class="section"' else
      if (table$Predictors[i] == "Marginal R2 / Conditional R2") ' class="r2"' else ''
    html <- c(html, paste0('<tr', cls, '>',
                           paste0('<td>', values, '</td>', collapse = ''), '</tr>'))
  }
  html <- c(html, '</tbody></table>',
            '<table style="margin-top:12px"><thead><tr><th>Model summary</th>',
            paste0('<th>', escape_html(unname(response_labels[included_responses])), '</th>', collapse = ""),
            '</tr></thead><tbody>')
  for (i in seq_len(nrow(summary))) {
    vals <- escape_html(unlist(summary[i, ], use.names = FALSE))
    vals[1] <- gsub("R2", "R<sup>2</sup>", vals[1], fixed = TRUE)
    html <- c(html, paste0('<tr>', paste0('<td>', vals, '</td>', collapse = ''), '</tr>'))
  }
  html <- c(html, '</tbody></table></body></html>')
  cat(html, file = file.path(outdir, paste0(filename, ".html")), sep = "\n")
  write.csv(table, file.path(outdir, paste0(filename, ".csv")),
            row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write.csv(summary, file.path(outdir, paste0(filename, "_summary.csv")),
            row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

# Save one combined table for all three signals.
write_paper_table(paper_table, responses, "paper_results_table")

# Also save one table per selected model, with short column names.
# Omit terms not present in that model; no candidate-model results are added.
model_file_labels <- c(pres_Pm = "sperm_whale",
                       del_clicks_all = "delphinid_clicks",
                       del_whistles_all = "delphinid_whistles")
for (j in seq_along(responses)) {
  response <- responses[j]
  columns <- c(1L, 1L + (j - 1L) * 3L + seq_len(3L))
  selected_terms <- coefficient_results$term[coefficient_results$response == response]
  keep_labels <- c(unname(term_labels[selected_terms]), summary_labels)
  individual <- paper_table[paper_table$Predictors %in% keep_labels, columns, drop = FALSE]
  names(individual) <- c("Predictors", "OR", "95% CI", "p")
  write_paper_table(individual, response,
                    paste0("paper_results_table_", model_file_labels[[response]]))
}

# ------------------------------------------------------------
# 13. Partial-effect panels: all models, including non-significant effects
# ------------------------------------------------------------
# Lines/points: original selected GLMM estimates at reference conditions.
# Ribbons/error bars: pointwise model-based Wald confidence intervals.
# One panel per covariate, one common legend, automatic panel y ranges.
# Do not infer an overall association test from overlap of these bands.

# Panel order.
predictor_labels <- c(
  bathy = "Bathymetry (m)",
  depth = "Glider depth (m)",
  MPA = "MPA status",
  diel = "Time of day"
)

# Colours and legend order.
colours <- c(
  "Delphinid echolocation clicks" = "#7D2D3F",
  "Delphinid whistles" = "#C97878",
  "Sperm whale" = "#FCA636FF"
)

legend_order <- c(
  "Delphinid echolocation clicks",
  "Delphinid whistles",
  "Sperm whale"
)

prediction_data$signal <- factor(
  as.character(prediction_data$signal),
  levels = legend_order
)

# Legend key: coloured line across a shaded rectangle.
draw_key_line_ribbon <- function(data, params, size) {
  grid::grobTree(
    grid::rectGrob(
      width = .9, height = .8,
      gp = grid::gpar(
        fill = scales::alpha(data$colour[1], .15),
        col = NA
      )
    ),
    grid::segmentsGrob(
      x0 = .05, x1 = .95, y0 = .5, y1 = .5,
      gp = grid::gpar(col = data$colour[1], lwd = 2)
    )
  )
}

make_panel <- function(v) {
  d <- prediction_data[prediction_data$predictor == v, ]
  
  p <- ggplot(
    d,
    aes(x, probability, colour = signal, fill = signal, group = signal)
  ) +
    # Border behind the prediction layers.
    annotate(
      "rect",
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf,
      fill = NA, colour = "black", linewidth = .5
    )
  
  if (v %in% c("depth", "bathy")) {
    p <- p +
      geom_ribbon(
        aes(ymin = lower, ymax = upper),
        alpha = .15, colour = NA, show.legend = FALSE
      ) +
      geom_line(linewidth = 1, show.legend = FALSE)
  } else {
    dodge <- position_dodge(width = .4)
    
    p <- p +
      geom_point(
        position = dodge, size = 2, show.legend = FALSE
      ) +
      # Error bars drawn over points and border.
      geom_errorbar(
        aes(ymin = lower, ymax = upper),
        position = dodge, width = .10,
        linewidth = .7, show.legend = FALSE
      ) +
      scale_x_continuous(
        breaks = 1:2,
        labels = if (v == "MPA") c("Outside", "Inside") else c("Day", "Night"),
        expand = expansion(add = .4)
      )
  }
  
  if (v == "bathy") {
    p <- p + scale_x_continuous(
      breaks = seq(0, 2500, by = 500),
      limits = c(0, NA),
      expand = expansion(mult = c(0, .05))
    )
  } else if (v == "depth") {
    p <- p + scale_x_continuous(
      breaks = seq(0, 800, by = 200),
      limits = c(0, NA),
      expand = expansion(mult = c(0, .05))
    )
  }
  legend_data <- d[!duplicated(d$signal), ]
  
  p +
    geom_point(
      data = legend_data,
      aes(x, probability, colour = signal),
      inherit.aes = FALSE,
      alpha = 0,
      show.legend = TRUE,
      key_glyph = draw_key_line_ribbon
    ) +
    scale_colour_manual(
      values = colours,
      limits = legend_order,
      breaks = legend_order,
      drop = FALSE,
      name = "Signal type"
    ) +
    scale_fill_manual(
      values = colours,
      limits = legend_order,
      drop = FALSE,
      name = NULL
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = .1),
      limits = c(0, NA),
      breaks = scales::breaks_pretty(n = 5),
      expand = expansion(mult = c(0, .08))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = predictor_labels[[v]],
      y = "Predicted detection probability"
    ) +
    guides(
      fill = "none",
      colour = guide_legend(
        title.position = "left",
        nrow = 1,
        override.aes = list(alpha = 1)
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = .3),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.width = grid::unit(.45, "cm"),
      legend.key.height = grid::unit(.4, "cm"),
      axis.text = element_text(size = 14, colour = "gray20"),
      axis.title.x = element_text(
        size = 15,
        margin = margin(t = 12)
      ),
      axis.title.y = element_text(
        size = 15,
        margin = margin(r = 12)
      ),
      plot.margin = margin(8, 10, 8, 8)
    )
}

# Build panels in the specified order.
panels <- lapply(names(predictor_labels), make_panel)

# Export individual panels.
for (i in seq_along(panels)) {
  for (extension in c("png", "pdf")) {
    ggsave(
      file.path(
        outdir,
        paste0("partial_effect_", names(predictor_labels)[i], ".", extension)
      ),
      panels[[i]],
      width = 18, height = 14, units = "cm",
      dpi = 400, bg = "white"
    )
  }
}

# Top left: bathymetry. Top right: glider depth.
# Bottom left: MPA. Bottom right: time of day.
# Shared legend, no panel letters or captions.
final_plot <- wrap_plots(
  panels,
  ncol = 2,
  byrow = TRUE,
  guides = "collect"
)

final_plot <- final_plot & theme(legend.position = "bottom")

#print(final_plot)

# Export combined figure.
for (extension in c("png", "pdf")) {
  ggsave(
    file.path(outdir, paste0("all_models_partial_effects.", extension)),
    final_plot,
    width = 29, height = 23, units = "cm",
    dpi = 400, bg = "white"
  )
}


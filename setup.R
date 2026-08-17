
#--------------------------------------------------------------------------------
# setup.R - Setup project
#--------------------------------------------------------------------------------

# Project name
project <- "ms-glider-cescab"

#--------------------------------------------------------------------------------
# 1. Install/load required packages
#--------------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")

pacman::p_load(

  # spatial data and GIS
  "sf", "terra", "lwgeom", "ncdf4",

  # data manipulation and import
  "data.table", "tidyverse", "lubridate",
  "fuzzyjoin", "R.matlab", "irr",

  # plotting and visualisation
  "ggspatial", "scico", "viridis", "tidyterra",
  "RColorBrewer", "egg", "scales", "png", "grid", "suncalc",

  install = FALSE # change to TRUE
)

#--------------------------------------------------------------------------------
# 2. Set main data path
#--------------------------------------------------------------------------------

# Change this path to the local directory where input/output data will be stored
main_dir <- file.path("C:/Users/greta/SML Dropbox/gitdata", project)

if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

#--------------------------------------------------------------------------------
# 3. Create input and output folders
#--------------------------------------------------------------------------------

input_dir <- file.path(main_dir, "input")
if (!dir.exists(input_dir)) dir.create(input_dir, recursive = TRUE)

output_dir <- file.path(main_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#--------------------------------------------------------------------------------
# 4. Shared cartography folder
#--------------------------------------------------------------------------------

carto_dir <- file.path(input_dir, "carto")

#--------------------------------------------------------------------------------
# 5. Define study area
#--------------------------------------------------------------------------------

# coordinate reference system
my_crs <- st_crs(4326) # WGS84

# CESCAB study area
study_area <- st_as_sfc(
  st_bbox(
    c(xmin = 2, xmax = 3.8, ymin = 37.5, ymax = 39.5),
    crs = my_crs
  )
)

study_area <- st_as_sf(study_area)

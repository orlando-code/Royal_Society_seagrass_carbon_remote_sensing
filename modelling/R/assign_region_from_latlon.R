#!/usr/bin/env Rscript

# Assign Region categories to points based on latitude/longitude.
#
# This script provides a simple, reproducible mapping from (lon, lat) to
# the Region factor used elsewhere in the project. It is intended for:
# - Backfilling Region for new core points
# - Optionally assigning Region to prediction-grid cells (if desired)
#
# NOTE: The current implementation uses approximate bounding boxes for
# broad European marine regions. You may refine these polygons later
# (e.g. using official shapefiles) without changing the interface.

rm(list = ls())
setwd(here::here())

source("modelling/R/helpers.R")
load_packages(c("here", "dplyr"))


load_packages(c("tidyverse"))


# region_shapes is now an sf data frame with regions: Baltic Sea, Mediterranean Sea, Black Sea, 
# North European Atlantic, and South European Atlantic with their respective geometries.

#' Assign Region based on lon/lat using simple geographic rules.
#'
#' @param dat Data frame with numeric columns longitude and latitude.
#' @return Same data frame with a character column region (if not present) or
#'         with missing values in region filled where possible.
assign_region_from_latlon <- function(dat) {
  if (!all(c("longitude", "latitude") %in% names(dat))) {
    stop("assign_region_from_latlon: dat must have 'longitude' and 'latitude' columns.")
  }
  dat <- dplyr::mutate(dat,
    longitude = as.numeric(.data$longitude),
    latitude  = as.numeric(.data$latitude)
  )

  # Preserve existing region labels if present; otherwise create an empty column
  if (!"region" %in% names(dat)) {
    dat$region <- NA_character_
  }

  # Helper for conditional assignment only where region is NA
  assign_where_na <- function(idx, value) {
    dat$region[idx & is.na(dat$region)] <<- value
  }

  lon <- dat$longitude
  lat <- dat$latitude

  # VERY APPROXIMATE BOUNDARIES (consistent with existing naming conventions):
  # These are intentionally conservative and can be refined later.

  # Mediterranean Sea
  assign_where_na(
    lon >= -6 & lon <= 36 & lat >= 30 & lat <= 46,
    "Mediterranean Sea"
  )

  # Baltic Sea
  assign_where_na(
    lon >= 10 & lon <= 30 & lat >= 53 & lat <= 66,
    "Baltic Sea"
  )

  # Black Sea
  assign_where_na(
    lon >= 27 & lon <= 42 & lat >= 40 & lat <= 47,
    "Black Sea"
  )

  # North European Atlantic (roughly NE Atlantic shelf seas)
  assign_where_na(
    lon >= -20 & lon <= 10 & lat >= 45 & lat <= 66,
    "North European Atlantic"
  )

  # South European Atlantic (Iberian and nearby)
  assign_where_na(
    lon >= -20 & lon <= 10 & lat >= 30 & lat < 45,
    "South European Atlantic"
  )

  # Any remaining NA region can stay NA (or be assigned "Unknown")
  dat$region[is.na(dat$region)] <- dat$region[is.na(dat$region)]  # no-op placeholder

  dat
}

# Example usage (uncomment to run interactively):
# dat <- readr::read_rds("data/all_extracted_new.rds")
# dat <- assign_region_from_latlon(dat)
# readr::write_rds(dat, "data/all_extracted_new_with_region.rds")


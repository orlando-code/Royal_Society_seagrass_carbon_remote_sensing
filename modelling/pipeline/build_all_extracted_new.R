#!/usr/bin/env Rscript

# Build (or rebuild) data/all_extracted_new.rds for the paper pipeline.
# - Uses nearest-neighbour extraction (this was shown to lead to models with highest R²).
# - Attaches key metadata columns.
# - Makes all column names lower-case before saving.
# - Clips remote-sensed covariates to zero to remove any processing errors.

# - If data/all_extracted_new.rds exists and FORCE_REBUILD_ALLEXTRACTED is FALSE:
#   it will do nothing (just print a message)
#   ELSE: it will re-create the file.

setwd(here::here())
load_packages(c("here", "dplyr", "ggplot2", "maps"))
source("modelling/R/extract_covariates_from_rasters.R")
source("modelling/R/helpers.R")

# save here
all_extracted_path <- "data/all_extracted_new.rds"

# Global flag: set TRUE to force re-extraction even if file exists
if (!exists("FORCE_REBUILD_ALLEXTRACTED")) {
  FORCE_REBUILD_ALLEXTRACTED <- TRUE
}

if (file.exists(all_extracted_path) && !isTRUE(FORCE_REBUILD_ALLEXTRACTED)) {
  cat("data/all_extracted_new.rds already exists. Skipping rebuild (set FORCE_REBUILD_ALLEXTRACTED <- TRUE to force).\n")
  quit(save = "no")
}

cat("======================================================================\n")
cat("BUILDING data/all_extracted_new.rds FOR PAPER PIPELINE\n")
cat("======================================================================\n\n")

cat("Loading original core-level data...\n")
dat <- read_rds("data/Updated_df/data/For_modeling_df_v3.rds")
cat("Loaded", nrow(dat), "observations\n\n")

cat("Extracting covariates at point locations with method = 'nearest'...\n")
cat("(This may take a few minutes)\n\n")

all_extracted <- extract_covariates_at_points(
  points = dat[, c("latitude", "longitude")],
  use_closest = TRUE,
  method = "nearest"
)

cat("\nExtraction complete.\n\n")

# Attach metadata
meta_cols <- c(
  "latitude", "longitude",
  "number_id_final_version",
  "seagrass_species",
  "Region",
  "median_carbon_density_100cm"
)

missing_meta <- setdiff(meta_cols, names(dat))
if (length(missing_meta) > 0) {
  stop(
    "Missing metadata columns in For_modeling_df_v3.rds: ",
    paste(missing_meta, collapse = ", ")
  )
}

# Attach covariates to survey site metadata
all_extracted_new <- dplyr::select(dat, dplyr::all_of(meta_cols)) %>%
  dplyr::bind_cols(all_extracted[, setdiff(names(all_extracted), c("latitude", "longitude"))])

# Clip to zero any values in remote-sensed covariates with negative values (sign of processing error)
all_extracted_new <- process_rs_covariates(all_extracted_new)

# Report on the dataset condition
cat("Rows:", nrow(all_extracted_new), "  Columns:", ncol(all_extracted_new), "\n\n")

# Basic NA diagnostics
covariate_cols <- grep("_closest$", names(all_extracted_new), value = TRUE)
covariate_cols <- setdiff(covariate_cols, meta_cols)
na_rows <- all_extracted_new[rowSums(is.na(all_extracted_new[, covariate_cols, drop = FALSE])) > 0, ]
cat("Rows with at least one NA in covariates:", nrow(na_rows), "\n\n")

# Count NA per point
if (length(covariate_cols) > 0) {
  all_extracted_new$n_na_per_point <- rowSums(is.na(all_extracted_new[, covariate_cols, drop = FALSE]))
}

# Lower-case all column names
names(all_extracted_new) <- tolower(names(all_extracted_new))

cat("Column names after lower-casing:\n")
cat("  ", paste(names(all_extracted_new)[seq_len(min(20, ncol(all_extracted_new)))], collapse = ", "), "\n\n")

# Save to file
write_rds(all_extracted_new, all_extracted_path)
cat("Saved new data to:", all_extracted_path, "\n")
cat("======================================================================\n")
cat("BUILD COMPLETE\n")
cat("======================================================================\n\n")

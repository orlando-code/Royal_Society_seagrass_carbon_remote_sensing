# Predict seagrass carbon density at sample coordinates using the pre-trained model.
#
# First install the renv package and run `renv::restore()` to install the
# necessary packages.

sys.source("modelling/R/init_repo.R", envir = .GlobalEnv)
project_root <- seagrass_init_repo(
  packages = c("dplyr", "readxl", "readr"),
  source_files = c(
    "modelling/R/helpers.R",
    "modelling/R/ml.R",
    "modelling/R/extract_covariates_from_rasters.R"
  ),
  include_helpers = FALSE,
  require_core_inputs = FALSE,
  check_renv = FALSE
)
sys.source("review/env_training_comparison.R", envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
dat_fp <- file.path("data", "review", "SMEEF Donor Meadow Coordinates.xlsx")
model_fp <- file.path("data", "review", "GPR_final.rds")
output_fp <- file.path("output", "review", "donor_meadow_carbon_density_predictions.csv")
env_comparison_output_dir <- file.path("output", "review", "donor_meadow_env_comparison")
train_data_fp <- file.path(project_root, "data", "all_extracted_new.rds")

# make sure output directory exists
dir.create(dirname(output_fp), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load sample coordinates from file 
# -----------------------------------------------------------------------------
if (grepl("\\.xlsx$", dat_fp, ignore.case = TRUE)) {
  df_raw <- readxl::read_excel(dat_fp, sheet = 1)
} else if (grepl("\\.csv$", dat_fp, ignore.case = TRUE)) {
  df_raw <- readr::read_csv(dat_fp, show_col_types = FALSE)
} else {
  stop("File extension not recognized. Please provide a .csv or .xlsx file.")
}

required_columns <- c("longitude", "latitude", "seagrass_species")
missing_cols <- setdiff(required_columns, names(df_raw))
if (length(missing_cols) > 0) {
  stop(
    "Input data is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

points <- df_raw %>%
  dplyr::select(dplyr::all_of(required_columns))

# -----------------------------------------------------------------------------
# Load model and extract environmental covariates at point locations
# -----------------------------------------------------------------------------
model <- readRDS(model_fp)
predictor_vars <- model$predictor_vars

raster_covars <- setdiff(
  predictor_vars,
  c("seagrass_species", "longitude", "latitude")
)
missing_rasters <- setdiff(tolower(raster_covars), raster_covariates)
if (length(missing_rasters) > 0) {
  stop(
    "Model requires raster covariate(s) not found under data/env_rasters directory: ",
    paste(missing_rasters, collapse = ", ")
  )
}

cat("\n\nExtracting environmental covariates at", nrow(points), "point(s)...\n")
pred_data <- extract_covariates_at_points(
  points = points,
  covariates = raster_covars,
  use_closest = TRUE
)
pred_data <- process_rs_covariates(pred_data)

na_covars <- raster_covars[
  vapply(raster_covars, function(v) any(is.na(pred_data[[v]]), na.rm = TRUE), logical(1))
]
if (length(na_covars) > 0) {
  warning(
    "Missing covariate values at some points: ",
    paste(na_covars, collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Compare environmental vectors to training data
# -----------------------------------------------------------------------------
train_data <- readRDS(train_data_fp)
env_cmp <- compare_prediction_env_to_training(
  training_data = train_data,
  prediction_data = pred_data,
  predictor_vars = predictor_vars,
  output_dir = env_comparison_output_dir,
  make_plots = FALSE
)

# -----------------------------------------------------------------------------
# Predict carbon density
# -----------------------------------------------------------------------------
cat("\n\nPredicting carbon density at", nrow(points), "point(s)...\n")
pred <- predict_model(model, pred_data, se = TRUE)

results <- dplyr::bind_cols(
  df_raw,
  tibble::tibble(
    predicted_carbon_density = pred$mean,
    predicted_se = pred$se
  ),
  env_cmp$point_summary %>%
    dplyr::select(
      env_similarity_score,
      env_flag_low_similarity,
      env_flag_any_predictor_outside_range,
      env_n_predictors_outside_range,
      env_pct_predictors_outside_range,
      dplyr::any_of("env_flag_novel_species"),
      env_flag_extrapolation
    )
)

# -----------------------------------------------------------------------------
# Append to original dataframe and write results to file
# -----------------------------------------------------------------------------
df_raw <- df_raw %>%
  dplyr::left_join(results, by = c("longitude", "latitude", "seagrass_species"))

dir.create("output", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(df_raw, output_fp)
cat("\n\nWrote predictions to", output_fp, "\n")


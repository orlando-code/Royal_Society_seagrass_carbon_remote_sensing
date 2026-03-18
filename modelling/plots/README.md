# Plot scripts

Figure-generation scripts called from `**run_paper.R**` (and can be run standalone if globals and inputs are set). Outputs go to `**output/<cv_regime>/cv_pipeline/**`, `**output/<cv_regime>/covariate_selection/**`, `**output/<cv_regime>/predictions/**`, and shared `**output/supplement/**` as below (`<cv_regime>` is derived from `cv_type` in `run_paper.R`).

---

## spatial_categorical_effect_all_models.R

**Called when:** Step 2b.

For each model (GPR, GAM, XGB), runs CV with several predictor configurations: env-only, +latitude, +longitude, +lat+lon, +region (if present), +lat+lon+region. Uses the same CV type and folds as the main pipeline. Produces a bar chart of mean R² by configuration and model, saved in `**output/<cv_regime>/cv_pipeline/**`. Useful to see how much spatial and regional terms add for each model.

---

## partial_dependence_all_models.R

**Called when:** Step 6.

Builds **iml::Predictor** wrappers for each final model (from `**output/<cv_regime>/final_models/<model>_final.rds**`), then computes partial dependence with **iml::FeatureEffect** (method `"pdp"**) for the top predictors (optionally limited by **`pruned_variables_to_include_.csv`** if present). Saves faceted PDPs as **`output/<cv_regime>/covariate_selection/pdp_.png`** for XGB, GAM, and GPR.

---

## spatial_prediction_maps.R

**Called when:** Step 7.

Loads `**output/<cv_regime>/final_models/<model>_final.rds**` for each model, builds a prediction grid from rasters (cached in `**output/cache/**`), applies the same bathymetry filter as elsewhere, and predicts on the grid. Saves a mean-prediction map per model (`**xgb_prediction_map.png**`, `**gam_prediction_map.png**`, `**gpr_prediction_map.png**`) and standard-error maps where available (`**gpr_se_map.png**`, `**gam_se_map.png**`) in `**output/<cv_regime>/predictions/**`. Uses an inline `**plot_prediction_map()**` helper (raster + world polygon + optional observation points); no dependency on the retired **plot_gpr_spatial_maps.R**.

---

## supplement.R

**Called when:** Step 8.

Produces all supplement figures: region outlines (MEOW-based), target histograms, correlation bar chart, prediction grid and similarity scores, and similarity histogram. Uses `**create_prediction_grid_from_rasters`** (grid cached in `**output/cache/**`). Saves everything to `**output/supplement/**`. Optional helpers **plot_raster_stack** and **plot_env_pairs** are commented out in the pipeline version; see the FLAG in the script to re-enable them for ad-hoc figures.
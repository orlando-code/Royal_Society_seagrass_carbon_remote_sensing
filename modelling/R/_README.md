# R helpers and core logic

Shared R code used by the pipeline and plot scripts. Startup is centralized via
`init_repo.R`; `helpers.R` sources `ml.R`.

---

## helpers.R

Central utility and modelling-help library. It:

- **Sources ml.R** (so loading helpers gives you scaling, encoding, fit_*, predict_model, etc.).
- **Package handling**: **load_packages()**, **install_packages()** for required/optional packages.
- **Spatial CV caching**: **get_cached_spatial_folds()** – builds or loads spatial (blockCV) or random folds, keyed by data hash, block size, n_folds, exclude_regions; cache path is **`output/cache/<cache_tag>_<hash>_folds.rds`**.
- **Data prep and response transform**: **transform_response()**, **inverse_response_transform()** for log-scale fitting.
- **Variable selection**: **prune_by_correlation()**, **select_top_env_then_species()**, **get_per_model_vars()**, **load_model_vars()**, **read_model_vars()**, **combine_pruned_model_variables()** – used by covariate pruning and downstream steps to get per-model predictor sets from **pruned_model_variables_perm.csv** / **pruned_model_variables_shap.csv**.
- **CV and tuning helpers**: **permutation_importance_cv()**, **compute_shap_importance()**, **run_cv()**, **prepare_data_for_model()**, **calculate_metrics()**, and model-specific tune helpers (**tune_gpr_cv**, **tune_xgboost**, etc.) used by **hyperparameter_tuning_pipeline.R** and **fit_final_models.R**.
- **Process RS covariates**: **process_rs_covariates()** – clips selected remote-sensing columns to ≥ 0.

Many pipeline scripts depend on **get_cached_spatial_folds**, **load_model_vars**, **permutation_importance_cv**, **compute_shap_importance**, and **ml::prepare_data_for_model** / **calculate_metrics**.

---

## ml.R

Machine learning and prediction layer:

- **Scaling and encoding**: **compute_scale_params()**, **apply_scaling()**, **ml::prepare_predictors_train()**, **prepare_predictors_train_numeric_only()**, **prepare_predictors_new()**, **apply_categorical_encoding()** – used for XGB, GAM, GPR (and RF where used).
- **Response transform**: **transform_response()**, **inverse_response_transform()** (used in concert with helpers).
- **Model fitting**: **fit_xgboost()**, **fit_gam()**, **fit_gpr()**, **fit_rf()** – take train (and optionally test) data and predictor names; **fit_gpr** supports **prediction_grid** and returns predictions + SE.
- **Unified prediction**: **predict_model()** – takes a saved model list (as from **output/<cv_regime>/final_models/*_final.rds**), new data, and optional **se = TRUE**; dispatches by model type (GPR returns SE when requested, others return **se = NULL**). **infer_model_type()** identifies GPR/XGB/GAM/RF from the object.
- **Data prep for CV**: **prepare_data_for_model()** – returns train/test (and predictor_vars) in the form expected by **fit_***.

Used by **cv_pipeline.R**, **hyperparameter_tuning_pipeline.R**, **fit_final_models.R**, **permutation_importance_final.R**, **shap_importance_final.R**, and **spatial_prediction_maps.R** (for **predict_model**).

---

## extract_covariates_from_rasters.R

Raster and grid handling:

- **NetCDF / raster config**: **build_covariate_config_from_dir()**, **inspect_nc_file()** – discovery of variables and lat/lon dimensions. Defines **raster_covariates** (and related config) used across the project.
- **Extraction at points**: **extract_covariates_at_points()** – extracts raster values at a set of lon/lat points (nearest or bilinear), with optional **use_closest** for filling coastal NAs.
- **Prediction grid**: **create_prediction_grid_from_rasters()** – builds a regular lon/lat grid, runs extraction for all **raster_covariates**, and optionally caches the result (cache path is passed in by the caller; pipeline uses **`output/cache/prediction_grid_cache_*.rds`**).

Sourced by **run_multiseed_pixel_grouped.R** (to set **raster_covariates**), **build_all_extracted_new.R**, **covariate_pruning_pipeline.R**, **cv_pipeline.R**, **spatial_prediction_maps.R**, **supplement.R**, and others that need the grid or the covariate list.

---

## assign_region_from_latlon.R

Maps longitude/latitude to a small set of regions (e.g. Mediterranean, North European Atlantic, Black Sea) by intersecting points with grouped MEOW ecoregion polygons.

- **meow_region_shapes()** – loads and groups the MEOW shapefile (cached in closure).
- **assign_region_from_latlon()** – adds a **region** column to a data frame that has **longitude** and **latitude**.

Used for **exclude_regions** filtering and for supplement/summary maps. Depends on **data/MEOW/meow_ecos.shp**.

---

## plot_config.R (in `modelling/plots/`)

Plot styling and labelling live in **`modelling/plots/plot_config.R`** (not under `modelling/R/`). It defines:

- **REGION_COLOURS**, **SPECIES_COLOURS**, **MODEL_COLOURS** – named colour vectors for consistent use in figures.
- **VAR_LABELS** – human-friendly names for covariates (used in axis/legend labels).
- **label_vars()** – maps variable names to display labels (falls back to the original name if not in **VAR_LABELS**).

Sourced by **`helpers.R`** (via `seagrass_init_repo` defaults), **`cv_pipeline.R`**, **`supplement.R`**, **`model_comparison.R`**, and plot scripts under **`modelling/plots/`** that need consistent colours or variable labels.

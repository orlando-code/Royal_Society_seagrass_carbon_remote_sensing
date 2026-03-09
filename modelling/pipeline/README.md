# Pipeline scripts

These scripts are run in sequence by **`run_paper.R`**. Each one expects globals (e.g. `target_var`, `exclude_regions`, `model_list`, `n_folds`, `cv_type`) to be set by the driver; they read from **`data/all_extracted_new.rds`** and write under **`output/`**.

---

## build_all_extracted_new.R

**Called when:** Step 0 (only if `data/all_extracted_new.rds` does not exist).

Builds the main extracted dataset from your raw core-level data and rasters. Uses nearest-neighbour extraction, lowercases column names, clips remote-sensed covariates to zero, and attaches metadata. Skips if the file already exists unless you set `FORCE_REBUILD_ALLEXTRACTED <- TRUE`. Depends on **`extract_covariates_from_rasters.R`** for the raster config and extraction helpers.

---

## covariate_pruning_pipeline.R

**Called when:** Step 1.

Runs variable selection for the paper: optional correlation filter, then per-model permutation importance (and optionally SHAP via `iml`) to get a ranked set of predictors per model. Writes importance CSVs and the combined pruned sets (**`pruned_model_variables_perm.csv`**, **`pruned_model_variables_shap.csv`**) to **`output/covariate_selection/`**. Uses globals: `use_correlation_filter`, `correlation_filter_threshold`, `permutation_max_vars`, `permutation_coverage`, `use_shap_per_model`, `model_list`.

---

## cv_pipeline.R

**Called when:** Steps 2 and 4.

Loads per-model predictor sets from the pruning CSVs (SHAP or permutation, depending on `use_shap_per_model`), builds spatial (or random) CV folds (using **`get_cached_spatial_folds`**; cache goes to **`output/cache/`**), and runs GPR, GAM, and XGB for each block size in `cv_blocksize_scan` (or a single `cv_blocksize`). Writes CV results, predictions, and comparison plots to **`output/cv_pipeline/`**. Does not run hyperparameter tuning (that’s step 3).

---

## hyperparameter_tuning_pipeline.R

**Called when:** Step 3.

Tunes GPR (kernel, nugget), XGBoost (nrounds, max_depth, learning_rate, subsample, colsample_bytree), and GAM (k_spatial) using the same folds and per-model predictor sets as the rest of the pipeline. Saves **`best_config_gpr.rds`**, **`best_config_gam.rds`**, **`best_config_xgb.rds`** in **`output/cv_pipeline/`**. These are read by **fit_final_models.R** and by the importance scripts (3b, 3c).

---

## permutation_importance_final.R

**Called when:** Step 3b.

Computes permutation feature importance for each model **after** tuning: uses the best config and best covariate set per model, runs **`permutation_importance_cv`** with the same spatial/random folds, and saves **`importance_perm_<model>.csv`** and **`importance_perm_<model>.png`** in **`output/cv_pipeline/`**. Requires the tuning RDS files and the pruned variable CSVs from step 1.

---

## shap_importance_final.R

**Called when:** Step 3c.

Computes SHAP (Shapley) feature importance for each model after tuning: fits each model with best config and best vars on the full data, then uses **`iml::Shapley`** to get mean |φ| per variable. Writes **`importance_shap_<model>.csv`** and **`importance_shap_<model>.png`** to **`output/cv_pipeline/`**. Uses the same best-config and pruned-variable inputs as the permutation script.

---

## fit_final_models.R

**Called when:** Step 5.

Loads the best config per model from **`output/cv_pipeline/best_config_*.rds`** and the per-model predictor sets from **`output/covariate_selection/`** (SHAP or permutation). Fits XGB, GAM, and GPR on **all** training data with those settings and saves **`XGB_final.rds`**, **`GAM_final.rds`**, **`GPR_final.rds`** to **`output/final_models/`**. Each RDS contains the fitted model, `predictor_vars`, `scale_params`, `encoding`, `hyperparams`, and `train_metrics`. These objects are used by the prediction maps and partial dependence scripts.

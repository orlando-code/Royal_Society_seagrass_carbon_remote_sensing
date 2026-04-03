# Pipeline scripts

These scripts are orchestrated by `modelling/run_multiseed_pixel_grouped.R` and configured via `modelling/config/pipeline_config.R`. Most scripts expect globals (for example `target_var`, `exclude_regions`, `model_list`, `n_folds`, `cv_type`) set by the driver. They read from `data/all_extracted_new.rds` and write under `output/`.

---

## build_all_extracted_new.R

**Called when:** Step 0 (only if `data/all_extracted_new.rds` is missing).

Builds the main extracted dataset from raw core-level data and rasters. Uses nearest-neighbour extraction, lowercases column names, clips remote-sensed covariates to zero, and attaches metadata. This script supports a `FORCE_REBUILD_ALLEXTRACTED` override for rebuild behavior; check the script defaults before relying on skip/rebuild semantics. Depends on `extract_covariates_from_rasters.R`.

---

## covariate_pruning_pipeline.R

**Called when:** Step 1.

Runs variable selection for the paper: optional correlation filter, then per-model permutation importance (and optionally SHAP via `iml`) to get a ranked set of predictors per model. Writes importance CSVs and the combined pruned sets (**`pruned_model_variables_perm.csv`**, **`pruned_model_variables_shap.csv`**) to **`output/<cv_regime>/covariate_selection/`**. Uses globals: `use_correlation_filter`, `correlation_filter_threshold`, `permutation_max_vars`, `permutation_coverage`, `use_shap_per_model`, `model_list`.

---

## cv_pipeline.R

**Called when:** Steps 2 and 4.

Loads per-model predictor sets from the pruning CSVs (SHAP or permutation, depending on `use_shap_per_model`), builds spatial (or random) CV folds (using **`get_cached_spatial_folds`**; cache goes to **`output/cache/`**), and runs GPR, GAM, and XGB for each block size in `cv_blocksize_scan` (or a single `cv_blocksize`). Writes CV results, predictions, and comparison plots to **`output/<cv_regime>/cv_pipeline/`**. Does not run hyperparameter tuning (that’s step 3).

---

## hyperparameter_tuning_pipeline.R

**Called when:** Warm-start Step 2 (optional).

Tunes GPR, XGBoost, and GAM using consistent folds and per-model predictor sets. Saves `best_config_gpr.rds`, `best_config_gam.rds`, `best_config_xgb.rds`, and `best_config_lr.rds` in `output/<cv_regime>/cv_pipeline/`. In robust workflows, multiseed tuning scripts under `modelling/multiseed/` are the primary path.

---

## permutation_importance_final.R

**Status:** Deprecated stub.

This file intentionally stops with a deprecation message. Robust importance/pruning is handled by `modelling/multiseed/robust_shap_covariate_pruning.R`.

---

## shap_importance_final.R

**Status:** Deprecated stub.

This file intentionally stops with a deprecation message. Robust SHAP pruning is handled by `modelling/multiseed/robust_shap_covariate_pruning.R`.

---

## fit_final_models.R

**Called when:** Step 5.

Loads the best config per model from **`output/<cv_regime>/cv_pipeline/best_config_*.rds`** and the per-model predictor sets from **`output/<cv_regime>/covariate_selection/`** (SHAP or permutation). Fits XGB, GAM, and GPR on **all** training data with those settings and saves **`XGB_final.rds`**, **`GAM_final.rds`**, **`GPR_final.rds`** to **`output/<cv_regime>/final_models/`**. Each RDS contains the fitted model, `predictor_vars`, `scale_params`, `encoding`, `hyperparams`, and `train_metrics`. These objects are used by the prediction maps and partial dependence scripts.

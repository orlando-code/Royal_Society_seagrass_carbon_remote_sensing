# Seagrass carbon mapping

Predict seagrass carbon density from remote sensing and environmental data (gap-filling within the study region). The workflow covers exploration (CV, variable selection), model choice (GPR hyperparameter tuning, categorical and spatial investigations), prediction (GPR maps, uncertainty, partial dependence), and supplement figures (region outlines, spatial collages).

**Core question**: Estimate carbon stocks of seagrass beds from environmental and remote-sensing covariates. The model is validated for gap-filling near training data; extrapolation beyond the study region should be contextualised by covariate shift / applicability domain.

---

## Repository structure

- **`data/`** – Input data and caches  
  - `all_extracted_new.rds`: core-level extracted covariates (built by `build_all_extracted_new.R` if missing).  
  - `prediction_grid_cache_*.rds`: prediction grids (cached by `gpr_predictions.R`).  
  - `VRE/data/env_rasters/`, `MEOW/`: raster and region data.

- **`modelling/`** – Analysis and figure-generation scripts  
  - **Driver**: `run_paper_figures.R` (at top level) runs the full pipeline (see below).  
  - **`modelling/R/`** – Shared utilities: `helpers.R`, `extract_covariates_from_rasters.R`, `plot_config.R`, `gpr_funs.R`, `assign_region_from_latlon.R`, `interpolation_funs.R`.  
  - **`modelling/pipeline/`** – Core pipeline: `build_all_extracted_new.R`, `cv_pipeline.R`, `model_permutation_pruning.R`, `gpr_hyperparameter_tuning.R`.  
  - **`modelling/plots/`** – Plot-only: `cv_pipeline_plots.R`, `model_permutation_pruning_plots.R`, `gpr_feature_importance_plot.R`, `gpr_region_collages.R`, `plot_region_outlines.R`.  
  - **`modelling/gpr/`** – GPR investigations and predictions: `gpr_categorical_investigation.R`, `gpr_spatial_investigation.R`, `gpr_parameter_search_summary.R`, `gpr_predictions.R`.

- **`figures/`** – All outputs (single canonical location)  
  - **`figures/cv_pipeline_output/`** – CV results, permutation pruning CSV, GPR best config, tuning and importance plots, parameter search summary.  
  - **`figures/interpolation/`** – GPR prediction maps, PDPs, permutation importance, prediction grid cache, region collage.  
  - **`figures/supplement/`** – e.g. `region_shapes.png`.  
  - **`figures/deprecated/`** – Old outputs and covariate-selection artefacts (kept for reference, not used by the pipeline).

- **`modelling/deprecated/`** – Scripts and docs no longer part of the main pipeline (diagnostics, alternative pruning, one-off investigations). Nothing here is run by `run_paper_figures.R`.

---

## Reproducing paper and supplement figures

From the project root in R:

```r
setwd("/path/to/seagrass")   # or use here::here()
source("modelling/run_paper_figures.R")
```

Or from the shell:

```bash
cd /path/to/seagrass
Rscript modelling/run_paper_figures.R
```

The driver runs, in order:

1. **Data** – Build `data/all_extracted_new.rds` if missing.  
2. **CV pipeline** – Spatial folds, model comparison, results → `figures/cv_pipeline_output/`.  
3. **Model permutation pruning** – Per-model predictor selection → `figures/cv_pipeline_output/cv_model_permutation_pruning_*.csv` (canonical); copies to `figures/covariate_selection/pruned_variables_to_include.csv` (GAM) and `pruned_variables_to_include_gpr.csv` (GPR).  
4. **GPR hyperparameter tuning** – Best kernel and config → `figures/cv_pipeline_output/gpr_best_config.rds`.  
5. **GPR investigations** – Categorical (species/region) and spatial (lat/lon, region) comparisons, then parameter search summary plot.  
6. **Plot scripts** – CV figures, permutation importance (combined + GPR-only for paper).  
7. **Best GPR feature importance** – Single plot from the fitted best model.  
8. **GPR predictions** – Maps, PDPs, optional importance; saves `gpr_prediction_grid.rds` for collages.  
9. **Region collage** – Spatial prediction panels per region → `figures/interpolation/gpr_region_collage.png`.  
10. **Supplement** – Region outlines → `figures/supplement/region_shapes.png`.

**Variable selection (step 2)**: The canonical output is `figures/cv_pipeline_output/cv_model_permutation_pruning_*.csv`. It has one row per (model, variable) with a `keep` column: variables are kept until cumulative permutation importance reaches `permutation_keep_frac` (default 0.95). So each model gets its own best feature set; the *number* of variables kept can be small (e.g. 2 for GPR) when importance is concentrated in a few predictors. Do not use `figures/deprecated/covariate_selection/` — that folder is legacy; the active copies are under `figures/covariate_selection/` and are written from the permutation CSV after step 2.

**Cached data**: The pipeline reuses cache where possible (e.g. `spatial_folds_cache.rds`, `gpr_tuning_results_cache.rds`, `gpr_prediction_grid.rds`, prediction grid RDS files in `data/`). Delete the relevant file(s) in `figures/cv_pipeline_output/` or `data/` to force recomputation.

**Run steps individually**: You can run scripts in the same order as above. Prerequisites:  
- Steps 2–3 require `data/all_extracted_new.rds`.  
- Step 3 requires Step 2 (spatial folds).  
- Steps 4–5 require Step 3 (pruning CSV).  
- Steps 6–7 require Steps 2 and 4 (CV results and GPR config).  
- Steps 8–9 require Step 4 (and Step 7 for the grid used by the collage).

---

## Main outputs (paths)

| Output | Path |
|--------|------|
| CV results (tables) | `figures/cv_pipeline_output/cv_results_*.csv`, `cv_predictions.csv` |
| Pruning (per-model predictor sets) | `figures/cv_pipeline_output/cv_model_permutation_pruning_*.csv` (canonical); `figures/covariate_selection/pruned_variables_to_include*.csv` (GAM and GPR copies) |
| GPR best config | `figures/cv_pipeline_output/gpr_best_config.rds` |
| CV and tuning plots | `figures/cv_pipeline_output/cv_*.png`, `gpr_hyperparameter_tuning_*.png`, `gpr_parameter_search_summary.png` |
| Importance plots | `figures/cv_pipeline_output/model_permutation_pruning_importance*.png`, `gpr_feature_importance_best_model.png` |
| GPR maps and PDPs | `figures/interpolation/gpr_*.png`, `gpr_prediction_grid.rds` |
| Region collage | `figures/interpolation/gpr_region_collage.png` |
| Supplement map | `figures/supplement/region_shapes.png` |

---

## Models and methods (short)

- **CV**: Random split, distance-buffered, and spatial block (e.g. 1 km) strategies; multiple models (GAM, GPR, RF, XGBoost, etc.).  
- **Variable selection**: Model-wise permutation importance (spatial CV); retain predictors up to a cumulative importance fraction (e.g. 95%).  
- **GPR**: Best kernel from tuning (1 km spatial CV); optional categorical (species, region) and spatial (lon/lat) terms; maps and uncertainty from the fitted model.

Earlier notes on depth component (median only), variable reduction, and regional structure are in the first lines of this README and in `modelling/deprecated/docs/` for reference.

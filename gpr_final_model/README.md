# Final GPR model for seagrass carbon density

This directory contains the optimal Gaussian Process Regression (GPR) model for predicting median seagrass carbon density from environmental and remote-sensing covariates.

## Prerequisites

- Run the main pipeline (via `run_paper_figures.R`) at least once so that the following exist:
  - `data/all_extracted_new.rds` (core-level extracted covariates)
  - `figures/cv_pipeline_output/gpr_best_config.rds` (best kernel from tuning)
  - `figures/cv_pipeline_output/cv_model_permutation_pruning_spatial_block_1000m.csv` (or equivalent) for the GPR predictor set

- R packages: `here`, `readr`, `dplyr`, `GauPro`, `ggplot2`, `viridisLite`

## Quick start

From the **project root**:

```r
setwd("/path/to/seagrass")   # or setwd(here::here())
source("gpr_final_model/run_final_gpr.R")
```

Or from the shell:

```bash
cd /path/to/seagrass
Rscript gpr_final_model/run_final_gpr.R
```

The script will:

1. Load data and the best GPR configuration (nugget, kernel, predictor set).
2. Fit the GPR on the full dataset and save a deployable object to `gpr_final_model/gpr_final_model.rds`.
3. Run a random train–test split (default 80/20), report R² and RMSE on the test set.
4. Build an **N×N spatial grid** over the dataset’s longitude–latitude extent. At each grid point, environmental covariates are extracted from the same raster dataset as the training data. The GPR predicts mean and standard error on this grid to produce a **spatial prediction map**. Default N = 100 per axis (10,000 cells); set `n_grid_side <- 5000` for 25M cells. Outputs are saved as RDS and optional PNG maps.

## Using the fitted model

### Load the saved model

```r
gpr <- readRDS("gpr_final_model/gpr_final_model.rds")
```

The object `gpr` is a list with:

- `model`: fitted `GauPro` model (from `GauPro::gpkm()`).
- `scale_params`: list with `x_center` and `x_scale` (training summary); required to scale new data before prediction.
- `predictor_vars`: character vector of predictor names (order matters for prediction).
- `value_var`: response name (`"median_carbon_density_100cm"`).
- `kernel`: kernel type used (e.g. `"matern52"`).
- `min.nug`: minimum nugget to use.

### Predict on new data

New data must be a data frame with the same predictor columns as in `gpr$predictor_vars`. Use the provided wrapper so that scaling is applied correctly:

```r
source("modelling/R/gpr_funs.R")   # for gpr_apply_scale
preds <- predict_gpr(gpr, newdata = my_new_data, se = TRUE)
# preds$mean: predicted median carbon density
# preds$se:   standard error of the prediction (if se = TRUE)
```

Implementation of `predict_gpr` (also in `run_final_gpr.R`):

```r
predict_gpr <- function(gpr, newdata, se = TRUE) {
  X <- newdata[, gpr$predictor_vars, drop = FALSE]
  X_scaled <- gpr_apply_scale(X, gpr$scale_params)
  out <- gpr$model$pred(X_scaled, se.fit = se, return_df = TRUE)
  list(mean = as.numeric(out$mean), se = if (se) as.numeric(out$se) else NULL)
}
```

### Fitting on your own data

To re-fit the same model specification on a different dataset:

1. Prepare a data frame with columns: `gpr$value_var` and the same columns as in `gpr$predictor_vars`.
2. Call `fit_gaussian_process_regression()` from `modelling/R/gpr_funs.R` with that data, the same `predictor_vars` and `kernel`, and a prediction grid (e.g. a subset of your data or a dummy grid) to obtain a new `gpr`-style list.

## Outputs

| Output | Description |
|--------|-------------|
| `gpr_final_model/gpr_final_model.rds` | Fitted model object (model, scale_params, predictor_vars, etc.) for loading and calling `predict_gpr()`. |
| `gpr_final_model/train_test_metrics.csv` | R² and RMSE from the random train–test split. |
| `gpr_final_model/grid_predictions.rds` | Spatial grid (longitude, latitude, env covariates, `gpr_mean`, `gpr_se`). |
| `gpr_final_model/grid_plot_mean.png`, `grid_plot_se.png` | Spatial maps of predicted mean and standard error. Optional. |
| `gpr_final_model/spatial_grid_n*.rds` | Cached N×N grid with extracted covariates (reused if present). |

## Options

At the top of `run_final_gpr.R` you can set:

- `train_frac`: fraction of data used for training in the train–test split (default 0.8).
- `n_grid_side`: number of points per axis for the N×N spatial grid (default 100 → 10k cells; 5000 → 25M cells).
- `save_grid_plots`: whether to save the spatial prediction maps (mean and SE).

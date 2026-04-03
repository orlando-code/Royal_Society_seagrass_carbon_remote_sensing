# Seagrass carbon mapping

Predict seagrass carbon density from remote sensing and environmental data (gap-filling within the study region). The workflow covers exploration (CV, variable selection), model choice (GPR, GAM, XGB with hyperparameter tuning), importance (permutation and SHAP), prediction (maps, uncertainty where possible, partial dependence), and supplement figures.

**Core question**: Estimate species-specific carbon stocks of seagrass beds from environmental and remote-sensing covariates. The models are validated for gap-filling near training data. 'Near' is by default set as at least 5km from the training core data to represent local (e.g. further up the coast) but be removed from at least some of the training environmental covariates.

---

## Quick start

From the project root:

```r
# Robust multiseed pipeline (recommended publication workflow)
source("modelling/run_multiseed_pixel_grouped.R")
```

Everything goes under `output/`: cache, covariate selection, robust tuning/evaluation, final models, prediction maps, diagnostics, and supplement. Shared caches (folds and prediction grids) live in `output/cache/`.

---

## Config (`pipeline_config.R`)

Pipeline settings are centralized in `modelling/config/pipeline_config.R` (`get_pipeline_config()`), then consumed by `modelling/run_multiseed_pixel_grouped.R`.

### Plotting / reporting

- `dpi`: Output figure resolution.
- `show_titles`: If `TRUE`, include plot titles/subtitles (useful for interactive runs); set `FALSE` for cleaner figure panels.

### Pipeline toggles

- `do_warm_start`: optional baseline pruning+tuning initialization.
- `do_shap_refined_tuning`: re-tune after robust SHAP pruning.
- `do_sensitivity`: run sensitivity suite and associated plots.
- `do_diagnostics`: run train/test fraction and variance diagnostics.
- `do_fit_final_models`, `do_supplement`: final model fitting and supplement outputs.

### Target and transforms

- `target_var`: Response column in `data/all_extracted_new.rds` (default `median_carbon_density_100cm`).
- `log_transform_target`: If `TRUE`, models are fit on `log(y)` and predictions are back-transformed for metrics.

### Spatial filtering

- `exclude_regions`: Region names to remove from all modelling stages; use `character(0)` to include all regions:

```r
exclude_regions <- c("Black Sea")   # exclude Black Sea
exclude_regions <- character(0)     # include all regions
```

This applies to covariate pruning, CV, tuning, final fits, and prediction maps.

### Covariate pruning / selection

- `use_correlation_filter`: If `TRUE`, drops highly correlated covariates before model-specific selection.
- `correlation_filter_threshold`: Absolute correlation threshold for pruning.
- `permutation_max_vars`: Maximum number of covariates retained after permutation selection.
- `n_permutations`: Replicates per variable when computing permutation importance.
- `permutation_coverage`: Cumulative importance coverage target.
- `use_shap_per_model`: Prefer per-model SHAP-selected covariate sets where available.

### Models

- `model_list`: models to run (default `c("GPR", "GAM", "XGB", "LR")`).

### Cross-validation design

- `n_folds`: Number of CV folds.
- `cv_type`: One of `"random"`, `"location_grouped"`, `"pixel_grouped"`, or `"spatial"` (block CV).
  - Use `"spatial"` to test extrapolation to new geography away from sampled cores; use `"pixel_grouped"` to avoid leakage from coarse raster pixels (recommended for tuning/pruning/importance).
- `cv_blocksize`: Spatial block size (metres) for single-block steps.
- `cv_blocksize_scan`: Vector of block sizes (metres) for CV comparison scans.

---

## Pipeline order (`run_multiseed_pixel_grouped.R`)


| Step | What it does | Outputs |
| ---- | ------------ | ------- |
| `-1` | Load config and create directories | `output/<cv_regime>/`, `output/cache/`, run metadata |
| `0` | Build `data/all_extracted_new.rds` if missing | `data/all_extracted_new.rds` |
| `1-2` | Optional warm-start baseline pruning+tuning | baseline pruned vars and `best_config_*.rds` |
| `3` | Robust multiseed hyperparameter tuning | robust tuning dirs under `output/<cv_regime>/cv_pipeline/` |
| `4` | Robust SHAP covariate pruning + SHAP plots | `output/<cv_regime>/covariate_selection/robust_pixel_grouped/` |
| `5` | Optional robust re-tuning with SHAP-selected covariates | updated robust configs |
| `6` | Robust held-out evaluation on disjoint eval seeds | run-specific evaluation folder in `cv_pipeline/` |
| `7` | Sensitivity analyses and plots (optional) | `.../sensitivity_suite/` |
| `8` | Train/test fraction and target-variance diagnostics (optional) | diagnostics CSV/plots |
| `9` | Fit final models from robust configs and robust covariates | `output/<cv_regime>/final_models/` |
| `10` | Benchmark testing against species-mean baseline | `output/<cv_regime>/analysis/` |
| `11-13` | Spatial maps, partial dependence, supplement | `predictions/`, `covariate_selection/`, `output/supplement/` |


---

## Directory structure

```text
seagrass_mapping/
├── data/                         # Input data and build artefacts
│   ├── all_extracted_new.rds     # Main extracted dataset (built by pipeline step 0)
│   ├── env_rasters/              # NetCDF rasters (download from Zenodo repository: see below)
│   ├── ICES_ecoregions/          # Shapefiles from which to assign new points without regions (download online: see below)
│   └── MEOW/                     # MEOW shapefile for region assignment (download online: see below)
├── modelling/
│   ├── run_multiseed_pixel_grouped.R  # Main driver for publication workflow
│   ├── config/                     # Central pipeline config
│   ├── multiseed/                  # Robust multiseed tuning/pruning/eval
│   ├── analysis/                   # Model comparison, sensitivity, diagnostics
│   ├── pipeline/                 # Data build, pruning, CV, tuning, importance, final fits
│   ├── plots/                    # Figures: PDPs, prediction maps, supplement
│   ├── R/                        # Shared R helpers, ML, raster extraction, plot config
├── output/                       # All pipeline outputs (replaces old figures/)
│   ├── cache/                    # Shared cached spatial folds and prediction grids
│   ├── <cv_regime>/             # Regime-specific outputs (based on `cv_type`)
│   │   ├── covariate_selection/ # Pruning results, importance, PDPs
│   │   ├── cv_pipeline/        # CV results, tuning configs, importance plots
│   │   ├── final_models/      # XGB_final.rds, GAM_final.rds, GPR_final.rds
│   │   ├── predictions/       # Prediction maps (+ SE where available)
│   └── supplement/              # Shared supplement outputs
└── README.md
```

`<cv_regime>` is controlled by `cv_regime_name` in `modelling/config/pipeline_config.R`:
- `"random"` -> `output/random`
- `"location_grouped"` -> `output/location_grouped`
- `"pixel_grouped"` -> `output/pixel_grouped`
- `"spatial"` -> `output/spatial_<cv_blocksize>m`

See `modelling/pipeline/README.md`, `modelling/plots/README.md`, `modelling/multiseed/README.md`, and `modelling/R/README.md` for per-directory details.

---

## Data flow

This is a short overview: see the README files in the relevant directories for more information.

- `pipeline_config.R` defines defaults and seed policy.
- `run_multiseed_pixel_grouped.R` orchestrates build -> robust selection -> robust evaluation -> final outputs.
- Effective run settings are written to `output/<cv_regime>/run_metadata/pipeline_config_effective.rds`.
- Final fitted models are saved to `output/<cv_regime>/final_models/` and reused by map/PDP scripts.

---

## Caching

The pipeline reuses caches where possible. All cache files live under `output/cache/` (spatial fold `.rds` files and prediction grid caches). Delete the relevant file to force recomputation of:

- `data/all_extracted_new.rds` - rebuild extracted data (step 0).
- `output/cache/*_folds.rds` - recompute spatial (or random) folds for CV/tuning/importance.
- `output/cache/prediction_grid_cache_*.rds` - rebuild prediction grids.
- `output/<cv_regime>/cv_pipeline/best_config_*.rds` - re-run hyperparameter tuning.

---

## Models and methods

- **Models**: GAM, GPR, XGBoost. All use the same per-model covariate set from pruning (permutation or SHAP). Where necessary, categorical variables (e.g. seagrass species, region) are encoded as integers.
- **CV**: pixel-grouped multiseed CV is the default publication workflow; core controls are in `pipeline_config.R`.
- **Variable selection**: Correlation filter plus per-model permutation importance (and optionally SHAP); top vars per model are written to `pruned_model_variables_perm.csv` / `pruned_model_variables_shap.csv`.

## Cross-Validation Approaches

The pipeline evaluates predictive performance under multiple fold construction strategies (see `modelling/pipeline/cv_pipeline.R`). This matters because the data are spatially clustered and some observations share identical `(longitude, latitude)` raster-derived covariates.

- `random_split`: Naive split across rows; train/test may contain points from the same location (and duplicate covariates can leak), so results are a highly optimistic baseline for out-of-sample prediction.
- `location_grouped_random`: All rows that share the exact same `(longitude, latitude)` are assigned to the same fold; this prevents leakage from duplicate coordinates but not from distinct locations that fall in the same coarse raster pixel.
- `pixel_grouped_random`: All rows with identical raster-derived covariate vectors are assigned to the same fold. This is strictly more conservative than `location_grouped_random`: distinct locations that map to the same raster pixel (common at ~4 km resolution) are also grouped, so the model is never evaluated on inputs it has seen verbatim during training. Used as the primary fold strategy for hyperparameter tuning, covariate pruning, and importance estimation.
- `spatial_block_<m>m`: The study area is partitioned into spatial blocks of size `<m>` metres; folds hold out whole blocks, testing extrapolation to new geography away from sampled cores and indicating how well covariates capture spatial structure at that scale.
- `region_stratified_<m>m`: Spatial blocks are built within each ecoregion `region`, and then folds are combined so each fold retains representation from multiple regions; this reduces “leave-one-region-out” artifacts caused by how blocks fall across coastal regions.

All fold types are cached/reused via `output/cache/` to make re-runs fast.

## Environmental data

First, ensure that a `data` directory exists in the main repository (see [Directory structure](#directory-structure)).

The NetCDF raster files containing environmental covariates (from remote sensing and re-analysis products) are not stored in this repository. Download them from the data archive associated with the paper (e.g. the Zenodo record referenced in the manuscript) and place all `.nc` files under `data/env_rasters/`. The pipeline will auto-discover these covariates at runtime using `raster_covariates` from `modelling/R/extract_covariates_from_rasters.R`.

---

## Reproduce results

1. Place required inputs under `data/` (`all_extracted_new.rds` can be built by step 0, rasters and region shapefiles must be downloaded as described below).
2. Run from the repository root:
   - `Rscript modelling/run_multiseed_pixel_grouped.R`
3. Record and archive:
   - `output/<cv_regime>/run_metadata/pipeline_config_effective.rds`
   - run-specific robust evaluation folder under `output/<cv_regime>/cv_pipeline/`.

---

## Assumptions and limitations

- This workflow is designed for gap-filling near sampled conditions; extrapolation to novel environments may degrade.
- Reported performance depends on fold construction and seed policy; robust multiseed evaluation is used to reduce split-variance artifacts.
- Data and large generated outputs are not versioned in git; reproducibility requires consistent external data placement and environment setup.

## Regions data

The following directories must be downloaded, unzipped (if necessary), and copied under the 'data' repository into directories titled `ICES_ecoregions` and `MEOW` respectively.

### ICES Ecoregions

The International Council for the exploration of the Sea provides ecoregion shapefiles for the ocean around Europe. This data can be accessed at [this page](https://gis.ices.dk/geonetwork/srv/api/records/4745e824-a612-4a1f-bc56-b540772166eb) via [this link](https://gis.ices.dk/shapefiles/ICES_ecoregions.zip).

Persistent identifier: [https://gis.ices.dk/geonetwork/srv/metadata/4745e824-a612-4a1f-bc56-b540772166eb](https://gis.ices.dk/geonetwork/srv/metadata/4745e824-a612-4a1f-bc56-b540772166eb)

### MEOW Ecoregions

The Marine Ecoregions of the World (MEOW) document global ecoregions. These are available via UNEP via [this link](https://data-gis.unep-wcmc.org/portal/home/item.html?id=80567b4443f4457b822f645a2f0d70cf#:~:text=Description-,Download%20Dataset,-This%20dataset%20combines).

Persistent identifier: [https://data-gis.unep-wcmc.org/server/rest/services/Hosted/WCMC036_MEOW_PPOW_2007_2012/FeatureServer](https://data-gis.unep-wcmc.org/server/rest/services/Hosted/WCMC036_MEOW_PPOW_2007_2012/FeatureServer)
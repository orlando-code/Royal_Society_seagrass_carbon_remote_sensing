# Seagrass carbon mapping

Predict seagrass carbon density from remote sensing and environmental data (gap-filling within the study region). The workflow covers exploration (CV, variable selection), model choice (GPR, GAM, XGB with hyperparameter tuning), importance (permutation and SHAP), prediction (maps, uncertainty where possible, partial dependence), and supplement figures.

**Core question**: Estimate species-specific carbon stocks of seagrass beds from environmental and remote-sensing covariates. The models are validated for gap-filling near training data. 'Near' is by default set as at least 5km from the training core data to represent local (e.g. further up the coast) but be removed from at least some of the training environmental covariates.

---

## Quick start

From the project root:

```r
# Full paper pipeline (all figures and outputs)
source("modelling/run_paper.R")
```

Everything goes under `**output/**`: cache, covariate selection, CV/tuning results, final models, prediction maps, and supplement. Caches (spatial folds, prediction grids) live in `**output/cache/**` so heavy steps can be skipped on re-runs.

---

## Config (run_paper.R)

All pipeline settings live at the top of `modelling/run_paper.R` (Step **-1**). These are the main knobs you’ll most commonly change.

### Plotting / reporting

- `**dpi`**: Output figure resolution (default `150`).
- `**show_titles`**: If `TRUE`, include plot titles/subtitles (useful for interactive runs); set `FALSE` for cleaner figure panels.

### What to run

- `**do_cv_on_defaults**`: If `TRUE` (default), runs the *default* CV pass (**Step 2**) and the default spatial/categorical investigation (**Step 2b**). Set to `FALSE` to skip those and only run cross-fold validation post model tuning.

### Target and transforms

- `**target_var`**: The response column in `data/all_extracted_new.rds` to model (default `median_carbon_density_100cm`).
- `**log_transform_target*`*: If `TRUE`, models are fit on \log(y) and predictions are back-transformed for metrics. This often stabilises variance for right-skewed carbon density.

### Spatial filtering

- `**exclude_regions**`: Character vector of region names to remove from *all* modelling stages (pruning, CV, tuning, final fits, prediction maps). Use `character(0)` to include everything ie:

```r
exclude_regions <- c("Black Sea")   # exclude Black Sea
exclude_regions <- character(0)     # include all regions
```

This applies to covariate pruning, CV, tuning, final fits, and prediction maps.

### Covariate pruning / selection

- `**use_correlation_filter**`: If `TRUE`, drops highly correlated covariates before model-specific selection.
- `**correlation_filter_threshold**`: Absolute correlation threshold for pruning (e.g. `0.8`).
- `**permutation_max_vars**`: Maximum number of covariates retained after permutation selection (e.g. `15L`).
- `**n_permutations**`: Replicates per variable when computing permutation importance (increase for more stable rankings).
- `**permutation_coverage**`: Cumulative importance coverage target (e.g. keep enough variables to explain `0.99` of total importance).
- `**use_shap_per_model**`: If `TRUE`, prefers per-model SHAP-selected covariate sets where available; otherwise uses permutation-selected sets.

### Models

- `**model_list**`: Which models to run. Defaults to `c("GPR", "GAM", "XGB")`.

### Cross-validation design

- `**n_folds**`: Number of CV folds (default `5L`).
- `**cv_type**`: One of `"random"`, `"location_grouped"`, `"pixel_grouped"`, or `"spatial"` (block CV).
  - Use `"spatial"` to test extrapolation to new geography away from sampled cores; use `"pixel_grouped"` to avoid leakage from coarse raster pixels (recommended for tuning/pruning/importance).
- `**cv_blocksize**`: Spatial block size (metres) for tuning and other “single block size” steps (default `5000L`).
- `**cv_blocksize_scan**`: Vector of block sizes (metres) to scan in the CV comparison plots (Step 2). Set to `NULL` or `integer(0)` to avoid scanning and use only `cv_blocksize`.

---

## Pipeline order (run_paper.R)


| Step   | What it does                                                                                       | Outputs                                                      |
| ------ | -------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| **-1** | Config (exclude_regions, model_list, pruning flags, n_folds, cv_type, etc.) and create output dirs | `output/` (and shared subdirs under `output/cache/`)       |
| **0**  | Build `data/all_extracted_new.rds` from raw rasters if missing                                     | `data/all_extracted_new.rds`                                 |
| **1**  | Covariate pruning: correlation filter + per-model permutation (and optionally SHAP) importance     | `output/<cv_regime>/covariate_selection/`                   |
| **2**  | CV pipeline: spatial (or random) folds, run GPR/GAM/XGB, compare block sizes                       | `output/<cv_regime>/cv_pipeline/`                           |
| **2b** | Spatial/categorical effect: env-only vs +lat vs +lon vs +region for each model                     | `output/<cv_regime>/cv_pipeline/`                           |
| **3**  | Hyperparameter tuning for GPR, GAM, XGB (best config per model)                                    | `output/<cv_regime>/cv_pipeline/best_config_*.rds`       |
| **3b** | Permutation importance with tuned models and best vars                                             | `output/<cv_regime>/cv_pipeline/importance_perm_*.csv/.png` |
| **3c** | SHAP importance with tuned models and best vars                                                    | `output/<cv_regime>/cv_pipeline/importance_shap_*.csv/.png` |
| **4**  | Cross-validation with tuned models and pruned covariates (second pass of `cv_pipeline.R`)          | `output/<cv_regime>/cv_pipeline/`                           |
| **4b** | Spatial/categorical effect: env-only vs +lat vs +lon vs +region for each model                     | `output/<cv_regime>/cv_pipeline/`                           |
| **5**  | Fit and save final models on all data (XGB, GAM, GPR)                                              | `output/<cv_regime>/final_models/*_final.rds`             |
| **6**  | Partial dependence plots for each final model                                                      | `output/<cv_regime>/covariate_selection/pdp_*.png`        |
| **7**  | Spatial prediction maps (mean + SE for GPR) for each model                                         | `output/<cv_regime>/predictions/*_prediction_map.png`, `{gpr,gam}_se_map.png` |
| **8**  | Supplement: region outlines, target histograms, correlation, similarity                            | `output/supplement/`                                       |


---

## Directory structure

```
seagrass/
├── data/                         # Input data and build artefacts
│   ├── all_extracted_new.rds     # Main extracted dataset (built by pipeline step 0)
│   ├── env_rasters/              # NetCDF rasters (download from Zenodo repository: see below)
│   ├── ICES_ecoregions/          # Shapefiles from which to assign new points without regions (download online: see below)
│   └── MEOW/                     # MEOW shapefile for region assignment (download online: see below)
├── modelling/
│   ├── run_paper.R       # Main driver – run this for the full pipeline
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

`<cv_regime>` is derived from `cv_type` in `modelling/run_paper.R`:
- `"random"` -> `output/random`
- `"location_grouped"` -> `output/location_grouped`
- `"pixel_grouped"` -> `output/pixel_grouped`
- `"spatial"` -> `output/spatial_<cv_blocksize>m`

See `**modelling/pipeline/README.md**`, `**modelling/plots/README.md**`, and `**modelling/R/README.md**` for a file-by-file description of each directory.

---

## Data flow

This is a short overview: see the README files in the relevant directories for more information.

Flowchart of repository workflow

---

## Caching

The pipeline reuses caches where possible. All cache files live under `**output/cache/**` (spatial fold `.rds` files and prediction grid caches). Delete the relevant file to force recomputation of:

- `**data/all_extracted_new.rds**` – rebuild extracted data (step 0).
- `**output/cache/*_folds.rds**` – recompute spatial (or random) folds for CV/tuning/importance.
- `**output/cache/prediction_grid_cache_*.rds**` – rebuild prediction grid used in supplement and spatial prediction maps.
- `**output/<cv_regime>/cv_pipeline/best_config_*.rds**` – re-run hyperparameter tuning (step 3).

---

## Models and methods

- **Models**: GAM, GPR, XGBoost. All use the same per-model covariate set from pruning (permutation or SHAP). Where necessary, categorical variables (e.g. seagrass species, region) are encoded as integers.
- **CV**: Spatial block CV (with configurable block size(s)) or random split; `**n_folds`** and `**cv_type`** are set in `**run_paper.R`**.
- **Variable selection**: Correlation filter plus per-model permutation importance (and optionally SHAP); top vars per model are written to `**pruned_model_variables_perm.csv`** / `**pruned_model_variables_shap.csv`**.

## Cross-Validation Approaches

The pipeline evaluates predictive performance under multiple fold construction strategies (see `modelling/pipeline/cv_pipeline.R`). This matters because the data are spatially clustered and some observations share identical `(longitude, latitude)` raster-derived covariates.

- `random_split`: Naive split across rows; train/test may contain points from the same location (and duplicate covariates can leak), so results are a highly optimistic baseline for out-of-sample prediction.
- `location_grouped_random`: All rows that share the exact same `(longitude, latitude)` are assigned to the same fold; this prevents leakage from duplicate coordinates but not from distinct locations that fall in the same coarse raster pixel.
- `pixel_grouped_random`: All rows with identical raster-derived covariate vectors are assigned to the same fold. This is strictly more conservative than `location_grouped_random`: distinct locations that map to the same raster pixel (common at ~4 km resolution) are also grouped, so the model is never evaluated on inputs it has seen verbatim during training. Used as the primary fold strategy for hyperparameter tuning, covariate pruning, and importance estimation.
- `spatial_block_<m>m`: The study area is partitioned into spatial blocks of size `<m>` metres; folds hold out whole blocks, testing extrapolation to new geography away from sampled cores and indicating how well covariates capture spatial structure at that scale.
- `region_stratified_<m>m`: Spatial blocks are built within each ecoregion `region`, and then folds are combined so each fold retains representation from multiple regions; this reduces “leave-one-region-out” artifacts caused by how blocks fall across coastal regions.

All fold types are cached/reused via `**output/cache/`** to make re-runs fast.

## Environmental data

First, ensure that a `data` directory exists in the main repository (see [Directory structure](#directory-structure)).

The NetCDF raster files containing environmental covariates (from remote sensing and re-analysis products) are not stored in this repository. Download them from the data archive associated with the paper (e.g. the Zenodo record referenced in the manuscript) and place all `.nc` files under `data/env_rasters/`. The pipeline will auto-discover these covariates at runtime using `raster_covariates` from `modelling/R/extract_covariates_from_rasters.R`.

## Regions data

The following directories must be downloaded, unzipped (if necessary), and copied under the 'data' repository into directories titled `ICES_ecoregions` and `MEOW` respectively.

### ICES Ecoregions

The International Council for the exploration of the Sea provides ecoregion shapefiles for the ocean around Europe. This data can be accessed at [this page](https://gis.ices.dk/geonetwork/srv/api/records/4745e824-a612-4a1f-bc56-b540772166eb) via [this link](https://gis.ices.dk/shapefiles/ICES_ecoregions.zip).

Persistent identifier: [https://gis.ices.dk/geonetwork/srv/metadata/4745e824-a612-4a1f-bc56-b540772166eb](https://gis.ices.dk/geonetwork/srv/metadata/4745e824-a612-4a1f-bc56-b540772166eb)

### MEOW Ecoregions

The Marine Ecoregions of the World (MEOW) document global ecoregions. These are available via UNEP via [this link](https://data-gis.unep-wcmc.org/portal/home/item.html?id=80567b4443f4457b822f645a2f0d70cf#:~:text=Description-,Download%20Dataset,-This%20dataset%20combines).

Persistent identifier: [https://data-gis.unep-wcmc.org/server/rest/services/Hosted/WCMC036_MEOW_PPOW_2007_2012/FeatureServer](https://data-gis.unep-wcmc.org/server/rest/services/Hosted/WCMC036_MEOW_PPOW_2007_2012/FeatureServer)
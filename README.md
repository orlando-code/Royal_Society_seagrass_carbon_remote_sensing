# Seagrass carbon mapping

This repository predicts carbon density of sediment cores on seagrass beds from remote sensing and environmental data. The workflow covers exploration (CV, variable selection); model choice – Generalised Additive Models (GAMs), Linear Regressors (LRs), Gaussian Process Regressors (GPRs), and XGBoost regressors (XGB) – with hyperparameter tuning; SHAP feature importance, prediction (maps, uncertainty, partial dependence), and supplement figures to accompany the paper, *Oceanographic drivers of carbon storage in European seagrass beds*, Gallo and Timmerman et al. (2026).

---

## Quick start

From the project root:

```r
# 0) Instantiate environment (required)
renv::restore(prompt = FALSE)

# 1) Optional (expensive): robust tuning seed sweep
source("modelling/analysis/run_tuning_seed_sweep.R")

# 2) Main robust multiseed pipeline (recommended publication workflow)
source("modelling/run_multiseed_pixel_grouped.R")
```

Generated files are stored under `output/`. The robust multiseed driver writes a **run-scoped** folder named from the selection seeds, e.g. `output/pixel_grouped_48-52-53-70-73/` (see `robust_fold_seed_list` / `seed_registry` in `pipeline_config.R`). That folder holds covariate selection, robust tuning/evaluation, diagnostics, and run metadata. Shared regime outputs (e.g. `output/pixel_grouped/covariate_selection/`) still use `cv_output_dir` from `cv_regime_name`. Cached artefacts shared across runs (e.g. prediction grids, fold caches) go in `output/cache/`.

---

## Config (`pipeline_config.R`)

Pipeline settings are centralized in `modelling/pipeline_config.R` (`get_pipeline_config()`), then consumed by `modelling/run_multiseed_pixel_grouped.R` and subsidiary scripts.

### Plotting / reporting

- `dpi` – Output figure resolution.
- `show_titles` – If `TRUE`, include plot titles/subtitles (useful for interactive runs); set `FALSE` for paper-ready figure panels.

### Pipeline toggles

- `do_shap_refined_tuning` – re-tune hyperparameters after robust SHAP pruning (multiseed driver Step 3).
- `do_tuning_seed_sweep_refined_tuning` – second robust tuning pass after SHAP pruning inside **standalone** `run_tuning_seed_sweep.R` subsets; set `TRUE` to mirror multiseed Step 3.
- `do_sensitivity` – run sensitivity suite and associated plots.
- `do_diagnostics` – run train/test fraction and variance diagnostics.
- `do_fit_final_models`, `do_supplement` – final model fitting and additional supplementary figures.
- `multiseed_run_output_id` – retained for compatibility; the multiseed driver’s primary output path is `output/pixel_grouped_<robust_seeds>/` derived from `robust_fold_seed_list`, not this id.
- `use_robust_seeds_from_tuning_sweep` – if `TRUE`, load `robust_fold_seed_list` (and `eval_fold_seed_list` if stored) from `output/tuning_seed_sweep_runs/chosen_seeds_latest.rds` after running `modelling/analysis/run_tuning_seed_sweep.R`. The tuning-seed sweep itself writes under `output/tuning_seed_sweep_runs/sweep_<run_id>/` (see that script’s header). Reusing a sweep folder requires both the same **planned subset manifest** and the same **effective sweep config** on disk; otherwise a new `sweep_<id>` is created.

### Target and transforms

- `target_var` – Response column in `data/all_extracted_new.rds` (default `median_carbon_density_100cm`).
- `log_transform_target` – If `TRUE`, models are fit on `log(y)`while predictions are back-transformed for reporting and metrics.

### Spatial filtering

- `exclude_regions` – Region names to remove from all modelling stages. Use `character(0)` to include all regions:

```r
exclude_regions <- c("Black Sea")   # exclude Black Sea
exclude_regions <- character(0)     # include all regions
```

This applies to covariate pruning, CV, tuning, final fits, and prediction maps.

### Covariate pruning / selection

- `use_correlation_filter` – If `TRUE`, drops highly correlated covariates (default `r > 0.8`) before model-specific selection.
- `correlation_filter_threshold` – Absolute correlation threshold for pruning (default `r > 0.8`).
- `permutation_max_vars` – Maximum number of covariates retained after permutation selection (default `15`).
- `n_permutations` – Replicates per variable when computing permutation importance (default `1` since very computationally intensive).
- `permutation_coverage` – Cumulative importance coverage target (default `0.99`).
- `use_shap_per_model` – Prefer per-model SHAP-selected covariate sets where available (default `TRUE`).

### Models

- `model_list` – models to run (default `c("GPR", "GAM", "XGB", "LR")`).

### Cross-validation design

- `n_folds` – Number of CV folds (default `5L` as chosen by sensitivity investigation).
- `cv_type` – One of `"random"`, `"location_grouped"`, `"pixel_grouped"`, or `"spatial"` (block CV).
  - Use `"spatial"` to test extrapolation to new geography away from sampled cores. `"pixel_grouped"` is the grouping reported in the paper.  It avoids leakage between train and test sets from coarse raster pixels.
- `cv_blocksize` – Spatial block size (metres) for single-block steps. Only relevant when `cv_type` is `"spatial".`
- `cv_blocksize_scan` – Vector of block sizes (metres) for CV comparison scans. Only relevant when `cv_type` is `"spatial".`

---

## Pipeline order (`run_multiseed_pixel_grouped.R`)


![Flowchart](report/flowchart.png)

Primary outputs for a run live under **`output/pixel_grouped_<robust_fold_seeds>/`** (hyphen-separated seed string). Shared baseline paths still use **`output/<cv_regime>/`** from `cv_regime_name` / `cv_output_dir`.

| Step    | What it does                                                   | Typical outputs (run folder unless noted)                                        |
| ------- | -------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| `-1`    | Load config; assign `run_output_dir`; write run metadata       | `output/pixel_grouped_<seeds>/_run_metadata/`, `output/<cv_regime>/`, `output/cache/` |
| `0`     | Build `data/all_extracted_new.rds` if missing                | `data/all_extracted_new.rds`                                                     |
| `1`     | Robust hyperparameter tuning (initial covariates)            | `.../cv_pipeline/robust_pixel_grouped_tuning/` (`best_config_*_robust.rds`, …)   |
| `2`     | Robust SHAP pruning + Step 2b SHAP importance plots            | `.../covariate_selection/robust_pixel_grouped/`                                  |
| `3`     | Optional re-tuning on SHAP-pruned covariates (`do_shap_refined_tuning`) | updates `.../cv_pipeline/robust_pixel_grouped_tuning/`                           |
| `4`     | Robust evaluation on `eval_fold_seed_list`                   | `across_seeds_summary.csv` and related CSVs under `output/pixel_grouped_<seeds>/` (run root) |
| `5`     | R² sensitivity suite + plots (optional, `do_sensitivity`)    | `.../sensitivity_suite/`                                                         |
| `6`     | Train/test fraction diagnostics (optional, `do_diagnostics`)   | diagnostics under run folder                                                     |
| `7`     | Fit final models (`modelling/R/fit_final_models.R`)          | `output/<cv_regime>/final_models/`                                               |
| `8`     | Model comparison vs baselines                                  | `output/<cv_regime>/analysis/`                                                   |
| `9`     | Spatial prediction maps                                        | `output/<cv_regime>/predictions/`                                                |
| `10`    | Partial dependence (final models)                              | `output/<cv_regime>/covariate_selection/` (PDP figures)                          |
| `11`    | Supplement figures                                             | `output/supplement/`                                                             |


---

## Directory structure

```text
seagrass_mapping/
├── data/                         # Input data and build artefacts
│   ├── all_extracted_new.rds     # Main extracted dataset (built by pipeline step 0)
│   ├── covariate_rasters/              # NetCDF rasters (download from Zenodo repository: see below)
│   ├── ICES_ecoregions/          # Shapefiles from which to assign new points without regions (download online: see below)
│   └── MEOW/                     # MEOW shapefile for region assignment (download online: see below)
├── modelling/
│   ├── run_multiseed_pixel_grouped.R  # Main driver for publication workflow
│   ├── pipeline_config.R         # Central pipeline + seed registry
│   ├── multiseed/                # Robust multiseed tuning/pruning/eval
│   ├── analysis/                 # Model comparison, sensitivity, tuning seed sweep, diagnostics
│   ├── pipeline/                 # Legacy/deprecation notes + historical stubs
│   ├── plots/                    # Figures: PDPs, prediction maps, supplement, plot_config.R
│   ├── R/                        # Shared R helpers, ML, raster extraction
├── output/                       # All pipeline outputs (replaces old figures/)
│   ├── cache/                    # Shared cached spatial folds and prediction grids
│   ├── tuning_seed_sweep_runs/   # Standalone sweep: sweep_<id>/, chosen_seeds_latest.rds
│   ├── pixel_grouped_<seeds>/    # Multiseed run folder (robust_fold_seed_list in the name)
│   │   ├── covariate_selection/  # Run-scoped robust SHAP outputs
│   │   ├── cv_pipeline/          # Run-scoped robust_pixel_grouped_tuning/
│   │   ├── _run_metadata/        # pipeline_config_effective.rds for that run
│   │   └── …                     # evaluation, sensitivity_suite, diagnostics, etc.
│   ├── <cv_regime>/              # Shared regime tree (`cv_output_dir`, e.g. output/pixel_grouped/)
│   │   ├── covariate_selection/  # Baseline pruning, PDPs from driver Step 10
│   │   ├── cv_pipeline/          # Legacy/baseline tuning paths
│   │   ├── final_models/         # XGB_final.rds, GAM_final.rds, GPR_final.rds
│   │   ├── predictions/          # Prediction maps (+ SE where available)
│   │   └── analysis/             # Model comparison outputs
│   └── supplement/               # Shared supplement outputs
├── report/                       # Figures for report and repository
└── README.md
```

`<cv_regime>` is controlled by `cv_regime_name` in `modelling/pipeline_config.R`:

- `"pixel_grouped"` -> `output/pixel_grouped` (default)
- `"random"` -> `output/random` (naïve random split)
- `"location_grouped"` -> `output/location_grouped` (group by unique latitude/longitude pairs)
- `"spatial"` -> `output/spatial_<cv_blocksize>m` (group by spatial block. Not recommended due to existing spatial clustering and small dataset)

See `modelling/pipeline/_README.md`, `modelling/plots/_README.md`, `modelling/multiseed/_README.md`, `modelling/analysis/_README.md`, and `modelling/R/_README.md` for per-directory details.

---

## Data flow

This is a short overview: see the README files in the relevant directories for more information.

- `pipeline_config.R` defines defaults and seed policy.
- `run_multiseed_pixel_grouped.R` orchestrates build → robust selection (tune → SHAP → optional re-tune) → robust evaluation → sensitivity/diagnostics → final outputs.
- Effective run settings are written to `output/pixel_grouped_<seeds>/_run_metadata/pipeline_config_effective.rds` for each multiseed run.
- Final fitted models are saved to `output/<cv_regime>/final_models/` and reused by map/PDP scripts.

---

## Caching

The pipeline reuses caches where possible. All cache files live under `output/cache/` (spatial fold `.rds` files and prediction grid caches). Delete the relevant file to force recomputation of:

- `data/all_extracted_new.rds` – rebuild extracted data (step 0).
- `output/cache/*_folds.rds` – recompute spatial (or random) folds for CV/tuning/importance.
- `output/cache/prediction_grid_cache_*.rds` – rebuild prediction grids.
- `output/<cv_regime>/cv_pipeline/best_config_*.rds` – re-run hyperparameter tuning.

---

## Models and methods

- **Models** – GAM, GPR, XGBoost. All use the same per-model covariate set from pruning (permutation or SHAP). Where necessary, categorical variables (e.g. seagrass species, region) are encoded as integers.
- **CV** – pixel-grouped multiseed CV is the default publication workflow; core controls are in `pipeline_config.R`.
- **Variable selection** – Correlation filter plus per-model permutation importance (and optionally SHAP); top vars per model are written to `pruned_model_variables_perm.csv` / `pruned_model_variables_shap.csv`.

## Cross-Validation Approaches

The pipeline evaluates predictive performance under multiple fold construction strategies (legacy implementation archived under `deprecated/cv_pipeline.R`). This matters because the data are spatially clustered and some observations share identical `(longitude, latitude)` raster-derived covariates.

- `random_split` – Naive split across rows; train/test may contain points from the same location (and duplicate covariates can leak), so results are a highly optimistic baseline for out-of-sample prediction.
- `location_grouped_random` – All rows that share the exact same `(longitude, latitude)` are assigned to the same fold; this prevents leakage from duplicate coordinates but not from distinct locations that fall in the same coarse raster pixel.
- `pixel_grouped_random` – All rows with identical raster-derived covariate vectors are assigned to the same fold. This is strictly more conservative than `location_grouped_random`: distinct locations that map to the same raster pixel (common at ~4 km resolution) are also grouped, so the model is never evaluated on inputs it has seen verbatim during training. Used as the primary fold strategy for hyperparameter tuning, covariate pruning, and importance estimation.
- `spatial_block_<m>m` – The study area is partitioned into spatial blocks of size `<m>` metres; folds hold out whole blocks, testing extrapolation to new geography away from sampled cores and indicating how well covariates capture spatial structure at that scale.
- `region_stratified_<m>m` – Spatial blocks are built within each ecoregion `region`, and then folds are combined so each fold retains representation from multiple regions; this reduces “leave-one-region-out” artifacts caused by how blocks fall across coastal regions.

All fold types are cached/reused via `output/cache/` to speed up re-runs.

---

## Reproduce results

1. Place required inputs under `data/` (`all_extracted_new.rds` can be built by step 0, rasters and region shapefiles must be downloaded as described below).
2. Run from the repository root:
  - `Rscript -e "renv::restore(prompt = FALSE)"`
  - optional: `Rscript modelling/analysis/run_tuning_seed_sweep.R`
  - `Rscript modelling/run_multiseed_pixel_grouped.R`
3. Record and archive:
  - `output/pixel_grouped_<seeds>/_run_metadata/pipeline_config_effective.rds`
  - the same `output/pixel_grouped_<seeds>/` tree (evaluation CSVs, `sensitivity_suite/` if run, etc.).
4. Once models have been saved and predictions generated, the `prediction_maps.ipynb` Jupyter notebook can be used to generate the predictive maps of carbon stocks. This requires downloading extra datasets (see below).
Seed policy is documented in `modelling/_SEED_REGISTRY.md`.

### Python environment for `prediction_maps.ipynb`

Use a clean virtual environment and install the notebook dependencies from `requirements.txt` via a `bash` terminal:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Environmental data

The NetCDF raster files containing environmental covariates (from remote sensing and re-analysis products) are stored via Zenodo (persistent identifier: [https://doi.org/10.5281/zenodo.19329403](https://doi.org/10.5281/zenodo.19329403)). Download and place all `.nc` files under `data/covariate_rasters/`. The pipeline will auto-discover these covariates at runtime via the `build_covariate_config_from_dir` function in `modelling/R/extract_covariates_from_rasters.R`.

## Regions data

The following directories must be downloaded, unzipped (if necessary), and copied under the 'data' repository into directories titled `ICES_ecoregions` and `MEOW` respectively.

### ICES Ecoregions

The International Council for the exploration of the Sea provides ecoregion shapefiles for the ocean around Europe. This data can be accessed at [this page](https://gis.ices.dk/geonetwork/srv/api/records/4745e824-a612-4a1f-bc56-b540772166eb) via [this link](https://gis.ices.dk/shapefiles/ICES_ecoregions.zip).

Persistent identifier: [https://gis.ices.dk/geonetwork/srv/metadata/4745e824-a612-4a1f-bc56-b540772166eb](https://gis.ices.dk/geonetwork/srv/metadata/4745e824-a612-4a1f-bc56-b540772166eb)

### MEOW Ecoregions

The Marine Ecoregions of the World (MEOW) document global ecoregions. These are available via UNEP via [this link](https://data-gis.unep-wcmc.org/portal/home/item.html?id=80567b4443f4457b822f645a2f0d70cf#:~:text=Description-,Download%20Dataset,-This%20dataset%20combines).

Persistent identifier: [https://data-gis.unep-wcmc.org/server/rest/services/Hosted/WCMC036_MEOW_PPOW_2007_2012/FeatureServer](https://data-gis.unep-wcmc.org/server/rest/services/Hosted/WCMC036_MEOW_PPOW_2007_2012/FeatureServer)

## Prediction Map datasets

### Seagrass Essential Ocean Variable (EOV)

The seagrass cover Essential Ocean Variable (EOV) was obtained from the EMODnet product catalogue from the EMODnet Seabed Habitats project via [this link](https://doi.org/10.34892/x6r3-d211). This allowed the calculation of national carbon stocks in seagrass beds from model-predicted carbon density. The species '`eunis_id`' codes were used to map the reported species to the classes in the core data dataset (see Supplementary Seagrass Essential Ocean Variable (EOV)).

### Exclusive Economic Zones (EEZs)

Exclusive Economic Zones (EEZs) were obtained from the Marine Regions website via [this link](https://www.marineregions.org/downloads.php). We used the `World EEZ v12 (2023-10-25, 122MB)` file in GeoPackage format. EEZs were used to attribute seagrass bed carbon to national territories.

---

## Assumptions and limitations

- This workflow is designed for gap-filling near sampled conditions; extrapolation to novel environments may degrade (as shown when datasets `cv_type = "spatial`).
- Reported performance depends on fold construction and seed policy; robust multiseed evaluation is used to reduce split-variance artifacts.


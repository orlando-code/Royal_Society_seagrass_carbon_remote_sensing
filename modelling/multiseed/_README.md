# Multi-seed robustness scripts

This directory contains the seeded `pixel_grouped` multi-seed workflow used by
the robust pipeline driver script.

Entry point: `modelling/run_multiseed_pixel_grouped.R` (see also `modelling/pipeline_config.R` for seed lists and toggles).

For **`modelling/run_multiseed_pixel_grouped.R`**, `run_output_dir` is **`output/pixel_grouped_<robust_fold_seeds>/`** (seeds joined with `-`, from `robust_fold_seed_list` / registry). Robust tuning, SHAP pruning outputs, evaluation summaries, and (when enabled) **`sensitivity_suite/`** and diagnostics live under that folder (e.g. `covariate_selection/robust_pixel_grouped/`, `cv_pipeline/robust_pixel_grouped_tuning/`). Legacy or other drivers may still use paths under **`output/<cv_regime>/cv_pipeline/`**; see the driver you run for the exact layout.

- `robust_shap_covariate_pruning.R`: robust SHAP-based variable selection across multiple fold seeds.
- `robust_hyperparameter_tuning.R`: robust hyperparameter selection across multiple fold seeds.
- `robust_evaluation.R`: robust evaluation across held-out fold seeds, with across-seed summaries.
- `train_test_fraction_diagnostic.R`: fold-size and target-variance diagnostics joined to robust evaluation results.

Legacy one-off seed-check scripts in `modelling/pipeline/` were removed in favour
of this consolidated directory.

# Multi-seed robustness scripts

This directory contains the seeded `pixel_grouped` multi-seed workflow used by
the robust pipeline driver script.

Entry point: `modelling/run_multiseed_pixel_grouped.R` (see also `modelling/config/pipeline_config.R` for seed lists and toggles).

- `robust_shap_covariate_pruning.R`: robust SHAP-based variable selection across multiple fold seeds.
- `robust_hyperparameter_tuning.R`: robust hyperparameter selection across multiple fold seeds.
- `robust_evaluation.R`: robust evaluation across held-out fold seeds, with across-seed summaries.
- `train_test_fraction_diagnostic.R`: fold-size and target-variance diagnostics joined to robust evaluation results.

Legacy one-off seed-check scripts in `modelling/pipeline/` were removed in favour
of this consolidated directory.

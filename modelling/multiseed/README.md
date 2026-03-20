# Multi-seed robustness scripts

This directory contains the `pixel_grouped_random` multi-seed workflow used by
`modelling/run_robust_pixel_grouped_random.R`.

- `robust_shap_covariate_pruning.R`: robust SHAP-based variable selection across multiple fold seeds.
- `robust_hyperparameter_tuning.R`: robust hyperparameter selection across multiple fold seeds.
- `robust_evaluation.R`: robust evaluation across held-out fold seeds, with across-seed summaries.
- `train_test_fraction_diagnostic.R`: fold-size and target-variance diagnostics joined to robust evaluation results.

Legacy one-off seed-check scripts in `modelling/pipeline/` were removed in favour
of this consolidated directory.

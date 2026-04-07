# Analysis scripts

Diagnostic and synthesis scripts for evaluating robustness, fold behavior, and
communication-ready summaries. Most are run standalone with `Rscript` and use
`modelling/pipeline_config.R` defaults (or `.GlobalEnv` overrides) plus outputs
already produced under `output/<cv_regime>/cv_pipeline/`.

---

## tuning_seed_sweep.R

**Called when:** Ad-hoc robust tuning sensitivity runs (often in `tmux`), optionally before the main driver.

Re-runs robust tuning + robust SHAP pruning + robust evaluation across different counts of robust seeds (for example `2/5/10/15/20`). Writes run-scoped outputs to `output/<cv_regime>/cv_pipeline/tuning_seed_sweep_runs/sweep_<run_id>/` (summaries, manifests, plots, and subset work folders) and updates `chosen_seeds_latest.rds` for optional reuse in `run_multiseed_pixel_grouped.R`.

---

## sensitivity_suite.R

**Called when:** Step 10 (if enabled) from `modelling/run_multiseed_pixel_grouped.R`, or standalone diagnostics.

Runs a fast R2 sensitivity suite with fixed tuned model setups (no re-tuning/re-selection): varies fold count and fold seeds, quantifies train-size and fold-composition effects, and optionally performs the expensive tuning-seed-count sweep block. Writes outputs under `<run_output_dir>/sensitivity_suite/` including per-fold, pooled-by-seed, summary, and seed-convergence/plateau CSVs.

---

## model_comparison.R

**Called when:** Standalone benchmark analysis after tuning/pruning artifacts exist.

Compares each fitted model against a species-mean baseline under CV, reporting fold-wise and pooled performance deltas. Reads tuned covariates and best configs, supports strict alignment to a specific robust evaluation run via `model_comparison_eval_run_dir`, and writes outputs to `output/<cv_regime>/analysis/<run_tag>/`.

---

## fold_sensitivity_check.R

**Called when:** Standalone diagnostic to justify fold design choices.

Checks how performance changes with `n_folds` and fold seed across configured CV types while holding tuned covariates/hyperparameters fixed. Produces by-fold and across-seed summaries plus RMSE/R2 plots in `output/<cv_regime>/cv_pipeline/fold_sensitivity_check/`.

---

## compile_model_performance_summaries.R

**Called when:** Standalone reporting synthesis step after CV and robust evaluation outputs exist.

Builds communication-ready summary tables by combining post-tuning CV summaries with negative-R2 diagnostics, and (when present) robust deployment summaries from the newest robust evaluation directory. Writes `model_performance_across_methods_summary.csv` and `deployment_robustness_pixel_grouped.csv` to `output/<cv_regime>/cv_pipeline/`.

---

## diagnose_r2_components_by_fold.R

**Called when:** Standalone troubleshooting for unstable/negative fold-level R2.

Decomposes fold-level R2 into `SS_res` and `SS_tot` components using `cv_predictions.csv`, quantifies sensitivity to low fold variance, and generates diagnostic scatter plots. Writes `r2_components_by_fold.csv`, `r2_components_summary.csv`, and `r2_components_diagnostics.png` to `output/<cv_regime>/cv_pipeline/`.

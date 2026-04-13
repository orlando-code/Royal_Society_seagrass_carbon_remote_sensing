# Seed Registry and Source of Truth

This project keeps seed policy centralized in `modelling/pipeline_config.R` under `seed_registry`.

## Canonical fields

- `seed_registry$active_robust_fold_seed_list`
  - active robust selection/tuning seeds used by default.
- `seed_registry$paper_robust_fold_seed_list`
  - frozen paper/report seeds; opt-in with `use_paper_seed_registry = TRUE`.
- `seed_registry$eval_fold_seed_list`
  - held-out evaluation seeds, kept disjoint from the tuning pool.

Derived fields used by scripts:

- `robust_fold_seed_list`
- `eval_fold_seed_list`

These are computed by `get_pipeline_config()` from `seed_registry` so downstream scripts can keep using existing names.

## Pipeline handoff from tuning sweep

Optional sweep -> robust run flow:

1. `Rscript modelling/analysis/tuning_seed_sweep.R`
2. `Rscript modelling/run_multiseed_pixel_grouped.R`

The sweep writes:

- `output/tuning_seed_sweep_runs/sweep_<run_id>/chosen_seeds_for_pipeline.rds`
- `output/tuning_seed_sweep_runs/chosen_seeds_latest.rds` (copy for the main driver)

If `use_robust_seeds_from_tuning_sweep = TRUE`, the robust pipeline reads
`output/tuning_seed_sweep_runs/chosen_seeds_latest.rds` and temporarily overrides robust/eval seeds for that run.

Sweep outputs (subset work, summaries, `_run_metadata/pipeline_config_effective.rds`) live under **`output/tuning_seed_sweep_runs/sweep_<run_id>/`**. Reusing an existing `sweep_<id>` folder requires both the same planned subset manifest and the same effective sweep configuration on disk; otherwise the script allocates the next free sweep id.

## Recommendations

- Keep `use_paper_seed_registry = FALSE` for exploratory reruns.
- Set `use_paper_seed_registry = TRUE` only for paper-grade reproduction.
- Keep all edits to defaults in `seed_registry`, not scattered across scripts.

# Purpose:
#   Re-run robust tuning + robust SHAP pruning + robust evaluation for different
#   robust seed counts (e.g., 2/5/10/15/20), then write outputs under a dedicated
#   run folder:
#   output/tuning_seed_sweep_runs/sweep_<run_id>/
#     - sensitivity_tuning_seed_sweep_summary.csv, manifests, plots, subset_work/<subset_id>/...
#   Sweep directory reuse (same run id) requires both an identical planned subset manifest *and*
#   an identical effective sweep config on disk (_run_metadata/pipeline_config_effective.rds);
#   otherwise a new sweep_<run_id> folder is allocated.
#   Writes chosen_seeds_latest.rds under output/tuning_seed_sweep_runs/ for optional use by
#   run_multiseed_pixel_grouped.R when use_robust_seeds_from_tuning_sweep = TRUE.


if (!exists("seagrass_init_repo", mode = "function", inherits = TRUE)) {
  init_path <- file.path("modelling", "R", "init_repo.R")
  if (!file.exists(init_path)) {
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(ff)) stop("Run from repo root or with: Rscript /path/to/run_tuning_seed_sweep.R", call. = FALSE)
    script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
    init_path <- normalizePath(file.path(dirname(script_path), "..", "R", "init_repo.R"), winslash = "/", mustWork = FALSE)
  }
  if (!file.exists(init_path)) stop("Missing bootstrap helper: modelling/R/init_repo.R", call. = FALSE)
  sys.source(init_path, envir = .GlobalEnv)
}
project_root <- seagrass_init_repo(
  packages = c("here", "dplyr", "readr", "ggplot2", "tidyr"),
  source_files = c("modelling/pipeline_config.R"),
  include_helpers = TRUE,
  require_core_inputs = TRUE,
  check_renv = TRUE
)

cfg <- get_pipeline_config()
apply_pipeline_defaults(
  cfg,
  c(
    "cv_regime_name", "cv_type_label", "tuning_seed_sweep_counts", "tuning_seed_pool",
    "eval_fold_seed_list", "tuning_seed_sampling",
    "do_tuning_seed_sweep_refined_tuning", "robust_pruned_importance_type",
    "tuning_seed_sweep_repeats", "tuning_seed_sweep_random_seed",
    "tuning_seed_sweep_skip_existing", "tuning_seed_sweep_unique_subsets",
    "tuning_seed_sweep_parallel_jobs", "model_list",
    "tuning_seed_sweep_force_recompute",
    "tuning_seed_sweep_run_id"
  ),
  envir = .GlobalEnv
)

cv_regime_name <- get("cv_regime_name", envir = .GlobalEnv)
cv_type_label <- get("cv_type_label", envir = .GlobalEnv)
tuning_seed_sweep_counts <- as.integer(get("tuning_seed_sweep_counts", envir = .GlobalEnv))
tuning_seed_pool <- as.integer(get("tuning_seed_pool", envir = .GlobalEnv))
eval_fold_seed_list <- as.integer(get("eval_fold_seed_list", envir = .GlobalEnv))
tuning_seed_sampling <- match.arg(
  get("tuning_seed_sampling", envir = .GlobalEnv),
  choices = c("prefix", "random")
)
do_tuning_seed_sweep_refined_tuning <- isTRUE(get("do_tuning_seed_sweep_refined_tuning", envir = .GlobalEnv))
robust_pruned_importance_type <- get("robust_pruned_importance_type", envir = .GlobalEnv)
tuning_seed_sweep_repeats <- as.integer(get("tuning_seed_sweep_repeats", envir = .GlobalEnv))
tuning_seed_sweep_random_seed <- as.integer(get("tuning_seed_sweep_random_seed", envir = .GlobalEnv))
tuning_seed_sweep_skip_existing <- isTRUE(get("tuning_seed_sweep_skip_existing", envir = .GlobalEnv))
tuning_seed_sweep_unique_subsets <- isTRUE(get("tuning_seed_sweep_unique_subsets", envir = .GlobalEnv))
tuning_seed_sweep_parallel_jobs <- as.integer(get("tuning_seed_sweep_parallel_jobs", envir = .GlobalEnv))
tuning_seed_sweep_force_recompute <- isTRUE(get("tuning_seed_sweep_force_recompute", envir = .GlobalEnv))
model_list <- intersect(get("model_list", envir = .GlobalEnv), c("GPR", "GAM", "XGB", "LR"))

tuning_sweep_runs_root <- file.path(project_root, "output", "tuning_seed_sweep_runs")
dir.create(tuning_sweep_runs_root, recursive = TRUE, showWarnings = FALSE)

# Keep only feasible counts and in ascending unique order.
tuning_seed_sweep_counts <- sort(unique(tuning_seed_sweep_counts[tuning_seed_sweep_counts >= 1L]))
tuning_seed_sweep_counts <- tuning_seed_sweep_counts[tuning_seed_sweep_counts <= length(tuning_seed_pool)]
if (length(tuning_seed_sweep_counts) == 0L) {
  stop("No valid tuning_seed_sweep_counts after filtering by tuning_seed_pool length.")
}

list_sweep_dirs <- function(root_dir) {
  d <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
  d[grepl("^sweep_", basename(d))]
}

read_registry_file <- function(path) {
  req <- c("n_tuning_seeds", "sweep_repeat_id", "tuning_seed_list")
  if (!file.exists(path)) return(NULL)
  x <- readr::read_csv(path, show_col_types = FALSE)
  if (!all(req %in% names(x))) return(NULL)
  x[, req, drop = FALSE]
}

load_subset_registry <- function() {
  # Prefer the newest existing sweep registry as canonical subset memory.
  sweep_dirs <- list_sweep_dirs(tuning_sweep_runs_root)
  if (length(sweep_dirs) > 0L) {
    reg_candidates <- file.path(sweep_dirs, "subset_registry.csv")
    reg_candidates <- reg_candidates[file.exists(reg_candidates)]
    if (length(reg_candidates) > 0L) {
      reg_candidates <- reg_candidates[order(file.info(reg_candidates)$mtime, decreasing = TRUE)]
      x <- read_registry_file(reg_candidates[[1L]])
      if (!is.null(x)) return(x)
    }
  }
}

# Save and restore globals modified by this sweep.
had_robust <- exists("robust_fold_seed_list", envir = .GlobalEnv, inherits = FALSE)
old_robust <- if (had_robust) get("robust_fold_seed_list", envir = .GlobalEnv, inherits = FALSE) else NULL
had_eval <- exists("eval_fold_seed_list", envir = .GlobalEnv, inherits = FALSE)
old_eval <- if (had_eval) get("eval_fold_seed_list", envir = .GlobalEnv, inherits = FALSE) else NULL
had_cv_type <- exists("cv_type", envir = .GlobalEnv, inherits = FALSE)
old_cv_type <- if (had_cv_type) get("cv_type", envir = .GlobalEnv, inherits = FALSE) else NULL
had_pruned_type <- exists("robust_pruned_importance_type", envir = .GlobalEnv, inherits = FALSE)
old_pruned_type <- if (had_pruned_type) get("robust_pruned_importance_type", envir = .GlobalEnv, inherits = FALSE) else NULL
had_run_output_dir <- exists("run_output_dir", envir = .GlobalEnv, inherits = FALSE)
old_run_output_dir <- if (had_run_output_dir) get("run_output_dir", envir = .GlobalEnv, inherits = FALSE) else NULL
on.exit({
  if (had_robust) {
    assign("robust_fold_seed_list", old_robust, envir = .GlobalEnv)
  } else if (exists("robust_fold_seed_list", envir = .GlobalEnv, inherits = FALSE)) {
    rm("robust_fold_seed_list", envir = .GlobalEnv)
  }

  if (had_eval) {
    assign("eval_fold_seed_list", old_eval, envir = .GlobalEnv)
  } else if (exists("eval_fold_seed_list", envir = .GlobalEnv, inherits = FALSE)) {
    rm("eval_fold_seed_list", envir = .GlobalEnv)
  }

  if (had_cv_type) {
    assign("cv_type", old_cv_type, envir = .GlobalEnv)
  } else if (exists("cv_type", envir = .GlobalEnv, inherits = FALSE)) {
    rm("cv_type", envir = .GlobalEnv)
  }

  if (had_pruned_type) {
    assign("robust_pruned_importance_type", old_pruned_type, envir = .GlobalEnv)
  } else if (exists("robust_pruned_importance_type", envir = .GlobalEnv, inherits = FALSE)) {
    rm("robust_pruned_importance_type", envir = .GlobalEnv)
  }

  if (had_run_output_dir) {
    assign("run_output_dir", old_run_output_dir, envir = .GlobalEnv)
  } else if (exists("run_output_dir", envir = .GlobalEnv, inherits = FALSE)) {
    rm("run_output_dir", envir = .GlobalEnv)
  }
}, add = TRUE)

make_subset_plan <- function(n_sel, repeats, seed_pool, sampling, enforce_unique = TRUE,
                             exclude_keys = character(0)) {
  repeats <- as.integer(max(1L, repeats))
  if (identical(sampling, "prefix")) {
    return(list(seed_pool[seq_len(n_sel)]))
  }

  out <- list()
  seen <- unique(as.character(exclude_keys))
  max_attempts <- max(20L, repeats * 20L)
  attempts <- 0L
  while (length(out) < repeats && attempts < max_attempts) {
    attempts <- attempts + 1L
    s <- sort(as.integer(sample(seed_pool, n_sel, replace = FALSE)))
    key <- paste(s, collapse = "-")
    if (!enforce_unique || !(key %in% seen)) {
      out[[length(out) + 1L]] <- s
      seen <- c(seen, key)
    }
  }
  if (length(out) < repeats) {
    warning(
      "Requested ", repeats, " unique subsets for n=", n_sel,
      " but generated ", length(out), ". Consider larger tuning_seed_pool or fewer repeats."
    )
  }
  out
}

#' One-time seed: copy canonical subsets from an old manifest_run.csv into the registry.
registry_from_manifest_run <- function() {
  sweep_dirs <- list_sweep_dirs(tuning_sweep_runs_root)
  mr_cands <- c(
    file.path(sweep_dirs, "sensitivity_tuning_seed_sweep_manifest_run.csv"),
    file.path(
      project_root, "output", cv_regime_name, "cv_pipeline", "sensitivity_suite",
      "sensitivity_tuning_seed_sweep_manifest_run.csv"
    )
  )
  mr <- mr_cands[file.exists(mr_cands)]
  if (length(mr) == 0L) {
    return(NULL)
  }
  mr <- mr[[1L]]
  x <- readr::read_csv(mr, show_col_types = FALSE)
  req <- c("n_tuning_seeds", "sweep_repeat_id", "tuning_seed_list")
  if (!all(req %in% names(x))) {
    return(NULL)
  }
  x <- x[, req, drop = FALSE]
  x <- dplyr::distinct(x, n_tuning_seeds, sweep_repeat_id, tuning_seed_list, .keep_all = TRUE)
  x <- x[order(x$n_tuning_seeds, x$sweep_repeat_id), , drop = FALSE]
  x
}

merge_registry <- function(old_reg, new_rows) {
  if (is.null(new_rows) || nrow(new_rows) == 0L) {
    return(old_reg)
  }
  if (is.null(old_reg) || nrow(old_reg) == 0L) {
    return(new_rows)
  }
  dplyr::bind_rows(old_reg, new_rows) %>%
    dplyr::distinct(n_tuning_seeds, sweep_repeat_id, .keep_all = TRUE) %>%
    dplyr::arrange(n_tuning_seeds, sweep_repeat_id)
}

subset_registry <- load_subset_registry()
if (is.null(subset_registry)) {
  subset_registry <- registry_from_manifest_run()
  if (!is.null(subset_registry)) {
    cat(
      "Loaded subset registry from existing manifest_run.\n",
      "Future runs can reuse these subsets when tuning_seed_sweep_repeats is reduced.\n",
      sep = ""
    )
  }
}

task_rows <- list()
registry_new_rows <- list()

for (n_sel in tuning_seed_sweep_counts) {
  need <- as.integer(tuning_seed_sweep_repeats)
  reg_n <- if (!is.null(subset_registry)) {
    subset_registry[subset_registry$n_tuning_seeds == n_sel, , drop = FALSE]
  } else {
    data.frame()
  }
  if (nrow(reg_n) > 0L) {
    reg_n <- reg_n[order(reg_n$sweep_repeat_id), , drop = FALSE]
  }

  if (nrow(reg_n) >= need) {
    take <- reg_n[seq_len(need), , drop = FALSE]
    for (i in seq_len(nrow(take))) {
      task_rows[[length(task_rows) + 1L]] <- data.frame(
        n_tuning_seeds = n_sel,
        sweep_repeat_id = as.integer(take$sweep_repeat_id[i]),
        tuning_seed_list = as.character(take$tuning_seed_list[i]),
        stringsAsFactors = FALSE
      )
    }
    next
  }

  if (nrow(reg_n) == 0L) {
    if (identical(tuning_seed_sampling, "random")) {
      set.seed(tuning_seed_sweep_random_seed + n_sel * 100000L)
    }
    subsets_for_n <- make_subset_plan(
      n_sel = n_sel,
      repeats = need,
      seed_pool = tuning_seed_pool,
      sampling = tuning_seed_sampling,
      enforce_unique = tuning_seed_sweep_unique_subsets,
      exclude_keys = character(0)
    )
    for (rep_idx in seq_along(subsets_for_n)) {
      s <- as.integer(subsets_for_n[[rep_idx]])
      tl <- paste(s, collapse = "-")
      task_rows[[length(task_rows) + 1L]] <- data.frame(
        n_tuning_seeds = n_sel,
        sweep_repeat_id = rep_idx,
        tuning_seed_list = tl,
        stringsAsFactors = FALSE
      )
      registry_new_rows[[length(registry_new_rows) + 1L]] <- data.frame(
        n_tuning_seeds = n_sel,
        sweep_repeat_id = rep_idx,
        tuning_seed_list = tl,
        stringsAsFactors = FALSE
      )
    }
    next
  }

  # Top up: keep existing registry rows, add new random subsets with higher sweep_repeat_id
  max_id <- max(reg_n$sweep_repeat_id, na.rm = TRUE)
  for (i in seq_len(nrow(reg_n))) {
    task_rows[[length(task_rows) + 1L]] <- data.frame(
      n_tuning_seeds = n_sel,
      sweep_repeat_id = as.integer(reg_n$sweep_repeat_id[i]),
      tuning_seed_list = as.character(reg_n$tuning_seed_list[i]),
      stringsAsFactors = FALSE
    )
  }
  n_more <- need - nrow(reg_n)
  exclude_keys <- as.character(reg_n$tuning_seed_list)
  if (identical(tuning_seed_sampling, "random")) {
    set.seed(tuning_seed_sweep_random_seed + n_sel * 100000L + max_id)
  }
  extra <- make_subset_plan(
    n_sel = n_sel,
    repeats = n_more,
    seed_pool = tuning_seed_pool,
    sampling = tuning_seed_sampling,
    enforce_unique = tuning_seed_sweep_unique_subsets,
    exclude_keys = exclude_keys
  )
  for (j in seq_along(extra)) {
    s <- as.integer(extra[[j]])
    tl <- paste(s, collapse = "-")
    rid <- max_id + j
    task_rows[[length(task_rows) + 1L]] <- data.frame(
      n_tuning_seeds = n_sel,
      sweep_repeat_id = rid,
      tuning_seed_list = tl,
      stringsAsFactors = FALSE
    )
    registry_new_rows[[length(registry_new_rows) + 1L]] <- data.frame(
      n_tuning_seeds = n_sel,
      sweep_repeat_id = rid,
      tuning_seed_list = tl,
      stringsAsFactors = FALSE
    )
  }
}

task_df <- dplyr::bind_rows(task_rows)

if (nrow(task_df) == 0L) stop("No sweep tasks were generated.")
task_df$sweep_subset_id <- paste0(
  "n", task_df$n_tuning_seeds,
  "_r", task_df$sweep_repeat_id,
  "_", gsub("-", "_", task_df$tuning_seed_list)
)
task_df$sweep_sampling <- tuning_seed_sampling

manifest_core <- function(df) {
  # Match on canonical planning keys; derived columns (subset IDs) can change
  # across script versions while representing the same planned sweep.
  req <- c("n_tuning_seeds", "sweep_repeat_id", "tuning_seed_list")
  out <- df[, req, drop = FALSE]
  out$n_tuning_seeds <- as.integer(out$n_tuning_seeds)
  out$sweep_repeat_id <- as.integer(out$sweep_repeat_id)
  out$tuning_seed_list <- as.character(out$tuning_seed_list)
  out <- out[order(out$n_tuning_seeds, out$sweep_repeat_id, out$tuning_seed_list), , drop = FALSE]
  rownames(out) <- NULL
  out
}

same_manifest_plan <- function(a, b) {
  aa <- as.data.frame(manifest_core(a), stringsAsFactors = FALSE)
  bb <- as.data.frame(manifest_core(b), stringsAsFactors = FALSE)
  rownames(aa) <- NULL
  rownames(bb) <- NULL
  identical(aa, bb)
}

find_matching_sweep_dir <- function(task_df) {
  target <- manifest_core(task_df)
  sweep_dirs <- list_sweep_dirs(tuning_sweep_runs_root)
  if (length(sweep_dirs) == 0L) return(NULL)
  for (d in sweep_dirs) {
    p <- file.path(d, "sensitivity_tuning_seed_sweep_manifest_planned.csv")
    if (!file.exists(p)) next
    old <- tryCatch(readr::read_csv(p, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(old)) next
    req <- c("n_tuning_seeds", "sweep_repeat_id", "tuning_seed_list")
    if (!all(req %in% names(old))) next
    if (same_manifest_plan(old, target)) return(d)
  }
  NULL
}

next_sequential_sweep_id <- function(root_dir, width = 3L) {
  sweep_dirs <- list_sweep_dirs(root_dir)
  ids <- integer()
  if (length(sweep_dirs) > 0L) {
    nm <- basename(sweep_dirs)
    hit <- regexec("^sweep_([0-9]+)$", nm)
    mt <- regmatches(nm, hit)
    ids <- as.integer(vapply(mt, function(x) if (length(x) == 2L) x[[2L]] else NA_character_, character(1)))
    ids <- ids[!is.na(ids)]
  }
  next_id <- if (length(ids) == 0L) 0L else max(ids) + 1L
  sprintf(paste0("%0", width, "d"), next_id)
}

# Effective sweep settings (excluding run folder id) used to decide whether an existing
# sweep_* directory from a manifest match or explicit tuning_seed_sweep_run_id is still valid.
pending_effective_sweep_cfg <- list(
  cv_regime_name = cv_regime_name,
  cv_type_label = cv_type_label,
  tuning_seed_sweep_counts = tuning_seed_sweep_counts,
  tuning_seed_pool = tuning_seed_pool,
  eval_fold_seed_list = eval_fold_seed_list,
  tuning_seed_sampling = tuning_seed_sampling,
  do_tuning_seed_sweep_refined_tuning = do_tuning_seed_sweep_refined_tuning,
  robust_pruned_importance_type = robust_pruned_importance_type,
  tuning_seed_sweep_repeats = tuning_seed_sweep_repeats,
  tuning_seed_sweep_random_seed = tuning_seed_sweep_random_seed,
  tuning_seed_sweep_skip_existing = tuning_seed_sweep_skip_existing,
  tuning_seed_sweep_unique_subsets = tuning_seed_sweep_unique_subsets,
  tuning_seed_sweep_parallel_jobs = tuning_seed_sweep_parallel_jobs,
  model_list = model_list,
  tuning_seed_sweep_force_recompute = tuning_seed_sweep_force_recompute
)

tuning_sweep_dir_saved_config_compatible <- function(sweep_dir, pending_cfg) {
  fp <- file.path(sweep_dir, "_run_metadata", "pipeline_config_effective.rds")
  if (!file.exists(fp)) {
    return(TRUE)
  }
  old <- tryCatch(readRDS(fp), error = function(e) NULL)
  if (is.null(old) || !is.list(old)) {
    return(FALSE)
  }
  ex <- "tuning_seed_sweep_run_id"
  a <- seagrass_strip_config_keys(pending_cfg, exclude_keys = ex)
  b <- seagrass_strip_config_keys(old, exclude_keys = ex)
  isTRUE(all.equal(a, b, check.attributes = FALSE))
}

tuning_seed_sweep_run_id <- get0("tuning_seed_sweep_run_id", envir = .GlobalEnv, ifnotfound = NULL)
explicit_out_dir <- NULL
if (!is.null(tuning_seed_sweep_run_id) && length(tuning_seed_sweep_run_id) == 1L && nzchar(as.character(tuning_seed_sweep_run_id))) {
  explicit_out_dir <- file.path(tuning_sweep_runs_root, paste0("sweep_", as.character(tuning_seed_sweep_run_id)))
}

matching_out_dir <- find_matching_sweep_dir(task_df)
if (!is.null(matching_out_dir) && !tuning_sweep_dir_saved_config_compatible(matching_out_dir, pending_effective_sweep_cfg)) {
  cat(
    "Found identical planned manifest under:\n  ", matching_out_dir,
    "\n  but saved effective sweep config differs from the current pipeline; allocating a new sweep directory.\n",
    sep = ""
  )
  matching_out_dir <- NULL
}

if (!is.null(explicit_out_dir) && dir.exists(explicit_out_dir) &&
    !tuning_sweep_dir_saved_config_compatible(explicit_out_dir, pending_effective_sweep_cfg)) {
  cat(
    "Explicit tuning_seed_sweep_run_id points to:\n  ", explicit_out_dir,
    "\n  but saved effective sweep config there differs from the current pipeline; allocating a new sweep directory (run id will be chosen automatically).\n",
    "  To reuse that path, align pipeline_config with that folder's _run_metadata/pipeline_config_effective.dput or remove the folder.\n",
    sep = ""
  )
  explicit_out_dir <- NULL
}

if (!is.null(matching_out_dir)) {
  sweep_out_dir <- matching_out_dir
  tuning_seed_sweep_run_id <- sub("^sweep_", "", basename(sweep_out_dir))
  cat("Found identical planned manifest and matching effective sweep config. Reusing sweep directory:\n  ", sweep_out_dir, "\n", sep = "")
} else if (!is.null(explicit_out_dir)) {
  sweep_out_dir <- explicit_out_dir
  if (!dir.exists(sweep_out_dir)) dir.create(sweep_out_dir, recursive = TRUE, showWarnings = FALSE)
  tuning_seed_sweep_run_id <- sub("^sweep_", "", basename(sweep_out_dir))
  cat(
    "No identical planned manifest found; using explicit sweep output directory:\n  ",
    sweep_out_dir, "\n",
    sep = ""
  )
} else {
  tuning_seed_sweep_run_id <- next_sequential_sweep_id(tuning_sweep_runs_root, width = 3L)
  sweep_out_dir <- file.path(tuning_sweep_runs_root, paste0("sweep_", tuning_seed_sweep_run_id))
  dir.create(sweep_out_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Starting new sweep directory:\n  ", sweep_out_dir, "\n", sep = "")
}
registry_csv <- file.path(sweep_out_dir, "subset_registry.csv")
manifest_planned_csv <- file.path(sweep_out_dir, "sensitivity_tuning_seed_sweep_manifest_planned.csv")
run_metadata_dir <- file.path(sweep_out_dir, "_run_metadata")
dir.create(run_metadata_dir, recursive = TRUE, showWarnings = FALSE)

effective_sweep_cfg <- utils::modifyList(
  pending_effective_sweep_cfg,
  list(tuning_seed_sweep_run_id = tuning_seed_sweep_run_id)
)
seagrass_confirm_same_config(
  current_cfg = effective_sweep_cfg,
  search_root = file.path(project_root, "output"),
  exclude_keys = c("tuning_seed_sweep_run_id"),
  prompt_prefix = "tuning seed sweep",
  path_regex = "output/tuning_seed_sweep_runs/sweep_[^/]+/_run_metadata/pipeline_config_effective\\.rds$"
)
saveRDS(effective_sweep_cfg, file.path(run_metadata_dir, "pipeline_config_effective.rds"))
dput(effective_sweep_cfg, file = file.path(run_metadata_dir, "pipeline_config_effective.dput"))

writeLines(
  c(
    "Tuning seed-count sweep run.",
    paste("Run id:", tuning_seed_sweep_run_id),
    paste("Created:", format(Sys.time(), tz = "UTC"), "UTC"),
    "Subset evaluations, tuning, and SHAP outputs for each subset live under subset_work/<subset_id>/.",
    "Promoted seeds for the main multiseed pipeline: see chosen_seeds_for_pipeline.rds (and chosen_seeds_latest.rds in the parent tuning_seed_sweep_runs/ folder)."
  ),
  file.path(sweep_out_dir, "_README.txt")
)

cat("Running standalone tuning seed-count sweep\n")
cat("  sweep output dir:", sweep_out_dir, "\n")
cat("  cv_regime_name:", cv_regime_name, "\n")
cat("  tuning_seed_sweep_counts:", paste(tuning_seed_sweep_counts, collapse = ", "), "\n")
cat("  tuning_seed_pool size:", length(tuning_seed_pool), "\n")
cat("  eval_fold_seed_list:", paste(eval_fold_seed_list, collapse = ", "), "\n")
cat("  tuning_seed_sampling:", tuning_seed_sampling, "\n")
cat("  tuning_seed_sweep_repeats:", tuning_seed_sweep_repeats, "\n")
cat("  tuning_seed_sweep_random_seed:", tuning_seed_sweep_random_seed, "\n")
cat("  tuning_seed_sweep_skip_existing:", tuning_seed_sweep_skip_existing, "\n")
cat("  tuning_seed_sweep_unique_subsets:", tuning_seed_sweep_unique_subsets, "\n")
cat("  tuning_seed_sweep_parallel_jobs:", tuning_seed_sweep_parallel_jobs, "\n")
cat("  tuning_seed_sweep_force_recompute:", tuning_seed_sweep_force_recompute, "\n")
cat("  do_tuning_seed_sweep_refined_tuning:", do_tuning_seed_sweep_refined_tuning, "\n")

# Persist registry after selecting output directory.
if (length(registry_new_rows) > 0L) {
  subset_registry <- merge_registry(subset_registry, dplyr::bind_rows(registry_new_rows))
}
if (!is.null(subset_registry) && nrow(subset_registry) > 0L) {
  readr::write_csv(subset_registry, registry_csv)
  if (length(registry_new_rows) > 0L) {
    cat(
      "Wrote subset_registry.csv (", nrow(subset_registry), " row(s)).\n",
      sep = ""
    )
  }
}

readr::write_csv(task_df, manifest_planned_csv)
cat("Wrote planned subset manifest to:\n  ", manifest_planned_csv, "\n", sep = "")

collect_eval_table <- function(eval_summary_csv, seeds_str, n_sel, rep_idx, subset_id, sampling,
                               cov_summary_dir = NULL) {
  eval_tbl <- readr::read_csv(eval_summary_csv, show_col_types = FALSE)
  eval_tbl$n_tuning_seeds <- n_sel
  eval_tbl$tuning_seed_list <- seeds_str
  eval_tbl$sweep_repeat_id <- rep_idx
  eval_tbl$sweep_subset_id <- subset_id
  eval_tbl$sweep_sampling <- sampling

  cov_base <- as.character(cov_summary_dir)
  if (is.null(cov_summary_dir) || !nzchar(cov_base)) {
    stop("collect_eval_table requires cov_summary_dir (subset covariate directory).")
  }

  shap_summary_csv <- file.path(
    cov_base,
    paste0("shap_importance_summary_robust_pixel_grouped_seeds_", seeds_str, ".csv")
  )
  if (file.exists(shap_summary_csv)) {
    shap_tbl <- readr::read_csv(shap_summary_csv, show_col_types = FALSE) %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(
        mean_shap_cv = mean(shap_importance_cv, na.rm = TRUE),
        median_shap_cv = stats::median(shap_importance_cv, na.rm = TRUE),
        .groups = "drop"
      )
    eval_tbl <- eval_tbl %>% dplyr::left_join(shap_tbl, by = "model")
  }

  pruned_csv <- file.path(
    cov_base,
    paste0("pruned_model_variables_shap_robust_pixel_grouped_seeds_", seeds_str, ".csv")
  )
  if (file.exists(pruned_csv)) {
    nvars_tbl <- readr::read_csv(pruned_csv, show_col_types = FALSE) %>%
      dplyr::count(model, name = "n_selected_vars")
    eval_tbl <- eval_tbl %>% dplyr::left_join(nvars_tbl, by = "model")
  }
  eval_tbl
}

run_one_subset <- function(task_row) {
  n_sel <- as.integer(task_row$n_tuning_seeds[[1]])
  rep_idx <- as.integer(task_row$sweep_repeat_id[[1]])
  seeds_str <- as.character(task_row$tuning_seed_list[[1]])
  sel_seeds <- as.integer(strsplit(seeds_str, "-", fixed = TRUE)[[1]])
  subset_id <- paste0("n", n_sel, "_r", rep_idx, "_", gsub("-", "_", seeds_str))
  cat("  Sweep n_tuning_seeds =", n_sel, " | repeat=", rep_idx, " | robust seeds:", seeds_str, "\n")

  subset_work_root <- file.path(sweep_out_dir, "subset_work", subset_id)
  eval_dir <- file.path(subset_work_root, "evaluation")
  cov_dir_subset <- file.path(subset_work_root, "covariates")
  tuning_dir_subset <- file.path(subset_work_root, "tuning")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cov_dir_subset, recursive = TRUE, showWarnings = FALSE)
  dir.create(tuning_dir_subset, recursive = TRUE, showWarnings = FALSE)

  eval_summary_csv <- file.path(eval_dir, "across_seeds_summary.csv")
  tune_cfg_paths <- file.path(
    tuning_dir_subset,
    paste0("best_config_", c("gpr", "gam", "xgb", "lr"), "_robust.rds")
  )
  tune_ready <- all(file.exists(tune_cfg_paths[match(tolower(model_list), c("gpr", "gam", "xgb", "lr"))]))

  shap_pruned_csv <- file.path(
    cov_dir_subset,
    paste0("pruned_model_variables_shap_robust_pixel_grouped_seeds_", seeds_str, ".csv")
  )
  shap_summary_csv <- file.path(
    cov_dir_subset,
    paste0("shap_importance_summary_robust_pixel_grouped_seeds_", seeds_str, ".csv")
  )
  shap_ready <- file.exists(shap_pruned_csv) && file.exists(shap_summary_csv)
  eval_ready <- file.exists(eval_summary_csv)

  determine_stage_mode <- function(eval_ready, tune_ready, shap_ready) {
    if (!isTRUE(tuning_seed_sweep_skip_existing)) return("full")
    # Prefer reusing existing evaluation CSVs when present. Only re-run evaluation
    # when tuning_seed_sweep_force_recompute is TRUE (opt-in stale refresh).
    if (isTRUE(eval_ready) && !isTRUE(tuning_seed_sweep_force_recompute)) {
      return("reuse_eval")
    }
    if (isTRUE(tuning_seed_sweep_force_recompute)) {
      if (tune_ready && shap_ready) return("eval_only")
      if (!tune_ready && shap_ready) return("tune_only")
      if (tune_ready && !shap_ready) return("shap_then_refine")
      return("full")
    }
    if (eval_ready) return("reuse_eval")
    if (tune_ready && shap_ready) return("eval_only")
    if (!tune_ready && shap_ready) return("tune_only")
    if (tune_ready && !shap_ready) return("shap_then_refine")
    "full"
  }

  run_stage <- function(stage_mode) {
    if (identical(stage_mode, "full")) {
      source(file.path(project_root, "modelling/multiseed/robust_hyperparameter_tuning.R"))
      source(file.path(project_root, "modelling/multiseed/robust_shap_covariate_pruning.R"))
      if (isTRUE(do_tuning_seed_sweep_refined_tuning)) {
        source(file.path(project_root, "modelling/multiseed/robust_hyperparameter_tuning.R"))
      }
      source(file.path(project_root, "modelling/multiseed/robust_evaluation.R"))
      return(invisible(NULL))
    }
    if (identical(stage_mode, "tune_only")) {
      source(file.path(project_root, "modelling/multiseed/robust_hyperparameter_tuning.R"))
      source(file.path(project_root, "modelling/multiseed/robust_evaluation.R"))
      return(invisible(NULL))
    }
    if (identical(stage_mode, "shap_then_refine")) {
      source(file.path(project_root, "modelling/multiseed/robust_shap_covariate_pruning.R"))
      if (isTRUE(do_tuning_seed_sweep_refined_tuning)) {
        source(file.path(project_root, "modelling/multiseed/robust_hyperparameter_tuning.R"))
      }
      source(file.path(project_root, "modelling/multiseed/robust_evaluation.R"))
      return(invisible(NULL))
    }
    if (identical(stage_mode, "eval_only")) {
      source(file.path(project_root, "modelling/multiseed/robust_evaluation.R"))
      return(invisible(NULL))
    }
    if (!identical(stage_mode, "reuse_eval")) {
      stop("Unknown stage_mode: ", stage_mode)
    }
    invisible(NULL)
  }

  stage_mode <- determine_stage_mode(
    eval_ready = eval_ready,
    tune_ready = tune_ready,
    shap_ready = shap_ready
  )
  cat("    stage_mode:", stage_mode, "\n")

  if (!identical(stage_mode, "reuse_eval")) {
    assign("robust_fold_seed_list", sel_seeds, envir = .GlobalEnv)
    assign("eval_fold_seed_list", eval_fold_seed_list, envir = .GlobalEnv)
    assign("cv_type", "pixel_grouped", envir = .GlobalEnv)
    assign("robust_pruned_importance_type", robust_pruned_importance_type, envir = .GlobalEnv)
    assign("run_output_dir", eval_dir, envir = .GlobalEnv)

    run_stage(stage_mode)
  } else {
    cat("    Reusing cached evaluation artifacts for seeds ", seeds_str, "\n", sep = "")
  }

  if (!file.exists(eval_summary_csv)) {
    warning("Missing eval summary for subset ", subset_id, ": ", eval_summary_csv)
    return(NULL)
  }
  out <- collect_eval_table(
    eval_summary_csv = eval_summary_csv,
    seeds_str = seeds_str,
    n_sel = n_sel,
    rep_idx = rep_idx,
    subset_id = subset_id,
    sampling = tuning_seed_sampling,
    cov_summary_dir = cov_dir_subset
  )
  out$stage_mode <- stage_mode
  out
}

task_list <- split(task_df, seq_len(nrow(task_df)))
jobs <- max(1L, as.integer(tuning_seed_sweep_parallel_jobs))
if (jobs > 1L && .Platform$OS.type == "unix" && length(task_list) > 1L) {
  jobs <- min(jobs, length(task_list))
  cat("\nRunning subset jobs in parallel with", jobs, "workers (process-level isolation)\n")
  sweep_rows <- parallel::mclapply(task_list, run_one_subset, mc.cores = jobs)
} else {
  if (jobs > 1L && .Platform$OS.type != "unix") {
    cat("\nParallel workers requested but not available on this OS; running sequentially.\n")
  }
  sweep_rows <- lapply(task_list, run_one_subset)
}
sweep_rows <- sweep_rows[!vapply(sweep_rows, is.null, logical(1))]

sweep_status_df <- task_df
sweep_status_df$stage_mode <- NA_character_
sweep_status_df$completed <- FALSE
if (length(sweep_rows) > 0L) {
  done_tbl <- dplyr::bind_rows(sweep_rows) %>%
    dplyr::select(sweep_subset_id, stage_mode) %>%
    dplyr::distinct()
  sweep_status_df <- sweep_status_df %>%
    dplyr::left_join(done_tbl, by = "sweep_subset_id", suffix = c("", "_done")) %>%
    dplyr::mutate(
      stage_mode = dplyr::coalesce(stage_mode_done, stage_mode),
      completed = !is.na(stage_mode)
    ) %>%
    dplyr::select(-stage_mode_done)
}
manifest_run_csv <- file.path(sweep_out_dir, "sensitivity_tuning_seed_sweep_manifest_run.csv")
readr::write_csv(sweep_status_df, manifest_run_csv)
cat("Wrote run subset manifest to:\n  ", manifest_run_csv, "\n", sep = "")

sweep_df <- dplyr::bind_rows(sweep_rows)
if (nrow(sweep_df) > 0L) {
  sweep_csv <- file.path(sweep_out_dir, "sensitivity_tuning_seed_sweep_summary.csv")
  readr::write_csv(sweep_df, sweep_csv)
  cat("Wrote tuning seed sweep summary to:\n  ", sweep_csv, "\n", sep = "")
  assign("plot_tuning_seed_sweep_disable_autorun_on_source", TRUE, envir = .GlobalEnv)
  source(file.path(project_root, "modelling/plots/plot_tuning_seed_sweep.R"))
  if (exists("plot_tuning_seed_sweep_disable_autorun_on_source", envir = .GlobalEnv, inherits = FALSE)) {
    rm("plot_tuning_seed_sweep_disable_autorun_on_source", envir = .GlobalEnv)
  }
  plot_tuning_seed_sweep(sweep_out_dir)
  rep5 <- pick_representative_tuning_seed_subset(
    sweep_df,
    n_tuning_seeds = 5L,
    metric_col = "mean_mean_rmse",
    metric_sd_col = "sd_mean_rmse"
  )
  if (!is.null(rep5)) {
    rep_row <- data.frame(
      robust_fold_seed_list = paste(rep5$robust_fold_seed_list, collapse = "-"),
      tuning_seed_list = rep5$tuning_seed_list,
      n_tuning_seeds = rep5$n_tuning_seeds,
      mean_mean_rmse_across_models = rep5$mean_metric_across_models,
      mean_sd_rmse_across_models = rep5$mean_sd_metric_across_models,
      mean_pooled_rmse_mean_across_models = rep5$mean_metric_across_models,
      median_pooled_rmse_mean_across_subsets_at_n = NA_real_,
      selection_rule = rep5$selection_rule,
      stringsAsFactors = FALSE
    )
    rep_csv <- file.path(sweep_out_dir, "sensitivity_tuning_seed_sweep_representative_n5.csv")
    readr::write_csv(rep_row, rep_csv)
    cat(
      "Representative 5-seed subset (closest to median mean RMSE across subsets):\n  ",
      rep_row$robust_fold_seed_list[1], "\n",
      "Wrote:\n  ", rep_csv, "\n",
      "selection_rule: ", rep_row$selection_rule[1], "\n",
      "Set pipeline_config robust_fold_seed_list to c(",
      paste0(rep5$robust_fold_seed_list, "L", collapse = ", "),
      ") if promoting this sweep.\n",
      sep = ""
    )

    chosen_payload <- list(
      robust_fold_seed_list = as.integer(rep5$robust_fold_seed_list),
      eval_fold_seed_list = as.integer(eval_fold_seed_list),
      tuning_seed_list = rep5$tuning_seed_list,
      n_tuning_seeds = rep5$n_tuning_seeds,
      sweep_run_id = tuning_seed_sweep_run_id,
      sweep_run_dir = sweep_out_dir,
      written_at = Sys.time()
    )
    saveRDS(chosen_payload, file.path(sweep_out_dir, "chosen_seeds_for_pipeline.rds"))
    saveRDS(chosen_payload, file.path(tuning_sweep_runs_root, "chosen_seeds_latest.rds"))
    cat(
      "Wrote chosen seeds for main pipeline:\n  ",
      file.path(sweep_out_dir, "chosen_seeds_for_pipeline.rds"), "\n  ",
      file.path(tuning_sweep_runs_root, "chosen_seeds_latest.rds"),
      "\nSet use_robust_seeds_from_tuning_sweep = TRUE in pipeline_config, or copy integers into robust_fold_seed_list.\n",
      sep = ""
    )
  } else {
    message("No n_tuning_seeds == 5 rows; skipped representative subset file.")
  }
} else {
  warning("Tuning seed sweep produced no rows.")
}

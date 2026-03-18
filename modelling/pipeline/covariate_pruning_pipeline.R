# Covariate Pruning Pipeline
#
# Global variables (set by run_paper.R):
#   use_correlation_filter, correlation_filter_threshold
#   permutation_max_vars   – max vars retained after permutation/SHAP importance
#   permutation_coverage   – cumulative importance coverage threshold
#   model_list             – character vector of model names
#   use_shap_per_model     – if TRUE, also compute per-model SHAP via iml (default FALSE)
#
# Outputs -> output/<cv_regime>/covariate_selection/
#   importance_permutation_<model>.csv   + importance_permutation_combined.csv
#   importance_shap_<model>.csv         + importance_shap_combined.csv (if use_shap_per_model)
#   pruned_model_variables_perm.csv     (combined top vars per model, permutation)
#   pruned_model_variables_shap.csv     (combined top vars per model, SHAP)

source("modelling/R/helpers.R")
source("modelling/R/assign_region_from_latlon.R")
load_packages(c("dplyr", "readr", "blockCV", "sf"))

# ---------------------------------------------------------------------------
# Data and region exclusion
# ---------------------------------------------------------------------------
dat <- readr::read_rds("data/all_extracted_new.rds")

# Exclude regions if specified
exclude_regions <- get0("exclude_regions", envir = .GlobalEnv, ifnotfound = character(0))
if (length(exclude_regions) > 0L) {
  if (!"region" %in% names(dat)) dat <- assign_region_from_latlon(dat)
  dat <- dat[is.na(dat$region) | !dat$region %in% exclude_regions, ]
  cat("Excluded region(s):", paste(exclude_regions, collapse = ", "), "\n")
}

# Fetch global variables with defaults as fallback
target_var           <- get0("target_var",           envir = .GlobalEnv, ifnotfound = "median_carbon_density_100cm")
log_transform_target <- get0("log_transform_target", envir = .GlobalEnv, ifnotfound = TRUE)
log_transform_target <- if (length(log_transform_target) != 1) TRUE else isTRUE(log_transform_target)
model_list           <- get0("model_list",           envir = .GlobalEnv, ifnotfound = c("GPR", "GAM", "XGB"))
# XXX faster debugging: reduce n_folds, n_permutations, permutation_max_vars
n_folds              <- get0("n_folds",               envir = .GlobalEnv, ifnotfound = 5L)
cv_type              <- get0("cv_type",               envir = .GlobalEnv, ifnotfound = "spatial")
cv_blocksize        <- get0("cv_blocksize",         envir = .GlobalEnv, ifnotfound = 5000L)
permutation_coverage <- get0("permutation_coverage",  envir = .GlobalEnv, ifnotfound = 0.99)
permutation_max_vars <- get0("permutation_max_vars",  envir = .GlobalEnv, ifnotfound = 15L)
n_permutations       <- get0("n_permutations",       envir = .GlobalEnv, ifnotfound = 1L) # TODO: increase for paper figures
use_shap_per_model   <- isTRUE(get0("use_shap_per_model", ifnotfound = FALSE))

# All predictor candidates (always defined, even if correlation filter is skipped)
cor_predictor_vars <- setdiff(
  colnames(dat),
  c("latitude", "longitude", "number_id_final_version",
    "seagrass_species", "region", target_var)
)

# ---------------------------------------------------------------------------
# Step i: Correlation filter
# ---------------------------------------------------------------------------
if (isTRUE(get0("use_correlation_filter"))) {
  cat("Pruning by correlation (threshold:", correlation_filter_threshold, ")...\n")
  cor_predictor_vars <- prune_by_correlation(dat, cor_predictor_vars, target_var,
    cor_threshold = correlation_filter_threshold)
  cat("Variable(s) after correlation filter:", length(cor_predictor_vars), "\n")
}


# ---------------------------------------------------------------------------
# Step ii: Per-model importance (permutation + optionally SHAP)
#          Both methods use env vars (cor_predictor_vars) + species as factor when present. Species are not pruned from models.
# ---------------------------------------------------------------------------

# Include species as a factor in the importance step (models fit with species by default)
importance_predictor_vars <- cor_predictor_vars
if ("seagrass_species" %in% names(dat)) {
  importance_predictor_vars <- c(importance_predictor_vars, "seagrass_species")
  cat("Including seagrass_species as an unprunable factor in importance ranking.\n")
}


imp_perm       <- NULL   # list: model -> data.frame(variable, rmse_increase)
imp_shap       <- NULL   # list: model -> data.frame(variable, shap_importance)
perm_shared_vars <- importance_predictor_vars  # fallback to cor_predictor_vars if permutation step is skipped

if (length(importance_predictor_vars) > permutation_max_vars) {
  supported_models <- intersect(model_list, c("RF", "GAM", "XGB", "GPR"))
  if (length(supported_models) == 0L)
    stop("No supported models for importance ranking.")

  # -- Permutation importance: folds from cv_type (run_paper.R) ----------
  cv_fold_info <- make_cv_folds(
    dat, cor_predictor_vars, n_folds, cv_type,
    cv_blocksize = cv_blocksize, exclude_regions = exclude_regions,
    cache_tag = "cov_pruning_perm"
  )
  fold_indices <- cv_fold_info$fold_indices
  cat(sprintf(
    "\nPERMUTATION IMPORTANCE (%d var(s), %d fold(s), %s):\n",
    length(importance_predictor_vars),
    n_folds,
    cv_fold_info$method_name
  ))

  imp_perm_all <- do.call(rbind, lapply(supported_models, function(model_name) {
    cat("\n  - ", model_name, " ... ", sep = "")
    out <- tryCatch({
      df <- permutation_importance_cv(
        dat, target_var, importance_predictor_vars, model_name,
        fold_indices, n_permutations = n_permutations, verbose = FALSE, log_response = log_transform_target
      )
      df$model <- model_name
      df
    }, error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n")
      NULL
    })
    if (is.null(out)) {
      data.frame(variable = character(), rmse_increase = numeric(), model = character(), stringsAsFactors = FALSE)
    } else {
      out
    }
  }))
  if (nrow(imp_perm_all) > 0L)
    imp_perm_all <- imp_perm_all[!is.na(imp_perm_all$variable) & imp_perm_all$variable != "", , drop = FALSE]

  imp_perm <- lapply(supported_models, function(model_name)
    imp_perm_all[imp_perm_all$model == model_name, c("variable", "rmse_increase"), drop = FALSE]
  )
  names(imp_perm) <- supported_models

  for (model_name in supported_models)
    print_importance(imp_perm[[model_name]], "rmse_increase", "Permutation", model_name,
                     max_show = permutation_max_vars)

  # -- SHAP importance (evaluated on the same importance_predictor_vars) --------
  if (use_shap_per_model) {
    cat(sprintf(
      "\nSHAP IMPORTANCE (iml::Shapley, same %d vars as permutation):\n",
      length(importance_predictor_vars)
    ))
    imp_shap <- lapply(supported_models, function(model_name) {
      cat("  - ", model_name, " ...\n", sep = "")
      shap_df <- tryCatch(
        compute_shap_importance(dat, target_var, importance_predictor_vars, model_name, log_response = log_transform_target),
        error = function(e) {
          cat("    ! SHAP error for", model_name, ":", conditionMessage(e), "\n")
          NULL
        }
      )
      if (!is.null(shap_df))
        print_importance(shap_df, "shap_importance", "SHAP", model_name,
                         max_show = permutation_max_vars)
      shap_df
    })
    names(imp_shap) <- supported_models
    imp_shap <- imp_shap[!vapply(imp_shap, is.null, logical(1))]
  }
}

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
covariate_dir <- file.path(get0("cv_output_dir", envir = .GlobalEnv, ifnotfound = "output"), "covariate_selection")
dir.create(covariate_dir, recursive = TRUE, showWarnings = FALSE)

# Per-model importance CSVs (permutation and SHAP)
if (!is.null(imp_perm)) {
  for (imp_perm_names in names(imp_perm)) {
    write.csv(imp_perm[[imp_perm_names]], file.path(covariate_dir, paste0("importance_permutation_", tolower(imp_perm_names), ".csv")), row.names = FALSE)
  }
  p_combined <- combine_pruned_model_variables(covariate_dir, "permutation")
  if (!is.null(p_combined)) cat("Saved importance_perm_<model>.csv +", basename(p_combined), "\n")
}
if (!is.null(imp_shap) && length(imp_shap) > 0L) {
  for (imp_shap_names in names(imp_shap)) {
    write.csv(imp_shap[[imp_shap_names]], file.path(covariate_dir, paste0("importance_shap_", tolower(imp_shap_names), ".csv")), row.names = FALSE)
  }
  s_combined <- combine_pruned_model_variables(covariate_dir, "shap")
  if (!is.null(s_combined)) cat("Saved importance_shap_<model>.csv +", basename(s_combined), "\n")
}

# Select top env vars only (species vars are never pruned and are appended)
select_top_env_then_species <- function(df, value_col, max_vars, coverage) {
  if (is.null(df) || nrow(df) == 0) return(character(0))
  species_idx <- is_species_var(df$variable)
  species_vars <- unique(df$variable[species_idx])
  env_df <- df[!species_idx, , drop = FALSE]
  if (nrow(env_df) == 0) return(species_vars)
  top_env <- select_top_vars(env_df, value_col, max_vars, coverage)
  unique(c(top_env, species_vars))
}

# Combined pruned sets (top env vars per model + all species vars; species is not pruned)
if (!is.null(imp_perm)) {
  perm_pruned <- do.call(rbind, lapply(names(imp_perm), function(imp_perm_names) {
    v <- select_top_env_then_species(imp_perm[[imp_perm_names]], "rmse_increase", permutation_max_vars, permutation_coverage)
    data.frame(model = imp_perm_names, variable = v, stringsAsFactors = FALSE)
  }))
  write.csv(perm_pruned, file.path(covariate_dir, "pruned_model_variables_perm.csv"), row.names = FALSE)
  cat("Saved pruned_model_variables_perm.csv\n")
}
if (!is.null(imp_shap) && length(imp_shap) > 0L) {
  shap_pruned <- do.call(rbind, lapply(names(imp_shap), function(imp_shap_names) {
    v <- select_top_env_then_species(imp_shap[[imp_shap_names]], "shap_importance", permutation_max_vars, permutation_coverage)
    data.frame(model = imp_shap_names, variable = v, stringsAsFactors = FALSE)
  }))
  write.csv(shap_pruned, file.path(covariate_dir, "pruned_model_variables_shap.csv"), row.names = FALSE)
  cat("Saved pruned_model_variables_shap.csv\n")
}

# Load dataframe for permutation and/or SHAP importance combined
if (file.exists(file.path(covariate_dir, "pruned_model_variables_perm.csv"))) {
  perm_pruned_df <- read.csv(file.path(covariate_dir, "pruned_model_variables_perm.csv"), stringsAsFactors = FALSE)
  cat("\nPERMUTATION IMPORTANCE OVERLAPS:\n\n")
  perm_overlaps <- compare_covariates_between_models(perm_pruned_df)
}
if (file.exists(file.path(covariate_dir, "pruned_model_variables_shap.csv"))) {
  shap_pruned_df <- read.csv(file.path(covariate_dir, "pruned_model_variables_shap.csv"), stringsAsFactors = FALSE)
  cat("\nSHAP IMPORTANCE OVERLAPS:\n\n")
  shap_overlaps <- compare_covariates_between_models(shap_pruned_df)
}
# Report on variable overlap between methods: Find variables shared by 3 models in both methods, 2 models, etc.
if (!is.null(perm_pruned_df) && !is.null(shap_pruned_df)) {
  cat("\nVARIABLE OVERLAP BETWEEN METHODS (variables shared by k models in BOTH methods):")
  if ("seagrass_species" %in% names(dat)) cat("\nN.B. Species variables are not pruned and are always included.\n\n")
  
  get_shared_vars <- function(overlaps, k) {
    tab <- overlaps[[as.character(k)]]
    if (!is.null(tab) && nrow(tab) > 0) setNames(tab$models, tab$variable) else list()
  }
  perm_k <- as.numeric(names(perm_overlaps))
  shap_k <- as.numeric(names(shap_overlaps))
  perm_k <- perm_k[is.finite(perm_k)]
  shap_k <- shap_k[is.finite(shap_k)]
  max_k  <- if (length(perm_k) && length(shap_k)) min(max(perm_k), max(shap_k)) else 0L
  for (k in if (max_k >= 2L) seq(max_k, 2L, by = -1L) else integer(0)) {
    perm_vars <- get_shared_vars(perm_overlaps, k)
    shap_vars <- get_shared_vars(shap_overlaps, k)
    shared <- intersect(names(perm_vars), names(shap_vars))
    if (length(shared)) {
      cat(sprintf("Variables shared by %d models in BOTH methods:\n", k))
      for (v in shared)
        cat(sprintf("  %s (PERM: %s; SHAP: %s)\n",
          v, paste(perm_vars[[v]], collapse=", "), paste(shap_vars[[v]], collapse=", ")))
      cat("\n")
    }
  }
}

cat("Covariate pruning complete.\n")

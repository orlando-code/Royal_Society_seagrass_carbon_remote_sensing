## Diagnose *why* R2 swings between CV folds.
##
## R2 decomposition (per fold):
##   SS_res = sum((y - yhat)^2)
##   SS_tot = sum((y - mean(y))^2)
##   R2 = 1 - SS_res / SS_tot
##
## Instability sources:
## - small SS_tot (low y-variance in a fold) -> R2 becomes very sensitive
## - large SS_res / poor predictions -> negative R2
##
## Inputs:
##   output/<cv_regime>/cv_pipeline/cv_predictions.csv
##
## Outputs:
##   output/<cv_regime>/cv_pipeline/r2_components_by_fold.csv
##   output/<cv_regime>/cv_pipeline/r2_components_summary.csv
##   output/<cv_regime>/cv_pipeline/r2_components_diagnostics.png
##
project_root <- here::here()
setwd(project_root)

source(file.path(project_root, "modelling/R/helpers.R"))
load_packages(c("here", "dplyr", "readr", "ggplot2", "tidyr"))

cv_regime_name <- get0("cv_regime_name", envir = .GlobalEnv, ifnotfound = "pixel_grouped")
target_method_pattern <- get0("target_method_pattern", envir = .GlobalEnv, ifnotfound = "pixel_grouped")
target_model_list <- get0("target_model_list", envir = .GlobalEnv, ifnotfound = c("XGB", "GPR", "GAM"))

pred_path <- file.path(project_root, "output", cv_regime_name, "cv_pipeline", "cv_predictions.csv")
stopifnot(file.exists(pred_path))

preds <- read.csv(pred_path, stringsAsFactors = FALSE)
required_cols <- c("method", "fold", "model", "observed", "predicted")
missing_cols <- setdiff(required_cols, names(preds))
if (length(missing_cols) > 0L) stop("Missing cols in cv_predictions.csv: ", paste(missing_cols, collapse = ", "))

preds <- preds[grepl(target_method_pattern, preds$method), , drop = FALSE]
preds <- preds[preds$model %in% target_model_list, , drop = FALSE]

cat("Loaded predictions:", nrow(preds), "rows from", pred_path, "\n")

components_by_fold <- preds %>%
  group_by(method, fold, model) %>%
  summarise(
    n_test = dplyr::n(),
    y_mean = mean(observed, na.rm = TRUE),
    y_sd = sd(observed, na.rm = TRUE),
    ss_res = sum((observed - predicted)^2, na.rm = TRUE),
    ss_tot = sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE),
    rmse = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
    r2 = ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_),
    r2_from_ratio = ifelse(ss_tot > 0, -(ss_res / ss_tot) + 1, NA_real_),
    ratio_ss_res_ss_tot = ifelse(ss_tot > 0, ss_res / ss_tot, NA_real_),
    .groups = "drop"
  )

out_dir <- file.path(project_root, "output", cv_regime_name, "cv_pipeline")
write.csv(components_by_fold, file.path(out_dir, "r2_components_by_fold.csv"), row.names = FALSE)

summary_by_model <- components_by_fold %>%
  group_by(model) %>%
  summarise(
    mean_r2 = mean(r2, na.rm = TRUE),
    sd_r2 = sd(r2, na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_y_sd = mean(y_sd, na.rm = TRUE),
    sd_y_sd = sd(y_sd, na.rm = TRUE),
    mean_ratio = mean(ratio_ss_res_ss_tot, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(summary_by_model, file.path(out_dir, "r2_components_summary.csv"), row.names = FALSE)

## Plots
df <- components_by_fold
p1 <- ggplot(df, aes(x = y_sd, y = r2, colour = model)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_minimal(base_size = 12) +
  labs(title = "R2 vs test-fold y spread (sd)", x = "sd(observed) per fold", y = "R2 per fold")

p2 <- ggplot(df, aes(x = ratio_ss_res_ss_tot, y = r2, colour = model)) +
  geom_point(size = 3, alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(title = "R2 vs SS_res/SS_tot ratio", x = "SS_res / SS_tot per fold", y = "R2 per fold")

p3 <- ggplot(df, aes(x = y_sd, y = rmse, colour = model)) +
  geom_point(size = 3, alpha = 0.85) +
  theme_minimal(base_size = 12) +
  labs(title = "RMSE vs test-fold y spread (sd)", x = "sd(observed) per fold", y = "RMSE per fold")

png_path <- file.path(out_dir, "r2_components_diagnostics.png")
p_combo <- patchwork::wrap_plots(p1, p2, p3, ncol = 2)
ggsave(png_path, p_combo, width = 11, height = 9, dpi = 200)

cat("Wrote:\n",
    file.path(out_dir, "r2_components_by_fold.csv"), "\n",
    file.path(out_dir, "r2_components_summary.csv"), "\n",
    file.path(out_dir, "r2_components_diagnostics.png"), "\n",
    sep = "")


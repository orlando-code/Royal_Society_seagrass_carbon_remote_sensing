# Model-specific, nested covariate selection for spatial CV
#
# For each model, performs greedy forward selection (or backward elimination) within each CV fold.
# Returns the set of covariates that maximize R2 for each model, for a given spatial split.
#
# Usage: source this script and call select_covariates_nested_cv(...)
# TODO: this doesn't have kriging methods
select_covariates_nested_cv <- function(core_data, predictor_vars, model_name, fold_indices, max_vars = 10, direction = "forward") {
    cat("\n[", model_name, "] Starting nested covariate selection with", length(predictor_vars), "candidates, max", max_vars, "variables.\n")
    # core_data: data.frame with all predictors and target
    # predictor_vars: character vector of candidate predictors
    # model_name: string, e.g. "GPR", "RF", etc.
    # fold_indices: vector of fold assignments (same length as nrow(core_data))
    # max_vars: maximum number of variables to select
    # direction: "forward" (default) or "backward"

    library(dplyr)
    library(purrr)

    n_folds <- length(unique(fold_indices))
    selected_vars <- c()
    remaining_vars <- predictor_vars
    best_r2 <- -Inf

    if (direction == "backward") {
        selected_vars <- predictor_vars
        remaining_vars <- c()
    }

    for (step in 1:max_vars) {
        candidates <- if (direction == "forward") remaining_vars else selected_vars
        if (length(candidates) == 0) break
        cat("[", model_name, "] Step", step, "- evaluating", length(candidates), "candidates (selected so far:", paste(selected_vars, collapse = ", "), ")\n")
        r2_by_var <- sapply(candidates, function(var) {
            vars_to_try <- if (direction == "forward") c(selected_vars, var) else setdiff(selected_vars, var)
            r2_folds <- c()
            for (fold in unique(fold_indices)) {
                train_idx <- which(fold_indices != fold)
                test_idx <- which(fold_indices == fold)
                train <- core_data[train_idx, c(vars_to_try, "median_carbon_density_100cm"), drop = FALSE]
                test <- core_data[test_idx, c(vars_to_try, "median_carbon_density_100cm"), drop = FALSE]
                # Fit model (simple version, replace with model-specific fit/predict)
                if (model_name == "RF") {
                    fit <- randomForest::randomForest(median_carbon_density_100cm ~ ., data = train)
                    pred <- predict(fit, test)
                } else if (model_name == "GAM") {
                    # Explicitly build formula with selected variables
                    formula_str <- paste("median_carbon_density_100cm ~", paste(vars_to_try, collapse = " + "))
                    fit <- mgcv::gam(as.formula(formula_str), data = train)
                    pred <- predict(fit, test)
                } else if (model_name == "GPR") {
                    # Use the robust utility from gpr_funs.R
                    if (!exists("fit_gpr_cv")) source("modelling/R/gpr_funs.R")
                    res <- fit_gpr_cv(train, test, value_var = "median_carbon_density_100cm", predictor_vars = vars_to_try)
                    pred <- res$predictions
                    if (!is.null(res$error) && !is.na(res$error)) {
                        cat("[GPR] Warning during fit:", res$error, "\n")
                    }
                } else {
                    return(NA)
                }
                r2 <- 1 - sum((test$median_carbon_density_100cm - pred)^2) / sum((test$median_carbon_density_100cm - mean(train$median_carbon_density_100cm))^2)
                r2_folds <- c(r2_folds, r2)
            }
            mean(r2_folds, na.rm = TRUE)
        })
        best_var <- names(which.max(r2_by_var))
        best_var_r2 <- max(r2_by_var, na.rm = TRUE)
        cat("[", model_name, "] Step", step, "- best candidate:", best_var, "(mean R2:", round(best_var_r2, 4), ")\n")
        if (best_var_r2 > best_r2) {
            best_r2 <- best_var_r2
            if (direction == "forward") {
                selected_vars <- c(selected_vars, best_var)
                remaining_vars <- setdiff(remaining_vars, best_var)
            } else {
                selected_vars <- setdiff(selected_vars, best_var)
            }
        } else {
            cat("[", model_name, "] No further improvement, stopping selection.\n")
            break
        }
    }
    cat("[", model_name, "] Final selected variables:", paste(selected_vars, collapse = ", "), "\n")
    selected_vars
}

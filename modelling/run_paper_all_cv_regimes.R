# Runs the full paper pipeline for multiple CV regimes.
#
# It reuses the default configuration in `modelling/run_paper.R`, but overrides
# `cv_type` and (for spatial) `cv_blocksize` by setting globals in `.GlobalEnv`
# before each run.
#
# Usage:
#   Rscript modelling/run_paper_all_cv_regimes.R

set.seed(42)

cv_runs <- list(
  list(name = "location_grouped",  cv_type = "location_grouped",  cv_blocksize = 1000L),
  list(name = "pixel_grouped",     cv_type = "pixel_grouped",     cv_blocksize = 1000L),
  list(name = "random",            cv_type = "random",            cv_blocksize = 1000L),
  list(name = "spatial_1000m",    cv_type = "spatial",          cv_blocksize = 1000L),
  list(name = "spatial_5000m",    cv_type = "spatial",          cv_blocksize = 5000L)
)

for (run in cv_runs) {
  cat("\n============================================================\n")
  cat("RUNNING CV REGIME:", run$name, "\n")
  cat("  cv_type     =", run$cv_type, "\n")
  cat("  cv_blocksize=", run$cv_blocksize, "\n")
  cat("============================================================\n")

  assign("cv_type", run$cv_type, envir = .GlobalEnv)
  assign("cv_blocksize", run$cv_blocksize, envir = .GlobalEnv)

  source("modelling/run_paper.R", local = FALSE)
}


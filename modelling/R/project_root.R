# Repository root discovery (base R only; safe to sys.source before helpers.R).
# Finds the directory that contains modelling/R/helpers.R.

seagrass_project_root <- function() {
  marker <- file.path("modelling", "R", "helpers.R")
  if (file.exists(marker)) {
    return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
  }
  walk_from <- function(start_dir) {
    d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
    for (.i in 1:40L) {
      if (file.exists(file.path(d, marker))) return(d)
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
    NA_character_
  }
  root <- walk_from(getwd())
  if (!is.na(root)) return(root)
  ca <- commandArgs(trailingOnly = FALSE)
  ff <- grep("^--file=", ca, value = TRUE)
  if (length(ff)) {
    sp <- sub("^--file=", "", ff[[1]])
    if (nzchar(sp)) {
      sd <- tryCatch(
        normalizePath(dirname(sp), winslash = "/", mustWork = FALSE),
        error = function(e) NA_character_
      )
      if (!is.na(sd) && nzchar(sd)) {
        root <- walk_from(sd)
        if (!is.na(root)) return(root)
      }
    }
  }
  if (requireNamespace("here", quietly = TRUE)) {
    hr <- tryCatch(
      normalizePath(here::here(), winslash = "/", mustWork = FALSE),
      error = function(e) NA_character_
    )
    if (!is.na(hr) && nzchar(hr)) {
      root <- walk_from(hr)
      if (!is.na(root)) return(root)
    }
  }
  stop(
    "Cannot find repository root (expected ", marker, ").\n",
    "  cd into the clone, or run: Rscript /path/to/any/script/inside/this/repo.R",
    call. = FALSE
  )
}

seagrass_setwd_project_root <- function() {
  root <- seagrass_project_root()
  setwd(root)
  invisible(root)
}

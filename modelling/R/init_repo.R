# Repository bootstrap helpers for reproducible script startup.

seagrass_source_project_root <- function() {
  if (exists("seagrass_setwd_project_root", mode = "function", inherits = TRUE)) {
    return(invisible(NULL))
  }
  pr <- file.path("modelling", "R", "project_root.R")
  if (!file.exists(pr)) {
    find_bootstrap <- function(start_dir) {
      d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
      for (.i in 1:40L) {
        candidate <- file.path(d, "modelling", "R", "project_root.R")
        if (file.exists(candidate)) return(candidate)
        parent <- dirname(d)
        if (identical(parent, d)) break
        d <- parent
      }
      NA_character_
    }
    starts <- c(getwd())
    ff <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(ff)) {
      script_path <- normalizePath(sub("^--file=", "", ff[[1]]), winslash = "/", mustWork = FALSE)
      starts <- c(starts, dirname(script_path))
    }
    for (st in unique(starts)) {
      candidate <- find_bootstrap(st)
      if (!is.na(candidate)) {
        pr <- candidate
        break
      }
    }
  }
  if (!file.exists(pr)) {
    stop(
      "Missing project bootstrap: modelling/R/project_root.R.\n",
      "Run from the repository root or call scripts with absolute paths.",
      call. = FALSE
    )
  }
  sys.source(pr, envir = .GlobalEnv)
  invisible(NULL)
}

load_packages <- function(packages,
                          optional = NULL,
                          auto_install = identical(
                            tolower(Sys.getenv("SEAGRASS_AUTO_INSTALL_PKGS", "false")),
                            "true"
                          )) {
  loaded <- character()
  missing_required <- character()
  for (pkg in unique(as.character(packages))) {
    if (isNamespaceLoaded(pkg)) next
    is_opt <- pkg %in% optional
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (isTRUE(auto_install)) {
        tryCatch(
          install.packages(pkg, quiet = TRUE),
          error = function(e) if (!is_opt) {
            stop("Failed to install '", pkg, "': ", conditionMessage(e), call. = FALSE)
          }
        )
      } else if (!is_opt) {
        missing_required <- c(missing_required, pkg)
        next
      } else {
        warning("Optional package '", pkg, "' is not installed.")
        next
      }
    }
    if (suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))) {
      loaded <- c(loaded, pkg)
    } else if (is_opt) {
      warning("Optional package '", pkg, "' could not be loaded.")
    } else {
      stop("Required package '", pkg, "' could not be loaded.", call. = FALSE)
    }
  }
  if (length(missing_required) > 0L) {
    missing_required <- unique(missing_required)
    stop(
      "Missing required package(s): ", paste(missing_required, collapse = ", "), "\n",
      "Use renv::restore() for reproducible setup.\n",
      "Or install explicitly: install.packages(c(",
      paste(sprintf("\"%s\"", missing_required), collapse = ", "), ")).\n",
      "To allow runtime installs, set SEAGRASS_AUTO_INSTALL_PKGS=true.",
      call. = FALSE
    )
  }
  if (length(loaded) > 0L) {
    cat("Loaded package(s):", paste(unique(loaded), collapse = ", "), "\n")
  }
  invisible(unique(loaded))
}

seagrass_check_renv <- function(project_root, strict = identical(
  tolower(Sys.getenv("SEAGRASS_STRICT_RENV", "true")),
  "true"
)) {
  lockfile <- file.path(project_root, "renv.lock")
  if (!file.exists(lockfile)) return(invisible(FALSE))
  using_renv_lib <- any(grepl("/renv/library/", normalizePath(.libPaths(), winslash = "/", mustWork = FALSE), fixed = TRUE))
  if (isTRUE(using_renv_lib)) return(invisible(TRUE))
  msg <- paste0(
    "renv.lock is present but renv library is not active.\n",
    "Run once from repo root:\n",
    "  Rscript -e \"renv::restore(prompt = FALSE)\""
  )
  if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg)
  invisible(FALSE)
}

build_core_data <- function(project_root) {
  all_extracted_path <- file.path(project_root, "data", "all_extracted_new.rds")
  if (file.exists(all_extracted_path)) return(invisible(TRUE))

  build_script <- file.path(project_root, "modelling", "R", "build_all_extracted_new.R")
  if (!file.exists(build_script)) return(invisible(FALSE))

  message(
    "Core input missing: ", all_extracted_path, "\n",
    "Attempting to build it now via modelling/R/build_all_extracted_new.R ..."
  )
  build_ok <- tryCatch(
    {
      # The build script sources extract_covariates_from_rasters.R and writes all_extracted_new.rds.
      sys.source(build_script, envir = .GlobalEnv)
      TRUE
    },
    error = function(e) {
      warning(
        "Automatic build attempt for data/all_extracted_new.rds failed: ",
        conditionMessage(e)
      )
      FALSE
    }
  )
  isTRUE(build_ok) && file.exists(all_extracted_path)
}

seagrass_require_core_inputs <- function(project_root) {
  core_data_path <- file.path(project_root, "data", "all_extracted_new.rds")
  raster_dir <- file.path(project_root, "data", "covariate_rasters")
  if (!file.exists(core_data_path)) {
    build_core_data(project_root)
  }
  missing <- character()
  if (!file.exists(core_data_path)) missing <- c(missing, core_data_path)
  if (!dir.exists(raster_dir)) missing <- c(missing, raster_dir)
  if (length(missing) > 0L) {
    stop(
      "Missing required core input(s):\n  ",
      paste(missing, collapse = "\n  "),
      "\nThese are required for pipeline execution.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

seagrass_init_repo <- function(packages = NULL,
                               source_files = character(),
                               include_helpers = TRUE,
                               require_core_inputs = FALSE,
                               check_renv = TRUE) {
  seagrass_source_project_root()
  project_root <- get("seagrass_setwd_project_root", mode = "function")()
  if (isTRUE(check_renv)) seagrass_check_renv(project_root)
  if (isTRUE(include_helpers)) {
    source(file.path(project_root, "modelling", "R", "helpers.R"))
  }
  if (length(source_files) > 0L) {
    for (f in source_files) source(file.path(project_root, f))
  }
  if (length(packages) > 0L) load_packages(packages)
  if (isTRUE(require_core_inputs)) seagrass_require_core_inputs(project_root)
  invisible(project_root)
}

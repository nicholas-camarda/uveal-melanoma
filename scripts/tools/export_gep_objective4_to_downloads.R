#!/usr/bin/env Rscript

# To run:
# Rscript scripts/tools/export_gep_objective4_to_downloads.R

#' Resolve the runtime analysis root for this project
#'
#' @return Character scalar absolute path to the runtime analysis root.
resolve_runtime_analysis_root <- function() {
    output_dir <- Sys.getenv("OUTPUT_DIR", unset = "")
    if (nzchar(trimws(output_dir))) {
        return(normalizePath(path.expand(output_dir), winslash = "/", mustWork = FALSE))
    }

    runtime_root <- Sys.getenv("OCULAR_RUNTIME_ROOT", unset = "")
    if (!nzchar(trimws(runtime_root))) {
        runtime_root <- file.path(path.expand("~/Workspaces"), "uveal-melanoma", "runtime")
    }

    file.path(normalizePath(path.expand(runtime_root), winslash = "/", mustWork = FALSE), "Analysis")
}

#' Return the cohorts that should be exported
#'
#' @return Named character vector mapping short cohort labels to runtime cohort directories.
get_objective4_cohort_map <- function() {
    c(
        full = "uveal_full",
        restricted = "uveal_restricted",
        gksrs = "gksrs"
    )
}

#' Determine the destination subdirectory for an exported artifact
#'
#' @param source_path Character scalar source file path.
#' @return `png` for figures, `diagnostics` for diagnostics workbooks, `excel` for the
#'   other workbooks.
get_objective4_artifact_category <- function(source_path) {
    extension <- tolower(tools::file_ext(source_path))
    if (identical(extension, "png")) {
        return("png")
    }
    if (identical(extension, "xlsx")) {
        return(if (grepl("_diagnostics\\.xlsx$", basename(source_path))) "diagnostics" else "excel")
    }
    stop(sprintf("Unsupported Objective 4 artifact type: %s", source_path), call. = FALSE)
}

#' List the shareable Objective 4 files of one cohort
#'
#' The package is the primary workbooks, the diagnostics workbooks and the figures,
#' the same files the publish artifact registry copies. The diagnostics workbooks list,
#' by patient id, every patient left out of or flagged in an analysis (`Record_Issues`),
#' so the package goes only to people who already have access to the study data.
#'
#' @param objective4_root Character scalar path to the cohort's `04_GEP_Validation` folder.
#' @return Character vector of absolute file paths, sorted.
list_objective4_shareable_files <- function(objective4_root) {
    files <- list.files(objective4_root, pattern = "\\.(xlsx|png)$", recursive = TRUE, full.names = TRUE)
    # Excel owner-lock files (~$name.xlsx) and hidden files are never results.
    sort(files[!grepl("^(~\\$|\\.)", basename(files))])
}

#' Copy shareable Objective 4 artifacts into a destination folder
#'
#' @param source_paths Character vector of source file paths.
#' @param export_root Character scalar target directory.
#' @param cohort_label Character scalar cohort subdirectory name.
#' @param dry_run Logical; when `TRUE`, report intended copies without writing.
#' @return Data frame with `cohort`, `category`, `source_path`, `destination_path`
#'   and `status` (`copied` or `would_copy`). A failed copy stops the export.
copy_objective4_artifacts <- function(source_paths, export_root, cohort_label, dry_run = FALSE) {
    categories <- vapply(source_paths, get_objective4_artifact_category, character(1))
    destination_paths <- file.path(export_root, cohort_label, categories, basename(source_paths))
    if (anyDuplicated(destination_paths) > 0L) {
        stop(sprintf("Objective 4 export has duplicate file names in cohort '%s'.", cohort_label), call. = FALSE)
    }

    if (!isTRUE(dry_run)) {
        for (index in seq_along(source_paths)) {
            dir.create(dirname(destination_paths[[index]]), recursive = TRUE, showWarnings = FALSE)
            if (!file.copy(source_paths[[index]], destination_paths[[index]], overwrite = FALSE)) {
                stop(sprintf("Failed to copy %s", source_paths[[index]]), call. = FALSE)
            }
        }
    }

    data.frame(
        cohort = rep(cohort_label, length(source_paths)),
        category = unname(categories),
        source_path = source_paths,
        destination_path = destination_paths,
        status = rep(if (isTRUE(dry_run)) "would_copy" else "copied", length(source_paths)),
        stringsAsFactors = FALSE
    )
}

#' Export the Objective 4 paper pack to Downloads
#'
#' Copies the primary workbooks, diagnostics workbooks (patient ids) and figures of the
#' full, restricted and GKSRS-only cohorts. A missing Objective 4 folder or an empty folder stops the
#' export, so a pack is never partial.
#'
#' @param destination_dir Optional destination directory. Defaults to a fresh
#'   folder under `~/Downloads`.
#' @param dry_run Logical; when `TRUE`, the script only reports what would be
#'   copied.
#' @return Invisibly returns the copy manifest.
export_gep_objective4_to_downloads <- function(
    destination_dir = NULL,
    dry_run = identical(Sys.getenv("GEP_OBJECTIVE4_DRY_RUN", unset = "0"), "1")
) {
    project_slug <- "uveal-melanoma"
    runtime_analysis_root <- resolve_runtime_analysis_root()

    if (is.null(destination_dir) || !nzchar(trimws(destination_dir))) {
        timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
        destination_dir <- file.path(
            path.expand("~/Downloads"),
            sprintf("%s_gep_objective4_%s", project_slug, timestamp)
        )
    } else {
        destination_dir <- normalizePath(path.expand(destination_dir), winslash = "/", mustWork = FALSE)
    }

    cohort_map <- get_objective4_cohort_map()
    manifest_rows <- lapply(names(cohort_map), function(cohort_label) {
        objective4_root <- file.path(runtime_analysis_root, cohort_map[[cohort_label]], "04_GEP_Validation")
        if (!dir.exists(objective4_root)) {
            stop(sprintf("Objective 4 folder not found for the %s cohort: %s", cohort_label, objective4_root), call. = FALSE)
        }
        source_paths <- list_objective4_shareable_files(objective4_root)
        if (length(source_paths) == 0L) {
            stop(sprintf("No Objective 4 workbooks or figures found for the %s cohort: %s", cohort_label, objective4_root), call. = FALSE)
        }
        copy_objective4_artifacts(
            source_paths = source_paths,
            export_root = destination_dir,
            cohort_label = cohort_label,
            dry_run = dry_run
        )
    })
    manifest <- do.call(rbind, manifest_rows)

    cat(sprintf("Objective 4 export target: %s\n", destination_dir))
    cat(sprintf("Files %s: %d\n", if (isTRUE(dry_run)) "that would be copied" else "copied", nrow(manifest)))
    print(manifest[, c("cohort", "category", "destination_path")], row.names = FALSE)

    invisible(manifest)
}

if (sys.nframe() == 0L) {
    destination_dir <- Sys.getenv("GEP_OBJECTIVE4_EXPORT_DIR", unset = "")
    if (!nzchar(trimws(destination_dir))) {
        destination_dir <- NULL
    }

    export_gep_objective4_to_downloads(destination_dir = destination_dir)
}

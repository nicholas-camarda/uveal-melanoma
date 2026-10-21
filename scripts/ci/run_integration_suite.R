#!/usr/bin/env Rscript

#' Build the child-process command for the actual-data integration suite
#'
#' The integration helper (`tests/integration/helper-bootstrap.R`) reads its
#' inputs only from two explicit environment variables. This function names
#' them for the child process alone, so the caller's shell is never modified.
#'
#' @param raw_data_dir Directory holding the raw source spreadsheet.
#' @param processed_data_dir Directory holding the analytic cohort RDS files.
#' @param repo_root Repository root containing `scripts/ci/run_testthat.R`;
#'   defaults to the checkout found from the working directory.
#' @return A list with `command`, `args` and `env` for `system2()`. `system2()`
#'   composes a shell command line, so paths are shell-quoted (the Project
#'   Vault path contains spaces).
build_integration_suite_command <- function(
    raw_data_dir,
    processed_data_dir,
    repo_root = here::here()
) {
    list(
        command = file.path(R.home("bin"), "Rscript"),
        args = shQuote(c(
            file.path(repo_root, "scripts", "ci", "run_testthat.R"),
            file.path(repo_root, "tests", "integration")
        )),
        env = c(
            paste0("OCULAR_INTEGRATION_RAW_DATA_DIR=", shQuote(raw_data_dir)),
            paste0(
                "OCULAR_INTEGRATION_PROCESSED_DATA_DIR=",
                shQuote(processed_data_dir)
            )
        )
    )
}

#' Run the actual-data integration suite from the canonical data locations
#'
#' Takes RAW_DATA_DIR and PROCESSED_DATA_DIR from the project path
#' configuration (the paths `assert_required_input_paths()` checks), fails
#' immediately if the raw input is missing, then runs `tests/integration`
#' through the failure-propagating `run_testthat.R`. The helper itself is
#' unchanged and still requires both variables, so no fallback exists there.
#'
#' @param raw_data_dir Raw-data directory; defaults to the configured
#'   RAW_DATA_DIR.
#' @param processed_data_dir Processed-data directory; defaults to the
#'   configured PROCESSED_DATA_DIR.
#' @param repo_root Repository root used to locate the runner and tests.
#' @return Invisibly returns `TRUE` when every integration test passes;
#'   otherwise throws an error.
run_integration_suite <- function(
    raw_data_dir = NULL,
    processed_data_dir = NULL,
    repo_root = here::here()
) {
    if (is.null(raw_data_dir) || is.null(processed_data_dir)) {
        config_env <- new.env(parent = globalenv())
        for (config_file in c("project_paths.R", "data_processing_policy.R")) {
            sys.source(
                file.path(repo_root, "scripts", "config", config_file),
                envir = config_env
            )
        }
        # Same check the pipeline performs before reading the raw spreadsheet.
        config_env$assert_required_input_paths()
        if (is.null(raw_data_dir)) {
            raw_data_dir <- config_env$RAW_DATA_DIR
        }
        if (is.null(processed_data_dir)) {
            processed_data_dir <- config_env$PROCESSED_DATA_DIR
        }
    }

    child <- build_integration_suite_command(
        raw_data_dir,
        processed_data_dir,
        repo_root
    )
    status <- system2(child$command, child$args, env = child$env)
    if (!identical(status, 0L)) {
        stop("Actual-data integration suite failed.", call. = FALSE)
    }
    invisible(TRUE)
}

if (sys.nframe() == 0L) {
    run_integration_suite()
}

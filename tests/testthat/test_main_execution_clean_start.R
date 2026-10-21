#' Run `main_execution()` with stubbed analyses
#'
#' @param cross_cohort_figure Function that replaces
#'   `create_gep_observed_vs_imported_panel_report()` for the run.
#' @return List with `result` (the `main_execution()` value) and `analyzed_datasets`.
run_main_execution_with_stubs <- function(cross_cohort_figure) {
    temp_processed_dir <- tempfile("processed-clean-start-")
    dir.create(temp_processed_dir, recursive = TRUE, showWarnings = FALSE)

    old_processed_data_dir <- PROCESSED_DATA_DIR
    old_run_objective_0 <- run_objective_0
    old_run_my_analysis <- run_my_analysis
    old_merge_baseline_tables_with_data <- merge_baseline_tables_with_data
    old_cross_cohort_figure <- create_gep_observed_vs_imported_panel_report

    analyzed_datasets <- character()

    on.exit({
        assign("PROCESSED_DATA_DIR", old_processed_data_dir, envir = .GlobalEnv)
        assign("run_objective_0", old_run_objective_0, envir = .GlobalEnv)
        assign("run_my_analysis", old_run_my_analysis, envir = .GlobalEnv)
        assign("merge_baseline_tables_with_data", old_merge_baseline_tables_with_data, envir = .GlobalEnv)
        assign("create_gep_observed_vs_imported_panel_report", old_cross_cohort_figure, envir = .GlobalEnv)
        unlink(temp_processed_dir, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("PROCESSED_DATA_DIR", temp_processed_dir, envir = .GlobalEnv)
    assign("run_objective_0", function() {
        saveRDS(tibble::tibble(id = 1), file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))
        saveRDS(tibble::tibble(id = 2), file.path(PROCESSED_DATA_DIR, "uveal_melanoma_restricted_cohort.rds"))
        list(
            success = TRUE,
            validated_cohorts = c("uveal_melanoma_full_cohort", "uveal_melanoma_restricted_cohort"),
            validation_errors = character(),
            created_datasets = c("uveal_melanoma_full_cohort", "uveal_melanoma_restricted_cohort")
        )
    }, envir = .GlobalEnv)
    assign("run_my_analysis", function(dataset_name, objectives_to_run = c(0, 1, 2, 3, 4)) {
        analyzed_datasets <<- c(analyzed_datasets, dataset_name)
        list(
            fatal_issues = character(),
            warning_issues = character(),
            run_state = "success",
            had_errors = FALSE,
            had_warnings = FALSE
        )
    }, envir = .GlobalEnv)
    assign("merge_baseline_tables_with_data", function(full_data, restricted_data, gksrs_only_data = NULL) {
        invisible(NULL)
    }, envir = .GlobalEnv)
    assign("create_gep_observed_vs_imported_panel_report", cross_cohort_figure, envir = .GlobalEnv)

    result <- main_execution()
    list(result = result, analyzed_datasets = analyzed_datasets)
}

test_that("main_execution discovers datasets created during objective 0 preflight", {
    run <- NULL
    expect_no_error(run <- run_main_execution_with_stubs(function(...) invisible(NULL)))
    expect_setequal(
        run$analyzed_datasets,
        c("uveal_melanoma_full_cohort", "uveal_melanoma_restricted_cohort")
    )
})

test_that("a failing cross-cohort figure fails the run state instead of aborting main_execution", {
    run <- NULL
    expect_no_error(run <- run_main_execution_with_stubs(function(...) stop("workbook missing", call. = FALSE)))
    expect_identical(run$result$run_state, "failed")
    expect_true(any(grepl("^objective4_cross_cohort_figure:workbook missing$", run$result$fatal_issues)))
})

test_that("project paths use distinct analysis and repository slugs", {
    expect_equal(PROJECT_SLUG, "uveal_melanoma")
    expect_equal(REPOSITORY_SLUG, "uveal-melanoma")
})

test_that("project paths follow the canonical workspace and Project Vault contract", {
    expect_equal(
        normalizePath(PROJECT_ROOT, winslash = "/", mustWork = FALSE),
        normalizePath(here::here(), winslash = "/", mustWork = FALSE)
    )
    expect_equal(
        normalizePath(DEFAULT_RUNTIME_ROOT, winslash = "/", mustWork = FALSE),
        file.path(path.expand("~/Workspaces"), "uveal-melanoma", "runtime")
    )
    expect_equal(
        normalizePath(DEFAULT_RAW_DATA_DIR, winslash = "/", mustWork = FALSE),
        file.path(
            path.expand("~/Library/CloudStorage/OneDrive-Personal"),
            "Project Vault",
            "Research",
            "uveal-melanoma",
            "Original Files"
        )
    )
    expect_equal(
        normalizePath(DEFAULT_PUBLISH_ROOT, winslash = "/", mustWork = FALSE),
        file.path(
            path.expand("~/Library/CloudStorage/OneDrive-Personal"),
            "Project Vault",
            "Research",
            "uveal-melanoma",
            "outputs"
        )
    )
})

test_that("standalone tools use the canonical workspace runtime", {
    configured_runtime_root <- RUNTIME_ROOT
    withr::local_envvar(c(
        OUTPUT_DIR = NA,
        OCULAR_RUNTIME_ROOT = NA,
        OCULAR_RUNTIME_PARENT_DIR = NA
    ))

    browse_env <- new.env(parent = baseenv())
    sys.source(
        here::here("scripts", "tools", "browse_diagnostics.R"),
        envir = browse_env
    )
    expect_equal(
        browse_env$get_default_output_dir(),
        file.path(path.expand("~/Workspaces"), "uveal-melanoma", "runtime", "Analysis")
    )
    expect_identical(RUNTIME_ROOT, configured_runtime_root)

    export_env <- new.env(parent = baseenv())
    sys.source(
        here::here("scripts", "tools", "export_gep_objective4_to_downloads.R"),
        envir = export_env
    )
    expect_equal(
        export_env$resolve_runtime_analysis_root("uveal_melanoma"),
        file.path(path.expand("~/Workspaces"), "uveal-melanoma", "runtime", "Analysis")
    )
})

test_that("Objective 4 exporter ships the same files as the publish registry and fails on an empty folder", {
    export_env <- new.env(parent = baseenv())
    sys.source(
        here::here("scripts", "tools", "export_gep_objective4_to_downloads.R"),
        envir = export_env
    )

    root <- withr::local_tempdir()
    files <- c(
        "full_cohort_gep_validation_5yr.xlsx",
        "full_cohort_gep_validation_mfs_curves.png",
        "full_cohort_gep_validation_diagnostics.xlsx",
        "d_no_gep/full_cohort_no_gep_5yr.xlsx",
        "d_no_gep/full_cohort_no_gep_group_risk.png",
        "d_no_gep/full_cohort_no_gep_diagnostics.xlsx"
    )
    for (file in files) {
        dir.create(dirname(file.path(root, file)), recursive = TRUE, showWarnings = FALSE)
        writeLines("x", file.path(root, file))
    }

    shared <- export_env$list_objective4_shareable_files(root)
    shared_relative <- sub(paste0(root, "/"), "", shared, fixed = TRUE)
    registry_relative <- files[vapply(
        files,
        function(file) is_publishable_relative_artifact(file.path("04_GEP_Validation", file), "cohort"),
        logical(1)
    )]
    expect_setequal(shared_relative, registry_relative)
    # Diagnostics workbooks (Record_Issues, by patient id) travel with the results.
    expect_equal(sum(grepl("_diagnostics\\.xlsx$", shared)), 2L)

    destination <- file.path(withr::local_tempdir(), "pack")
    manifest <- export_env$copy_objective4_artifacts(shared, destination, "full")
    expect_true(all(file.exists(manifest$destination_path)))
    expect_setequal(manifest$category, c("excel", "diagnostics", "png"))
    dry <- export_env$copy_objective4_artifacts(shared, file.path(destination, "dry"), "full", dry_run = TRUE)
    expect_true(all(dry$status == "would_copy"))
    expect_false(dir.exists(file.path(destination, "dry")))

    expect_error(export_env$get_objective4_artifact_category("notes.md"), "Unsupported Objective 4 artifact type")
    duplicate_dir <- withr::local_tempdir()
    dir.create(file.path(duplicate_dir, "a"))
    dir.create(file.path(duplicate_dir, "b"))
    writeLines("x", file.path(duplicate_dir, "a", "same.png"))
    writeLines("x", file.path(duplicate_dir, "b", "same.png"))
    expect_error(
        export_env$copy_objective4_artifacts(
            export_env$list_objective4_shareable_files(duplicate_dir), withr::local_tempdir(), "full"
        ),
        "duplicate file names"
    )
})

test_that("portable path tools do not embed a maintainer home directory", {
    tool_paths <- c(
        here::here("scripts", "tools", "browse_diagnostics.R"),
        here::here("scripts", "tools", "export_gep_objective4_to_downloads.R")
    )

    expect_true(all(file.exists(tool_paths)))
    tool_text <- unlist(lapply(tool_paths, readLines, warn = FALSE), use.names = FALSE)
    expect_false(any(grepl("/Users/ncamarda", tool_text, fixed = TRUE)))
})

test_that("portable test runner propagates testthat failures", {
    runner <- here::here("scripts", "ci", "run_testthat.R")
    expect_true(file.exists(runner))

    passing_dir <- tempfile("intentional-testthat-pass-")
    failing_dir <- tempfile("intentional-testthat-failure-")
    dir.create(passing_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(failing_dir, recursive = TRUE, showWarnings = FALSE)
    withr::defer(unlink(passing_dir, recursive = TRUE, force = TRUE), envir = parent.frame())
    withr::defer(unlink(failing_dir, recursive = TRUE, force = TRUE), envir = parent.frame())
    writeLines(
        "testthat::test_that('intentional pass', testthat::expect_true(TRUE))",
        file.path(passing_dir, "test-intentional-pass.R")
    )
    writeLines(
        "testthat::test_that('intentional failure', testthat::expect_true(FALSE))",
        file.path(failing_dir, "test-intentional-failure.R")
    )

    passing_status <- system2(
        file.path(R.home("bin"), "Rscript"),
        c(shQuote(runner), shQuote(passing_dir)),
        stdout = FALSE,
        stderr = FALSE
    )
    failing_status <- suppressWarnings(
        system2(
            file.path(R.home("bin"), "Rscript"),
            c(shQuote(runner), shQuote(failing_dir)),
            stdout = FALSE,
            stderr = FALSE
        )
    )
    expect_identical(passing_status, 0L)
    expect_true(failing_status != 0L)
})

test_that("clean-checkout production and CI files are not excluded", {
    ignore_text <- readLines(here::here(".gitignore"), warn = FALSE)

    expect_true(file.exists(here::here("scripts", "tools", "export_gep_objective4_to_downloads.R")))
    expect_true(file.exists(here::here(".github", "workflows", "portable-tests.yml")))
    expect_true(file.exists(here::here(".lintr")))
    expect_false(any(trimws(ignore_text) == "scripts/tools/export_gep_objective4_to_downloads.R"))
    expect_false(any(trimws(ignore_text) == ".github/"))
    expect_false(any(trimws(ignore_text) == ".lintr"))
})

test_that("configured raw, runtime, and publish paths retain their storage roles", {
    expect_equal(RAW_DATA_DIR, file.path(EXPORT_ROOT, "Original Files"))
    expect_equal(DATA_DIR, EXPORT_ROOT)
    expect_equal(EXPORT_ANALYSIS_DIR, file.path(EXPORT_ROOT, "outputs"))
    expect_equal(PROCESSED_DATA_DIR, file.path(RUNTIME_ROOT, "Analytic Dataset"))
    expect_equal(OUTPUT_DIR, file.path(RUNTIME_ROOT, "Analysis"))
    expect_equal(TOOLS_OUTPUT_DIR, file.path(RUNTIME_ROOT, "tools_output"))
})

test_that("propensity sensitivity output is registered only for the restricted cohort", {
    output_root <- tempfile("cohort-output-routing-")
    old_output <- OUTPUT_DIR
    on.exit({
        assign("OUTPUT_DIR", old_output, envir = .GlobalEnv)
        unlink(output_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    assign("OUTPUT_DIR", output_root, envir = .GlobalEnv)

    full <- setup_cohort_outputs("uveal_melanoma_full_cohort", "uveal_full")
    restricted <- setup_cohort_outputs(
        OBJECTIVE1_PROPENSITY_DATASET,
        "uveal_restricted"
    )
    gksrs <- setup_cohort_outputs("uveal_melanoma_gksrs_only_cohort", "gksrs")

    expect_false("obj1_propensity_sensitivity" %in% names(full$output_dirs))
    expect_true("obj1_propensity_sensitivity" %in% names(restricted$output_dirs))
    expect_false("obj1_propensity_sensitivity" %in% names(gksrs$output_dirs))
    expect_false(dir.exists(file.path(
        full$cohort_base_dir,
        "01_Efficacy",
        "h_propensity_score_sensitivity"
    )))
    expect_true(dir.exists(restricted$output_dirs$obj1_propensity_sensitivity))
    expect_false(dir.exists(file.path(
        gksrs$cohort_base_dir,
        "01_Efficacy",
        "h_propensity_score_sensitivity"
    )))
})

test_that("initialize_runtime_dirs creates configured runtime directories", {
    tmp_root <- tempfile("runtime-init-")
    runtime_root <- file.path(tmp_root, "runtime")

    old_runtime_root <- RUNTIME_ROOT
    old_processed <- PROCESSED_DATA_DIR
    old_output <- OUTPUT_DIR
    old_logs <- LOGS_DIR
    old_tools <- TOOLS_OUTPUT_DIR
    old_test <- TEST_OUTPUT_DIR
    old_merged <- MERGED_TABLES_DIR

    on.exit({
        assign("RUNTIME_ROOT", old_runtime_root, envir = .GlobalEnv)
        assign("PROCESSED_DATA_DIR", old_processed, envir = .GlobalEnv)
        assign("OUTPUT_DIR", old_output, envir = .GlobalEnv)
        assign("LOGS_DIR", old_logs, envir = .GlobalEnv)
        assign("TOOLS_OUTPUT_DIR", old_tools, envir = .GlobalEnv)
        assign("TEST_OUTPUT_DIR", old_test, envir = .GlobalEnv)
        assign("MERGED_TABLES_DIR", old_merged, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("RUNTIME_ROOT", runtime_root, envir = .GlobalEnv)
    assign("PROCESSED_DATA_DIR", file.path(runtime_root, "Analytic Dataset"), envir = .GlobalEnv)
    assign("OUTPUT_DIR", file.path(runtime_root, "Analysis"), envir = .GlobalEnv)
    assign("LOGS_DIR", file.path(runtime_root, "logs"), envir = .GlobalEnv)
    assign("TOOLS_OUTPUT_DIR", file.path(runtime_root, "tools_output"), envir = .GlobalEnv)
    assign("TEST_OUTPUT_DIR", file.path(runtime_root, "test_output"), envir = .GlobalEnv)
    assign("MERGED_TABLES_DIR", file.path(runtime_root, "Analysis", "merged_tables"), envir = .GlobalEnv)

    initialize_runtime_dirs()

    expect_true(dir.exists(RUNTIME_ROOT))
    expect_true(dir.exists(PROCESSED_DATA_DIR))
    expect_true(dir.exists(OUTPUT_DIR))
    expect_true(dir.exists(LOGS_DIR))
    expect_true(dir.exists(TOOLS_OUTPUT_DIR))
    expect_true(dir.exists(TEST_OUTPUT_DIR))
    expect_true(dir.exists(MERGED_TABLES_DIR))
})

test_that("assert_required_input_paths fails fast for missing raw inputs", {
    tmp_root <- tempfile("runtime-raw-check-")
    dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

    old_raw <- RAW_DATA_DIR
    old_input <- INPUT_FILENAME

    on.exit({
        assign("RAW_DATA_DIR", old_raw, envir = .GlobalEnv)
        assign("INPUT_FILENAME", old_input, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("RAW_DATA_DIR", file.path(tmp_root, "missing_raw"), envir = .GlobalEnv)
    assign("INPUT_FILENAME", "missing_input.xlsx", envir = .GlobalEnv)

    expect_error(
        assert_required_input_paths(),
        regexp = "Required raw input path checks failed"
    )
})

test_that("get_export_snapshot_dir constructs dated snapshot path", {
    tmp_root <- tempfile("runtime-export-")
    export_root <- file.path(tmp_root, "export")
    export_analysis_root <- file.path(export_root, "Analysis")

    old_export_root <- EXPORT_ROOT
    old_export_analysis <- EXPORT_ANALYSIS_DIR

    on.exit({
        assign("EXPORT_ROOT", old_export_root, envir = .GlobalEnv)
        assign("EXPORT_ANALYSIS_DIR", old_export_analysis, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("EXPORT_ROOT", export_root, envir = .GlobalEnv)
    assign("EXPORT_ANALYSIS_DIR", export_analysis_root, envir = .GlobalEnv)

    snapshot_dir <- get_export_snapshot_dir("2026-03-25")
    expect_equal(
        snapshot_dir,
        file.path(EXPORT_ANALYSIS_DIR, "2026-03-25")
    )
})

test_that("resolve_config_path rejects relative runtime or export overrides", {
    expect_error(
        resolve_config_path("relative/runtime", "/tmp/runtime-default"),
        regexp = "must be absolute"
    )
})

test_that("artifact registry allowlists only explicit publishable outputs", {
    expect_true(is_publishable_relative_artifact("01_Efficacy/summary.xlsx", "cohort"))
    expect_true(is_publishable_relative_artifact("merged_tables/final_table.csv", "merged_tables"))
    expect_true(is_publishable_relative_artifact("merged_tables/objective4_cross_cohort_figures/three_panel.png", "merged_tables"))
    expect_true(is_publishable_relative_artifact("merged_tables/objective4_cross_cohort_figures/three_panel_report.md", "merged_tables"))
    expect_true(is_publishable_relative_artifact("04_GEP_Validation/full_cohort_gep_validation_5yr.xlsx", "cohort"))
    expect_true(is_publishable_relative_artifact("04_GEP_Validation/d_no_gep/full_cohort_no_gep_group_risk.png", "cohort"))
    # The GEP diagnostics workbooks carry the Record_Issues sheets and are part of the package.
    expect_true(is_publishable_relative_artifact("04_GEP_Validation/full_cohort_gep_validation_diagnostics.xlsx", "cohort"))
    expect_true(is_publishable_relative_artifact("04_GEP_Validation/d_no_gep/full_cohort_no_gep_diagnostics.xlsx", "cohort"))
    expect_false(is_publishable_relative_artifact("01_Efficacy/model.rds", "cohort"))
    expect_false(is_publishable_relative_artifact("01_Efficacy/summary_diagnostics.xlsx", "cohort"))
    expect_false(is_publishable_relative_artifact("cache/intermediate.csv", "cohort"))
    expect_true(is_publishable_relative_artifact(
        "01_Efficacy/h_propensity_score_sensitivity/restricted_cohort_propensity_overlap_summary.md",
        "cohort"
    ))
    expect_false(is_publishable_relative_artifact(
        "01_Efficacy/h_propensity_score_sensitivity/restricted_cohort_propensity_design_audit.rds",
        "cohort"
    ))
    expect_false(is_publishable_relative_artifact(
        "01_Efficacy/h_propensity_score_sensitivity/patient_weights.csv",
        "cohort"
    ))
})

test_that("propensity publication dry run includes six reader artifacts and excludes audit RDS", {
    tmp_root <- tempfile("propensity-publish-dryrun-")
    install_publish_contract_fixture(tmp_root)
    output_root <- file.path(tmp_root, "runtime", "Analysis")
    export_root <- file.path(tmp_root, "export")
    propensity_dir <- file.path(
        output_root,
        "uveal_restricted",
        "01_Efficacy",
        "h_propensity_score_sensitivity"
    )
    dir.create(propensity_dir, recursive = TRUE)
    reader_basenames <- paste0(
        "restricted_cohort_",
        unname(OBJECTIVE1_PROPENSITY_ARTIFACT_BASENAMES[names(
            OBJECTIVE1_PROPENSITY_ARTIFACT_BASENAMES
        ) != "audit"])
    )
    audit_basename <- paste0(
        "restricted_cohort_",
        OBJECTIVE1_PROPENSITY_ARTIFACT_BASENAMES[["audit"]]
    )
    for (basename in c(reader_basenames, audit_basename)) {
        writeLines("unit artifact", file.path(propensity_dir, basename))
    }

    old_output <- OUTPUT_DIR
    old_merged <- MERGED_TABLES_DIR
    old_export_root <- EXPORT_ROOT
    old_export_analysis <- EXPORT_ANALYSIS_DIR
    on.exit({
        assign("OUTPUT_DIR", old_output, envir = .GlobalEnv)
        assign("MERGED_TABLES_DIR", old_merged, envir = .GlobalEnv)
        assign("EXPORT_ROOT", old_export_root, envir = .GlobalEnv)
        assign("EXPORT_ANALYSIS_DIR", old_export_analysis, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    assign("OUTPUT_DIR", output_root, envir = .GlobalEnv)
    assign("MERGED_TABLES_DIR", file.path(output_root, "merged_tables"), envir = .GlobalEnv)
    assign("EXPORT_ROOT", export_root, envir = .GlobalEnv)
    assign("EXPORT_ANALYSIS_DIR", file.path(export_root, "outputs"), envir = .GlobalEnv)

    result <- publish_outputs(
        cohorts = OBJECTIVE1_PROPENSITY_DATASET,
        snapshot_id = "2026-08-04-propensity-unit",
        include_merged_tables = FALSE,
        dry_run = TRUE
    )
    would_copy <- result$manifest$source_path[result$manifest$status == "would_copy"]

    expect_equal(result$summary$would_copy, 9L)
    expect_true(all(reader_basenames %in% basename(would_copy)))
    expect_false(audit_basename %in% basename(would_copy))
})

test_that("publish_outputs dry run reports publishable outputs and excludes runtime artifacts", {
    tmp_root <- tempfile("runtime-publish-dryrun-")
    install_publish_contract_fixture(tmp_root)
    runtime_root <- file.path(tmp_root, "runtime")
    output_root <- file.path(runtime_root, "Analysis")
    export_root <- file.path(tmp_root, "export")
    export_analysis_root <- file.path(export_root, "Analysis")

    dir.create(file.path(output_root, "uveal_full", "01_Efficacy"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_root, "uveal_full", "cache"), recursive = TRUE, showWarnings = FALSE)

    writeLines("report", file.path(output_root, "uveal_full", "01_Efficacy", "summary.xlsx"))
    writeLines("model", file.path(output_root, "uveal_full", "01_Efficacy", "model.rds"))
    writeLines("diagnostics", file.path(output_root, "uveal_full", "01_Efficacy", "summary_diagnostics.xlsx"))
    writeLines("cache", file.path(output_root, "uveal_full", "cache", "cache.csv"))

    old_output <- OUTPUT_DIR
    old_merged <- MERGED_TABLES_DIR
    old_export_root <- EXPORT_ROOT
    old_export_analysis <- EXPORT_ANALYSIS_DIR

    on.exit({
        assign("OUTPUT_DIR", old_output, envir = .GlobalEnv)
        assign("MERGED_TABLES_DIR", old_merged, envir = .GlobalEnv)
        assign("EXPORT_ROOT", old_export_root, envir = .GlobalEnv)
        assign("EXPORT_ANALYSIS_DIR", old_export_analysis, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("OUTPUT_DIR", output_root, envir = .GlobalEnv)
    assign("MERGED_TABLES_DIR", file.path(output_root, "merged_tables"), envir = .GlobalEnv)
    assign("EXPORT_ROOT", export_root, envir = .GlobalEnv)
    assign("EXPORT_ANALYSIS_DIR", export_analysis_root, envir = .GlobalEnv)

    result <- publish_outputs(snapshot_id = "2026-03-25-unit", dry_run = TRUE)

    expect_true(result$dry_run)
    expect_equal(result$summary$would_copy, 8)
    expect_equal(result$summary$copied, 0)
    expect_true(any(result$manifest$status == "skipped_not_publishable"))
    expect_true(any(result$manifest$status == "optional_root_absent"))
    expect_equal(result$summary$missing, 0)
    expect_false(dir.exists(result$snapshot_dir))
})

test_that("publish CLI report gives dry-run scope and exact execute command", {
    report <- format_publish_outputs_cli_report(
        result = list(
            snapshot_dir = "/tmp/export/2026-07-28",
            dry_run = TRUE,
            summary = list(
                publishable_files = 2,
                analytic_data_files = 0,
                analysis_log_basename = "run_log_20260728_120000.txt",
                copied = 0,
                would_copy = 2,
                skipped = 1,
                missing = 0,
                failed = 0,
                snapshot_exists = FALSE
            )
        ),
        opts = list(
            snapshot_id = "2026-07-28",
            snapshot_id_supplied = TRUE,
            cohorts = c("uveal_melanoma_full_cohort", "gksrs"),
            include_merged_tables = FALSE
        )
    )

    expect_match(report, "DRY RUN: no files were copied", fixed = TRUE)
    expect_match(report, "Would copy: 2", fixed = TRUE)
    expect_match(report, "Excluded by registry: 1", fixed = TRUE)
    expect_match(
        report,
        "Rscript scripts/workflow/publish_outputs.R --execute --snapshot-id 2026-07-28 --cohorts uveal_melanoma_full_cohort,gksrs --no-merged-tables",
        fixed = TRUE
    )

    default_snapshot_command <- publish_execute_command(list(
        snapshot_id_supplied = FALSE,
        cohorts = NULL,
        include_merged_tables = TRUE
    ))
    expect_identical(
        default_snapshot_command,
        "Rscript scripts/workflow/publish_outputs.R --execute"
    )
})

test_that("default publish snapshots add letter suffixes on the same day", {
    tmp_root <- tempfile("runtime-publish-suffix-")
    export_root <- file.path(tmp_root, "export")
    export_analysis_root <- file.path(export_root, "Analysis")
    dir.create(file.path(export_analysis_root, "2026-08-04"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(export_analysis_root, "2026-08-04-a"), recursive = TRUE, showWarnings = FALSE)

    old_export_root <- EXPORT_ROOT
    old_export_analysis <- EXPORT_ANALYSIS_DIR
    on.exit({
        assign("EXPORT_ROOT", old_export_root, envir = .GlobalEnv)
        assign("EXPORT_ANALYSIS_DIR", old_export_analysis, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("EXPORT_ROOT", export_root, envir = .GlobalEnv)
    assign("EXPORT_ANALYSIS_DIR", export_analysis_root, envir = .GlobalEnv)

    expect_identical(
        next_available_publish_snapshot_id(default_snapshot_id = "2026-08-04"),
        "2026-08-04-b"
    )
})

test_that("publish CLI main prints the concise report without manifest rows", {
    tmp_root <- tempfile("runtime-publish-cli-")
    install_publish_contract_fixture(tmp_root)
    output_root <- file.path(tmp_root, "runtime", "Analysis")
    export_root <- file.path(tmp_root, "export")
    export_analysis_root <- file.path(export_root, "Analysis")
    dir.create(file.path(output_root, "uveal_full", "01_Efficacy"), recursive = TRUE, showWarnings = FALSE)
    writeLines("report", file.path(output_root, "uveal_full", "01_Efficacy", "summary.xlsx"))

    old_output <- OUTPUT_DIR
    old_merged <- MERGED_TABLES_DIR
    old_export_root <- EXPORT_ROOT
    old_export_analysis <- EXPORT_ANALYSIS_DIR
    on.exit({
        assign("OUTPUT_DIR", old_output, envir = .GlobalEnv)
        assign("MERGED_TABLES_DIR", old_merged, envir = .GlobalEnv)
        assign("EXPORT_ROOT", old_export_root, envir = .GlobalEnv)
        assign("EXPORT_ANALYSIS_DIR", old_export_analysis, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("OUTPUT_DIR", output_root, envir = .GlobalEnv)
    assign("MERGED_TABLES_DIR", file.path(output_root, "merged_tables"), envir = .GlobalEnv)
    assign("EXPORT_ROOT", export_root, envir = .GlobalEnv)
    assign("EXPORT_ANALYSIS_DIR", export_analysis_root, envir = .GlobalEnv)

    output <- capture_output(main(list(
        snapshot_id = "2026-07-28-cli",
        dry_run = TRUE,
        include_merged_tables = FALSE,
        cohorts = "uveal_full",
        help = FALSE
    )))

    expect_match(output, "DRY RUN: no files were copied", fixed = TRUE)
    expect_match(output, "Next step (performs the copy):", fixed = TRUE)
    expect_false(grepl("source_path", output, fixed = TRUE))
})

test_that("publish_outputs copies the Objective 4 workbooks, diagnostics workbooks and figures", {
    tmp_root <- tempfile("runtime-publish-obj4-")
    install_publish_contract_fixture(tmp_root)
    runtime_root <- file.path(tmp_root, "runtime")
    output_root <- file.path(runtime_root, "Analysis")
    export_root <- file.path(tmp_root, "export")
    export_analysis_root <- file.path(export_root, "Analysis")
    obj4_dir <- file.path(output_root, "uveal_full", "04_GEP_Validation")
    obj4_no_gep_dir <- file.path(obj4_dir, "d_no_gep")
    obj1_dir <- file.path(output_root, "uveal_full", "01_Efficacy")

    dir.create(obj4_no_gep_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(obj1_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("model diagnostics", file.path(obj1_dir, "full_cohort_overall_survival_diagnostics.xlsx"))
    writeLines("primary", file.path(obj4_dir, "full_cohort_gep_validation_5yr.xlsx"))
    writeLines("figure", file.path(obj4_dir, "full_cohort_gep_validation_mfs_curves.png"))
    writeLines("record identifiers", file.path(obj4_dir, "full_cohort_gep_validation_diagnostics.xlsx"))
    writeLines("primary", file.path(obj4_no_gep_dir, "full_cohort_no_gep_5yr.xlsx"))
    writeLines("record identifiers", file.path(obj4_no_gep_dir, "full_cohort_no_gep_diagnostics.xlsx"))

    old_output <- OUTPUT_DIR
    old_merged <- MERGED_TABLES_DIR
    old_export_root <- EXPORT_ROOT
    old_export_analysis <- EXPORT_ANALYSIS_DIR

    on.exit({
        assign("OUTPUT_DIR", old_output, envir = .GlobalEnv)
        assign("MERGED_TABLES_DIR", old_merged, envir = .GlobalEnv)
        assign("EXPORT_ROOT", old_export_root, envir = .GlobalEnv)
        assign("EXPORT_ANALYSIS_DIR", old_export_analysis, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("OUTPUT_DIR", output_root, envir = .GlobalEnv)
    assign("MERGED_TABLES_DIR", file.path(output_root, "merged_tables"), envir = .GlobalEnv)
    assign("EXPORT_ROOT", export_root, envir = .GlobalEnv)
    assign("EXPORT_ANALYSIS_DIR", export_analysis_root, envir = .GlobalEnv)

    result <- publish_outputs(snapshot_id = "2026-04-02-obj4", dry_run = TRUE)

    would_copy_sources <- result$manifest$source_path[result$manifest$status == "would_copy"]
    skipped_sources <- result$manifest$source_path[result$manifest$status == "skipped_not_publishable"]

    expect_true(any(grepl("full_cohort_gep_validation_5yr\\.xlsx$", would_copy_sources)))
    expect_true(any(grepl("full_cohort_gep_validation_mfs_curves\\.png$", would_copy_sources)))
    expect_true(any(grepl("full_cohort_no_gep_5yr\\.xlsx$", would_copy_sources)))
    expect_true(any(grepl("full_cohort_gep_validation_diagnostics\\.xlsx$", would_copy_sources)))
    expect_true(any(grepl("full_cohort_no_gep_diagnostics\\.xlsx$", would_copy_sources)))
    # Model diagnostics workbooks of the other objectives stay out of snapshots.
    expect_true(any(grepl("full_cohort_overall_survival_diagnostics\\.xlsx$", skipped_sources)))
})

test_that("publish_outputs creates a new snapshot and rejects existing snapshot overwrite", {
    tmp_root <- tempfile("runtime-publish-copy-")
    install_publish_contract_fixture(tmp_root)
    runtime_root <- file.path(tmp_root, "runtime")
    output_root <- file.path(runtime_root, "Analysis")
    export_root <- file.path(tmp_root, "export")
    export_analysis_root <- file.path(export_root, "Analysis")

    dir.create(file.path(output_root, "uveal_full", "01_Efficacy"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_root, "merged_tables"), recursive = TRUE, showWarnings = FALSE)
    writeLines("report", file.path(output_root, "uveal_full", "01_Efficacy", "summary.xlsx"))
    writeLines("merged", file.path(output_root, "merged_tables", "merged_baseline.csv"))

    old_output <- OUTPUT_DIR
    old_logs <- LOGS_DIR
    old_merged <- MERGED_TABLES_DIR
    old_export_root <- EXPORT_ROOT
    old_export_analysis <- EXPORT_ANALYSIS_DIR

    on.exit({
        assign("OUTPUT_DIR", old_output, envir = .GlobalEnv)
        assign("LOGS_DIR", old_logs, envir = .GlobalEnv)
        assign("MERGED_TABLES_DIR", old_merged, envir = .GlobalEnv)
        assign("EXPORT_ROOT", old_export_root, envir = .GlobalEnv)
        assign("EXPORT_ANALYSIS_DIR", old_export_analysis, envir = .GlobalEnv)
        unlink(tmp_root, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    assign("OUTPUT_DIR", output_root, envir = .GlobalEnv)
    assign("LOGS_DIR", file.path(runtime_root, "logs"), envir = .GlobalEnv)
    assign("MERGED_TABLES_DIR", file.path(output_root, "merged_tables"), envir = .GlobalEnv)
    assign("EXPORT_ROOT", export_root, envir = .GlobalEnv)
    assign("EXPORT_ANALYSIS_DIR", export_analysis_root, envir = .GlobalEnv)

    first_publish <- publish_outputs(snapshot_id = "2026-03-25-publish", dry_run = FALSE)
    expect_true(dir.exists(first_publish$snapshot_dir))
    expect_true(file.exists(file.path(first_publish$snapshot_dir, "uveal_full", "01_Efficacy", "summary.xlsx")))
    expect_true(file.exists(file.path(first_publish$snapshot_dir, "merged_tables", "merged_baseline.csv")))
    expect_true(file.exists(file.path(first_publish$snapshot_dir, basename(first_publish$latest_analysis_log))))
    expect_equal(length(list.files(file.path(first_publish$snapshot_dir, "analytic_data"))), 6L)
    expect_true(file.exists(file.path(first_publish$snapshot_dir, "publish_manifest.csv")))

    expect_error(
        publish_outputs(snapshot_id = "2026-03-25-publish", dry_run = FALSE),
        regexp = "Snapshot target already exists"
    )
})
test_that("configured raw input uses the versioned NDC workbook", {
    expect_identical(
        INPUT_FILENAME,
        "Ocular Melanoma Master Spreadsheet REVISED FOR STATS (08-11-26 NDC).xlsx"
    )
})

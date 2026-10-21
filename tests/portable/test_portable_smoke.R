test_that("portable smoke path writes GEP validation outputs and dry-run publish on synthetic data", {
    cohort_root <- file.path(OUTPUT_DIR, "uveal_restricted")
    output_dirs <- list(obj4_gep = file.path(cohort_root, "04_GEP_Validation"))

    for (dir_path in output_dirs) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
    withr::defer(unlink(cohort_root, recursive = TRUE), envir = parent.frame())

    log_basename <- paste0(
        "run_log_",
        format(Sys.time() + 2, "%Y%m%d_%H%M%S")
    )
    configured_log_path <- file.path(LOGS_DIR, paste0(log_basename, ".txt"))
    setup_logging(
        log_path = configured_log_path,
        level = "INFO",
        progress = FALSE,
        context_in_file = TRUE
    )
    withr::defer(
        setup_logging(log_path = NULL, level = "INFO", progress = FALSE),
        envir = parent.frame()
    )
    set_log_context(replace = TRUE)
    log_phase("MAIN EXECUTION PHASE")

    analytic_data <- create_synthetic_ci_dataset()
    dir.create(PROCESSED_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
    save_cohorts(list(uveal_melanoma_restricted_cohort = analytic_data))

    validation <- run_gep_definitive_validation(analytic_data, n_boot = 20L)
    outputs <- write_gep_validation_outputs(
        validation, analytic_data, "uveal_melanoma_restricted_cohort", output_dirs, "smoke_"
    )

    expect_true(all(file.exists(c(outputs$primary, outputs$diagnostics, outputs$figures))))

    logger::log_info(">>> ALL ANALYSES COMPLETED SUCCESSFULLY!")
    logger::log_info(">>> Datasets analyzed: 3")
    logger::log_info(">>> COMPLETED MAIN EXECUTION PHASE (Duration: 0.0 seconds)")

    text_log_path <- file.path(LOGS_DIR, "txt", paste0(log_basename, ".txt"))
    expect_true(file.exists(text_log_path))

    publish_result <- publish_outputs(
        cohorts = "uveal_melanoma_restricted_cohort",
        snapshot_id = "portable-smoke",
        include_merged_tables = FALSE,
        dry_run = TRUE
    )

    expect_true(publish_result$dry_run)
    expect_gt(publish_result$summary$would_copy, 0)
    expect_equal(
        publish_result$snapshot_dir,
        file.path(EXPORT_ANALYSIS_DIR, "portable-smoke")
    )
    expect_true(any(grepl("04_GEP_Validation", publish_result$manifest$source_path)))
    # The primary workbook and the diagnostics workbook (with the Record_Issues sheet,
    # by patient id) are both part of the published package.
    published_sources <- normalizePath(publish_result$manifest$source_path[publish_result$manifest$status == "would_copy"], winslash = "/", mustWork = FALSE)
    expect_true(normalizePath(outputs$primary, winslash = "/", mustWork = FALSE) %in% published_sources)
    expect_true(normalizePath(outputs$diagnostics, winslash = "/", mustWork = FALSE) %in% published_sources)
    normalized_manifest_sources <- normalizePath(
        publish_result$manifest$source_path,
        winslash = "/",
        mustWork = FALSE
    )
    normalized_text_log_path <- normalizePath(
        text_log_path,
        winslash = "/",
        mustWork = TRUE
    )
    log_manifest_index <- normalized_manifest_sources == normalized_text_log_path
    expect_true(any(log_manifest_index))
    expect_false(any(grepl("\\.jsonl$", publish_result$manifest$source_path)))
    expect_equal(
        normalizePath(
            dirname(publish_result$manifest$destination_path[log_manifest_index]),
            winslash = "/",
            mustWork = FALSE
        ),
        normalizePath(publish_result$snapshot_dir, winslash = "/", mustWork = FALSE)
    )
})

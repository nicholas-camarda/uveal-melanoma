build_objective4_mfs_sensitivity_output_dirs <- function(test_output_dir) {
    list(
        obj4_mfs = file.path(test_output_dir, "04_GEP_Validation", "a_mfs_validation"),
        obj4_mss = file.path(test_output_dir, "04_GEP_Validation", "b_mss_validation"),
        obj4_ph_diagnostics = file.path(test_output_dir, "04_GEP_Validation", "c_proportional_hazards_diagnostics")
    )
}

test_that("Objective 4 MFS sensitivity collector builds follow-up and repeat-radiation summaries", {
    test_data <- create_test_dataset() %>%
        dplyr::mutate(
            recurrence1_treatment_clean = factor(
                c("GKSRS", "Plaque", rep(NA_character_, 18)),
                levels = c("GKSRS", "Plaque")
            ),
            initial_gk = replace(initial_gk, 1, "Y"),
            initial_plaque = replace(initial_plaque, 1, "Y")
        )

    sensitivity <- collect_objective4_mfs_sensitivity_results(
        data = test_data,
        dataset_name = "uveal_melanoma_full_cohort"
    )

    expect_true(all(c(
        "mfs_followup_sensitivity",
        "mfs_treatment_mix_sensitivity",
        "mfs_repeat_radiation_sensitivity",
        "guardrail_notes"
    ) %in% names(sensitivity)))

    horizon_overall <- sensitivity$mfs_followup_sensitivity$horizon_overall
    expect_equal(sum(horizon_overall$n), 20)
    expect_equal(
        horizon_overall$n[horizon_overall$five_year_followup_view == "event_by_5yr"],
        8
    )
    expect_equal(
        horizon_overall$n[horizon_overall$five_year_followup_view == "followup_ge_5yr"],
        4
    )
    expect_equal(
        horizon_overall$n[horizon_overall$five_year_followup_view == "censored_pre_5yr"],
        8
    )

    repeat_quality <- sensitivity$mfs_repeat_radiation_sensitivity$data_quality_checks
    expect_equal(repeat_quality$repeat_radiation_exposure_n, 2)
    expect_equal(repeat_quality$both_initial_modalities_n, 1)

    pbt_only <- sensitivity$mfs_treatment_mix_sensitivity$pbt_only
    expect_true(nrow(pbt_only) > 0)
    expect_equal(sum(pbt_only$n), 10)
})

test_that("Objective 4 MFS sensitivity guardrail reports PBT rows inside gksrs-only cohorts", {
    test_data <- create_test_dataset() %>%
        dplyr::mutate(
            recurrence1_treatment_clean = factor(
                rep(NA_character_, 20),
                levels = c("GKSRS", "Plaque")
            )
        )

    sensitivity <- collect_objective4_mfs_sensitivity_results(
        data = test_data,
        dataset_name = "uveal_melanoma_gksrs_only_cohort"
    )

    expect_true(any(grepl(
        "PBT-labeled MFS-eligible row",
        sensitivity$guardrail_notes$note,
        fixed = TRUE
    )))
})

test_that("Objective 4 MFS sensitivity writer saves workbook and summary text", {
    test_output_dir <- file.path(TEST_OUTPUT_DIR, "objective4_mfs_sensitivity_writer")
    dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)
    withr::defer(unlink(test_output_dir, recursive = TRUE), envir = parent.frame())

    sensitivity <- collect_objective4_mfs_sensitivity_results(
        data = create_test_dataset(),
        dataset_name = "uveal_melanoma_full_cohort"
    )

    paths <- write_objective4_mfs_sensitivity_outputs(
        sensitivity_results = sensitivity,
        output_dir = test_output_dir,
        prefix = "writer_"
    )

    expect_true(file.exists(paths$workbook))
    expect_true(file.exists(paths$summary))

    workbook_sheets <- readxl::excel_sheets(paths$workbook)
    expect_true(all(c(
        "Followup_Operational",
        "Followup_5yr",
        "TxMix_ByClass",
        "Repeat_Comparison",
        "Guardrail_Notes"
    ) %in% workbook_sheets))

    summary_text <- paste(readLines(paths$summary, warn = FALSE), collapse = "\n")
    expect_match(summary_text, "^# Objective 4 MFS Sensitivity Summary", perl = TRUE)
    expect_match(summary_text, "`followup_ge_5yr` means")
    expect_match(summary_text, "## Event-Row Diagnostics", fixed = TRUE)
    expect_match(summary_text, "\\| View \\| Count \\|", perl = TRUE)
    expect_match(summary_text, "event row IDs:")
})

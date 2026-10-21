# Test file for Objective 4: GEP Analysis
# Tests the actual content and statistical results of the GEP analysis pipeline
# Run with: testthat::test_dir('tests/testthat')

# CRITICAL: Set test environment variables BEFORE sourcing load_all.R
# This prevents directory creation in the wrong location

# Load required packages first
library(dplyr)

# Override project constants to prevent test interference

# Load the project environment with ALL of the variables and functions
# You do not need to load libraries separately

# Source the helper file for test data creation
extract_risk_table_y_limits <- function(risk_table_plot) {
    y_scales <- Filter(function(scale) "y" %in% scale$aesthetics, risk_table_plot$scales$scales)
    y_scales[[length(y_scales)]]$limits
}

extract_displayed_risk_counts <- function(risk_table_plot, time_point = 0) {
    built_layer <- suppressWarnings(ggplot2::ggplot_build(risk_table_plot)$data[[1]])
    built_layer %>%
        dplyr::filter(.data$x == time_point) %>%
        dplyr::arrange(dplyr::desc(.data$y)) %>%
        dplyr::pull(.data$label) %>%
        as.integer()
}

# Test 4a: GEP MFS validation produces valid statistical results
# Test 4b: GEP MSS validation produces valid statistical results
test_that("survival helper can separate KM display groups from Cox model groups", {
    test_data <- tibble::tibble(
        tt_mets_months = c(12, 24, 18, 30),
        id = sprintf("T%03d", seq_along(tt_mets_months)),
        mets_event = c(0, 1, 0, 1),
        biopsy1_gep = factor(c("Class 1", "Class 1", GEP_GROUP_TEST_FAILED, GEP_GROUP_TEST_FAILED)),
        biopsy1_gep_model = factor(c("Class 1", "Class 1", GEP_GROUP_TEST_FAILED, GEP_GROUP_TEST_FAILED))
    )

    convergence_warnings <- character()
    result <- withCallingHandlers(
        analyze_time_to_event_outcomes(
            data = test_data,
            time_var = "tt_mets_months",
            event_var = "mets_event",
            group_var = "biopsy1_gep",
            model_group_var = "biopsy1_gep_model",
            confounders = NULL,
            ylab = "Metastasis-Free Survival Probability",
            analysis_type = "post_treatment_only",
            dataset_name = "unit_test_display_model_split",
            output_dirs = create_output_structure(withr::local_tempdir()),
            prefix = "unit_",
            route_key = "obj1_mets"
        ),
        warning = function(warning_condition) {
            convergence_warnings <<- c(
                convergence_warnings,
                conditionMessage(warning_condition)
            )
            invokeRestart("muffleWarning")
        }
    )

    expect_true(length(convergence_warnings) > 0)
    expect_true(all(grepl(
        "Ran out of iterations and did not converge",
        convergence_warnings,
        fixed = TRUE
    )))
    expect_true(any(grepl(GEP_GROUP_TEST_FAILED, names(result$fit$strata), fixed = TRUE)))
    expect_null(result$cox_model)
    expect_match(result$diagnostics$raw_model_output, "After sparse-level exclusions")
    expect_match(result$diagnostics$raw_model_output, "requires at least two non-missing")
})

test_that("risk table keeps displayed row labels aligned with their counts when legend order differs from fit order", {
    desired_order <- GEP_DETAIL_LEVELS
    baseline_counts <- c(5, 4, 3, 2, 1, 1, 6, 7)
    names(baseline_counts) <- desired_order

    test_data <- tibble::tibble(
        tt_mets_months = seq_len(sum(baseline_counts)) * 6,
        id = sprintf("T%03d", seq_along(tt_mets_months)),
        mets_event = rep(c(0, 1), length.out = sum(baseline_counts)),
        biopsy1_gep = factor(
            unlist(Map(rep, names(baseline_counts), baseline_counts)),
            levels = rev(desired_order)
        )
    )

    result <- suppressWarnings(analyze_time_to_event_outcomes(
        data = test_data,
        time_var = "tt_mets_months",
        event_var = "mets_event",
        group_var = "biopsy1_gep",
        model_group_var = "biopsy1_gep",
        confounders = NULL,
        ylab = "Metastasis-Free Survival Probability",
        analysis_type = "post_treatment_only",
        dataset_name = "unit_test_risk_table_alignment",
        legend_labels = desired_order,
        output_dirs = create_output_structure(withr::local_tempdir()),
        prefix = "unit_",
        route_key = "obj1_mets"
    ))

    expect_equal(rev(extract_risk_table_y_limits(result$plot$table)), desired_order)
    expect_equal(extract_displayed_risk_counts(result$plot$table, time_point = 0), unname(baseline_counts[desired_order]))
})

test_that("risk table preserves requested row order for a subset of present strata", {
    desired_order <- GEP_DETAIL_LEVELS
    baseline_counts <- c(5, 2, 7)
    names(baseline_counts) <- desired_order[c(1, 4, 8)]
    present_order <- desired_order[desired_order %in% names(baseline_counts)]

    test_data <- tibble::tibble(
        tt_mets_months = seq_len(sum(baseline_counts)) * 6,
        id = sprintf("T%03d", seq_along(tt_mets_months)),
        mets_event = rep(c(0, 1), length.out = sum(baseline_counts)),
        biopsy1_gep = factor(
            unlist(Map(rep, names(baseline_counts), baseline_counts)),
            levels = rev(desired_order)
        )
    )

    result <- suppressWarnings(analyze_time_to_event_outcomes(
        data = test_data,
        time_var = "tt_mets_months",
        event_var = "mets_event",
        group_var = "biopsy1_gep",
        model_group_var = "biopsy1_gep",
        confounders = NULL,
        ylab = "Metastasis-Free Survival Probability",
        analysis_type = "post_treatment_only",
        dataset_name = "unit_test_risk_table_alignment_subset",
        legend_labels = desired_order,
        output_dirs = create_output_structure(withr::local_tempdir()),
        prefix = "unit_",
        route_key = "obj1_mets"
    ))

    expect_equal(rev(extract_risk_table_y_limits(result$plot$table)), present_order)
    expect_equal(extract_displayed_risk_counts(result$plot$table, time_point = 0), unname(baseline_counts[present_order]))
})

test_that("Cox diagnostics keep group counts and event rates in Term_support", {
    test_output_dir <- file.path(TEST_OUTPUT_DIR, "objective4_group_counts_in_raw_output")
    dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

    test_data <- tibble::tibble(
        tt_mets_months = c(8, 16, 24, 32, 10, 20, 30, 40),
        id = sprintf("T%03d", seq_along(tt_mets_months)),
        mets_event = c(1, 0, 0, 1, 0, 1, 0, 1),
        biopsy1_gep_model = factor(
            c("Class 1", "Class 1", "Class 1", "Class 1", "Class 2", "Class 2", "Class 2", "Class 2"),
            levels = c("Class 1", "Class 2")
        )
    )

    result <- generate_regression_table(
        data = test_data,
        outcome_var = "mets_event",
        predictor_vars = "biopsy1_gep_model",
        confounders = NULL,
        model_type = "cox",
        effect_measure = "HR",
        analysis_name = "unit_test_group_counts",
        dataset_name = "unit_test_group_counts",
        output_dir = test_output_dir,
        prefix = "unit_",
        time_var = "tt_mets_months",
        event_var = "mets_event",
        treatment_var = "biopsy1_gep_model"
    )

    expect_true(file.exists(result$output_files$diagnostics_path))

    term_support <- readxl::read_excel(result$output_files$diagnostics_path, sheet = "Term_support")
    expect_true(all(c("n", "events", "non_events", "event_rate_percent", "is_reference") %in% names(term_support)))

    reference_row <- term_support[term_support$variable == "biopsy1_gep_model" & term_support$level == "Class 1", , drop = FALSE]
    expect_equal(nrow(reference_row), 1)
    expect_true(reference_row$is_reference[[1]])
    expect_equal(reference_row$n[[1]], 4)
    expect_equal(reference_row$events[[1]], 2)
    expect_equal(reference_row$non_events[[1]], 2)
    expect_equal(reference_row$event_rate_percent[[1]], 50)

    comparison_row <- term_support[term_support$variable == "biopsy1_gep_model" & term_support$level == "Class 2", , drop = FALSE]
    expect_equal(nrow(comparison_row), 1)
    expect_false(comparison_row$is_reference[[1]])
    expect_equal(comparison_row$n[[1]], 4)
    expect_equal(comparison_row$events[[1]], 2)
    expect_equal(comparison_row$non_events[[1]], 2)
    expect_equal(comparison_row$event_rate_percent[[1]], 50)
})

# Test 4c: Competing risk analysis produces valid statistical results
# Test 4d: Time-specific analysis produces valid statistical results
# Eligibility keys on the Class 1/Class 2 availability group and the endpoint fields.
test_that("Eligibility filters properly exclude invalid data", {
    test_data <- tibble::tibble(
        gep_group = c("Class 1", "Class 2", GEP_GROUP_CLASS_PRAME_MISSING, GEP_GROUP_TEST_FAILED, GEP_GROUP_NOT_TESTED),
        biopsy1_gep_mfs = c(0.8, 0.4, 0.6, 0.5, 0.7),
        biopsy1_gep_mss = c(0.9, 0.6, 0.7, 0.8, 0.8),
        mets_free_at_baseline = rep(TRUE, 5),
        tt_mets_months_analysis = c(24, 36, 18, 12, 48),
        mets_event_analysis = c(0L, 1L, 0L, 1L, 0L),
        tt_death_years = c(2, 3, 1.5, 1, 4),
        melanoma_death_event = c(0, 1, 0, 1, 0),
        competing_death_event = c(0, 0, 1, 0, 0),
        prame_status = c("Negative", "Positive", "Unknown", "Not Available", "Not Available")
    )

    flagged <- derive_gep_analysis_flags(test_data)

    expect_identical(flagged$mfs_analysis_eligible, c(TRUE, TRUE, FALSE, FALSE, FALSE))
    expect_identical(flagged$mss_analysis_eligible, c(TRUE, TRUE, FALSE, FALSE, FALSE))
    expect_identical(flagged$gep_validation_set, c("Eligible", "Eligible", "No GEP Data", "No GEP Data", "No GEP Data"))
})

test_that("Eligibility removes every non-definitive group from the cohort analyses", {
    cohort_names <- c(
        "uveal_melanoma_full_cohort",
        "uveal_melanoma_restricted_cohort",
        "uveal_melanoma_gksrs_only_cohort"
    )

    for (dataset_name in cohort_names) {
        actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, paste0(dataset_name, ".rds")))
        flagged <- derive_gep_analysis_flags(actual_data)
        mfs_data <- flagged %>% filter(mfs_analysis_eligible)
        mss_data <- flagged %>% filter(mss_analysis_eligible)

        expect_true(all(as.character(mfs_data$gep_group) %in% GEP_DEFINITIVE_GROUP_LEVELS),
            info = sprintf("MFS-eligible rows should retain only Class 1/Class 2 for %s", dataset_name)
        )
        expect_true(all(as.character(mss_data$gep_group) %in% GEP_DEFINITIVE_GROUP_LEVELS),
            info = sprintf("MSS-eligible rows should retain only Class 1/Class 2 for %s", dataset_name)
        )
        expect_true(all(as.character(mfs_data$biopsy1_gep) %in% GEP_VALIDATION_GROUP_LEVELS),
            info = sprintf("MFS-eligible rows should carry a Class x PRAME label for %s", dataset_name)
        )
    }
})

test_that("Canonical GEP variables retain original levels without cohort-wide collapse", {
    actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))
    expect_false(file.exists(file.path(PROCESSED_DATA_DIR, "other_map.rds")))
    expect_true(all(c("biopsy1_gep_raw", "biopsy1_gep", "gep_group", "prame_status", "gep12_prame_status", "location") %in% names(actual_data)))
    expect_false("Other" %in% levels(actual_data$location))
    expect_false("Other" %in% levels(actual_data$biopsy1_gep_raw))
    expect_false(any(as.character(actual_data$biopsy1_gep_raw) == "Other", na.rm = TRUE))
    expect_true(any(as.character(actual_data$location) == "Cilio-Choroidal"))
    expect_identical(levels(actual_data$biopsy1_gep), GEP_DETAIL_LEVELS)
    expect_identical(levels(actual_data$gep_group), GEP_GROUP_LEVELS)
    expect_true(any(as.character(actual_data$gep_group) == GEP_GROUP_NOT_TESTED))
})

test_that("MFS eligibility excludes metastasis present at baseline", {
    test_data <- tibble::tibble(
        gep_group = c("Class 1", "Class 1"),
        biopsy1_gep_mfs = c(0.80, 0.80),
        biopsy1_gep_mss = c(0.85, 0.85),
        mets_free_at_baseline = c(FALSE, TRUE),
        tt_mets_months_analysis = c(NA_real_, 36),
        mets_event_analysis = c(NA_integer_, 0L),
        tt_death_years = c(3, 3),
        melanoma_death_event = c(0L, 0L),
        competing_death_event = c(0L, 0L),
        prame_status = c("Positive", "Positive")
    )

    flagged <- derive_gep_analysis_flags(test_data)

    expect_false(flagged$mfs_analysis_eligible[1])
    expect_true(flagged$mfs_analysis_eligible[2])
    expect_true(all(flagged$mss_analysis_eligible))
})

test_that("Eligibility derivation fails when a required column is missing", {
    expect_error(
        derive_gep_analysis_flags(tibble::tibble(gep_group = "Class 1")),
        "GEP analysis flags require columns"
    )
})

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
        mets_event = c(0, 1, 0, 1),
        biopsy1_gep = factor(c("Class 1", "Class 1", "GEP Failed/Indeterminate", "GEP Failed/Indeterminate")),
        biopsy1_gep_model = factor(c("Class 1", "Class 1", "GEP Failed/Indeterminate", "GEP Failed/Indeterminate"))
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
            output_dirs = NULL,
            prefix = NULL
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
    expect_true(any(grepl("GEP Failed/Indeterminate", names(result$fit$strata))))
    expect_null(result$cox_model)
    expect_match(result$diagnostics$raw_model_output, "After sparse-level exclusions")
    expect_match(result$diagnostics$raw_model_output, "requires at least two non-missing")
})

test_that("risk table keeps displayed row labels aligned with their counts when legend order differs from fit order", {
    desired_order <- c(
        "Class 1 PRAME Negative",
        "Class 1 PRAME Positive",
        "Class 2 PRAME Negative",
        "Class 2 PRAME Positive",
        "GEP Failed/Indeterminate",
        "GEP Not Tested"
    )
    baseline_counts <- c(5, 4, 3, 2, 6, 7)
    names(baseline_counts) <- desired_order

    test_data <- tibble::tibble(
        tt_mets_months = seq_len(sum(baseline_counts)) * 6,
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
        output_dirs = NULL,
        prefix = NULL
    ))

    expect_equal(rev(extract_risk_table_y_limits(result$plot$table)), desired_order)
    expect_equal(extract_displayed_risk_counts(result$plot$table, time_point = 0), unname(baseline_counts[desired_order]))
})

test_that("risk table preserves requested row order for a subset of present strata", {
    desired_order <- c(
        "Class 1 PRAME Negative",
        "Class 1 PRAME Positive",
        "Class 2 PRAME Negative",
        "Class 2 PRAME Positive",
        "GEP Failed/Indeterminate",
        "GEP Not Tested"
    )
    baseline_counts <- c(
        "Class 1 PRAME Negative" = 5,
        "Class 2 PRAME Positive" = 2,
        "GEP Not Tested" = 7
    )
    present_order <- desired_order[desired_order %in% names(baseline_counts)]

    test_data <- tibble::tibble(
        tt_mets_months = seq_len(sum(baseline_counts)) * 6,
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
        output_dirs = NULL,
        prefix = NULL
    ))

    expect_equal(rev(extract_risk_table_y_limits(result$plot$table)), present_order)
    expect_equal(extract_displayed_risk_counts(result$plot$table, time_point = 0), unname(baseline_counts[present_order]))
})

test_that("four-group collapsed MFS KM retains both no-GEP strata as separate rows", {
    test_output_dir <- tempfile("objective4_four_group_km_")
    dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

    test_data <- tibble::tibble(
        id = seq_len(8),
        tt_mets_months = c(12, 24, 18, 36, 6, 30, 9, 42),
        mets_event = c(0, 1, 0, 1, 1, 0, 1, 0),
        mets_free_at_baseline = TRUE,
        tt_mets_months_analysis = c(12, 24, 18, 36, 6, 30, 9, 42),
        mets_event_analysis = c(0, 1, 0, 1, 1, 0, 1, 0),
        biopsy1_gep = c(
            "Class 1 PRAME Negative",
            "Class 1 PRAME Positive",
            "Class 2 PRAME Negative",
            "Class 2 PRAME Positive",
            "GEP Not Tested",
            "GEP Not Tested",
            "GEP Failed/Indeterminate",
            "GEP Failed/Indeterminate"
        ),
        gep_class_simple = c("Class 1", "Class 1", "Class 2", "Class 2", NA, NA, NA, NA)
    )

    result <- suppressWarnings(
        create_mfs_four_group_survival_curves(
            data = test_data,
            output_dir = test_output_dir,
            km_output_dir = test_output_dir,
            prefix = "unit_",
            dataset_name = "unit_test_four_group_km",
            return_plot = TRUE,
            save_plot = TRUE
        )
    )

    expect_equal(
        result$present_levels,
        c("Class 1", "Class 2", "GEP Not Tested", "GEP Failed/Indeterminate")
    )
    expect_true(all(
        c("GEP Not Tested", "GEP Failed/Indeterminate") %in%
            unique(as.character(stats::na.omit(result$plot_data$gep_km_simple)))
    ))
    expect_equal(
        rev(extract_risk_table_y_limits(result$plot$table)),
        c("Class 1", "Class 2", "GEP Not Tested", "GEP Failed/Indeterminate")
    )
    expect_equal(extract_displayed_risk_counts(result$plot$table, time_point = 0), c(2L, 2L, 2L, 2L))
    expect_true(file.exists(file.path(test_output_dir, "unit_mfs_four_group_gep_km.png")))
})

test_that("Cox diagnostics keep group counts and event rates in Raw_model_output", {
    test_output_dir <- file.path(TEST_OUTPUT_DIR, "objective4_group_counts_in_raw_output")
    dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

    test_data <- tibble::tibble(
        tt_mets_months = c(8, 16, 24, 32, 10, 20, 30, 40),
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

    raw_output <- readxl::read_excel(result$output_files$diagnostics_path, sheet = "Raw_model_output")
    expect_true(all(c(
        "group_n", "group_events", "group_non_events", "group_event_rate_pct",
        "reference_level", "reference_n", "reference_events", "reference_non_events",
        "reference_event_rate_pct"
    ) %in% names(raw_output)))
    expect_false("group_event_summary" %in% names(raw_output))

    reference_row <- raw_output[raw_output$row_type == "Reference Level" & raw_output$variable == "Class 1", , drop = FALSE]
    expect_equal(nrow(reference_row), 1)
    expect_equal(reference_row$group_n[[1]], 4)
    expect_equal(reference_row$group_events[[1]], 2)
    expect_equal(reference_row$group_non_events[[1]], 2)
    expect_equal(reference_row$group_event_rate_pct[[1]], 50)

    coefficient_row <- raw_output[raw_output$row_type == "Coefficient" & raw_output$variable == "Class 2", , drop = FALSE]
    expect_equal(nrow(coefficient_row), 1)
    expect_equal(coefficient_row$group_n[[1]], 4)
    expect_equal(coefficient_row$group_events[[1]], 2)
    expect_equal(coefficient_row$group_non_events[[1]], 2)
    expect_equal(coefficient_row$group_event_rate_pct[[1]], 50)
    expect_equal(coefficient_row$reference_level[[1]], "Class 1")
    expect_equal(coefficient_row$reference_n[[1]], 4)
    expect_equal(coefficient_row$reference_events[[1]], 2)
    expect_equal(coefficient_row$reference_non_events[[1]], 2)
    expect_equal(coefficient_row$reference_event_rate_pct[[1]], 50)
})

# Test 4c: Competing risk analysis produces valid statistical results
# Test 4d: Time-specific analysis produces valid statistical results
# Test that eligibility filters properly exclude invalid data
test_that("Eligibility filters properly exclude invalid data", {
    test_data <- tibble::tibble(
        biopsy1_gep = c(
            "Class 1 PRAME Negative",
            "Class 2 PRAME Positive",
            "GEP Failed/Indeterminate",
            "GEP Failed/Indeterminate",
            "GEP Not Tested"
        ),
        biopsy1_gep_raw = c(
            "Class_1A_PRAME_negative",
            "Class_2_PRAME_positive",
            "Class_1A_PRAME_not_reported",
            "Unknown",
            "No"
        ),
        gep_class_simple = c("Class 1", "Class 2", "Class 1", "GEP Failed/Indeterminate", "GEP Not Tested"),
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

    refreshed_data <- refresh_gep_analysis_flags(test_data)

    expect_true(refreshed_data$mfs_analysis_eligible[1],
        info = "Definitive Class 1 should be included in MFS analysis"
    )
    expect_true(refreshed_data$mfs_analysis_eligible[2],
        info = "Definitive Class 2 should be included in MFS analysis"
    )
    expect_false(refreshed_data$mfs_analysis_eligible[3],
        info = "Class 1 not reported should be excluded from MFS analysis"
    )
    expect_false(refreshed_data$mfs_analysis_eligible[4],
        info = "Unknown raw GEP rows should be excluded from MFS analysis"
    )
    expect_false(refreshed_data$mfs_analysis_eligible[5],
        info = "GEP not tested rows should be excluded from MFS analysis"
    )

    expect_true(refreshed_data$mss_analysis_eligible[1],
        info = "Definitive Class 1 should be included in MSS analysis"
    )
    expect_true(refreshed_data$mss_analysis_eligible[2],
        info = "Definitive Class 2 should be included in MSS analysis"
    )
    expect_false(refreshed_data$mss_analysis_eligible[3],
        info = "Class 1 not reported should be excluded from MSS analysis"
    )
    expect_false(refreshed_data$mss_analysis_eligible[4],
        info = "Unknown raw GEP rows should be excluded from MSS analysis"
    )
    expect_false(refreshed_data$mss_analysis_eligible[5],
        info = "GEP not tested rows should be excluded from MSS analysis"
    )
})

test_that("Objective 4 eligibility refresh removes failed rows from cohort analyses", {
    cohort_names <- c(
        "uveal_melanoma_full_cohort",
        "uveal_melanoma_restricted_cohort",
        "uveal_melanoma_gksrs_only_cohort"
    )

    for (dataset_name in cohort_names) {
        actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, paste0(dataset_name, ".rds")))
        refreshed_data <- refresh_gep_analysis_flags(actual_data)
        mfs_data <- refreshed_data %>% filter(mfs_analysis_eligible)
        mss_data <- refreshed_data %>% filter(mss_analysis_eligible)

        expect_false(any(as.character(mfs_data$biopsy1_gep) %in% c("GEP Failed/Indeterminate", "GEP Not Tested")),
            info = sprintf("MFS-eligible rows should exclude failed and not-tested labels for %s", dataset_name)
        )
        expect_false(any(as.character(mss_data$biopsy1_gep) %in% c("GEP Failed/Indeterminate", "GEP Not Tested")),
            info = sprintf("MSS-eligible rows should exclude failed and not-tested labels for %s", dataset_name)
        )
        expect_true(all(as.character(mfs_data$gep_class_simple) %in% GEP_DEFINITIVE_SIMPLE_LEVELS),
            info = sprintf("MFS-eligible rows should retain only definitive simple classes for %s", dataset_name)
        )
        expect_true(all(as.character(mss_data$gep_class_simple) %in% GEP_DEFINITIVE_SIMPLE_LEVELS),
            info = sprintf("MSS-eligible rows should retain only definitive simple classes for %s", dataset_name)
        )
    }
})

test_that("Canonical GEP variables retain original levels without cohort-wide collapse", {
    actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))
    expect_false(file.exists(file.path(PROCESSED_DATA_DIR, "other_map.rds")))
    expect_true(all(c("biopsy1_gep_raw", GEP_DISPLAY_VARIABLES, "location") %in% names(actual_data)))
    expect_false("Other" %in% levels(actual_data$location))
    expect_false("Other" %in% levels(actual_data$biopsy1_gep_raw))
    expect_false(any(as.character(actual_data$biopsy1_gep_raw) == "Other", na.rm = TRUE))
    expect_true(any(as.character(actual_data$location) == "Cilio-Choroidal"))
    expect_true(any(as.character(actual_data$biopsy1_gep) == "GEP Not Tested"))
})

test_that("Eligibility depends on definitive raw labels even when text raw retains definitive class", {
    test_data <- tibble::tibble(
        biopsy1_gep = factor(c("Class 1 PRAME Positive", "Class 1 PRAME Positive")),
        biopsy1_gep_raw = factor(c("Unknown", "Class_1A_PRAME_positive")),
        biopsy1_gep_text_raw = c("Class_1A_PRAME_positive", "Class_1A_PRAME_positive"),
        gep_class_simple = factor(c("Class 1", "Class 1")),
        biopsy1_gep_mfs = c(0.80, 0.80),
        biopsy1_gep_mss = c(0.85, 0.85),
        mets_free_at_baseline = c(TRUE, TRUE),
        tt_mets_months_analysis = c(36, 36),
        mets_event_analysis = c(0L, 0L),
        tt_death_years = c(3, 3),
        melanoma_death_event = c(0, 0),
        competing_death_event = c(0, 0),
        prame_status = factor(c("Positive", "Positive"))
    )

    refreshed_data <- refresh_gep_analysis_flags(test_data)

    expect_false(refreshed_data$mfs_analysis_eligible[1])
    expect_true(refreshed_data$mfs_analysis_eligible[2])
    expect_false(refreshed_data$mss_analysis_eligible[1])
    expect_true(refreshed_data$mss_analysis_eligible[2])
})

test_that("MFS eligibility excludes metastasis present at baseline", {
    test_data <- tibble::tibble(
        biopsy1_gep = c("Class 1 PRAME Positive", "Class 1 PRAME Positive"),
        biopsy1_gep_raw = c("Class_1A_PRAME_positive", "Class_1A_PRAME_positive"),
        gep_class_simple = c("Class 1", "Class 1"),
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

    refreshed_data <- refresh_gep_analysis_flags(test_data)

    expect_false(refreshed_data$mfs_analysis_eligible[1])
    expect_true(refreshed_data$mfs_analysis_eligible[2])
    expect_true(all(refreshed_data$mss_analysis_eligible))
})

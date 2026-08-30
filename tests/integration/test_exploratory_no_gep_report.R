# Tests for exploratory no-GEP reporting
library(dplyr)

test_that("shared endpoint datasets use corrected fields and eligibility", {
    prepared <- list(full_data = tibble::tibble(
        exploratory_gep_group = factor(c("Class 1", "Class 2", "GEP Not Tested", "GEP Failed/Indeterminate")),
        mets_free_at_baseline = c(TRUE, TRUE, TRUE, FALSE),
        tt_mets_months_analysis = c(12, 18, 72, 120),
        mets_event_analysis = c(0L, 1L, 1L, 1L),
        objective4_mfs_event_type = c(0L, 1L, 0L, 0L),
        tt_mets_months = c(99, 99, 99, 99),
        mets_event = c(1L, 1L, 0L, 0L),
        tt_death_months = c(12, 18, 24, 30),
        objective4_mss_event_type = c(0L, 1L, 2L, 1L)
    ))
    mfs <- prepare_exploratory_mfs_analysis_data(prepared)
    mss <- prepare_exploratory_mss_analysis_data(prepared)

    expect_equal(mfs$tt_mets_months_analysis, c(12, 18, 72))
    expect_equal(mfs$mets_event_analysis, c(0L, 1L, 1L))
    expect_false("objective4_mfs_event_type" %in% names(mfs))
    expect_true(all(mfs$mets_free_at_baseline))
    expect_true(all(mss$objective4_mss_event_type %in% c(0L, 1L, 2L)))
    expect_equal(mfs$tt_mets_months, c(99, 99, 99))
    expect_equal(mfs$mets_event, c(1L, 1L, 0L))
})

test_that("shared endpoint bundles retain one fit and global comparison result", {
    actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))
    prepared <- prepare_exploratory_no_gep_data(actual_data)
    mfs_data <- prepare_exploratory_mfs_analysis_data(prepared)
    mss_data <- prepare_exploratory_mss_analysis_data(prepared)
    mfs <- fit_exploratory_mfs_analysis(mfs_data)
    mss <- fit_exploratory_mss_cif(mss_data)

    expect_s3_class(mfs$fit, "survfit")
    expect_equal(nrow(mfs$data), nrow(mfs_data))
    expect_true(is.finite(mfs$global_test$p_value))
    expect_s3_class(mss$fit, "tidycuminc")
    expect_equal(nrow(mss$data), nrow(mss_data))
    expect_true(is.finite(mss$gray_test$p_value))
    expect_true(all(c("time", "outcome", "strata", "estimate", "n.risk", "n.censor") %in% names(mss$tidy)))
})

test_that("no-GEP MFS uses shared censor and risk-table machinery", {
    fixture <- tibble::tibble(
        exploratory_gep_group = factor(c("Class 1", "Class 2", "GEP Not Tested", "GEP Failed/Indeterminate")),
        mets_free_at_baseline = TRUE,
        tt_mets_months_analysis = c(12, 18, 24, 30),
        mets_event_analysis = c(0L, 1L, 0L, 1L),
        objective4_mfs_event_type = c(0L, 1L, 0L, 1L),
        tt_mets_months = c(99, 99, 99, 99),
        mets_event = c(1L, 1L, 1L, 1L)
    )
    result <- create_exploratory_mfs_km_plot(fixture, tempfile(fileext = ".png"), return_plot = TRUE)

    expect_true(sum(result$fit$n.censor) > 0)
    expect_true(!is.null(result$plot$table))
    expect_true(nrow(result$plot_data) > 0)
    expect_true(all(result$plot_data$tt_mets_months_analysis <= 30))
    expect_false(result$p_value_annotation)
})

test_that("standard GEP KM uses the shared incident-MFS representation", {
    fixture <- tibble::tibble(
        id = c("baseline", "late", "censored", "class2_event"),
        biopsy1_gep = c(
            "Class 1 PRAME Negative",
            "Class 1 PRAME Positive",
            "Class 1 PRAME Negative",
            "Class 2 PRAME Negative"
        ),
        gep_class_simple = factor(
            c("Class 1", "Class 1", "Class 1", "Class 2"),
            levels = c("Class 1", "Class 2", "GEP Failed/Indeterminate", "GEP Not Tested")
        ),
        exploratory_gep_group = factor(
            c("Class 1", "Class 1", "Class 1", "Class 2"),
            levels = c("Class 1", "Class 2", "GEP Failed/Indeterminate", "GEP Not Tested")
        ),
        mets_free_at_baseline = c(FALSE, TRUE, TRUE, TRUE),
        tt_mets_months_analysis = c(NA_real_, 72, 90, 18),
        mets_event_analysis = c(NA_integer_, 1L, 0L, 1L),
        objective4_mfs_event_type = c(NA_integer_, 0L, 0L, 1L),
        tt_mets_months = c(0, 12, 90, 18),
        mets_event = c(1L, 0L, 0L, 1L)
    )

    result <- create_mfs_collapsed_survival_curves(
        data = fixture,
        output_dir = tempfile(),
        prefix = "",
        subtitle_suffix = "test",
        output_filename = "unused.png",
        return_plot = TRUE,
        save_plot = FALSE
    )

    expect_false("baseline" %in% result$plot_data$id)
    expect_true("late" %in% result$plot_data$id)
    expect_equal(
        result$plot_data$tt_mets_months_analysis[result$plot_data$id == "late"],
        72
    )
    expect_equal(
        result$plot_data$mets_event_analysis[result$plot_data$id == "late"],
        1L
    )
    expect_true(any(result$fit$time >= 72 & result$fit$n.event == 1))
})

test_that("standard GEP and no-GEP KM preparation retain the same incident-MFS rows", {
    fixture <- tibble::tibble(
        id = c("baseline", "late", "censored", "class2_event"),
        biopsy1_gep = c(
            "Class 1 PRAME Negative",
            "Class 1 PRAME Positive",
            "Class 1 PRAME Negative",
            "Class 2 PRAME Negative"
        ),
        gep_class_simple = factor(
            c("Class 1", "Class 1", "Class 1", "Class 2"),
            levels = c("Class 1", "Class 2", "GEP Failed/Indeterminate", "GEP Not Tested")
        ),
        exploratory_gep_group = factor(
            c("Class 1", "Class 1", "Class 1", "Class 2"),
            levels = c("Class 1", "Class 2", "GEP Failed/Indeterminate", "GEP Not Tested")
        ),
        mets_free_at_baseline = c(FALSE, TRUE, TRUE, TRUE),
        tt_mets_months_analysis = c(NA_real_, 72, 90, 18),
        mets_event_analysis = c(NA_integer_, 1L, 0L, 1L),
        objective4_mfs_event_type = c(NA_integer_, 0L, 0L, 1L),
        tt_mets_months = c(0, 12, 90, 18),
        mets_event = c(1L, 0L, 0L, 1L)
    )

    standard <- prepare_incident_mfs_km_data(fixture)
    no_gep <- prepare_exploratory_mfs_analysis_data(list(full_data = fixture))

    expect_setequal(standard$id, no_gep$id)
    expect_false("baseline" %in% standard$id)
    expect_true("late" %in% standard$id)
    expect_equal(
        standard %>% dplyr::arrange(.data$id) %>% dplyr::select(id, tt_mets_months_analysis, mets_event_analysis),
        no_gep %>% dplyr::arrange(.data$id) %>% dplyr::select(id, tt_mets_months_analysis, mets_event_analysis),
        ignore_attr = TRUE
    )
})

test_that("poster GEP KM preparation excludes baseline metastasis and keeps late events", {
    fixture <- tibble::tibble(
        id = c("baseline", "late", "censored"),
        gep_class_simple = factor(c("Class 1", "Class 2", "Class 1"), levels = c("Class 1", "Class 2")),
        mets_free_at_baseline = c(FALSE, TRUE, TRUE),
        tt_mets_months_analysis = c(NA_real_, 72, 90),
        mets_event_analysis = c(NA_integer_, 1L, 0L),
        tt_mets_months = c(0, 12, 90),
        mets_event = c(1L, 0L, 0L)
    )

    prepared <- prepare_mfs_simple_binary_poster_km_data(fixture)

    expect_false("baseline" %in% prepared$id)
    expect_true("late" %in% prepared$id)
    expect_equal(prepared$tt_mets_months_analysis[prepared$id == "late"], 72)
    expect_equal(prepared$mets_event_analysis[prepared$id == "late"], 1L)
})

test_that("MSS CIF uses Aalen-Johansen coding and a single shared fit", {
    fixture <- tibble::tibble(
        exploratory_gep_group = factor(rep(c("Class 1", "Class 2", "GEP Not Tested", "GEP Failed/Indeterminate"), each = 4)),
        tt_death_months = c(12, 18, 24, 30, 10, 20, 35, 40, 8, 16, 28, 44, 14, 22, 32, 48),
        objective4_mss_event_type = c(0L, 1L, 2L, 0L, 0L, 1L, 0L, 2L, 0L, 1L, 2L, 0L, 0L, 1L, 0L, 2L)
    )
    fitted <- fit_exploratory_mss_cif(fixture)
    result <- create_exploratory_mss_cif_plot(
        fixture,
        tempfile(fileext = ".png"),
        analysis_fit = fitted,
        return_plot = TRUE
    )

    expect_true(inherits(fitted$fit, "tidycuminc"))
    expect_true(is.finite(fitted$gray_test$p_value))
    expect_true(all(fitted$data$.mss_outcome %in% c("censored", "melanoma_death", "other_death")))
    expect_true(all(c("time", "outcome", "strata", "estimate", "n.risk", "n.censor") %in% names(fitted$tidy)))
    expect_identical(result$fit, fitted$fit)
    expect_s3_class(result$plot, "ggplot")
})

test_that("unsupported MSS comparisons are explicit and non-fatal", {
    fixture <- tibble::tibble(
        exploratory_gep_group = factor(c("Class 1", "Class 1")),
        tt_death_months = c(12, 18),
        objective4_mss_event_type = c(0L, 0L)
    )
    fitted <- fit_exploratory_mss_cif(fixture)
    expect_true(is.na(fitted$gray_test$p_value))
    expect_true(fitted$gray_test$status %in% c("skipped", "no_event_of_interest", "fit_failed"))
})

test_that("exploratory no-GEP dataset preparation isolates reference and scoring cohorts", {
    actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))

    prepared <- prepare_exploratory_no_gep_data(actual_data)

    expect_identical(levels(prepared$full_data$exploratory_gep_group), levels(actual_data$gep_class_simple))
    expect_identical(levels(prepared$full_data$sex), levels(actual_data$sex))
    expect_identical(levels(prepared$full_data$location), levels(actual_data$location))
    expect_identical(levels(prepared$full_data$initial_t_stage_simple), levels(actual_data$initial_t_stage_simple))
    expect_identical(levels(prepared$full_data$internal_reflectivity), levels(actual_data$internal_reflectivity))
    expect_identical(levels(prepared$full_data$srf), levels(actual_data$srf))

    expect_true(all(prepared$definitive_reference$exploratory_gep_group %in% c("Class 1", "Class 2")))
    expect_true(all(prepared$no_gep_scoring$exploratory_gep_group %in% c("GEP Failed/Indeterminate", "GEP Not Tested")))
    expect_true(all(unique(prepared$no_gep_scoring$no_gep_group) %in% c("GEP Failed/Indeterminate", "GEP Not Tested")))
    expect_true(all(prepared$predictors %in% names(prepared$definitive_reference)))
    expect_true(any(prepared$full_data$ciliary_involvement == 1, na.rm = TRUE))

    cilio_rows <- prepared$full_data %>%
        dplyr::filter(.data$location == "Cilio-Choroidal")
    choroidal_rows <- prepared$full_data %>%
        dplyr::filter(.data$location == "Choroidal")

    expect_true(nrow(cilio_rows) > 0)
    expect_true(all(cilio_rows$ciliary_involvement == 1))
    expect_true(all(choroidal_rows$ciliary_involvement == 0))

    screening_predictors <- prepared$predictor_screening %>%
        dplyr::filter(.data$status == "retained") %>%
        dplyr::pull(.data$predictor)
    expect_setequal(prepared$predictors, screening_predictors)
})

test_that("exploratory no-GEP KM verification checks displayed counts and cohort counts", {
    actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))
    prepared <- prepare_exploratory_no_gep_data(actual_data)

    verification <- verify_exploratory_no_gep_km_fix(
        actual_data,
        prepared_data = prepared,
        expected_group_counts = prepared$group_snapshot
    )
    expect_equal(verification$observed_n, prepared$group_snapshot$expected_n)
    expect_equal(verification$expected_n, prepared$group_snapshot$expected_n)
    expect_true(all(verification$status == "matched"))
    expect_equal(as.character(stats::na.omit(verification$simple_km_display_order)), c("Class 1", "Class 2", "GEP Not Tested"))
    expect_equal(as.integer(stats::na.omit(verification$simple_km_displayed_n)), c(58L, 27L, 161L))

    modified_data <- actual_data
    class1_idx <- which(as.character(modified_data$exploratory_gep_group) == "Class 1")[1]
    modified_data$exploratory_gep_group[class1_idx] <- factor(
        "Class 2",
        levels = levels(modified_data$exploratory_gep_group)
    )

    expect_error(
        verify_exploratory_no_gep_km_fix(
            modified_data,
            expected_group_counts = prepared$group_snapshot
        ),
        "group counts do not match the prepared-dataset snapshot"
    )
})

test_that("exploratory no-GEP preparation fails fast when Objective 0 columns are absent", {
    broken_data <- tibble::tibble(
        gep_class_simple = factor(c("Class 1", "Class 2"))
    )

    expect_error(
        prepare_exploratory_no_gep_data(broken_data),
        "expects the Objective 0 prepared cohort"
    )
})

test_that("exploratory direct targets preserve MFS eligibility and MSS competing outcomes", {
    actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))
    fixture <- actual_data
    initially_prepared <- prepare_exploratory_no_gep_data(actual_data)

    baseline_mets_id <- initially_prepared$mfs_model_data$id[[1]]
    competing_id <- initially_prepared$mss_model_data$id[[1]]
    early_censor_id <- initially_prepared$mss_model_data$id[[2]]
    baseline_row <- match(baseline_mets_id, fixture$id)
    mss_rows_index <- match(c(competing_id, early_censor_id), fixture$id)

    fixture$mets_free_at_baseline[[baseline_row]] <- FALSE
    fixture$tt_mets_months_analysis[[baseline_row]] <- NA_real_
    fixture$mfs_event_5yr[[baseline_row]] <- 1L

    fixture$tt_death_months[mss_rows_index] <- 12
    fixture$melanoma_death_event[mss_rows_index] <- 0L
    fixture$competing_death_event[mss_rows_index] <- c(1L, 0L)
    fixture$mss_event_5yr[mss_rows_index] <- 0L

    prepared <- prepare_exploratory_no_gep_data(fixture)

    expect_false(baseline_mets_id %in% prepared$mfs_model_data$id)
    expect_true(all(prepared$mfs_model_data$mets_free_at_baseline))
    expect_true(all(is.finite(prepared$mfs_model_data$tt_mets_months_analysis)))

    mss_rows <- prepared$mss_model_data %>%
        dplyr::filter(.data$id %in% c(competing_id, early_censor_id)) %>%
        dplyr::arrange(match(.data$id, c(competing_id, early_censor_id)))
    mss_status <- derive_horizon_status(
        mss_rows$tt_death_months,
        mss_rows$objective4_mss_event_type,
        60
    )

    expect_identical(mss_rows$objective4_mss_event_type, c(2L, 0L))
    expect_identical(mss_status$horizon_event, c(0L, NA_integer_))
    expect_identical(mss_status$known_status, c(TRUE, FALSE))

    mss_weight_check <- derive_fold_ipcw_payload(
        training = tibble::tibble(
            time = c(12, 18, 24, 36, 60, 72),
            event_type = c(1L, 2L, 0L, 1L, 0L, 0L)
        ),
        assessment = tibble::tibble(
            time = mss_rows$tt_death_months,
            event_type = mss_rows$objective4_mss_event_type
        ),
        time_var = "time",
        event_type_var = "event_type",
        horizon_months = 60
    )
    expect_gt(mss_weight_check$assessment$ipcw_weight[[1]], 0)
    expect_identical(mss_weight_check$assessment$ipcw_weight[[2]], 0)
})

test_that("scoped OOF performance uses one keyed prediction set and fails closed", {
    oof <- tibble::tibble(
        stable_id = rep(as.character(1:6), 2),
        repeat_id = rep(1:2, each = 6),
        prediction = rep(c(0.9, 0.7, 0.6, 0.4, 0.3, 0.1), 2),
        horizon_event = rep(c(1L, 1L, 0L, 0L, NA_integer_, 0L), 2),
        ipcw_weight = rep(c(2, 1, 1, 2, 0, 1), 2),
        exploratory_gep_group = rep(
            c("Class 2", "GEP Failed/Indeterminate", "Class 1", "GEP Not Tested", "GEP Not Tested", "Class 1"),
            2
        )
    )

    scoped <- summarize_scoped_ipcw_oof_performance(
        oof,
        group_var = "exploratory_gep_group"
    )

    expect_setequal(unique(scoped$performance_scope), c("Overall", "No GEP"))
    expect_equal(nrow(scoped), 4)
    expect_true(all(scoped$evaluation_method == "outer-training-fold IPCW weighted OOF AUC/Brier/calibration"))
    expect_true(all(scoped$not_tested_n == 2L))
    expect_true(all(scoped$failed_indeterminate_n == 1L))
    expect_true(all(scoped$auc_status[scoped$performance_scope == "Overall"] == "ok"))
    expect_true(all(scoped$auc_status[scoped$performance_scope == "No GEP"] == "ok"))

    no_gep_one_class <- oof %>%
        dplyr::mutate(
            horizon_event = dplyr::if_else(
                .data$exploratory_gep_group %in% c("GEP Failed/Indeterminate", "GEP Not Tested") & .data$ipcw_weight > 0,
                0L,
                .data$horizon_event
            )
        )
    unsupported <- summarize_scoped_ipcw_oof_performance(
        no_gep_one_class,
        group_var = "exploratory_gep_group"
    )
    expect_true(all(is.na(unsupported$cv_auc[unsupported$performance_scope == "No GEP"])))
    expect_true(all(unsupported$auc_status[unsupported$performance_scope == "No GEP"] == "unsupported_no_weighted_cases"))

    absent_target_population <- summarize_scoped_ipcw_oof_performance(
        oof %>% dplyr::filter(!.data$exploratory_gep_group %in% c("GEP Failed/Indeterminate", "GEP Not Tested")),
        group_var = "exploratory_gep_group"
    )
    absent_rows <- absent_target_population %>% dplyr::filter(.data$performance_scope == "No GEP")
    expect_equal(nrow(absent_rows), 2)
    expect_true(all(absent_rows$auc_status == "unsupported_no_positive_weight"))
    expect_true(all(is.na(absent_rows$cv_auc)))
})

test_that("exploratory horizon summaries use censoring-aware event estimates", {
    prediction_data <- tibble::tibble(
        no_gep_group = c("GEP Failed/Indeterminate", "GEP Failed/Indeterminate", "GEP Failed/Indeterminate"),
        tt_mets_months = c(48, 24, 24),
        tt_mets_months_analysis = c(48, 24, 24),
        mets_event_analysis = c(1L, 0L, 0L),
        mets_event = c(1, 0, 0),
        mets_free_at_baseline = TRUE,
        objective4_mfs_event_type = c(1L, 0L, 0L),
        tt_death_months = c(48, 24, 24),
        melanoma_death_event = c(1, 0, 0),
        competing_death_event = c(0, 0, 0),
        mfs_event_5yr = c(1L, 0L, 0L),
        mss_event_5yr = c(1L, 0L, 0L),
        surrogate_class2_probability = c(0.5, 0.4, 0.6),
        predicted_mfs_5yr_risk = c(0.5, 0.4, 0.6),
        predicted_mss_5yr_risk = c(0.5, 0.4, 0.6)
    )

    summary_tbl <- summarize_no_gep_predictions(prediction_data)

    expect_equal(summary_tbl$mfs_observed_method[[1]], "kaplan_meier_at_horizon")
    expect_equal(summary_tbl$mss_observed_method[[1]], "aalen_johansen_cif_at_horizon")
    expect_true(summary_tbl$observed_mfs_5yr_event_rate[[1]] > mean(prediction_data$mfs_event_5yr))
    expect_true(summary_tbl$observed_mss_5yr_event_rate[[1]] > mean(prediction_data$mss_event_5yr))
})

test_that("exploratory pooled summaries tolerate bins with no melanoma failures", {
    prediction_data <- tibble::tibble(
        no_gep_group = c("GEP Failed/Indeterminate", "GEP Failed/Indeterminate", "GEP Not Tested", "GEP Not Tested"),
        surrogate_probability_bin = c("Low", "High", "Low", "High"),
        mfs_risk_bin = c("Low", "High", "Low", "High"),
        mss_risk_bin = c("Low", "High", "Low", "High"),
        surrogate_class2_probability = c(0.2, 0.8, 0.3, 0.7),
        predicted_mfs_5yr_risk = c(0.1, 0.6, 0.2, 0.5),
        predicted_mss_5yr_risk = c(0.05, 0.4, 0.1, 0.3),
        tt_mets_months = c(24, 48, 36, 60),
        tt_mets_months_analysis = c(24, 48, 36, 60),
        mets_event_analysis = c(0L, 1L, 0L, 1L),
        mets_event = c(0, 1, 0, 1),
        mets_free_at_baseline = TRUE,
        objective4_mfs_event_type = c(0L, 1L, 0L, 1L),
        tt_death_months = c(24, 48, 36, 60),
        melanoma_death_event = c(0, 1, 0, 0),
        competing_death_event = c(0, 0, 0, 0),
        mfs_event_5yr = c(0L, 1L, 0L, 1L),
        mss_event_5yr = c(0L, 1L, 0L, 0L)
    )

    pooled_summary <- summarize_pooled_no_gep_sensitivity(prediction_data)
    risk_strata_summary <- summarize_no_gep_risk_strata(prediction_data)

    expect_true(nrow(pooled_summary) > 0)
    expect_true(nrow(risk_strata_summary) > 0)
    expect_true(all(stats::na.omit(pooled_summary$mss_observed_method) == "aalen_johansen_cif_at_horizon"))
    expect_true(all(stats::na.omit(risk_strata_summary$MSS_Observed_Method) == "aalen_johansen_cif_at_horizon"))
    expect_equal(
        risk_strata_summary$Observed_MSS_5yr_Event_Rate[
            risk_strata_summary$No_GEP_Group == "GEP Not Tested" &
                risk_strata_summary$Analysis == "Direct_60mo_Melanoma_Death_Cumulative_Incidence_Risk"
        ],
        c(0, 0)
    )
})

test_that("exploratory pooled summaries omit unbinned predictions", {
    prediction_data <- tibble::tibble(
        no_gep_group = rep("GEP Not Tested", 3),
        surrogate_probability_bin = c("Low", "High", "High"),
        mfs_risk_bin = c("Low", "High", NA_character_),
        mss_risk_bin = c("Low", "High", "High"),
        surrogate_class2_probability = c(0.2, 0.8, 0.7),
        predicted_mfs_5yr_risk = c(0.1, 0.6, NA_real_),
        predicted_mss_5yr_risk = c(0.05, 0.4, 0.3),
        tt_mets_months = c(24, 48, NA_real_),
        tt_mets_months_analysis = c(24, 48, NA_real_),
        mets_event_analysis = c(0L, 1L, NA_integer_),
        mets_event = c(0, 1, NA_integer_),
        mets_free_at_baseline = c(TRUE, TRUE, FALSE),
        objective4_mfs_event_type = c(0L, 1L, NA_integer_),
        tt_death_months = c(24, 48, 36),
        melanoma_death_event = c(0, 0, 0),
        competing_death_event = c(0, 0, 0),
        mfs_event_5yr = c(0L, 1L, NA_integer_),
        mss_event_5yr = c(0L, 0L, 0L)
    )

    pooled_summary <- summarize_pooled_no_gep_sensitivity(prediction_data)
    direct_mfs_summary <- pooled_summary %>%
        dplyr::filter(.data$analysis == "Direct_MFS_5yr_Risk")

    expect_false(any(is.na(pooled_summary$bin)))
    expect_true(all(is.finite(pooled_summary$mean_predicted)))
    expect_equal(nrow(direct_mfs_summary), 2L)
    expect_equal(sum(direct_mfs_summary$n), 2L)
})

test_that("exploratory predictor screening ignores unused Other factor levels", {
    exploratory_data <- tibble::tibble(
        exploratory_gep_group = factor(
            c(rep("Class 1", 41), rep("Class 2", 40)),
            levels = c("Class 1", "Class 2", "GEP Failed/Indeterminate", "GEP Not Tested")
        ),
        mfs_event_5yr = c(rep(0L, 50), rep(1L, 31)),
        mss_event_5yr = c(rep(0L, 55), rep(1L, 26)),
        tt_mets_months_analysis = rep(72, 81),
        mets_free_at_baseline = TRUE,
        objective4_mfs_event_type = as.integer(c(rep(0L, 50), rep(1L, 31))),
        tt_death_months = rep(72, 81),
        objective4_mss_event_type = as.integer(c(rep(0L, 55), rep(1L, 26))),
        location = factor(
            c(rep("Choroidal", 74), rep("Cilio-Choroidal", 7)),
            levels = c("Choroidal", "Cilio-Choroidal", "Other")
        )
    )

    screening <- screen_exploratory_predictors(
        data = exploratory_data,
        candidate_predictors = "location",
        factor_predictors = "location",
        completeness_threshold = 0.9,
        min_level_count = 5
    )

    expect_equal(screening$status, "retained")
    expect_match(screening$reason, "passes completeness and sparse-level screening", fixed = TRUE)

    model_data <- build_exploratory_model_dataset(
        data = exploratory_data,
        predictors = "location",
        factor_predictors = "location",
        outcome_var = "class2_outcome",
        group_levels = c("Class 1", "Class 2")
    )

    expect_equal(levels(model_data$location), c("Choroidal", "Cilio-Choroidal"))
    expect_false("Other" %in% levels(model_data$location))
})

test_that("exploratory no-GEP report writes workbook, summary, and plots", {
    actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))
    test_output_dir <- file.path(TEST_OUTPUT_DIR, "exploratory_no_gep")
    withr::defer(unlink(test_output_dir, recursive = TRUE), teardown_env())

    results <- run_exploratory_no_gep_report(
        dataset_name = "uveal_melanoma_full_cohort",
        output_dir = test_output_dir,
        verify_km_fix = FALSE,
        data = actual_data
    )

    expect_true(file.exists(results$output_paths$workbook))
    expect_true(file.exists(results$output_paths$summary))
    expect_true(file.exists(results$output_paths$mfs_km))
    expect_true(file.exists(results$output_paths$mss_cif))
    expect_true(file.exists(results$output_paths$surrogate_density))

    workbook_sheets <- openxlsx::getSheetNames(results$output_paths$workbook)
    expect_equal(
        workbook_sheets,
        c(
            "Start_Here",
            "Key_Findings_5yr",
            "Risk_Ladder_5yr",
            "No_GEP_Subgroups",
            "Follow_Up_Context",
            "Model_Performance",
            "Parsimonious_Sensitivity",
            "Surrogate_Model_Coefficients",
            "Direct_MFS_Coefficients",
            "Direct_MSS_Coefficients",
            "Model_Calibration",
            "Predictor_Contribution",
            "Overlap_Diagnostics",
            "Baseline_Comparisons",
            "Data_Audit",
            "No_GEP_Predictions",
            "Sensitivity_Pooled_No_GEP",
            "KM_Corrected_MFS",
            "KM_Corrected_MSS"
        )
    )

    start_here_sheet <- openxlsx::read.xlsx(results$output_paths$workbook, sheet = "Start_Here")
    key_findings_sheet <- openxlsx::read.xlsx(results$output_paths$workbook, sheet = "Key_Findings_5yr")
    risk_ladder_sheet <- openxlsx::read.xlsx(results$output_paths$workbook, sheet = "Risk_Ladder_5yr")
    model_performance_sheet <- openxlsx::read.xlsx(results$output_paths$workbook, sheet = "Model_Performance")
    km_mfs_sheet <- openxlsx::read.xlsx(results$output_paths$workbook, sheet = "KM_Corrected_MFS")
    km_mss_sheet <- openxlsx::read.xlsx(results$output_paths$workbook, sheet = "KM_Corrected_MSS")

    expect_false(any(c("section", "item", "detail", "guide_text") %in% names(key_findings_sheet)))
    expect_false(any(c("section", "item", "detail", "guide_text") %in% names(risk_ladder_sheet)))
    expect_false(any(c("section", "item", "detail", "guide_text") %in% names(model_performance_sheet)))
    expect_equal(nrow(key_findings_sheet), 4)
    expect_equal(
        key_findings_sheet$group,
        c("Class 1", "GEP Not Tested", "GEP Failed/Indeterminate", "Class 2")
    )
    expect_true(
        "median_predicted_60mo_melanoma_death_cumulative_incidence_risk" %in%
            names(key_findings_sheet)
    )
    expect_false("median_predicted_5yr_mss_risk" %in% names(key_findings_sheet))
    expect_true(all(c("section", "label", "value") %in% names(start_here_sheet)))
    expect_true(all(c("group", "n", "observed_5yr_mfs_event_rate", "median_predicted_5yr_mfs_risk") %in% names(risk_ladder_sheet)))
    expect_true(all(c(
        "model", "model_method", "reported_risk_scale", "cv_auc",
        "cv_auc_stability_interval", "calibration_slope_stability_interval",
        "practical_read"
    ) %in% names(model_performance_sheet)))
    expect_true(all(c(
        "cumulative_incidence_conf_low",
        "cumulative_incidence_conf_high",
        "mss_probability_conf_low",
        "mss_probability_conf_high",
        "gray_test_global_curve_p_value",
        "gray_test_global_curve_test_status",
        "gray_test_global_curve_test_reason"
    ) %in% names(km_mss_sheet)))

    summary_text <- paste(readLines(results$output_paths$summary), collapse = "\n")
    expect_match(summary_text, "descriptive only", fixed = TRUE)
    expect_match(summary_text, "homogeneous intermediate-risk group", fixed = TRUE)
    expect_match(summary_text, "## Follow-Up Context", fixed = TRUE)
    expect_match(summary_text, "no-GEP scoring cohort", fixed = TRUE)
    expect_match(summary_text, "## Key Findings at 5 Years", fixed = TRUE)
    expect_match(summary_text, "95% repeated-partition stability interval", fixed = TRUE)
    expect_match(summary_text, "censoring weights are estimated in each outer training fold", fixed = TRUE)
    expect_match(summary_text, "60-month melanoma-death cumulative-incidence risk", fixed = TRUE)
    expect_match(
        summary_text,
        "Median predicted 60-month melanoma-death cumulative-incidence risk",
        fixed = TRUE
    )
    reader_facing_model_text <- paste(
        summary_text,
        paste(start_here_sheet$value, collapse = "\n"),
        paste(model_performance_sheet$model, collapse = "\n"),
        paste(model_performance_sheet$prediction_target, collapse = "\n"),
        collapse = "\n"
    )
    expect_false(grepl("predicted 5-year MSS risk", reader_facing_model_text, fixed = TRUE))
    expect_false(grepl("Direct 5-year MSS", reader_facing_model_text, fixed = TRUE))
    expect_false(grepl("melanoma-specific model", reader_facing_model_text, fixed = TRUE))
    expect_match(summary_text, "## Parsimonious Sensitivity Check", fixed = TRUE)
    expect_match(summary_text, "## Retained Baseline Predictors", fixed = TRUE)
    expect_match(summary_text, "Std. coef.", fixed = TRUE)
    expect_match(summary_text, "ranked highly in the penalized model", fixed = TRUE)
    expect_match(summary_text, "P\\(Class 2-like \\| baseline features\\)")
    expect_match(summary_text, "Cilio-Choroidal", fixed = TRUE)
    expect_match(summary_text, "0-1 probability scale", fixed = TRUE)
    expect_match(summary_text, "overlap diagnostic", ignore.case = TRUE)

    expect_s3_class(results$surrogate_model$model, "cv.glmnet")
    expect_s3_class(results$direct_models$mfs$model, "cv.glmnet")
    expect_s3_class(results$direct_models$mss$model, "cv.glmnet")
    expect_s3_class(results$parsimonious_models$mfs$model, "cv.glmnet")
    expect_s3_class(results$parsimonious_models$mss$model, "cv.glmnet")
    expect_true("calibration_status" %in% names(results$surrogate_model$metrics))
    expect_true("calibration_status" %in% names(results$direct_models$mfs$metrics))
    expect_true("calibration_status" %in% names(results$direct_models$mss$metrics))
    expect_equal(results$direct_models$mfs$metrics$model_mode_used[[1]], "ipcw_horizon_mfs")
    expect_equal(results$direct_models$mss$metrics$model_mode_used[[1]], "ipcw_horizon_competing_risk_mss")
    expect_equal(
        results$direct_models$mss$metrics$prediction_target[[1]],
        "60-month melanoma-death cumulative-incidence risk"
    )
    expect_setequal(
        unique(results$direct_models$mfs$scoped_oof_performance$performance_scope),
        c("Overall", "No GEP")
    )
    expect_setequal(
        unique(results$direct_models$mss$scoped_oof_performance$performance_scope),
        c("Overall", "No GEP")
    )
    expect_true(all(c(
        "performance_scope", "evaluation_method", "prediction_target",
        "metric_status", "positive_weight_n", "case_n", "control_n",
        "weighted_cases", "weighted_controls",
        "failed_indeterminate_n", "not_tested_n", "uncertainty_method"
    ) %in% names(model_performance_sheet)))
    expect_true(all(c(
        "log_rank_global_curve_p_value", "log_rank_global_curve_test_status",
        "log_rank_global_curve_test_reason"
    ) %in% names(km_mfs_sheet)))
    expect_true(all(c(
        "gray_test_global_curve_p_value", "gray_test_global_curve_test_status",
        "gray_test_global_curve_test_reason"
    ) %in% names(km_mss_sheet)))
    expect_equal(unique(km_mfs_sheet$log_rank_global_curve_p_value), results$mfs_analysis$global_test$p_value)
    expect_equal(unique(km_mss_sheet$gray_test_global_curve_p_value), results$mss_analysis$gray_test$p_value)
    expect_false("model_fallback_reason" %in% names(results$direct_models$mfs$metrics))
    expect_false("model_fallback_reason" %in% names(results$direct_models$mss$metrics))
    expect_false("raw_backtest" %in% names(results$direct_models$mfs))
    expect_false("raw_backtest" %in% names(results$direct_models$mss))
    expect_equal(
        nrow(results$direct_models$mfs$oof_predictions),
        nrow(prepare_exploratory_no_gep_data(actual_data)$mfs_model_data) * GEP_EXPLORATORY_CV_REPEATS
    )
    expect_equal(
        nrow(results$direct_models$mss$oof_predictions),
        nrow(prepare_exploratory_no_gep_data(actual_data)$mss_model_data) * GEP_EXPLORATORY_CV_REPEATS
    )
    prepared_for_oof <- prepare_exploratory_no_gep_data(actual_data)
    expect_setequal(
        unique(results$direct_models$mfs$oof_predictions$stable_id),
        as.character(prepared_for_oof$mfs_model_data$id)
    )
    expect_false(any(
        as.character(actual_data$id[!actual_data$mets_free_at_baseline]) %in%
            results$direct_models$mfs$oof_predictions$stable_id
    ))
    expect_identical(
        anyDuplicated(results$direct_models$mfs$oof_predictions[c("repeat_id", "stable_id")]),
        0L
    )
    expect_identical(
        anyDuplicated(results$direct_models$mss$oof_predictions[c("repeat_id", "stable_id")]),
        0L
    )
    expect_true(all(c("cv_auc_ci_lower", "cv_auc_ci_upper", "cv_repeats") %in% names(results$direct_models$mfs$metrics)))
    expect_true(all(c("cv_auc_ci_lower", "cv_auc_ci_upper", "cv_repeats") %in% names(results$direct_models$mss$metrics)))
    expect_true(all(unique(results$no_gep_predictions$no_gep_group) %in% c("GEP Failed/Indeterminate", "GEP Not Tested")))
    expect_true("start_here" %in% names(results))
    expect_true("key_findings_5yr" %in% names(results))
    expect_true("no_gep_subgroups" %in% names(results))
    expect_true("model_performance" %in% names(results))
    expect_true("surrogate_model_coefficients" %in% names(results))
    expect_true("model_calibration" %in% names(results))
    expect_true("predictor_contribution" %in% names(results))
    expect_true("risk_ladder" %in% names(results))
    expect_true("overlap_diagnostics" %in% names(results))
    expect_true("parsimonious_sensitivity" %in% names(results))
    expect_true(any(results$predictor_contribution$section == "model_contribution"))
    expect_true(
        "predicted_60mo_melanoma_death_cumulative_incidence_risk" %in%
            names(results$no_gep_predictions)
    )
    expect_false("predicted_mss_5yr_risk" %in% names(results$no_gep_predictions))
    expect_true(all(c("Group", "Interpretation_Note") %in% names(results$unified_no_gep_overview)))
    expect_true(all(c(
        "Model", "Model_Method", "Reported_Risk_Scale", "Top_Predictor_1",
        "Use_Case", "CV_AUC_Stability_Lower", "CV_AUC_Stability_Upper"
    ) %in% names(results$unified_no_gep_model_comparison)))
    expect_true(all(c("No_GEP_Group", "Analysis", "Bin") %in% names(results$unified_no_gep_risk_strata)))
    expect_true(all(c("Group", "Observed_MFS_Method", "Observed_MFS_5yr_Event_Rate", "Median_Predicted_MFS_5yr_Risk", "Reported_Risk_Scale") %in% names(results$unified_no_gep_risk_ladder)))
    expect_true(
        "Median_Predicted_60mo_Melanoma_Death_Cumulative_Incidence_Risk" %in%
            names(results$unified_no_gep_risk_ladder)
    )
    expect_false("Median_Predicted_MSS_5yr_Risk" %in% names(results$unified_no_gep_risk_ladder))
    expect_equal(
        as.character(results$risk_ladder$group),
        c("Class 1", "GEP Not Tested", "GEP Failed/Indeterminate", "Class 2")
    )
})

test_that("report-native no-GEP figures reconcile to source tables", {
    actual_data <- readRDS(file.path(PROCESSED_DATA_DIR, "uveal_melanoma_full_cohort.rds"))
    test_output_dir <- file.path(TEST_OUTPUT_DIR, "exploratory_no_gep_figures")
    withr::defer(unlink(test_output_dir, recursive = TRUE), teardown_env())
    results <- run_exploratory_no_gep_report(
        output_dir = test_output_dir,
        verify_km_fix = FALSE,
        data = actual_data
    )

    expect_true(file.exists(results$output_paths$subgroup_outcomes))
    expect_true(file.exists(results$output_paths$direct_model_contributions))
    expect_setequal(results$no_gep_subgroups$no_gep_group, c("GEP Failed/Indeterminate", "GEP Not Tested"))
    followup_by_group <- results$follow_up_context %>%
        dplyr::filter(.data$summary_scope == "no_gep_group")
    expect_true(setequal(followup_by_group$no_gep_group, results$no_gep_subgroups$no_gep_group))
    subgroup_plot <- create_exploratory_no_gep_subgroup_outcomes_plot(
        results$no_gep_subgroups,
        results$follow_up_context,
        tempfile(fileext = ".png"),
        return_plot = TRUE
    )
    actual_subgroup_source <- subgroup_plot$plot_data %>%
        dplyr::mutate(
            no_gep_group = as.character(.data$no_gep_group),
            measure = as.character(.data$measure)
        ) %>%
        dplyr::select("no_gep_group", "n", "median_followup_years", "measure", "risk") %>%
        tidyr::pivot_wider(names_from = "measure", values_from = "risk") %>%
        dplyr::rename(
            observed_5yr_mfs_event_rate = "Observed 5-year MFS",
            predicted_5yr_mfs_risk = "Predicted 5-year MFS",
            observed_5yr_mss_event_rate = "Observed 5-year MSS",
            predicted_60mo_mss_risk = "Predicted 60-month MSS"
        ) %>%
        dplyr::select(
            "no_gep_group", "n", "median_followup_years",
            "observed_5yr_mfs_event_rate", "predicted_5yr_mfs_risk",
            "observed_5yr_mss_event_rate", "predicted_60mo_mss_risk"
        ) %>%
        dplyr::arrange(.data$no_gep_group)
    expected_subgroup_source <- results$no_gep_subgroups %>%
        dplyr::transmute(
            no_gep_group = as.character(.data$no_gep_group),
            n = .data$n,
            observed_5yr_mfs_event_rate = .data$observed_5yr_mfs_event_rate,
            predicted_5yr_mfs_risk = .data$median_predicted_5yr_mfs_risk,
            observed_5yr_mss_event_rate = .data$observed_5yr_mss_event_rate,
            predicted_60mo_mss_risk = .data$median_predicted_60mo_melanoma_death_cumulative_incidence_risk
        ) %>%
        dplyr::left_join(
            followup_by_group %>% dplyr::select("no_gep_group", "median_followup_years"),
            by = "no_gep_group"
        ) %>%
        dplyr::select(
            "no_gep_group", "n", "median_followup_years",
            "observed_5yr_mfs_event_rate", "predicted_5yr_mfs_risk",
            "observed_5yr_mss_event_rate", "predicted_60mo_mss_risk"
        ) %>%
        dplyr::arrange(.data$no_gep_group)
    expect_equal(actual_subgroup_source, expected_subgroup_source, tolerance = 1e-12)
    contribution_rows <- results$predictor_contribution %>%
        dplyr::filter(
            .data$section == "model_contribution",
            .data$model %in% c("Direct 5-Year MFS Risk", "Direct 60-Month Melanoma-Death Cumulative-Incidence Risk")
        )
    expect_true(setequal(
        unique(contribution_rows$model),
        c("Direct 5-Year MFS Risk", "Direct 60-Month Melanoma-Death Cumulative-Incidence Risk")
    ))
    contributor_plot <- create_exploratory_no_gep_direct_model_contributions_plot(
        results$predictor_contribution,
        tempfile(fileext = ".png"),
        return_plot = TRUE
    )
    actual_contributor_source <- contributor_plot$plot_data %>%
        dplyr::mutate(
            model = as.character(.data$model),
            predictor = as.character(.data$predictor)
        ) %>%
        dplyr::select("model", "predictor", "standardized_abs_coefficient", "direction", "rank") %>%
        dplyr::arrange(.data$model, .data$rank, .data$predictor)
    expected_contributor_source <- contribution_rows %>%
        dplyr::mutate(
            model = as.character(.data$model),
            predictor = as.character(.data$predictor)
        ) %>%
        dplyr::select("model", "predictor", "standardized_abs_coefficient", "direction", "rank") %>%
        dplyr::arrange(.data$model, .data$rank, .data$predictor)
    expect_equal(actual_contributor_source, expected_contributor_source, tolerance = 1e-12)
    expect_false(any(contributor_plot$plot_data$model == "Surrogate Class 2 Probability"))
    ranks_by_model <- split(contribution_rows$rank, contribution_rows$model)
    expect_true(all(vapply(ranks_by_model, function(x) {
        identical(sort(as.integer(x)), seq_len(length(x)))
    }, logical(1))))
})

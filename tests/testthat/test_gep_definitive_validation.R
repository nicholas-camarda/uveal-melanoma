synthetic_validation_cohort <- function(n_per_group = 10L, seed = 11L) {
    withr::with_seed(seed, {
        groups <- rep(GEP_VALIDATION_GROUP_LEVELS, each = n_per_group)
        imported_survival <- rep(c(0.95, 0.80, 0.58, 0.45), each = n_per_group)
        event_time <- stats::rexp(length(groups), rate = -log(imported_survival) / 60)
        censor_time <- stats::runif(length(groups), 20, 120)
        data.frame(
            gep_group = ifelse(grepl("Class 1", groups), "Class 1", "Class 2"),
            prame_status = ifelse(grepl("positive", groups), "Positive", "Negative"),
            biopsy1_gep_mfs = imported_survival,
            biopsy1_gep_mss = pmin(0.99, imported_survival + 0.03),
            tt_mets_months_analysis = pmin(event_time, censor_time),
            mfs_event_type = as.integer(event_time <= censor_time),
            tt_death_months = pmin(event_time, censor_time),
            mss_event_type = as.integer(event_time <= censor_time),
            mss_event_type_unknown_cause_melanoma = as.integer(event_time <= censor_time),
            mfs_analysis_eligible = TRUE,
            mss_analysis_eligible = TRUE,
            gep_test_timing = GEP_TEST_TIMING_LEVELS[[2]],
            gep_test_days = 0
        )
    })
}

test_that("validation groups cover exactly the four Class x PRAME cells", {
    grouped <- assign_gep_validation_group(data.frame(
        gep_group = c("Class 1", "Class 1", "Class 2", "Class 2", "Class 2", GEP_GROUP_NOT_TESTED),
        prame_status = c("Negative", "Positive", "Negative", "Positive", "Unknown", "Not Available")
    ))
    expect_identical(
        as.character(grouped$validation_group),
        c(GEP_VALIDATION_GROUP_LEVELS, NA_character_, NA_character_)
    )
    expect_identical(levels(grouped$validation_group), GEP_VALIDATION_GROUP_LEVELS)
})

test_that("Kaplan-Meier horizon risk is censoring-aware and never extrapolated", {
    # Ten patients: one metastasis at 12 months, nine followed past 60 months.
    estimate <- estimate_km_risk_at_horizon(
        time_months = c(12, rep(70, 9)), status = c(1L, rep(0L, 9)), horizon_months = 60
    )
    expect_true(estimate$estimable)
    expect_equal(estimate$risk, 0.1, tolerance = 1e-8)
    expect_identical(estimate$events, 1L)
    expect_identical(estimate$event_free_through, 9L)
    expect_false(is.na(estimate$ci_lower))

    # A patient censored at 12 months is not a known 5-year survivor: with one
    # censored at 12 and one metastasis at 60, KM survival at 60 is 0 (not 0.5).
    censored_early <- estimate_km_risk_at_horizon(c(12, 60), c(0L, 1L), horizon_months = 60)
    expect_equal(censored_early$risk, 1, tolerance = 1e-8)
    expect_identical(censored_early$censored_before, 1L)

    # No follow-up to the horizon: no value is carried forward.
    short <- estimate_km_risk_at_horizon(c(10, 30, 40), c(1L, 0L, 0L), horizon_months = 60)
    expect_false(short$estimable)
    expect_true(is.na(short$risk))

    # No events by the horizon: point estimate only, no degenerate interval.
    no_events <- estimate_km_risk_at_horizon(rep(70, 5), rep(0L, 5), horizon_months = 60)
    expect_equal(no_events$risk, 0)
    expect_true(is.na(no_events$ci_lower) && is.na(no_events$ci_upper))
})

test_that("cumulative incidence keeps competing deaths in the risk set", {
    # No censoring before the horizon: CIF equals the share of melanoma deaths.
    estimate <- estimate_cif_risk_at_horizon(
        time_months = c(10, 20, 30, 40, 50, 55, 58, 59, 59.5, 70),
        status = c(1L, 1L, 2L, 2L, 0L, 0L, 0L, 0L, 0L, 0L),
        horizon_months = 60
    )
    expect_equal(estimate$risk, 0.2, tolerance = 1e-8)
    expect_identical(estimate$events, 2L)
    expect_identical(estimate$competing, 2L)

    competing_only <- estimate_cif_risk_at_horizon(c(10, 70, 80), c(2L, 0L, 0L), 60)
    expect_equal(competing_only$risk, 0)
    expect_true(is.na(competing_only$ci_lower))

    short <- estimate_cif_risk_at_horizon(c(10, 20), c(1L, 0L), 60)
    expect_false(short$estimable)
})

test_that("Harrell's C counts ties as one half and reports the tied-pair share", {
    cohort <- synthetic_validation_cohort()
    risk <- 1 - cohort$biopsy1_gep_mfs
    time <- cohort$tt_mets_months_analysis
    status <- cohort$mfs_event_type
    horizon <- 60

    truncated_time <- pmin(time, horizon)
    truncated_event <- as.integer(status == 1L & time <= horizon)
    concordant <- discordant <- tied <- 0
    for (i in which(truncated_event == 1L)) {
        later <- which(truncated_time > truncated_time[i])
        concordant <- concordant + sum(risk[i] > risk[later])
        discordant <- discordant + sum(risk[i] < risk[later])
        tied <- tied + sum(risk[i] == risk[later])
    }
    manual_c <- (concordant + 0.5 * tied) / (concordant + discordant + tied)

    result <- calculate_gep_harrell_c(time, status, risk, horizon, n_boot = 50L, seed = 1L)
    expect_true(result$estimable)
    expect_equal(result$c_index, manual_c, tolerance = 1e-8)
    expect_equal(result$tied_pair_share, tied / (concordant + discordant + tied), tolerance = 1e-8)
    expect_gt(result$tied_pair_share, 0)
    expect_identical(result$distinct_risk_values, 4L)
    # Dropping tied pairs (the previous behaviour) would overstate C.
    expect_lt(result$c_index, concordant / (concordant + discordant))
    expect_true(result$ci_lower <= result$c_index && result$c_index <= result$ci_upper)
})

test_that("discrimination is reported as not estimable instead of forced", {
    few_events <- calculate_gep_harrell_c(
        time_months = c(rep(70, 29), 10), status = c(rep(0L, 29), 1L),
        risk = rep(c(0.1, 0.4), 15), horizon_months = 60, n_boot = 10L, seed = 1L
    )
    expect_false(few_events$estimable)
    expect_match(few_events$reason, "Fewer than 3 events")

    no_variation <- calculate_gep_harrell_c(
        time_months = rep(c(10, 70), 15), status = rep(c(1L, 0L), 15),
        risk = rep(0.2, 30), horizon_months = 60, n_boot = 10L, seed = 1L
    )
    expect_false(no_variation$estimable)

    # Events after the horizon do not count toward the horizon C.
    late_events <- calculate_gep_harrell_c(
        time_months = rep(c(80, 90), 15), status = rep(1L, 30),
        risk = rep(c(0.1, 0.4), 15), horizon_months = 60, n_boot = 10L, seed = 1L
    )
    expect_false(late_events$estimable)
})

test_that("bootstrap intervals are reproducible and leave the global RNG untouched", {
    cohort <- synthetic_validation_cohort()
    risk <- 1 - cohort$biopsy1_gep_mfs
    set.seed(99)
    expected_next <- stats::runif(1)
    set.seed(99)
    first <- calculate_gep_harrell_c(cohort$tt_mets_months_analysis, cohort$mfs_event_type, risk, 60, 40L, 7L)
    actual_next <- stats::runif(1)
    second <- calculate_gep_harrell_c(cohort$tt_mets_months_analysis, cohort$mfs_event_type, risk, 60, 40L, 7L)
    expect_identical(first$ci_lower, second$ci_lower)
    expect_identical(first$ci_upper, second$ci_upper)
    expect_equal(actual_next, expected_next)
})

test_that("group table compares observed and imported risk without extrapolating", {
    cohort <- synthetic_validation_cohort()
    prepared <- prepare_gep_validation_endpoint(cohort, "MFS")
    table <- build_gep_group_table(prepared, "MFS", horizon_months = 60)
    expect_identical(table$group, c(GEP_VALIDATION_GROUP_LEVELS, "Overall (case-mix standardised)"))
    expect_equal(table$n[1:4], rep(10L, 4))
    expect_equal(table$imported_risk[1:4], 1 - c(0.95, 0.80, 0.58, 0.45), tolerance = 1e-8)
    expect_equal(
        table$difference_pct_points[1:4],
        100 * (table$observed_risk[1:4] - table$imported_risk[1:4]),
        tolerance = 1e-8
    )
    # The overall row is the group-size-weighted mean of the group rows.
    expect_equal(table$observed_risk[5], sum(table$n[1:4] * table$observed_risk[1:4]) / 40, tolerance = 1e-8)

    # A group that is not followed to the horizon yields NA and no overall value.
    truncated <- cohort
    in_group <- truncated$biopsy1_gep_mfs == 0.45
    truncated$tt_mets_months_analysis[in_group] <- pmin(truncated$tt_mets_months_analysis[in_group], 30)
    short_table <- build_gep_group_table(prepare_gep_validation_endpoint(truncated, "MFS"), "MFS", 60)
    expect_false(short_table$estimable[4])
    expect_true(is.na(short_table$observed_risk[4]))
    expect_match(short_table$note[4], "Not followed to 60 months")
    expect_true(is.na(short_table$observed_risk[5]))
})

test_that("competing-risk AUC uses the horizon element whether or not competing deaths occur", {
    cohort <- synthetic_validation_cohort(n_per_group = 25L, seed = 5L)
    status <- cohort$mss_event_type
    time <- cohort$tt_death_months
    risk <- 1 - cohort$biopsy1_gep_mss

    without_competing <- calculate_gep_competing_risk_auc(time, status, risk, 60, n_boot = 20L, seed = 2L)
    expect_true(without_competing$estimable)
    expect_gt(without_competing$auc, 0.5)

    competing_status <- status
    competing_status[which(status == 0L)[1:3]] <- 2L
    with_competing <- calculate_gep_competing_risk_auc(time, competing_status, risk, 60, n_boot = 20L, seed = 2L)
    expect_true(with_competing$estimable)
    expect_true(is.finite(with_competing$ci_lower))
    expect_equal(with_competing$valid_draws + with_competing$failed_draws, 20L)

    sparse <- calculate_gep_competing_risk_auc(
        time_months = c(rep(80, 28), 10, 20), status = c(rep(0L, 28), 1L, 1L),
        risk = rep(c(0.1, 0.4), 15), horizon_months = 60, n_boot = 10L, seed = 1L
    )
    expect_false(sparse$estimable)
    expect_match(sparse$reason, "Fewer than 3 melanoma-specific deaths")
})

test_that("proportional-hazards check reports one term with an explicit status", {
    prepared <- prepare_gep_validation_endpoint(synthetic_validation_cohort(n_per_group = 25L), "MFS")
    result <- check_gep_proportional_hazards(prepared)
    expect_identical(nrow(result), 1L)
    expect_true(result$status %in% c("no concern", "concern"))

    skipped <- check_gep_proportional_hazards(prepared[prepared$status == 0L, ])
    expect_identical(skipped$status, "not assessed")
    expect_match(skipped$reason, "Too few events")
})

test_that("the definitive validation module reports no 7- or 10-year horizons", {
    result <- run_gep_definitive_validation(synthetic_validation_cohort(n_per_group = 15L), n_boot = 20L)
    expect_identical(result$horizon_months, 60)
    module_text <- paste(readLines(here::here("scripts", "gep", "validation", "gep_definitive_validation.R")), collapse = "\n")
    expect_false(grepl("7yr|10yr|expected_mfs_7|expected_mss_7|\\^\\(7", module_text))
    expect_false(grepl("concordance.index", module_text, fixed = TRUE))
})

test_that("validation outputs reconcile flow counts with Record_Issues and keep identifiers out of the primary workbook", {
    cohort <- synthetic_validation_cohort(n_per_group = 12L)
    extra <- data.frame(
        gep_group = c(GEP_GROUP_TEST_FAILED, GEP_GROUP_TEST_FAILED, rep(GEP_GROUP_NOT_TESTED, 3)),
        prame_status = "Not Available",
        biopsy1_gep_mfs = NA_real_, biopsy1_gep_mss = NA_real_,
        tt_mets_months_analysis = 30, mfs_event_type = 0L, tt_death_months = 30, mss_event_type = 0L,
        mss_event_type_unknown_cause_melanoma = 0L,
        mfs_analysis_eligible = FALSE, mss_analysis_eligible = FALSE,
        gep_test_timing = GEP_TEST_TIMING_LEVELS[[5]], gep_test_days = NA_real_
    )
    cohort <- rbind(cohort, extra)
    cohort$id <- sprintf("P%03d", seq_len(nrow(cohort)))
    cohort$biopsy1_gep_raw <- c(
        ifelse(cohort$gep_group[1:48] == "Class 1", "Class_1A_PRAME_negative", "Class_2_PRAME_negative"),
        "Failed", "Class_2_PRAME_not_reported", NA, NA, NA
    )
    cohort$gep_group <- classify_gep_group(cohort$biopsy1_gep_raw)
    # One primary patient was sampled at recurrence (late) and one has no biopsy date.
    cohort$gep_test_days[3] <- 400
    cohort$gep_test_timing[3] <- GEP_TEST_TIMING_LATE
    cohort$gep_test_timing[4] <- GEP_TEST_TIMING_LEVELS[[5]]
    cohort$gep_test_days[4] <- NA_real_
    # One primary patient has a biopsy about six years before treatment (date review).
    cohort$gep_test_timing[5] <- GEP_TEST_TIMING_LEVELS[[1]]
    cohort$gep_test_days[5] <- -2197
    # One primary patient is metastatic at baseline: ineligible for MFS, still eligible for MSS.
    cohort$mfs_analysis_eligible[6] <- FALSE
    cohort$mets_free_at_baseline <- TRUE
    cohort$mets_free_at_baseline[6] <- FALSE

    results <- run_gep_definitive_validation(cohort, n_boot = 20L)
    out_dir <- withr::local_tempdir()
    paths <- write_gep_validation_outputs(
        results, cohort, "uveal_melanoma_full_cohort",
        output_dirs = list(obj4_mfs = file.path(out_dir, "04_GEP_Validation", "a_metastasis_free_survival")),
        prefix = "test_"
    )

    expect_true(all(file.exists(c(paths$primary, paths$diagnostics, paths$figures))))
    expect_identical(
        readxl::excel_sheets(paths$primary),
        c("Group_Table_MFS", "Group_Table_MSS", "Discrimination", "Sensitivity_Late_Biopsy", "Sensitivity_Cause_of_Death")
    )
    expect_identical(
        readxl::excel_sheets(paths$diagnostics),
        c("Specification_Coding", "Denominators_Flow", "Checks", "Warnings_Limitations", "Record_Issues")
    )
    # Identifiers live only in the diagnostics workbook.
    for (sheet in readxl::excel_sheets(paths$primary)) {
        expect_false("id" %in% names(readxl::read_excel(paths$primary, sheet = sheet)))
    }

    flow <- readxl::read_excel(paths$diagnostics, sheet = "Denominators_Flow")
    issues <- readxl::read_excel(paths$diagnostics, sheet = "Record_Issues")
    total <- flow$n[flow$step == "All patients in cohort"]
    primary <- flow$n[grepl("^Definitive Class 1/2", flow$step)]
    expect_identical(total, as.numeric(nrow(cohort)))
    # Every patient outside the primary set, plus primary patients flagged for biopsy timing.
    flagged <- flow$n[grepl("^Biopsy timing in the primary set: (More than 180|No biopsy)", flow$step)]
    expect_identical(sum(flagged), 2)
    review <- flow$n[grepl("date review", flow$step)]
    expect_identical(review, 1)
    not_eligible <- flow$n[grepl("^Primary-set patients not eligible", flow$step)]
    expect_identical(not_eligible, c(1, 0))
    expect_identical(nrow(issues), as.integer(total - primary + sum(flagged) + review + sum(not_eligible)))
    ineligible_row <- issues[startsWith(issues$step, "Primary-set patient not eligible"), ]
    expect_identical(ineligible_row$record_key, "P006")
    expect_match(ineligible_row$step, "MFS")
    expect_match(ineligible_row$current_value, "mets_free_at_baseline = FALSE", fixed = TRUE)
    outside <- issues[issues$step == "GEP availability group", ]
    expect_setequal(unique(outside$reason), c(GEP_GROUP_TEST_FAILED, GEP_GROUP_CLASS_PRAME_MISSING, GEP_GROUP_NOT_TESTED))
    expect_true(all(outside$field == "biopsy1_gep_raw"))
    timing <- issues[startsWith(issues$step, "Biopsy timing"), ]
    expect_setequal(timing$reason, c(GEP_TEST_TIMING_LATE, GEP_TEST_TIMING_LEVELS[[5]], GEP_TEST_REVIEW_BEFORE_REASON))
    expect_identical(timing$current_value[timing$reason == GEP_TEST_REVIEW_BEFORE_REASON], "-2197 days from treatment")
    expect_equal(sum(flow$n[startsWith(flow$step, "Biopsy timing in the primary set")]), primary)

    sensitivity <- readxl::read_excel(paths$primary, sheet = "Sensitivity_Late_Biopsy")
    expect_true(all(c("endpoint", "quantity", "primary", "excluding_late_biopsy", "difference", "unit") %in% names(sensitivity)))
    expect_identical(unname(results$late_biopsy_sensitivity$patients_excluded[["MFS"]]), 1L)
    expect_true(all(issues$key_source == "id"))
    expect_true(all(startsWith(issues$record_key, "P")))

    # Without an identifier column the key is the row position, and says so.
    anonymous <- build_gep_record_issues(cohort[, setdiff(names(cohort), "id")])
    expect_identical(nrow(anonymous), nrow(issues))
    expect_true(all(anonymous$key_source == "row position (no identifier column)"))

    # A cohort in which every patient is in the primary set yields an empty, well-formed sheet.
    complete <- cohort[cohort$gep_group %in% c("Class 1", "Class 2") & cohort$gep_test_timing == GEP_TEST_TIMING_LEVELS[[2]] &
        cohort$mfs_analysis_eligible, ]
    expect_identical(nrow(build_gep_record_issues(complete)), 0L)

    checks <- readxl::read_excel(paths$diagnostics, sheet = "Checks")
    expect_true(all(checks$status %in% c("no concern", "concern", "not assessed")))
    expect_false(anyNA(checks$status))
})

test_that("validation reporting fails clearly when the GEP group column is missing", {
    cohort <- synthetic_validation_cohort(n_per_group = 10L)
    results <- run_gep_definitive_validation(cohort, n_boot = 10L)
    cohort$gep_group <- NULL
    expect_error(
        write_gep_validation_outputs(
            results, cohort, "uveal_melanoma_full_cohort",
            output_dirs = list(obj4_mfs = file.path(withr::local_tempdir(), "a")), prefix = "t_"
        ),
        "requires columns: gep_group, biopsy1_gep_raw"
    )
})

test_that("curve figures draw groups without melanoma deaths as flat zero curves", {
    cohort <- synthetic_validation_cohort(n_per_group = 12L)
    cohort$mss_event_type[cohort$biopsy1_gep_mss > 0.95] <- 0L
    prepared <- prepare_gep_validation_endpoint(cohort, "MSS")
    table <- build_gep_group_table(prepared, "MSS", horizon_months = 60)
    plot <- plot_gep_validation_curves(prepared, table, "MSS", 60, "Synthetic cohort")
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
})

test_that("unknown-cause deaths are recounted as melanoma deaths in the sensitivity only", {
    endpoint_data <- data.frame(
        validation_group = factor(rep(GEP_VALIDATION_GROUP_LEVELS[1:2], each = 10), levels = GEP_VALIDATION_GROUP_LEVELS),
        time_months = rep(c(70, 20), 10),
        status = c(rep(0L, 10), rep(c(0L, 2L), 5)),
        stringsAsFactors = FALSE
    )
    # Three of the five other-cause deaths in the second group have an unknown cause.
    endpoint_data$status_unknown_melanoma <- endpoint_data$status
    endpoint_data$status_unknown_melanoma[c(12, 14, 16)] <- 1L

    comparison <- compare_unknown_cause_of_death_assumption(endpoint_data, "validation_group", 60)
    first <- comparison[comparison$group == GEP_VALIDATION_GROUP_LEVELS[1], ]
    second <- comparison[comparison$group == GEP_VALIDATION_GROUP_LEVELS[2], ]
    expect_identical(first$deaths_reclassified, 0L)
    expect_equal(first$difference_pct_points, 0)
    expect_identical(second$deaths_reclassified, 3L)
    expect_gt(second$unknown_cause_as_melanoma_risk, second$primary_risk)
    # The primary MSS status is untouched.
    expect_identical(sum(endpoint_data$status == 1L), 0L)
})

test_that("the MSS AUC treats competing deaths as controls (AUC_2), not as excluded patients (AUC_1)", {
    n <- 60L
    risk <- seq(0.05, 0.95, length.out = n)
    status <- rep(0L, n)
    time <- rep(100, n)
    # Melanoma deaths among high-risk patients; competing deaths among low-risk patients.
    status[c(55, 52, 49, 45, 40, 36)] <- 1L
    time[c(55, 52, 49, 45, 40, 36)] <- c(12, 20, 25, 30, 41, 50)
    status[c(3, 8, 14, 20, 26, 31)] <- 2L
    time[c(3, 8, 14, 20, 26, 31)] <- c(10, 15, 22, 28, 35, 44)
    fit <- timeROC::timeROC(T = time, delta = status, marker = risk, cause = 1, weighting = "marginal", times = 60, iid = FALSE)
    auc_1 <- fit$AUC_1[[length(fit$AUC_1)]]
    auc_2 <- fit$AUC_2[[length(fit$AUC_2)]]
    expect_gt(abs(auc_1 - auc_2), 1e-4)

    result <- calculate_gep_competing_risk_auc(time, status, risk, 60, n_boot = 10L, seed = 4L)
    expect_equal(result$auc, auc_2, tolerance = 1e-10)
    expect_gt(abs(result$auc - auc_1), 1e-4)
})

test_that("the retired survcomp package (which dropped tied pairs) is neither used nor locked", {
    lock <- jsonlite::fromJSON(here::here("renv.lock"), simplifyVector = FALSE)
    expect_false(any(c("survcomp", "survivalROC", "SuppDists", "rmeta", "bootstrap") %in% names(lock$Packages)))
    script_files <- list.files(here::here("scripts"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
    script_text <- vapply(script_files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1))
    expect_identical(basename(script_files)[grepl("survcomp", script_text, fixed = TRUE)], character(0))
})

test_that("an eligible patient without a validation group stops the run instead of vanishing", {
    cohort <- synthetic_validation_cohort()
    cohort$prame_status[1] <- "Unknown"
    expect_error(
        prepare_gep_validation_endpoint(cohort, "MFS"),
        "1 patient\\(s\\) eligible for the MFS endpoint have no validation group"
    )
})

test_that("the late-biopsy sensitivity reports point estimates without running a bootstrap", {
    cohort <- synthetic_validation_cohort(n_per_group = 15L)
    cohort$gep_test_timing[2] <- GEP_TEST_TIMING_LATE
    sensitivity <- run_gep_late_biopsy_sensitivity(
        prepare_gep_validation_endpoint(cohort, "MFS"), prepare_gep_validation_endpoint(cohort, "MSS"), 60
    )
    expect_identical(unname(sensitivity$patients_excluded), c(1L, 1L))
    expect_true(all(is.finite(sensitivity$table$primary[grepl("Harrell|AUC", sensitivity$table$quantity)])))
    point <- calculate_gep_harrell_c(cohort$tt_mets_months_analysis, cohort$mfs_event_type, 1 - cohort$biopsy1_gep_mfs, 60, n_boot = 0L, seed = 1L)
    expect_true(point$estimable)
    expect_true(is.na(point$ci_lower))
    expect_identical(point$valid_draws + point$failed_draws, 0L)
})

test_that("every bootstrap uses one interval rule with a valid-draw requirement", {
    draws <- c(seq(0.1, 0.9, length.out = 95), rep(NA_real_, 5))
    interval <- bootstrap_percentile_interval(draws, 100L)
    expect_identical(c(interval$valid_draws, interval$failed_draws), c(95L, 5L))
    expect_equal(interval$ci_lower, unname(stats::quantile(draws, 0.025, na.rm = TRUE)))

    too_many_failures <- bootstrap_percentile_interval(c(draws[1:80], rep(NA_real_, 20)), 100L)
    expect_true(is.na(too_many_failures$ci_lower))
    expect_identical(too_many_failures$failed_draws, 20L)
    expect_true(is.na(bootstrap_percentile_interval(numeric(0), 0L)$ci_upper))

    for (file in c("gep_definitive_validation.R", "gep_no_gep_analysis.R")) {
        module_text <- paste(readLines(here::here("scripts", "gep", "validation", file)), collapse = "\n")
        expect_false(grepl("0.9 * n_boot", module_text, fixed = TRUE), info = file)
    }
})

test_that("the cross-cohort figure refuses workbooks from an earlier run and missing workbooks", {
    output_root <- withr::local_tempdir()
    sources <- data.frame(
        cohort_dir = c("a_cohort", "b_cohort"), prefix = c("a_", "b_"),
        cohort_label = c("Cohort A", "Cohort B"), stringsAsFactors = FALSE
    )
    table <- build_gep_group_table(
        prepare_gep_validation_endpoint(synthetic_validation_cohort(), "MFS"), "MFS", 60
    )
    for (i in seq_len(nrow(sources))) {
        dir <- file.path(output_root, sources$cohort_dir[i], "04_GEP_Validation")
        dir.create(dir, recursive = TRUE)
        write_gep_workbook(list(Group_Table_MFS = table), file.path(dir, paste0(sources$prefix[i], "gep_validation_5yr.xlsx")))
    }
    rlang::local_bindings(OUTPUT_DIR = output_root, .env = .GlobalEnv)
    figure_dir <- file.path(output_root, "figures")

    path <- create_gep_observed_vs_imported_panel_report(Sys.time() - 3600, figure_dir, sources)
    expect_true(file.exists(path))
    expect_error(
        create_gep_observed_vs_imported_panel_report(Sys.time() + 3600, figure_dir, sources),
        "was not written in this run"
    )
    unlink(file.path(output_root, "b_cohort", "04_GEP_Validation"), recursive = TRUE)
    expect_error(
        create_gep_observed_vs_imported_panel_report(Sys.time() - 3600, figure_dir, sources),
        "Validation workbook not found"
    )
})

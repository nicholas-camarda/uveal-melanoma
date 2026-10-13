synthetic_validation_cohort <- function(n_per_group = 10L, seed = 11L) {
    withr::with_seed(seed, {
        groups <- rep(GEP_VALIDATION_GROUP_LEVELS, each = n_per_group)
        imported_survival <- rep(c(0.95, 0.80, 0.58, 0.45), each = n_per_group)
        event_time <- stats::rexp(length(groups), rate = -log(imported_survival) / 60)
        censor_time <- stats::runif(length(groups), 20, 120)
        data.frame(
            gep_class_simple = ifelse(grepl("Class 1", groups), "Class 1", "Class 2"),
            prame_status = ifelse(grepl("positive", groups), "Positive", "Negative"),
            biopsy1_gep_mfs = imported_survival,
            biopsy1_gep_mss = pmin(0.99, imported_survival + 0.03),
            tt_mets_months_analysis = pmin(event_time, censor_time),
            mfs_event_type = as.integer(event_time <= censor_time),
            tt_death_months = pmin(event_time, censor_time),
            mss_event_type = as.integer(event_time <= censor_time),
            mfs_analysis_eligible = TRUE,
            mss_analysis_eligible = TRUE
        )
    })
}

test_that("validation groups cover exactly the four Class x PRAME cells", {
    grouped <- assign_gep_validation_group(data.frame(
        gep_class_simple = c("Class 1", "Class 1", "Class 2", "Class 2", "Class 2", "GEP Not Tested"),
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

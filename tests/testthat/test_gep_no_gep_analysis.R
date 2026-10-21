synthetic_no_gep_cohort <- function(n = 200L, seed = 21L) {
    withr::with_seed(seed, {
        diameter <- round(stats::rnorm(n, 12, 3), 1)
        t_stage <- as.character(cut(diameter, c(-Inf, 9, 12, 15, Inf), labels = paste0("T", 1:4)))
        event_rate <- -log(1 - stats::plogis(-2 + 0.25 * (diameter - 12))) / 60
        event_time <- stats::rexp(n, rate = event_rate)
        censor_time <- stats::runif(n, 10, 150)
        time <- pmin(event_time, censor_time)
        status <- as.integer(event_time <= censor_time)
        data.frame(
            id = seq_len(n),
            gep_group = sample(c("Class 1", "Class 2", "Not tested", "Failed or unknown result"), n, TRUE, prob = c(0.2, 0.1, 0.65, 0.05)),
            age_at_diagnosis = stats::rnorm(n, 64, 10),
            sex = sample(c("Female", "Male"), n, TRUE),
            initial_tumor_height = round(stats::rnorm(n, 5, 1.5), 1),
            initial_tumor_diameter = diameter,
            location = c("Iris", "Iris", sample(c("Choroidal", "Ciliary Body"), n - 2L, TRUE, prob = c(0.8, 0.2))),
            initial_t_stage_simple = t_stage,
            internal_reflectivity = sample(c("Low", "Medium", "Unknown"), n, TRUE, prob = c(0.2, 0.2, 0.6)),
            srf = sample(c("No", "Yes"), n, TRUE),
            initial_vision = stats::rnorm(n, 0.4, 0.3),
            optic_nerve_involvement = sample(0:1, n, TRUE),
            treatment_group = sample(c("PBT", "GKSRS"), n, TRUE),
            treatment_year = sample(2008:2024, n, TRUE),
            mets_free_at_baseline = c(FALSE, rep(TRUE, n - 1L)),
            tt_mets_months_analysis = time, mfs_event_type = status,
            tt_death_months = time, mss_event_type = status,
            mss_event_type_unknown_cause_melanoma = status,
            stringsAsFactors = FALSE
        )
    })
}

test_that("no-GEP preparation recodes Unknown reflectivity to missing and requires unique identifiers", {
    cohort <- synthetic_no_gep_cohort()
    prepared <- prepare_no_gep_data(cohort)
    expect_false("Unknown" %in% levels(prepared$internal_reflectivity))
    expect_equal(sum(is.na(prepared$internal_reflectivity)), sum(cohort$internal_reflectivity == "Unknown"))
    expect_identical(levels(prepared$gep_group), GEP_GROUP_LEVELS)

    duplicated_ids <- cohort
    duplicated_ids$id[2] <- duplicated_ids$id[1]
    expect_error(prepare_no_gep_data(duplicated_ids), "unique, observed patient identifier")
    expect_error(prepare_no_gep_data(cohort[, setdiff(names(cohort), "srf")]), "requires columns: srf")

    unclassified <- cohort
    unclassified$gep_group[1] <- NA
    expect_error(prepare_no_gep_data(unclassified), "gep_group is missing for 1 patient")
})

test_that("predictor screening drops incomplete predictors and keeps sparse-level patients out of the model data", {
    prepared <- prepare_no_gep_data(synthetic_no_gep_cohort())
    sets <- list(
        MFS = prepare_no_gep_endpoint(prepared, "MFS"),
        MSS = prepare_no_gep_endpoint(prepared, "MSS"),
        Surrogate = prepared[prepared$gep_group %in% c("Class 1", "Class 2"), ]
    )
    screening <- screen_no_gep_predictors(sets)
    expect_identical(screening$status[screening$predictor == "internal_reflectivity"], "dropped")
    expect_match(screening$reason[screening$predictor == "internal_reflectivity"], "completeness below 90%")
    expect_identical(screening$status[screening$predictor == "initial_tumor_diameter"], "retained")

    # The patient with metastasis at baseline is not eligible for MFS but is for MSS.
    expect_identical(nrow(sets$MFS), nrow(sets$MSS) - 1L)

    predictors <- screening$predictor[screening$status == "retained"]
    built <- build_no_gep_model_data(sets$MSS, predictors, "MSS")
    # Two patients have a location level (Iris) with fewer than the minimum patients.
    expect_identical(built$counts$removed_sparse_level, 2L)
    expect_identical(
        nrow(sets$MSS) - built$counts$removed_sparse_level - built$counts$removed_missing_predictor,
        built$counts$model_n
    )
    expect_identical(nrow(built$data), built$counts$model_n)
    expect_identical(nrow(built$removed), built$counts$removed_sparse_level + built$counts$removed_missing_predictor)
    expect_true(all(built$removed$reason %in% c("Factor level with fewer than the minimum patients", "Missing baseline predictor")))
    expect_false(any(built$data$id %in% as.numeric(built$removed$id)))
})

test_that("out-of-fold performance compares models with no-fit comparators on identical patients", {
    n <- 80L
    outcome <- rep(c(0L, 1L), length.out = n)
    make_oof <- function(prediction) {
        do.call(rbind, lapply(1:2, function(r) {
            data.frame(
                stable_id = as.character(seq_len(n)), repeat_id = r, prediction = prediction, horizon_event = outcome,
                ipcw_weight = 1, gep_group = rep(c("Not tested", "Not tested", "Class 1"), length.out = n),
                initial_t_stage_simple = ifelse(outcome == 1L, "T4", "T1"),
                initial_tumor_diameter = 10 + 5 * outcome, stringsAsFactors = FALSE
            )
        }))
    }
    informative <- make_oof(ifelse(outcome == 1L, 0.8, 0.2))
    uninformative <- make_oof(rep(0.5, n))
    performance <- summarize_no_gep_model_performance(
        list(`Penalized, all retained predictors` = informative, `Penalized, parsimonious predictors` = uninformative),
        "MFS", n_boot = 30L, seed = 3L
    )

    no_gep <- performance[performance$performance_scope == "No GEP", ]
    expect_setequal(no_gep$model, c(
        "Penalized, all retained predictors", "Penalized, parsimonious predictors", "T stage alone", "Tumor diameter alone"
    ))
    expect_equal(no_gep$auc[no_gep$model == "T stage alone"], 1)
    expect_equal(no_gep$auc[no_gep$model == "Penalized, parsimonious predictors"], 0.5)
    expect_equal(no_gep$auc[no_gep$model == "Penalized, all retained predictors"], 1)
    # Unit weights: the no-information Brier score is p(1 - p) = 0.25 for a 50% prevalence.
    informative_row <- no_gep[no_gep$model == "Penalized, all retained predictors", ]
    expect_equal(informative_row$no_information_brier, 0.25, tolerance = 1e-8)
    expect_equal(informative_row$brier, 0.04, tolerance = 1e-8)
    expect_equal(informative_row$brier_skill, 1 - 0.04 / 0.25, tolerance = 1e-8)
    expect_true(is.na(no_gep$brier[no_gep$model == "T stage alone"]))
    # Same patients in every row of a scope.
    expect_length(unique(no_gep$patients), 1L)
    expect_identical(t_stage_to_score(c("T1", "T4")), c(1, 4))
})

test_that("bootstrap intervals are reproducible and leave the global RNG untouched", {
    n <- 60L
    oof <- data.frame(
        stable_id = as.character(seq_len(n)), repeat_id = 1L,
        prediction = seq(0.05, 0.95, length.out = n), horizon_event = rep(c(0L, 1L), length.out = n), ipcw_weight = 1
    )
    set.seed(5)
    expected_next <- stats::runif(1)
    set.seed(5)
    first <- summarize_no_gep_oof_metric(oof, "prediction", "auc", 40L, 9L)
    actual_next <- stats::runif(1)
    second <- summarize_no_gep_oof_metric(oof, "prediction", "auc", 40L, 9L)
    expect_identical(first$ci_lower, second$ci_lower)
    expect_identical(first$ci_upper, second$ci_upper)
    expect_equal(actual_next, expected_next)
    expect_true(first$ci_lower <= first$estimate && first$estimate <= first$ci_upper)
})

test_that("group outcomes never extrapolate and report reverse Kaplan-Meier follow-up", {
    endpoint_data <- data.frame(
        gep_group = factor(
            rep(c("Not tested", "Failed or unknown result"), c(20, 6)), levels = GEP_GROUP_LEVELS
        ),
        time_months = c(rep(c(80, 30), 10), c(12, 20, 24, 30, 36, 40)),
        status = c(rep(c(0L, 1L), 10), rep(0L, 6))
    )
    outcomes <- build_no_gep_group_outcomes(endpoint_data, "MFS", 60)
    expect_identical(outcomes$group, GEP_GROUP_LEVELS)
    tested <- outcomes[outcomes$group == "Not tested", ]
    failed <- outcomes[outcomes$group == "Failed or unknown result", ]
    expect_true(tested$estimable)
    expect_equal(tested$observed_risk, 0.5, tolerance = 1e-8)
    expect_false(failed$estimable)
    expect_true(is.na(failed$observed_risk))
    expect_match(failed$note, "Not followed to 60 months \\(maximum 40 months\\)")
    expect_identical(outcomes$n[outcomes$group == "Class 1"], 0L)
    expect_true(is.finite(tested$median_followup_months))
})

test_that("baseline table is descriptive and covers every group", {
    prepared <- prepare_no_gep_data(synthetic_no_gep_cohort())
    table <- build_no_gep_baseline_table(prepared)
    expect_setequal(setdiff(names(table), c("characteristic", "statistic")), c("All patients", GEP_GROUP_LEVELS))
    expect_identical(table[["All patients"]][table$characteristic == "Patients"], as.character(nrow(prepared)))
    expect_true(all(c(sprintf("Treated in %d or later", GEP_TESTING_ERA_START_YEAR), "Age at diagnosis, years") %in% table$characteristic))
})

test_that("no-GEP outputs reconcile their flow with Record_Issues and keep identifiers out of the primary workbook", {
    cohort <- synthetic_no_gep_cohort()
    results <- run_no_gep_analysis(cohort, n_boot = 20L)
    out_dir <- withr::local_tempdir()
    paths <- write_no_gep_outputs(results, prepare_no_gep_data(cohort), out_dir, "test_")

    expect_true(all(file.exists(c(paths$primary, paths$diagnostics, paths$figures))))
    expect_identical(
        readxl::excel_sheets(paths$primary),
        c("Group_Outcomes_MFS", "Group_Outcomes_MSS", "Baseline_By_Group", "Model_Performance", "AUC_Difference", "Sensitivity_Treatment_Era", "Sensitivity_Cause_of_Death")
    )
    expect_identical(
        readxl::excel_sheets(paths$diagnostics),
        c("Specification_Coding", "Denominators_Flow", "Predictor_Screening", "Checks", "Warnings_Limitations", "Record_Issues")
    )
    for (sheet in readxl::excel_sheets(paths$primary)) {
        expect_false("record_key" %in% names(readxl::read_excel(paths$primary, sheet = sheet)))
    }

    flow <- results$model_flow
    expect_identical(
        flow$candidates - flow$not_eligible_for_endpoint - flow$removed_sparse_level - flow$removed_missing_predictor,
        flow$model_n
    )
    issues <- readxl::read_excel(paths$diagnostics, sheet = "Record_Issues")
    expect_identical(
        nrow(issues),
        as.integer(sum(flow$not_eligible_for_endpoint + flow$removed_sparse_level + flow$removed_missing_predictor) + results$treatment_era_excluded_n)
    )
    # Every patient treated before the testing era is listed by id, and the flow counts them.
    # Patients not eligible for an endpoint show the fields behind the exclusion, not a placeholder.
    endpoint_issues <- issues[grepl("^Not eligible for the endpoint", issues$reason), ]
    expect_gt(nrow(endpoint_issues), 0L)
    expect_false(any(endpoint_issues$current_value == "see cohort record"))
    expect_true(any(grepl("mets_free_at_baseline = FALSE", endpoint_issues$current_value, fixed = TRUE)))
    era_issues <- issues[issues$analysis_set == "Treatment-era sensitivity", ]
    expect_identical(nrow(era_issues), as.integer(sum(cohort$treatment_year < GEP_TESTING_ERA_START_YEAR)))
    expect_setequal(era_issues$record_key, as.character(cohort$id[cohort$treatment_year < GEP_TESTING_ERA_START_YEAR]))
    diagnostics_flow <- readxl::read_excel(paths$diagnostics, sheet = "Denominators_Flow")
    expect_identical(
        diagnostics_flow$n[grepl("^Treated before", diagnostics_flow$step)],
        as.numeric(results$treatment_era_excluded_n)
    )
    checks <- readxl::read_excel(paths$diagnostics, sheet = "Checks")
    expect_true(all(checks$status %in% c("no concern", "concern", "not assessed")))
    expect_false(anyNA(checks$status))

    # Predictor screening, the treatment-era counts and the group follow-up are each reported once.
    expect_false(any(c("Predictor retained", sprintf("Patients treated before %d", GEP_TESTING_ERA_START_YEAR)) %in% checks$check))
    limitations <- readxl::read_excel(paths$diagnostics, sheet = "Warnings_Limitations")
    expect_identical(limitations$limitation, c("Selection into testing", "Internal validation only", "Interval construction"))

    # One denominator per endpoint: the group outcome counts sum to the model analysis set.
    mfs_group_n <- sum(results$group_outcomes$n[results$group_outcomes$endpoint == "MFS"])
    expect_identical(mfs_group_n, flow$candidates[flow$analysis_set == "MFS"] - flow$not_eligible_for_endpoint[flow$analysis_set == "MFS"])
})

test_that("the no-information Brier score weights patients by their censoring weights", {
    rows <- data.frame(
        horizon_event = c(1L, 0L, 0L, 1L, 0L, 0L),
        ipcw_weight = c(4, 1, 1, 1, 2, 1),
        prediction = 0.5
    )
    weighted_prevalence <- sum(rows$ipcw_weight * rows$horizon_event) / sum(rows$ipcw_weight)
    expected <- sum(rows$ipcw_weight * (rows$horizon_event - weighted_prevalence)^2) / sum(rows$ipcw_weight)
    unweighted <- mean(rows$horizon_event)
    expect_gt(abs(weighted_prevalence - unweighted), 0.05)
    expect_equal(no_gep_repeat_metric(rows, rows$prediction, "null_brier"), expected, tolerance = 1e-10)
})

test_that("reverse Kaplan-Meier follow-up matches a hand calculation", {
    # Eight patients; censoring (status 0) at 60, 60, 80, 90; events at 10, 20, 30, 40.
    endpoint_data <- data.frame(
        gep_group = factor(rep("Not tested", 8), levels = GEP_GROUP_LEVELS),
        time_months = c(10, 20, 30, 40, 60, 60, 80, 90),
        status = c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L)
    )
    outcomes <- build_no_gep_group_outcomes(endpoint_data, "MFS", 60)
    # Reverse KM: the four censorings are the "events". At 60 months 4 are at risk and 2 are
    # censored, so S = 0.5 exactly; at 80 months 2 are at risk and 1 is censored, so S = 0.25.
    # S equals 0.5 from 60 until 80, so by the survfit convention the median is the midpoint, 70.
    expect_equal(outcomes$median_followup_months[outcomes$group == "Not tested"], 70)
})

test_that("the treatment-era sensitivity restricts each group to the testing era and says when it cannot", {
    era <- GEP_TESTING_ERA_START_YEAR
    endpoint_data <- data.frame(
        gep_group = factor(
            c(rep("Not tested", 10), rep("Class 1", 6), rep("Class 2", 4)),
            levels = GEP_GROUP_LEVELS
        ),
        treatment_year = c(rep(era - 3, 5), rep(era + 2, 5), rep(era + 3, 6), rep(era - 5, 4)),
        time_months = c(rep(12, 5), rep(70, 5), 20, rep(70, 5), rep(70, 4)),
        status = c(rep(1L, 5), rep(0L, 5), 1L, rep(0L, 5), rep(0L, 4))
    )
    table <- build_no_gep_treatment_era_sensitivity(endpoint_data, "MFS", 60)
    row <- function(group) table[table$group == group, ]

    # Half the group had a metastasis at 12 months but was treated before the era.
    not_tested <- row("Not tested")
    expect_identical(c(not_tested$patients, not_tested$patients_in_era, not_tested$patients_before_era), c(10L, 5L, 5L))
    expect_equal(not_tested$primary_risk, 0.5)
    expect_equal(not_tested$risk_in_era, 0)
    expect_equal(not_tested$difference_pct_points, -50)
    expect_match(not_tested$note, "No events by the horizon in the era")

    # A group that is entirely in the era is unchanged.
    class_1 <- row("Class 1")
    expect_equal(class_1$risk_in_era, class_1$primary_risk)
    expect_identical(class_1$patients_before_era, 0L)
    expect_match(class_1$note, "Every patient was treated in the era")

    # A group with no 5-year estimate says so, not that it is unchanged.
    short <- build_no_gep_treatment_era_sensitivity(
        transform(endpoint_data[endpoint_data$gep_group == "Class 1", ], time_months = 30, status = 0L), "MFS", 60
    )
    expect_match(short$note[short$group == "Class 1"], "Not followed to 60 months in the era")

    # A group with nobody in the era has no estimate, and the table says so.
    class_2 <- row("Class 2")
    expect_identical(class_2$patients_in_era, 0L)
    expect_true(is.na(class_2$risk_in_era))
    expect_false(class_2$estimable_in_era)
    expect_match(class_2$note, "No patient treated in the era")

    # Empty groups stay in the table.
    expect_setequal(table$group, GEP_GROUP_LEVELS)
    expect_error(
        prepare_no_gep_data(transform(synthetic_no_gep_cohort(), treatment_year = replace(treatment_year, 1, NA))),
        "treatment_year is missing for 1 patient"
    )
})

test_that("the model AUC figure keeps intervals that reach below the usual axis range", {
    performance <- data.frame(
        endpoint = c("MFS", "MFS"), performance_scope = "No GEP",
        model = c("Penalized, all retained predictors", "T stage alone"),
        auc = c(0.45, 0.60), auc_ci_lower = c(0.22, 0.45), auc_ci_upper = c(0.70, 0.75),
        stringsAsFactors = FALSE
    )
    plot <- plot_no_gep_model_auc(performance)
    built <- ggplot2::ggplot_build(plot)
    bars <- built$data[[2]]
    expect_equal(nrow(bars), 2L)
    expect_false(anyNA(bars$xmin))
    expect_lte(built$layout$panel_params[[1]]$x.range[1], 0.22)
})

test_that("the model AUC figure marks comparators by the policy labels, not by their wording", {
    comparator_labels <- c("Stage only" = "initial_t_stage_simple", "Diameter only" = "initial_tumor_diameter")
    rlang::local_bindings(GEP_NO_GEP_COMPARATOR_PREDICTORS = comparator_labels, .env = .GlobalEnv)
    performance <- data.frame(
        endpoint = "MFS", performance_scope = "No GEP",
        model = c("Penalized, all retained predictors", "Stage only", "Diameter only"),
        auc = c(0.6, 0.7, 0.65), auc_ci_lower = c(0.5, 0.55, 0.5), auc_ci_upper = c(0.75, 0.8, 0.78),
        stringsAsFactors = FALSE
    )
    points <- ggplot2::ggplot_build(plot_no_gep_model_auc(performance))$data[[3]]
    # Shape 16 is the filled ridge-model dot and shape 5 the open comparator diamond.
    expect_identical(sort(points$shape), c(5, 5, 16))
})

test_that("the no-GEP dataset is defined once in the policy", {
    workflow_text <- paste(readLines(here::here("scripts", "workflow", "objective_4_gep_analysis.R"), warn = FALSE), collapse = "\n")
    expect_true(grepl("GEP_NO_GEP_DATASET", workflow_text, fixed = TRUE))
    expect_false(grepl("uveal_melanoma_full_cohort\"", workflow_text, fixed = TRUE))
    expect_identical(GEP_NO_GEP_DATASET, PUBLISH_ANALYTIC_DATASET_IDS[[1]])
})

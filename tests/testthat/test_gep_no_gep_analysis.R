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
        built$counts$eligible - built$counts$removed_sparse_level - built$counts$removed_missing_predictor,
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
    expect_true(all(c("Treated in 2013 or later", "Age at diagnosis, years") %in% table$characteristic))
})

test_that("no-GEP outputs reconcile their flow with Record_Issues and keep identifiers out of the primary workbook", {
    cohort <- synthetic_no_gep_cohort()
    results <- run_no_gep_analysis(cohort, n_boot = 20L)
    out_dir <- withr::local_tempdir()
    paths <- write_no_gep_outputs(results, prepare_no_gep_data(cohort), out_dir, "test_")

    expect_true(all(file.exists(c(paths$primary, paths$diagnostics, paths$figures))))
    expect_identical(
        readxl::excel_sheets(paths$primary),
        c("Group_Outcomes_MFS", "Group_Outcomes_MSS", "Baseline_By_Group", "Model_Performance")
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
        as.integer(sum(flow$not_eligible_for_endpoint + flow$removed_sparse_level + flow$removed_missing_predictor))
    )
    checks <- readxl::read_excel(paths$diagnostics, sheet = "Checks")
    expect_true(all(checks$status %in% c("no concern", "concern", "not assessed")))
    expect_false(anyNA(checks$status))

    # One denominator per endpoint: the group outcome counts sum to the model analysis set.
    mfs_group_n <- sum(results$group_outcomes$n[results$group_outcomes$endpoint == "MFS"])
    expect_identical(mfs_group_n, flow$candidates[flow$analysis_set == "MFS"] - flow$not_eligible_for_endpoint[flow$analysis_set == "MFS"])
})

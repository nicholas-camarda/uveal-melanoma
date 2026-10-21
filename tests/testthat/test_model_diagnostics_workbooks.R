synthetic_model_cohort <- function(n = 120L, seed = 7L) {
    withr::with_seed(seed, {
        data <- data.frame(
            id = sprintf("P%03d", seq_len(n)),
            treatment_group = factor(sample(c("PBT", "GKSRS"), n, TRUE), levels = c("PBT", "GKSRS")),
            age = round(stats::rnorm(n, 62, 10)),
            location = factor(sample(c("Choroidal", "Ciliary", "Rare"), n, TRUE, prob = c(0.7, 0.27, 0.03))),
            time = stats::rexp(n, 0.02) + 1,
            event = stats::rbinom(n, 1, 0.35),
            stringsAsFactors = FALSE
        )
        data$age[c(5, 6)] <- NA
        data
    })
}

fit_synthetic_model <- function(kind, output_dir, data = synthetic_model_cohort(), min_level_count = THRESHOLD_RARITY) {
    exclusion <- apply_sparse_level_exclusions(
        data, variables = c("treatment_group", "location"), analysis_name = "unit_model",
        min_level_count = min_level_count, id_col = "id"
    )
    result <- suppressWarnings(generate_regression_table(
        data = exclusion$data, outcome_var = "event", predictor_vars = "treatment_group",
        confounders = c("age", "location"), model_type = kind,
        effect_measure = if (kind == "cox") "HR" else "OR", analysis_name = paste0("unit_", kind),
        dataset_name = "unit_data", output_dir = output_dir, prefix = "t_",
        time_var = if (kind == "cox") "time" else NULL, event_var = if (kind == "cox") "event" else NULL,
        treatment_var = "treatment_group",
        sparse_level_diagnostics = exclusion$sparse_level_diagnostics, filter_stats = exclusion$filter_stats
    ))
    list(result = result, exclusion = exclusion, path = file.path(output_dir, paste0("t_unit_", kind, "_diagnostics.xlsx")))
}

test_that("sparse-level exclusions list every removed patient by id and stop without an id", {
    data <- synthetic_model_cohort()
    exclusion <- apply_sparse_level_exclusions(data, variables = "location", analysis_name = "unit", id_col = "id")
    records <- exclusion$filter_stats$removed_records

    expect_equal(sort(records$record_key), sort(data$id[data$location == "Rare"]))
    expect_true(all(records$step == MODEL_RECORD_STEP_SPARSE_LEVEL))
    expect_true(all(records$field == "location" & records$current_value == "Rare"))
    expect_equal(length(unique(records$record_key)), exclusion$removed_row_count)

    expect_error(
        apply_sparse_level_exclusions(data[, setdiff(names(data), "id")], variables = "location", id_col = NULL),
        "without a patient identifier"
    )
    # With nothing to remove the identifier is not needed.
    expect_equal(
        nrow(apply_sparse_level_exclusions(data[data$location != "Rare", setdiff(names(data), "id")], variables = "location")$filter_stats$removed_records),
        0L
    )
})

test_that("a pre-fit missing-data filter is listed with the missing variables and checked against the filter", {
    data <- synthetic_model_cohort()
    kept <- data[!is.na(data$age), ]
    stats <- append_prefit_filter_records(
        list(removed_records = build_model_record_table()), data, kept, columns = c("age", "event"),
        id_col = "id", reason = "Missing a required variable"
    )
    expect_equal(stats$removed_records$record_key, data$id[is.na(data$age)])
    expect_true(all(stats$removed_records$field == "age" & stats$removed_records$step == MODEL_RECORD_STEP_PREFIT_FILTER))
    # A filter that removed different rows than the columns imply stops instead of listing the wrong patients.
    expect_error(
        append_prefit_filter_records(list(removed_records = build_model_record_table()), data, data, "age", "id", "x"),
        "removed 0 row"
    )
})

test_that("Cox and logistic workbooks share one sheet family and reconcile every excluded patient", {
    output_dir <- withr::local_tempdir()
    for (kind in c("cox", "logistic")) {
        fit <- fit_synthetic_model(kind, output_dir)
        sheets <- readxl::excel_sheets(fit$path)
        expect_equal(sheets, c("Specification_Coding", "Denominators_Flow", "Checks", "Term_support", "Record_Issues"))

        records <- readxl::read_xlsx(fit$path, sheet = "Record_Issues")
        data <- synthetic_model_cohort()
        expect_setequal(records$record_key[records$step == MODEL_RECORD_STEP_SPARSE_LEVEL], data$id[data$location == "Rare"])
        expect_setequal(records$record_key[records$step == MODEL_RECORD_STEP_COMPLETE_CASE], data$id[is.na(data$age)])
        expect_false(anyNA(records$record_key))

        flow <- readxl::read_xlsx(fit$path, sheet = "Denominators_Flow")
        expect_equal(
            flow$excluded_from_previous_n[flow$stage == "Model input after pre-fit exclusions"],
            length(unique(records$record_key[records$step == MODEL_RECORD_STEP_SPARSE_LEVEL]))
        )

        term_support <- readxl::read_xlsx(fit$path, sheet = "Term_support")
        # A continuous covariate is one row, never one row per observed value.
        expect_equal(sum(term_support$variable == "age"), 1L)
        expect_equal(term_support$level[term_support$variable == "age"], "Continuous")

        checks <- readxl::read_xlsx(fit$path, sheet = "Checks")
        expect_true(all(checks$status %in% c("no concern", "concern", "not assessed")))
        expect_true("Events per estimated coefficient" %in% checks$check)
    }
})

test_that("each count and attribute appears in one sheet of a model workbook", {
    fit <- fit_synthetic_model("logistic", withr::local_tempdir())
    sheets <- lapply(
        stats::setNames(readxl::excel_sheets(fit$path), readxl::excel_sheets(fit$path)),
        function(sheet) names(readxl::read_xlsx(fit$path, sheet = sheet))
    )
    shared_keys <- c("variable", "level", "reason", "n")
    repeated <- table(unlist(lapply(sheets, setdiff, shared_keys)))
    expect_true(all(repeated == 1L), info = paste(names(repeated)[repeated > 1], collapse = ", "))
    expect_false(any(grepl("estimate|ci_lower|ci_upper|p_value|hazard_ratio|odds_ratio", unlist(sheets))))
})

test_that("a thin model is flagged by the events-per-coefficient check", {
    data <- synthetic_model_cohort()
    data$event <- 0L
    data$event[1:6] <- 1L
    fit <- fit_synthetic_model("logistic", withr::local_tempdir(), data = data)
    checks <- readxl::read_xlsx(fit$path, sheet = "Checks")
    expect_equal(checks$status[checks$check == "Events per estimated coefficient"], "concern")
})

test_that("a Record_Issues count that does not match the workbook counts stops the build", {
    records <- build_model_record_table(
        record_key = c("P1", "P2"), key_source = "id", step = MODEL_RECORD_STEP_SPARSE_LEVEL,
        reason = "r", field = "location", current_value = "Rare"
    )
    expect_silent(reconcile_model_record_counts(records, 2, 0, "unit"))
    expect_error(reconcile_model_record_counts(records, 3, 0, "unit"), "removed before the fit but 2 are listed")
    expect_error(reconcile_model_record_counts(records, 2, 1, "unit"), "dropped by the model fit but 0 are listed")
})

test_that("a diagnostics write failure fails the run instead of being logged", {
    expect_error(write_diagnostics_workbook(NULL, tempfile(fileext = ".xlsx")), "No diagnostics were provided")
    fit <- fit_synthetic_model("logistic", withr::local_tempdir())
    expect_error(
        suppressWarnings(write_diagnostics_workbook(
            fit$result$diagnostics, file.path(withr::local_tempdir(), "missing_folder", "d.xlsx")
        )),
        "workbook was not written"
    )
})

test_that("an error inside the Cox regression layer fails the survival analysis instead of becoming a skipped workbook", {
    data <- synthetic_model_cohort(n = 200L)
    data <- data[!is.na(data$age) & data$location != "Rare", ]
    data$location <- droplevels(data$location)
    rlang::local_bindings(
        generate_regression_table = function(...) stop("regression layer failure", call. = FALSE),
        .env = .GlobalEnv
    )
    expect_error(
        suppressWarnings(analyze_time_to_event_outcomes(
            data = data, time_var = "time", event_var = "event", group_var = "treatment_group",
            confounders = "location", analysis_type = "all_patients", dataset_name = "unit_test",
            ylab = "Unit Test Survival", legend_labels = c("PBT", "GKSRS"),
            output_dirs = build_peer_review_objective1_output_dirs(withr::local_tempdir()),
            prefix = "unit_", route_key = "obj1_os"
        )),
        "regression layer failure"
    )
})

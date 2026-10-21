test_that("reviewer-response adjusted models use continuous age rather than dichotomized age", {
    expect_true("age_at_diagnosis" %in% confounders)
    expect_false("age_at_diagnosis_general_pop_median" %in% confounders)
})

test_that("visual acuity sensitivity consumes Objective 0 treatment year directly", {
    source_text <- paste(
        readLines(here::here("scripts", "analysis", "vision_safety_analysis.R"), warn = FALSE),
        collapse = "\n"
    )

    expect_match(source_text, '"Treatment year", "treatment_year", "treatment_year"', fixed = TRUE)
    expect_false(grepl("add_visual_acuity_treatment_year", source_text, fixed = TRUE))
    expect_false(grepl("treatment_year_centered", source_text, fixed = TRUE))
})

test_that("tracked methods document the restricted overlap-weighted sensitivity", {
    methods_text <- paste(
        readLines(here::here("docs", "STATISTICAL_METHODS.md"), warn = FALSE),
        collapse = "\n"
    )
    required_phrases <- c(
        "restricted cohort", "overlap weight", "overlap population",
        "age_at_diagnosis", "sex", "location", "initial_tumor_height",
        "initial_tumor_diameter", "srf", "treatment_year", "robust sandwich",
        "local recurrence", "metastatic progression", "overall survival",
        "progression-free survival", "exploratory sensitivity"
    )
    for (phrase in required_phrases) {
        expect_match(methods_text, phrase, ignore.case = TRUE)
    }
    expect_match(methods_text, "adjusted Cox models remain the primary analyses", ignore.case = TRUE)
    expect_match(methods_text, "death before.*censored", ignore.case = TRUE)
    expect_match(methods_text, "unweighted descriptive", ignore.case = TRUE)
    expect_match(methods_text, "does not estimate.*Fine-Gray", ignore.case = TRUE)
    expect_match(methods_text, "Objective 0.*derives.*once", ignore.case = TRUE)
    expect_false(grepl("primary propensity", methods_text, ignore.case = TRUE))
    expect_false(grepl("propensity[^\n]*(full cohort|stabilized IPTW|matching was used|ATE for the full)", methods_text, ignore.case = TRUE))
})

test_that("Objective 1 returns Cox-led local recurrence and metastasis time-to-event analyses", {
    pipeline <- get_objective1_pipeline()

    expect_true("recurrence_time_to_event" %in% names(pipeline$results))
    expect_true("mets_time_to_event" %in% names(pipeline$results))
    expect_s3_class(pipeline$results$recurrence_time_to_event$cox_model, "coxph")
    expect_s3_class(pipeline$results$mets_time_to_event$cox_model, "coxph")

    recurrence_summary <- file.path(
        pipeline$output_dirs$obj1_recurrence_cox,
        "test_local_recurrence_free_probability_effect_summary.xlsx"
    )
    metastasis_summary <- file.path(
        pipeline$output_dirs$obj1_mets_cox,
        "test_metastasis_free_survival_probability_effect_summary.xlsx"
    )
    expect_true(file.exists(recurrence_summary))
    expect_true(file.exists(metastasis_summary))

    if (file.exists(recurrence_summary) && file.exists(metastasis_summary)) {
        recurrence_rows <- readxl::read_xlsx(recurrence_summary)
        metastasis_rows <- readxl::read_xlsx(metastasis_summary)
        expect_true(any(recurrence_rows$effect_measure == "HR"))
        expect_true(any(metastasis_rows$effect_measure == "HR"))
        expect_false(any(recurrence_rows$effect_measure == "OR"))
        expect_false(any(metastasis_rows$effect_measure == "OR"))
    }
})

test_that("Objective 1 recurrence and metastasis descriptive summaries are not labeled co-primary", {
    pipeline <- get_objective1_pipeline()

    recurrence_summary_path <- file.path(pipeline$output_dirs$obj1_recurrence_event_support, "test_recurrence1_event_support_summary.xlsx")
    mets_summary_path <- file.path(pipeline$output_dirs$obj1_mets_event_support, "test_mets_progression_event_support_summary.xlsx")

    expect_true(file.exists(recurrence_summary_path))
    expect_true(file.exists(mets_summary_path))

    for (summary_path in c(recurrence_summary_path, mets_summary_path)) {
        if (!file.exists(summary_path)) {
            next
        }
        expect_true(all(c(
            "descriptive_event_counts",
            "cumulative_incidence",
            "competing_risk_support",
            "estimand_notes"
        ) %in% readxl::excel_sheets(summary_path)))

        estimand_notes <- readxl::read_xlsx(summary_path, sheet = "estimand_notes")
        expect_true(all(estimand_notes$role %in% c("descriptive_support", "supportive_time_to_event_context")))
        expect_false(any(estimand_notes$role == "co-primary"))
    }
})

test_that("Objective 1 writes five-year capped OS and PFS Cox sensitivity summaries", {
    pipeline <- get_objective1_pipeline()

    os_summary <- file.path(
        pipeline$output_dirs$obj1_os_sensitivity,
        "test_overall_survival_probability_5yr_capped_effect_summary.xlsx"
    )
    pfs_summary <- file.path(
        pipeline$output_dirs$obj1_pfs_sensitivity,
        "test_progression_free_survival_probability_5yr_capped_effect_summary.xlsx"
    )

    expect_true(file.exists(os_summary))
    expect_true(file.exists(pfs_summary))

    if (file.exists(os_summary) && file.exists(pfs_summary)) {
        os_rows <- readxl::read_xlsx(os_summary)
        pfs_rows <- readxl::read_xlsx(pfs_summary)
        expect_true(any(grepl("5-year capped", os_rows$model_label, fixed = TRUE)))
        expect_true(any(grepl("5-year capped", pfs_rows$model_label, fixed = TRUE)))
        expect_true(all(os_rows$n_patients >= os_rows$n_events, na.rm = TRUE))
        expect_true(all(pfs_rows$n_patients >= pfs_rows$n_events, na.rm = TRUE))
    }

    ph_or_skip_patterns <- c(
        "test_overall_survival_probability_5yr_capped_ph_diagnostics\\.xlsx$",
        "test_overall_survival_probability_5yr_capped_skipped.*\\.(xlsx|txt)$",
        "test_progression_free_survival_probability_5yr_capped_ph_diagnostics\\.xlsx$",
        "test_progression_free_survival_probability_5yr_capped_skipped.*\\.(xlsx|txt)$"
    )
    for (pattern in ph_or_skip_patterns[c(1, 2)]) {
        os_hits <- list.files(pipeline$output_dirs$obj1_os_sensitivity, pattern = pattern, full.names = TRUE)
        if (length(os_hits) > 0) break
    }
    for (pattern in ph_or_skip_patterns[c(3, 4)]) {
        pfs_hits <- list.files(pipeline$output_dirs$obj1_pfs_sensitivity, pattern = pattern, full.names = TRUE)
        if (length(pfs_hits) > 0) break
    }
    expect_true(length(os_hits) > 0, info = "Capped OS model must write PH diagnostics or an explicit skip artifact.")
    expect_true(length(pfs_hits) > 0, info = "Capped PFS model must write PH diagnostics or an explicit skip artifact.")
})

test_that("Objective 1 KM figures cap display at SURVIVAL_XAXIS_MAX_MONTHS while Cox models keep full follow-up", {
    pipeline <- get_objective1_pipeline()
    data <- pipeline$input_data

    recurrence_plot <- pipeline$results$recurrence_time_to_event$plot$plot
    pfs_plot <- pipeline$results$pfs_analysis$plot$plot
    recurrence_x_range <- ggplot2::ggplot_build(recurrence_plot)$layout$panel_params[[1]]$x.range
    pfs_x_range <- ggplot2::ggplot_build(pfs_plot)$layout$panel_params[[1]]$x.range
    axis_cap_tolerance <- SURVIVAL_XAXIS_MAX_MONTHS * 0.05

    expect_lte(max(recurrence_x_range), SURVIVAL_XAXIS_MAX_MONTHS + axis_cap_tolerance)
    expect_lte(max(pfs_x_range), SURVIVAL_XAXIS_MAX_MONTHS + axis_cap_tolerance)
    expect_gt(max(data$tt_recurrence_months, na.rm = TRUE), SURVIVAL_XAXIS_MAX_MONTHS)
    expect_gt(max(data$tt_pfs_months, na.rm = TRUE), SURVIVAL_XAXIS_MAX_MONTHS)
})

test_that("each route resolves to its own configured directory and nothing else", {
    output_dirs <- create_output_structure(withr::local_tempdir())

    expect_identical(resolve_route_output_dir(output_dirs, "obj1_mets"), output_dirs$obj1_mets)
    expect_identical(resolve_route_output_dir(output_dirs, "obj1_mets", "cox"), output_dirs$obj1_mets_cox)
    expect_identical(resolve_route_output_dir(output_dirs, "obj1_os", "km"), output_dirs$obj1_os_km)
    expect_false(identical(
        resolve_route_output_dir(output_dirs, "obj1_mets", "cox"),
        resolve_route_output_dir(output_dirs, "obj1_os", "cox")
    ))
    expect_error(resolve_route_output_dir(output_dirs, "not_a_route"), "Unknown output route")
    expect_error(resolve_route_output_dir(output_dirs, NULL), "Unknown output route")
    # An artifact the route's profile does not define has no directory and is an error.
    expect_error(resolve_route_output_dir(output_dirs, "obj2_nvg", "km"), "No output directory is configured for `obj2_nvg_km`")
})

test_that("a route whose directory is not configured fails instead of using a default location", {
    invalid_output_dirs <- list(
        NULL,
        list(baseline_characteristics = tempfile("baseline")),
        list(obj1_os = NULL),
        list(obj1_os = ""),
        list(obj1_os = character()),
        list(obj1_os = NA_character_)
    )

    for (output_dirs in invalid_output_dirs) {
        expect_error(
            resolve_route_output_dir(output_dirs, "obj1_os"),
            "No output directory is configured for `obj1_os`",
            fixed = TRUE
        )
    }
    expect_error(get_output_dir(list(), "obj1_root"), "No output directory is configured for `obj1_root`", fixed = TRUE)
})

test_that("the survival analysis requires explicit directories and a route", {
    data <- create_test_dataset()
    arguments <- list(
        data = data, time_var = "tt_death_months", event_var = "death_event",
        ylab = "Overall Survival Probability", dataset_name = "unit", prefix = "unit_"
    )
    expect_error(do.call(analyze_time_to_event_outcomes, c(arguments, list(route_key = "obj1_os"))), "output_dirs")
    expect_error(
        do.call(analyze_time_to_event_outcomes, c(arguments, list(output_dirs = NULL, route_key = "obj1_os"))),
        "No output directory is configured"
    )
    expect_error(
        do.call(analyze_time_to_event_outcomes, c(arguments, list(output_dirs = create_output_structure(withr::local_tempdir()), route_key = "obj9_x"))),
        "Unknown output route"
    )
})

test_that("the output structure creates exactly the routes in use, with no Objective 4 survival folders", {
    root <- withr::local_tempdir()
    output_dirs <- create_output_structure(root)

    legacy <- grep("^obj4_(mfs|mss|ph)", names(output_dirs), value = TRUE)
    expect_identical(legacy, character(0))
    expect_identical(output_dirs$obj4_gep, file.path(root, "04_GEP_Validation"))
    expect_identical(output_dirs$obj4_no_gep, file.path(root, "04_GEP_Validation", "d_no_gep"))
    created <- list.dirs(file.path(root, "04_GEP_Validation"), recursive = TRUE, full.names = FALSE)
    expect_setequal(created[nzchar(created)], "d_no_gep")

    # Every endpoint route registers one numbered folder per artifact type of its profile.
    for (route_key in names(ENDPOINT_PROFILE_BY_ROUTE_PREFIX)) {
        profile <- ENDPOINT_PROFILES[[ENDPOINT_PROFILE_BY_ROUTE_PREFIX[[route_key]]]]
        expect_true(all(paste0(route_key, "_", profile) %in% names(output_dirs)), info = route_key)
        expect_true(all(dir.exists(unlist(output_dirs[paste0(route_key, "_", profile)]))), info = route_key)
    }
})

test_that("output routing code has no fallback directory", {
    for (file in c("scripts/utils/output_utilities.R", "scripts/analysis/survival_outcomes.R", "scripts/analysis/rmst_visualization.R")) {
        text <- paste(readLines(here::here(file), warn = FALSE), collapse = "\n")
        expect_false(grepl("getwd()", text, fixed = TRUE), info = file)
        expect_false(grepl('"test_output"', text, fixed = TRUE), info = file)
        expect_false(grepl("determine_survival_output_dir|resolve_endpoint_output_dir|obj4_mfs|obj4_mss", text), info = file)
    }
})

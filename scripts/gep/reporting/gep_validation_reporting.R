# Reporting for the definitive-GEP 5-year validation: workbooks and figures.
# Statistical work lives in scripts/gep/validation/gep_definitive_validation.R.

#' Count patients through the GEP validation flow
#'
#' @param data Cohort data frame with `gep_group` and the endpoint eligibility flags.
#' @param results Output of `run_gep_definitive_validation()`.
#' @return Data frame with `step`, `n`, and `note`; the first row is all patients and
#'   the GEP availability groups sum to it.
summarize_gep_validation_flow <- function(data, results) {
    group_counts <- table(factor(data$gep_group, levels = GEP_GROUP_LEVELS))
    primary <- data$gep_group %in% GEP_DEFINITIVE_GROUP_LEVELS
    primary_n <- sum(primary)
    flow <- data.frame(
        step = c(
            "All patients in cohort",
            paste("GEP group:", GEP_GROUP_LEVELS),
            "Definitive Class 1/2 with PRAME reported (primary validation set)",
            "Eligible for the MFS endpoint", "Analysed for MFS",
            "Eligible for the MSS endpoint", "Analysed for MSS",
            "Primary-set patients not eligible for the MFS endpoint",
            "Primary-set patients not eligible for the MSS endpoint",
            paste("Biopsy timing in the primary set:", GEP_TEST_TIMING_LEVELS),
            sprintf("Biopsy more than %d days before treatment in the primary set (date review)", GEP_TEST_REVIEW_BEFORE_DAYS)
        ),
        n = c(
            nrow(data), as.integer(group_counts), as.integer(primary_n),
            sum(data$mfs_analysis_eligible %in% TRUE), nrow(results$mfs$prepared),
            sum(data$mss_analysis_eligible %in% TRUE), nrow(results$mss$prepared),
            sum(primary & !(data$mfs_analysis_eligible %in% TRUE)), sum(primary & !(data$mss_analysis_eligible %in% TRUE)),
            as.integer(table(factor(data$gep_test_timing[primary], levels = GEP_TEST_TIMING_LEVELS))),
            sum(primary & data$gep_test_days < -GEP_TEST_REVIEW_BEFORE_DAYS, na.rm = TRUE)
        ),
        note = c(
            "", rep("", length(GEP_GROUP_LEVELS)),
            "Class 1 + Class 2; patients with a known class but PRAME not reported keep no imputed prediction and are excluded",
            "Definitive call, valid imported prediction, metastasis-free at baseline",
            "Equals eligible; an eligible patient without a validation group, follow-up time or status stops the run",
            "Definitive call, valid imported prediction, melanoma-death follow-up available",
            "Equals eligible; an eligible patient without a validation group, follow-up time or status stops the run",
            "Listed in Record_Issues with the field values behind the exclusion",
            "Listed in Record_Issues with the field values behind the exclusion",
            c("", "", "", "Flagged in Record_Issues; excluded in the late-biopsy sensitivity analysis", "Flagged in Record_Issues"),
            "Counted in the timing row above; flagged in Record_Issues and kept in the analysis"
        ),
        stringsAsFactors = FALSE
    )
    flow
}

#' List patients outside the primary validation set
#'
#' One row per patient not in the Class 1/Class 2 primary set, with the step and
#' reason, the raw laboratory label, and nothing else; plus the primary-set patients who are
#' not eligible for an endpoint (with the fields that decide it) or whose biopsy timing needs a look. When the cohort has an
#' identifier column the key is that identifier; otherwise it is the row position and
#' `key_source` says so. Contains identifiers: keep it in the diagnostics workbook,
#' which is part of the results package.
#'
#' @param data Cohort data frame with `gep_group`, `biopsy1_gep_raw`, `gep_test_timing`
#'   and `gep_test_days`.
#' @return Data frame with `record_key`, `key_source`, `step`, `reason`, `field`,
#'   `current_value`.
build_gep_record_issues <- function(data) {
    id_col <- pick_sparse_level_id_col(data)
    record_key <- if (is.null(id_col)) as.character(seq_len(nrow(data))) else as.character(data[[id_col]])
    key_source <- if (is.null(id_col)) "row position (no identifier column)" else id_col
    outside <- !(data$gep_group %in% GEP_DEFINITIVE_GROUP_LEVELS)
    n_outside <- sum(outside)
    group_rows <- data.frame(
        record_key = record_key[outside],
        key_source = rep(key_source, n_outside),
        step = rep("GEP availability group", n_outside),
        reason = as.character(data$gep_group[outside]),
        field = rep("biopsy1_gep_raw", n_outside),
        current_value = ifelse(is.na(data$biopsy1_gep_raw[outside]), "(no result recorded)", as.character(data$biopsy1_gep_raw[outside])),
        stringsAsFactors = FALSE
    )
    # Primary-set patients whose biopsy timing needs a look (kept in the analysis).
    review_before <- !is.na(data$gep_test_days) & data$gep_test_days < -GEP_TEST_REVIEW_BEFORE_DAYS
    flagged <- !outside & (as.character(data$gep_test_timing) %in% c(GEP_TEST_TIMING_LATE, GEP_TEST_TIMING_LEVELS[[5]]) | review_before)
    n_flagged <- sum(flagged)
    timing_rows <- data.frame(
        record_key = record_key[flagged],
        key_source = rep(key_source, n_flagged),
        step = rep("Biopsy timing (patient kept in the primary set)", n_flagged),
        reason = ifelse(review_before[flagged], GEP_TEST_REVIEW_BEFORE_REASON, as.character(data$gep_test_timing[flagged])),
        field = rep("gep_test_days", n_flagged),
        current_value = ifelse(is.na(data$gep_test_days[flagged]), "missing", sprintf("%.0f days from treatment", data$gep_test_days[flagged])),
        stringsAsFactors = FALSE
    )
    # Primary-set patients who cannot enter an endpoint, with the fields that decide it.
    not_eligible_rows <- function(endpoint, eligible_column, field_columns) {
        idx <- !outside & !(data[[eligible_column]] %in% TRUE)
        fields <- intersect(field_columns, names(data))
        data.frame(
            record_key = record_key[idx],
            key_source = rep(key_source, sum(idx)),
            step = rep(sprintf("Primary-set patient not eligible for the %s endpoint", endpoint), sum(idx)),
            reason = rep("Metastatic at baseline, or follow-up time or event status missing", sum(idx)),
            field = rep(paste(fields, collapse = ", "), sum(idx)),
            current_value = format_record_values(data[idx, , drop = FALSE], fields),
            stringsAsFactors = FALSE
        )
    }
    endpoint_rows <- rbind(
        not_eligible_rows("MFS", "mfs_analysis_eligible", c("mfs_analysis_eligible", "mets_free_at_baseline", "tt_mets_months_analysis", "mfs_event_type")),
        not_eligible_rows("MSS", "mss_analysis_eligible", c("mss_analysis_eligible", "tt_death_months", "mss_event_type"))
    )
    rbind(group_rows, endpoint_rows, timing_rows)
}

#' Evaluate the quality checks that qualify the validation estimates
#'
#' Every check carries a result and a status (`no concern`, `concern`, or
#' `not assessed`) so a check that could not run never looks as if it passed.
#'
#' @param results Output of `run_gep_definitive_validation()`.
#' @return Data frame with `endpoint`, `check`, `scope`, `result`, `status`, `reason`.
build_gep_validation_checks <- function(results) {
    rows <- list()
    add <- function(endpoint, check, scope, result, status, reason) {
        rows[[length(rows) + 1L]] <<- data.frame(
            endpoint = endpoint, check = check, scope = scope, result = result, status = status, reason = reason,
            stringsAsFactors = FALSE
        )
    }
    for (endpoint in c("MFS", "MSS")) {
        part <- results[[tolower(endpoint)]]
        for (i in seq_len(nrow(part$group_table))) {
            g <- part$group_table[i, ]
            if (grepl("^Overall", g$group)) next
            add(endpoint, "Followed to the horizon", g$group,
                sprintf("maximum follow-up %.0f months", g$max_followup_months),
                if (isTRUE(g$estimable)) "no concern" else "concern",
                if (isTRUE(g$estimable)) "Observed risk estimated at the horizon" else g$note)
            add(endpoint, "Patients known event-free at the horizon", g$group,
                as.character(g$event_free_through_horizon),
                if (g$event_free_through_horizon >= GEP_MIN_AT_RISK_AT_HORIZON) "no concern" else "concern",
                sprintf("Fewer than %d patients known event-free at the horizon makes the observed risk imprecise", GEP_MIN_AT_RISK_AT_HORIZON))
            add(endpoint, "Events by the horizon", g$group, as.character(g$events_by_horizon),
                if (g$events_by_horizon > 0) "no concern" else "concern",
                "No events: the observed risk is 0 and its interval cannot be estimated")
        }
        d <- part$discrimination
        if (identical(endpoint, "MFS")) {
            add(endpoint, "Tied pairs in Harrell's C", "All patients",
                if (isTRUE(d$estimable)) sprintf("%.1f%% of comparable pairs tied on imported risk", 100 * d$tied_pair_share) else "not computed",
                if (!isTRUE(d$estimable)) "not assessed" else if (d$tied_pair_share > 0) "concern" else "no concern",
                "The imported risk has four values; tied pairs count one half, so C is a ranking of four groups")
        }
        add(endpoint, "Discrimination estimable", "All patients",
            if (isTRUE(d$estimable)) "estimated" else d$reason,
            if (isTRUE(d$estimable)) "no concern" else "not assessed",
            "Minimum events and variation in imported risk required")
        if (isTRUE(d$estimable)) {
            add(endpoint, "Bootstrap draws", "All patients",
                sprintf("%d valid, %d failed", d$valid_draws, d$failed_draws),
                if (d$failed_draws == 0L) "no concern" else "concern",
                "Failed draws lack an event or imported-risk contrast")
            if (d$events < 2 * GEP_MIN_DISCRIMINATION_EVENTS) {
                add(endpoint, "Events supporting discrimination", "All patients", as.character(d$events), "concern",
                    "Very few events: the interval is wide and the estimate unstable")
            }
        }
        ph <- part$proportional_hazards
        add(endpoint, "Proportional hazards (Schoenfeld, imported risk)", "All patients",
            if (is.na(ph$p_value)) "not computed" else sprintf("p = %.2f", ph$p_value), ph$status, ph$reason)
    }
    late_n <- results$late_biopsy_sensitivity$patients_excluded[["MFS"]]
    add("MFS", "GEP sample taken more than 180 days after treatment", "Primary set",
        sprintf("%d patients", late_n), if (late_n == 0L) "no concern" else "concern",
        "Such a sample (for example at local recurrence) describes a different tumor state; see Sensitivity_Late_Biopsy")
    describe_check_reasons(dplyr::bind_rows(rows))
}

#' Keep the reason text only for checks that raise a concern or were not assessed
#'
#' @param checks Data frame with `status` and `reason`.
#' @return `checks` with `reason` set to "No issue found" where the status is `no concern`.
describe_check_reasons <- function(checks) {
    checks$reason <- ifelse(checks$status == "no concern", "No issue found", checks$reason)
    checks
}

#' Describe the methods behind each reported quantity
#'
#' @param horizon_months Numeric horizon in months.
#' @return Data frame with method, estimand, effect measure, units, direction and a
#'   plain-language interpretation for each reported quantity.
build_gep_validation_specification <- function(horizon_months) {
    years <- horizon_months / 12
    data.frame(
        quantity = c(
            "Observed MFS risk", "Observed MSS risk", "Imported risk",
            "Difference", "Overall row", "Harrell's C (MFS)", "Competing-risk AUC (MSS)",
            "Proportional-hazards check", "Late-biopsy sensitivity", "Unknown-cause-of-death sensitivity"
        ),
        method = c(
            sprintf("Kaplan-Meier risk (1 - survival) at %d months, log-log 95%% interval", horizon_months),
            sprintf("Aalen-Johansen cumulative incidence at %d months, non-melanoma death competing, log-log interval", horizon_months),
            "1 minus the imported 5-year survival probability for the group (laboratory class-level estimate, not a patient-level prediction)",
            "Observed minus imported risk",
            "Group risks weighted by group size",
            sprintf("Concordance of imported risk with follow-up censored at %d months; ties count one half; patient-level bootstrap interval", horizon_months),
            sprintf("Cumulative/dynamic time-dependent AUC at %d months, inverse-probability-of-censoring weights, non-melanoma deaths as controls; patient-level bootstrap interval", horizon_months),
            "Schoenfeld residual test of imported risk in a Cox model (cause-specific for MSS)",
            sprintf("Group tables and discrimination recomputed without patients whose first biopsy was more than %d days after treatment", GEP_TEST_TIMING_LATE_DAYS),
            "Observed melanoma-specific risk recomputed counting deaths with an unknown recorded cause (Unknown, Unknown_Cancer) as melanoma deaths"
        ),
        estimand = c(
            sprintf("Probability of metastasis by %.0f years (death before metastasis is censored)", years),
            sprintf("Probability of melanoma death by %.0f years in the presence of competing death", years),
            "Published risk for the Class x PRAME group",
            "Observed risk minus imported risk", "Group-size-weighted observed risk",
            "Probability that a patient with an earlier event has the higher imported risk",
            "Probability that a patient with melanoma death by the horizon has a higher imported risk than a patient without",
            "Whether the effect of imported risk is constant over time",
            "Whether the results depend on patients sampled long after treatment",
            "Whether the results depend on how deaths of unknown cause are counted"
        ),
        units = c("proportion", "proportion", "proportion", "percentage points", "proportion", "0.5 = chance, 1 = perfect", "0.5 = chance, 1 = perfect", "p-value", "as the quantity", "proportion"),
        direction = c(
            "Higher = more metastasis", "Higher = more melanoma deaths", "Higher = worse", "Positive = observed worse than imported",
            "Higher = worse", "Higher imported risk = earlier event", "Higher imported risk = melanoma death", "p < 0.05 flagged as concern",
            "Difference = value without late samples minus primary value", "Difference = sensitivity risk minus primary risk"
        ),
        interpretation = c(
            "Compare with the imported risk in the same row; wide intervals reflect few patients and events",
            "Compare with the imported risk in the same row; wide intervals reflect few patients and events",
            "Four distinct values, so discrimination reflects ranking of four groups",
            "If the imported risk lies inside the observed 95% interval, the data are compatible with it",
            "No interval is reported; see group rows",
            "Because the predictor has four values, C describes how well the groups are ordered, not individual-level prediction",
            "Based on few melanoma deaths; see events and interval",
            "A concern would mean the imported risk's effect changes with follow-up time",
            "Small differences mean the late samples do not drive the results (see the sheet for each group)",
            "The primary analysis counts these deaths as competing; this sheet shows the opposite extreme"
        ),
        stringsAsFactors = FALSE
    )
}

#' Compile the limitations that apply to this cohort's estimates
#'
#' @param checks Data frame from `build_gep_validation_checks()`.
#' @param dataset_name Character cohort identifier.
#' @return Data frame with `limitation` and `detail`.
build_gep_validation_limitations <- function(checks, dataset_name) {
    sparse_groups <- unique(checks$scope[checks$check == "Patients known event-free at the horizon" & checks$status == "concern"])
    no_event_groups <- unique(checks$scope[checks$check == "Events by the horizon" & checks$status == "concern"])
    data.frame(
        limitation = c(
            "Imported predictions are group-level", "No extrapolation", "Sparse groups", "Groups without events", "Overlapping cohorts"
        ),
        detail = c(
            "Each patient carries one of four published class-level risks. Agreement with the imported values is not a patient-level calibration.",
            "Only the 5-year horizon is reported. Risks are not carried past a group's last follow-up; a group not followed to 60 months has no observed risk.",
            if (length(sparse_groups) > 0) paste("Fewer than", GEP_MIN_AT_RISK_AT_HORIZON, "patients known event-free at 60 months:", paste(sparse_groups, collapse = "; ")) else "None",
            if (length(no_event_groups) > 0) paste("No events by 60 months, so the interval cannot be estimated:", paste(no_event_groups, collapse = "; ")) else "None",
            sprintf("%s is one of three overlapping analytic views of the same patients; the views are not independent replications.", dataset_name)
        ),
        stringsAsFactors = FALSE
    )
}

#' Plot cumulative risk curves by group with the imported 5-year risk
#'
#' One step curve per Class x PRAME group (Kaplan-Meier risk for MFS, Aalen-Johansen
#' cumulative incidence for MSS), drawn to the horizon, with the imported
#' horizon risk as an open diamond.
#'
#' @param prepared Data frame from `prepare_gep_validation_endpoint()`.
#' @param group_table Group table from `build_gep_group_table()`.
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @param horizon_months Numeric horizon in months.
#' @param cohort_label Character cohort name for the subtitle.
#' @return A `ggplot` object.
plot_gep_validation_curves <- function(prepared, group_table, endpoint, horizon_months, cohort_label) {
    curves <- lapply(levels(prepared$validation_group), function(group_name) {
        g <- prepared[prepared$validation_group == group_name, , drop = FALSE]
        if (nrow(g) == 0L) return(NULL)
        if (identical(endpoint, "MFS")) {
            fit <- survival::survfit(survival::Surv(time_months, status == 1L) ~ 1, data = g)
            steps <- data.frame(time = c(0, fit$time), risk = c(0, 1 - fit$surv))
        } else if (!any(g$status == 1L)) {
            # cuminc() needs at least one event; with none the cumulative incidence is 0.
            steps <- data.frame(time = 0, risk = 0)
        } else {
            fit <- cmprsk::cuminc(ftime = g$time_months, fstatus = g$status, cencode = 0)
            curve <- fit[[grep(" 1$", names(fit))[1]]]
            steps <- data.frame(time = c(0, curve$time), risk = c(0, curve$est))
        }
        steps <- steps[steps$time <= horizon_months, , drop = FALSE]
        last_time <- min(max(g$time_months), horizon_months)
        steps <- rbind(steps, data.frame(time = last_time, risk = steps$risk[nrow(steps)]))
        event_free <- group_table$event_free_through_horizon[match(group_name, group_table$group)]
        steps$group <- sprintf("%s (n = %d; %d event-free at %.0f y)", group_name, nrow(g), event_free, horizon_months / 12)
        steps$group_name <- group_name
        steps
    })
    curve_data <- dplyr::bind_rows(curves)
    curve_data$time_years <- curve_data$time / 12
    labels <- unique(curve_data[, c("group", "group_name")])
    palette <- get_gep_class_palette(labels$group_name)
    colors <- stats::setNames(unname(palette[labels$group_name]), labels$group)
    imported <- group_table[!grepl("^Overall", group_table$group), , drop = FALSE]
    imported$label <- vapply(imported$group, function(x) labels$group[match(x, labels$group_name)], character(1))
    ylab <- if (identical(endpoint, "MFS")) "Cumulative risk of metastasis" else "Cumulative incidence of melanoma death"
    ggplot2::ggplot(curve_data, ggplot2::aes(x = .data$time_years, y = .data$risk, color = .data$group)) +
        ggplot2::geom_step(linewidth = 1.1) +
        ggplot2::geom_point(
            data = imported, ggplot2::aes(x = horizon_months / 12, y = .data$imported_risk, color = .data$label),
            shape = 5, size = 3.5, stroke = 1.2, inherit.aes = FALSE
        ) +
        ggplot2::scale_color_manual(values = colors, name = NULL) +
        ggplot2::scale_x_continuous(limits = c(0, horizon_months / 12 + 0.15), breaks = 0:(horizon_months / 12)) +
        ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
        ggplot2::labs(
            x = "Years since treatment", y = ylab,
            title = sprintf("%s by imported GEP group", ylab), subtitle = cohort_label,
            caption = "Open diamonds: imported 5-year risk. Curves end at each group's last follow-up (no extrapolation)."
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(legend.position = "bottom", legend.direction = "vertical", panel.grid.minor = ggplot2::element_blank())
}

#' Plot observed against imported 5-year risk for each group
#'
#' @param group_tables Named list of group tables (one per panel, for example per
#'   cohort or endpoint) as produced by `build_gep_group_table()`; panels keep list order.
#' @param title Character plot title.
#' @return A `ggplot` object.
plot_gep_observed_vs_imported <- function(group_tables, title = "Observed vs imported 5-year risk") {
    plot_data <- dplyr::bind_rows(group_tables, .id = "panel") %>%
        dplyr::filter(!grepl("^Overall", .data$group)) %>%
        dplyr::mutate(
            panel = factor(.data$panel, levels = names(group_tables)),
            group = factor(.data$group, levels = rev(GEP_VALIDATION_GROUP_LEVELS)),
            no_events = .data$events_by_horizon == 0
        )
    colors <- get_gep_class_palette(levels(plot_data$group))
    ggplot2::ggplot(plot_data, ggplot2::aes(y = .data$group)) +
        ggplot2::geom_errorbar(
            ggplot2::aes(xmin = .data$observed_ci_lower, xmax = .data$observed_ci_upper, color = .data$group),
            width = 0.2, linewidth = 0.9, orientation = "y", na.rm = TRUE
        ) +
        ggplot2::geom_point(ggplot2::aes(x = .data$observed_risk, color = .data$group, shape = .data$no_events), size = 3.4, stroke = 1.2, na.rm = TRUE) +
        ggplot2::geom_point(ggplot2::aes(x = .data$imported_risk), shape = 5, size = 3.6, stroke = 1.1, color = "grey20") +
        ggplot2::facet_wrap(~panel) +
        ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
        ggplot2::scale_color_manual(values = colors, guide = "none") +
        ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
        ggplot2::labs(
            x = "5-year risk", y = NULL, title = title,
            caption = wrap_plot_text("Filled points: observed risk with 95% interval. Open circles: no events occurred, so no interval can be estimated. Open diamonds: imported risk.", 100)
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title.position = "plot", plot.caption.position = "plot")
}

#' Write the validation workbooks and figures for one cohort
#'
#' Primary workbook (publishable): group tables and discrimination. Diagnostics
#' workbook (part of the results package; `Record_Issues` lists patients by id): specification,
#' denominators, checks, limitations, and record issues.
#'
#' @param results Output of `run_gep_definitive_validation()`.
#' @param data Cohort data frame.
#' @param dataset_name Character cohort identifier.
#' @param output_dirs List of output directories (uses `obj4_mfs` to locate the GEP root).
#' @param prefix Character file prefix.
#' @return Invisibly, a list of written file paths.
write_gep_validation_outputs <- function(results, data, dataset_name, output_dirs, prefix) {
    required <- c("gep_group", "biopsy1_gep_raw", "gep_test_timing", "gep_test_days", "mfs_analysis_eligible", "mss_analysis_eligible")
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0L) {
        stop(sprintf("GEP validation reporting requires columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
    if (anyNA(data$gep_group)) {
        stop(sprintf("gep_group is missing for %d patient(s); every patient must be classified.", sum(is.na(data$gep_group))), call. = FALSE)
    }
    gep_root <- dirname(output_dirs$obj4_mfs)
    ensure_output_dir(gep_root)
    cohort_label <- switch(
        dataset_name,
        uveal_melanoma_full_cohort = "Full cohort",
        uveal_melanoma_restricted_cohort = "Restricted cohort",
        uveal_melanoma_gksrs_only_cohort = "GKSRS-only cohort",
        dataset_name
    )
    horizon <- results$horizon_months

    discrimination <- data.frame(
        endpoint = c("MFS", "MSS"),
        metric = c("Harrell's C (ties = 0.5)", "Competing-risk AUC"),
        estimate = c(results$mfs$discrimination$c_index, results$mss$discrimination$auc),
        ci_lower = c(results$mfs$discrimination$ci_lower, results$mss$discrimination$ci_lower),
        ci_upper = c(results$mfs$discrimination$ci_upper, results$mss$discrimination$ci_upper),
        patients = c(nrow(results$mfs$prepared), nrow(results$mss$prepared)),
        events_by_horizon = c(results$mfs$discrimination$events, results$mss$discrimination$events),
        tied_pair_share = c(results$mfs$discrimination$tied_pair_share, NA_real_),
        status = c(
            if (results$mfs$discrimination$estimable) "estimated" else results$mfs$discrimination$reason,
            if (results$mss$discrimination$estimable) "estimated" else results$mss$discrimination$reason
        ),
        stringsAsFactors = FALSE
    )
    checks <- build_gep_validation_checks(results)

    primary_path <- file.path(gep_root, paste0(prefix, "gep_validation_5yr.xlsx"))
    write_gep_workbook(
        list(
            Group_Table_MFS = results$mfs$group_table,
            Group_Table_MSS = results$mss$group_table,
            Discrimination = discrimination,
            Sensitivity_Late_Biopsy = results$late_biopsy_sensitivity$table,
            Sensitivity_Cause_of_Death = results$cause_of_death_sensitivity
        ),
        primary_path
    )
    diagnostics_path <- file.path(gep_root, paste0(prefix, "gep_validation_diagnostics.xlsx"))
    write_gep_workbook(
        list(
            Specification_Coding = build_gep_validation_specification(horizon),
            Denominators_Flow = summarize_gep_validation_flow(data, results),
            Checks = checks,
            Warnings_Limitations = build_gep_validation_limitations(checks, dataset_name),
            Record_Issues = build_gep_record_issues(data)
        ),
        diagnostics_path
    )

    figure_paths <- character()
    for (endpoint in c("MFS", "MSS")) {
        part <- results[[tolower(endpoint)]]
        curve_path <- file.path(gep_root, paste0(prefix, "gep_validation_", tolower(endpoint), "_curves.png"))
        ggplot2::ggsave(
            curve_path, plot_gep_validation_curves(part$prepared, part$group_table, endpoint, horizon, cohort_label),
            width = 8, height = 7, dpi = 300, bg = "white"
        )
        figure_paths <- c(figure_paths, curve_path)
    }
    comparison_path <- file.path(gep_root, paste0(prefix, "gep_validation_observed_vs_imported.png"))
    ggplot2::ggsave(
        comparison_path,
        plot_gep_observed_vs_imported(
            list(`Metastasis (MFS)` = results$mfs$group_table, `Melanoma death (MSS)` = results$mss$group_table),
            title = sprintf("%s: observed vs imported 5-year risk", cohort_label)
        ),
        width = 10, height = 4.5, dpi = 300, bg = "white"
    )
    invisible(list(primary = primary_path, diagnostics = diagnostics_path, figures = c(figure_paths, comparison_path)))
}

#' Cross-cohort observed-vs-imported MFS figure
#'
#' Reads the group tables from each cohort's primary validation workbook and draws
#' one panel per cohort. A workbook written before `run_started` belongs to an earlier
#' run, so the figure stops instead of mixing runs.
#'
#' @param run_started `POSIXct` time the current run began.
#' @param output_dir Character directory to write the PNG to.
#' @param cohort_sources Data frame with `cohort_dir`, `prefix`, `cohort_label`.
#' @return Invisibly, the written file path.
create_gep_observed_vs_imported_panel_report <- function(run_started,
                                                         output_dir = file.path(MERGED_TABLES_DIR, "objective4_cross_cohort_figures"),
                                                         cohort_sources = data.frame(
                                                             cohort_dir = c("uveal_full", "uveal_restricted", "gksrs"),
                                                             prefix = c("full_cohort_", "restricted_cohort_", "gksrs_only_cohort_"),
                                                             cohort_label = c("Full cohort", "Restricted cohort", "GKSRS-only cohort"),
                                                             stringsAsFactors = FALSE
                                                         )) {
    tables <- lapply(seq_len(nrow(cohort_sources)), function(i) {
        path <- file.path(OUTPUT_DIR, cohort_sources$cohort_dir[i], "04_GEP_Validation",
                          paste0(cohort_sources$prefix[i], "gep_validation_5yr.xlsx"))
        if (!file.exists(path)) {
            stop(sprintf("Validation workbook not found: %s", path), call. = FALSE)
        }
        if (file.mtime(path) < run_started) {
            stop(sprintf("Validation workbook was not written in this run: %s", path), call. = FALSE)
        }
        readxl::read_excel(path, sheet = "Group_Table_MFS")
    })
    names(tables) <- cohort_sources$cohort_label
    ensure_output_dir(output_dir)
    figure_path <- file.path(output_dir, "objective4_gep_observed_vs_imported_mfs.png")
    ggplot2::ggsave(
        figure_path,
        plot_gep_observed_vs_imported(tables, title = "Observed vs imported 5-year metastasis risk (MFS), by cohort"),
        width = 12, height = 4.5, dpi = 300, bg = "white"
    )
    invisible(figure_path)
}

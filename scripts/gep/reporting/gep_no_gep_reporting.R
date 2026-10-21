# Reporting for the no-GEP analysis: workbooks and figures.
# Statistical work lives in scripts/gep/validation/gep_no_gep_analysis.R.

#' Patient flow through the no-GEP analysis
#'
#' @param data Prepared cohort with `gep_group`.
#' @param results Output of `run_no_gep_analysis()`.
#' @return Data frame with `step` and `n`; GEP groups sum to the cohort, and each
#'   analysis set reconciles as candidates minus exclusions equals model patients.
summarize_no_gep_flow <- function(data, results) {
    group_counts <- table(factor(data$gep_group, levels = GEP_GROUP_LEVELS))
    flow <- results$model_flow
    set_rows <- do.call(rbind, lapply(seq_len(nrow(flow)), function(i) {
        f <- flow[i, ]
        data.frame(
            step = c(
                sprintf("%s: candidates", f$analysis_set),
                sprintf("%s: not eligible for the endpoint", f$analysis_set),
                sprintf("%s: removed, factor level with fewer than %d patients", f$analysis_set, GEP_NO_GEP_MIN_LEVEL_COUNT),
                sprintf("%s: removed, missing baseline predictor", f$analysis_set),
                sprintf("%s: patients in the model", f$analysis_set)
            ),
            n = c(f$candidates, f$not_eligible_for_endpoint, f$removed_sparse_level, f$removed_missing_predictor, f$model_n),
            stringsAsFactors = FALSE
        )
    }))
    rbind(
        data.frame(step = "All patients in cohort", n = nrow(data), stringsAsFactors = FALSE),
        data.frame(step = paste("GEP group:", GEP_GROUP_LEVELS), n = as.integer(group_counts), stringsAsFactors = FALSE),
        set_rows,
        data.frame(
            step = sprintf("Treated before %d (left out of the treatment-era sensitivity only)", GEP_TESTING_ERA_START_YEAR),
            n = results$treatment_era_excluded_n, stringsAsFactors = FALSE
        )
    )
}

#' Checks that qualify the no-GEP estimates
#'
#' Each check has a result and a status (`no concern`, `concern`, `not assessed`).
#'
#' @param results Output of `run_no_gep_analysis()`.
#' @return Data frame with `endpoint`, `check`, `scope`, `result`, `status`, `reason`.
build_no_gep_checks <- function(results) {
    rows <- list()
    add <- function(endpoint, check, scope, result, status, reason) {
        rows[[length(rows) + 1L]] <<- data.frame(
            endpoint = endpoint, check = check, scope = scope, result = result, status = status, reason = reason,
            stringsAsFactors = FALSE
        )
    }
    outcomes <- results$group_outcomes
    for (i in seq_len(nrow(outcomes))) {
        o <- outcomes[i, ]
        if (o$n == 0L) next
        add(o$endpoint, "Followed to the horizon", o$group,
            sprintf("maximum follow-up %.0f months", o$max_followup_months),
            if (isTRUE(o$estimable)) "no concern" else "concern",
            if (isTRUE(o$estimable)) "Observed risk estimated at the horizon" else o$note)
        add(o$endpoint, "Patients known event-free at the horizon", o$group, as.character(o$event_free_through_horizon),
            if (o$event_free_through_horizon >= GEP_MIN_AT_RISK_AT_HORIZON) "no concern" else "concern",
            sprintf("Fewer than %d known event-free patients makes the observed risk imprecise", GEP_MIN_AT_RISK_AT_HORIZON))
    }
    era <- results$treatment_era_sensitivity
    era <- era[era$endpoint == "MFS", , drop = FALSE]
    for (i in seq_len(nrow(era))) {
        e <- era[i, ]
        if (e$patients == 0L) next
        add("All", sprintf("Patients treated before %d", GEP_TESTING_ERA_START_YEAR), e$group,
            sprintf("%d of %d patients", e$patients_before_era, e$patients),
            if (e$patients_before_era == 0L) "no concern" else "concern",
            "Testing appears in this cohort only from that year on, so the group's risk may mix eras; see Sensitivity_Treatment_Era")
    }
    perf <- results$performance
    for (i in seq_len(nrow(perf))) {
        p <- perf[i, ]
        if (!is.finite(p$auc)) {
            add(p$endpoint, "AUC estimable", paste(p$performance_scope, p$model, sep = ": "), "not estimable", "not assessed",
                "Insufficient weighted cases or controls")
            next
        }
        add(p$endpoint, "AUC interval excludes chance (0.5)", paste(p$performance_scope, p$model, sep = ": "),
            sprintf("AUC %.2f (%.2f-%.2f)", p$auc, p$auc_ci_lower, p$auc_ci_upper),
            if (is.finite(p$auc_ci_lower) && p$auc_ci_lower > 0.5) "no concern" else "concern",
            "Lower bound of the bootstrap interval at or below 0.5 means the data do not separate the model from chance")
        if (is.finite(p$brier_skill)) {
            add(p$endpoint, "Brier score better than no-information reference", paste(p$performance_scope, p$model, sep = ": "),
                sprintf("skill %.3f", p$brier_skill), if (p$brier_skill >= GEP_NO_GEP_MIN_BRIER_SKILL) "no concern" else "concern",
                sprintf("Brier skill below %.2f means the predicted risks add little squared-error improvement over the overall prevalence", GEP_NO_GEP_MIN_BRIER_SKILL))
        }
        if (p$bootstrap_failed_draws > 0) {
            add(p$endpoint, "Bootstrap draws", paste(p$performance_scope, p$model, sep = ": "),
                sprintf("%d valid, %d failed", p$bootstrap_valid_draws, p$bootstrap_failed_draws), "concern",
                "Failed draws lack cases or controls in the resample")
        }
    }
    for (endpoint in names(results$design_columns)) {
        cases <- perf$cases_by_horizon[perf$endpoint == endpoint & perf$performance_scope == "All patients" & grepl("all retained", perf$model)]
        columns <- results$design_columns[[endpoint]][["full"]]
        add(endpoint, "Events per model parameter (full model)", "All patients",
            sprintf("%d events / %d design columns = %.1f", cases, columns, cases / columns),
            if (cases / columns >= GEP_NO_GEP_MIN_EVENTS_PER_PARAMETER) "no concern" else "concern",
            sprintf("Fewer than %d events per design column makes an unpenalized fit unreliable; ridge shrinkage reduces but does not remove the problem", GEP_NO_GEP_MIN_EVENTS_PER_PARAMETER))
    }
    dropped <- results$predictor_screening[results$predictor_screening$status == "dropped", , drop = FALSE]
    for (i in seq_len(nrow(dropped))) {
        add("All", "Predictor retained", dropped$predictor[i], "dropped", "concern", dropped$reason[i])
    }
    describe_check_reasons(dplyr::bind_rows(rows))
}

#' Methods and interpretation for each reported quantity
#'
#' @param horizon_months Numeric horizon in months.
#' @return Data frame with method, estimand, units, direction and interpretation.
build_no_gep_specification <- function(horizon_months) {
    base <- data.frame(
        quantity = c("Observed risk by group", "Median follow-up", "Out-of-fold AUC", "Comparators", "Brier skill", "Class 2-like score AUC", "AUC difference vs T stage", "Unknown-cause-of-death sensitivity"),
        method = c(
            sprintf("Kaplan-Meier risk (MFS) or Aalen-Johansen cumulative incidence (MSS) at %d months; no value is carried past a group's last follow-up", horizon_months),
            "Reverse Kaplan-Meier from treatment (censoring treated as the event)",
            sprintf("Repeated nested cross-validation of a ridge model; censoring weights estimated inside each training fold; cumulative/dynamic AUC at %d months; patient-level bootstrap with one random partition repeat per draw", horizon_months),
            "T stage (1 to 4) and tumor diameter used directly as the risk score, with no fitting, on the same out-of-fold patients and weights",
            "1 minus model Brier score divided by the Brier score of the overall weighted prevalence",
            "Ridge logistic model of Class 2 versus Class 1 among definitive patients; out-of-fold AUC with bootstrap interval",
            "Model AUC minus T-stage-alone AUC on the same resampled patients and the same partition repeat in each bootstrap draw",
            "Observed melanoma-specific risk per group with deaths of unknown recorded cause counted as melanoma deaths instead of competing deaths"
        ),
        estimand = c(
            sprintf("Probability of the event by %.0f years", horizon_months / 12), "Typical follow-up of patients still under observation",
            "Probability that a patient with the event by the horizon has a higher score than one without",
            "Same as the model AUC", "Proportional reduction in squared error versus no information", "Probability that a Class 2 patient scores higher than a Class 1 patient",
            "Change in ranking ability relative to a single baseline measurement",
            "Probability of melanoma death by 5 years under the alternative cause-of-death assumption"
        ),
        units = c("proportion", "months", "0.5 = chance, 1 = perfect", "0.5 = chance, 1 = perfect", "0 = no better than prevalence", "0.5 = chance, 1 = perfect", "AUC difference (positive = model ranks better)", "proportion"),
        direction = c("Higher = more events", "Longer = more complete follow-up", "Higher = better ranking", "Higher = better ranking", "Higher = better", "Higher = better ranking", "Positive = model better than T stage", "Higher = more melanoma deaths"),
        interpretation = c(
            "Compare groups only where the estimate is present; groups not followed to the horizon have no estimate",
            "Short follow-up means a group's 5-year risk cannot be assessed",
            "Interval includes 0.5 when the data do not separate the model from chance",
            "A fitted model that does not exceed these has not added information beyond one baseline measurement",
            "Zero or negative: predicted risks do not improve on the overall prevalence",
            "Baseline clinical features and molecular class are related only weakly if the interval includes 0.5",
            "An interval that includes 0 means the model is not shown to rank better than T stage alone",
            "Compare with the primary risk in the same row"
        ),
        stringsAsFactors = FALSE
    )
    era <- data.frame(
        quantity = "Treatment-era sensitivity",
        method = sprintf("Observed risk per group recomputed among patients treated in %d or later (the era in which testing appears in this cohort), with a log-log 95%% interval", GEP_TESTING_ERA_START_YEAR),
        estimand = sprintf("Probability of the event by %.0f years among patients treated in the testing era", horizon_months / 12),
        units = "proportion; difference in percentage points",
        direction = "Difference = era risk minus primary risk",
        interpretation = "Small differences mean the group comparison does not depend on treatment era; a group with every patient in the era is unchanged",
        stringsAsFactors = FALSE
    )
    rbind(base, era)
}

#' Limitations for the no-GEP results
#'
#' @param results Output of `run_no_gep_analysis()`.
#' @return Data frame with `limitation` and `detail`.
build_no_gep_limitations <- function(results) {
    outcomes <- results$group_outcomes
    not_estimable <- unique(outcomes$group[!outcomes$estimable & outcomes$n > 0L])
    data.frame(
        limitation = c("Groups without a 5-year estimate", "Selection into testing", "Internal validation only", "Interval construction", "Predictors dropped"),
        detail = c(
            if (length(not_estimable) > 0) paste("Not followed to 60 months, so no observed risk:", paste(not_estimable, collapse = "; ")) else "None",
            "Testing coverage depends on treatment era and tumor size (see Baseline_By_Group and Sensitivity_Treatment_Era); group differences describe who was tested, not the effect of testing.",
            "Performance is out-of-fold within this cohort; there is no external validation.",
            "Bootstrap intervals resample patients and out-of-fold predictions; the model is not refit on each bootstrap sample.",
            paste(results$predictor_screening$predictor[results$predictor_screening$status == "dropped"], collapse = "; ")
        ),
        stringsAsFactors = FALSE
    )
}

#' Wrap text for plot titles and captions
#'
#' @param text Character scalar.
#' @param width Integer maximum line width in characters.
#' @return Text with line breaks inserted.
wrap_plot_text <- function(text, width = 95L) {
    paste(strwrap(text, width = width), collapse = "\n")
}

#' Plot observed 5-year risk for every GEP availability group
#'
#' @param group_outcomes Data frame from `build_no_gep_group_outcomes()` for both endpoints.
#' @return A `ggplot` object.
plot_no_gep_group_risk <- function(group_outcomes) {
    plot_data <- group_outcomes[group_outcomes$n >= GEP_MIN_REPORTED_GROUP_SIZE, , drop = FALSE]
    plot_data$no_events <- plot_data$events_by_horizon == 0
    plot_data$label <- factor(
        sprintf("%s (n = %d)", plot_data$group, plot_data$n),
        levels = unique(sprintf("%s (n = %d)", plot_data$group, plot_data$n)[order(match(plot_data$group, rev(GEP_GROUP_LEVELS)))])
    )
    plot_data$support <- ifelse(
        plot_data$estimable,
        sprintf("%s%d known event-free", ifelse(plot_data$no_events, "0 events; ", ""), plot_data$event_free_through_horizon),
        sprintf("no 5-year estimate; %d event-free", plot_data$event_free_through_horizon)
    )
    colors <- get_gep_class_palette(unique(as.character(plot_data$group)))
    ggplot2::ggplot(plot_data, ggplot2::aes(y = .data$label)) +
        ggplot2::geom_errorbar(
            ggplot2::aes(xmin = .data$observed_ci_lower, xmax = .data$observed_ci_upper, color = .data$group),
            width = 0.25, linewidth = 0.9, orientation = "y", na.rm = TRUE
        ) +
        ggplot2::geom_point(ggplot2::aes(x = .data$observed_risk, color = .data$group, shape = .data$no_events), size = 3.4, stroke = 1.2, na.rm = TRUE) +
        ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
        ggplot2::geom_text(ggplot2::aes(x = 1.04, label = .data$support), hjust = 0, size = 3.4, color = "grey30") +
        ggplot2::facet_wrap(~endpoint, scales = "free_y", labeller = ggplot2::as_labeller(c(MFS = "Metastasis", MSS = "Melanoma death"))) +
        ggplot2::scale_color_manual(values = stats::setNames(unname(colors), names(colors)), guide = "none") +
        ggplot2::scale_x_continuous(limits = c(0, 2.1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = scales::percent_format(accuracy = 1)) +
        ggplot2::labs(
            x = "Observed 5-year risk (95% interval)", y = NULL,
            title = "5-year risk by GEP availability group",
            subtitle = wrap_plot_text("Kaplan-Meier (metastasis) and Aalen-Johansen (melanoma death). No interval is shown when no events occurred.", 110),
            caption = wrap_plot_text("Known event-free: patients observed event-free through 60 months. Open circles: no events occurred. Groups not followed to 60 months have no estimate; groups with fewer than 3 patients are not shown.", 110)
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_line(color = "grey92"), plot.title.position = "plot", plot.caption.position = "plot")
}

#' Plot out-of-fold AUC for the models and comparators
#'
#' @param performance Data frame of model performance.
#' @return A `ggplot` object.
plot_no_gep_model_auc <- function(performance) {
    plot_data <- performance[
        (performance$performance_scope %in% c("No GEP", "Definitive Class 1/2 patients")) & is.finite(performance$auc), ,
        drop = FALSE
    ]
    plot_data$kind <- ifelse(grepl("alone", plot_data$model), "No-fit comparator", "Ridge model")
    # The axis only zooms: setting scale limits would delete intervals reaching below the
    # lower bound, so the bound follows the lowest interval instead.
    lowest <- suppressWarnings(min(plot_data$auc_ci_lower, na.rm = TRUE))
    axis_lower <- if (is.finite(lowest)) min(0.3, floor(lowest * 10) / 10) else 0.3
    plot_data$model <- factor(plot_data$model, levels = rev(unique(plot_data$model)))
    plot_data$endpoint <- factor(plot_data$endpoint, levels = c("MFS", "MSS", "Class 2 vs Class 1"),
                                 labels = c("Metastasis (no-GEP patients)", "Melanoma death (no-GEP patients)", "Class 2 vs Class 1 (definitive patients)"))
    ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$auc, y = .data$model, shape = .data$kind)) +
        ggplot2::geom_vline(xintercept = 0.5, linetype = "longdash", color = "grey50", linewidth = 0.8) +
        ggplot2::geom_errorbar(
            ggplot2::aes(xmin = .data$auc_ci_lower, xmax = .data$auc_ci_upper),
            width = 0.2, linewidth = 0.9, color = "grey25", orientation = "y", na.rm = TRUE
        ) +
        ggplot2::geom_point(size = 3.6, color = "grey10") +
        ggplot2::facet_wrap(~endpoint, ncol = 1, scales = "free_y") +
        ggplot2::scale_shape_manual(values = c("Ridge model" = 16, "No-fit comparator" = 5), name = NULL) +
        ggplot2::scale_x_continuous(breaks = seq(axis_lower, 1, by = 0.1)) +
        ggplot2::coord_cartesian(xlim = c(axis_lower, 1)) +
        ggplot2::labs(
            x = "Out-of-fold AUC at 5 years (95% bootstrap interval)", y = NULL,
            title = "Baseline features vs chance and single measurements",
            subtitle = wrap_plot_text("Dashed line: chance (0.5). Open diamonds: T stage or tumor diameter used directly, with no fitting.", 80),
            caption = wrap_plot_text("AUC: probability that a patient with the outcome scores higher than one without.", 80)
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), plot.title.position = "plot", plot.caption.position = "plot")
}

#' Blank the results of groups too small to describe without identifying a patient
#'
#' @param table Data frame with a patient-count column.
#' @param count_col Character name of the patient-count column.
#' @param keep_cols Character vector of columns kept as they are (labels and the count).
#' @return `table` with all other columns set to `NA` and a `note` for groups smaller
#'   than `GEP_MIN_REPORTED_GROUP_SIZE`.
suppress_small_groups <- function(table, count_col, keep_cols) {
    small <- table[[count_col]] < GEP_MIN_REPORTED_GROUP_SIZE
    if (!"note" %in% names(table)) {
        table$note <- ""
    }
    for (column in setdiff(names(table), c(keep_cols, "note"))) {
        table[[column]][small] <- NA
    }
    table$note[small] <- sprintf("Fewer than %d patients: results not shown", GEP_MIN_REPORTED_GROUP_SIZE)
    table
}

#' Write the no-GEP workbooks and figures
#'
#' Primary workbook (publishable): group outcomes, baseline table, performance.
#' Diagnostics workbook (part of the results package; `Record_Issues` lists patients by id).
#'
#' @param results Output of `run_no_gep_analysis()`.
#' @param data Prepared cohort with `gep_group`.
#' @param output_dir Character output directory.
#' @param prefix Character file prefix.
#' @return Invisibly, a list of written paths.
write_no_gep_outputs <- function(results, data, output_dir, prefix) {
    ensure_output_dir(output_dir)
    checks <- build_no_gep_checks(results)
    primary_path <- file.path(output_dir, paste0(prefix, "no_gep_5yr.xlsx"))
    write_gep_workbook(
        list(
            Group_Outcomes_MFS = suppress_small_groups(results$group_outcomes[results$group_outcomes$endpoint == "MFS", ], "n", c("endpoint", "group", "n")),
            Group_Outcomes_MSS = suppress_small_groups(results$group_outcomes[results$group_outcomes$endpoint == "MSS", ], "n", c("endpoint", "group", "n")),
            Baseline_By_Group = results$baseline,
            Model_Performance = results$performance,
            AUC_Difference = results$auc_differences,
            Sensitivity_Treatment_Era = suppress_small_groups(results$treatment_era_sensitivity, "patients_in_era", c("endpoint", "group", "patients", "patients_in_era", "patients_before_era")),
            Sensitivity_Cause_of_Death = suppress_small_groups(results$cause_of_death_sensitivity, "patients", c("group", "patients"))
        ),
        primary_path
    )
    id_col <- pick_sparse_level_id_col(data)
    removed <- results$removed_records
    removed$key_source <- rep(id_col %||% "row position", nrow(removed))
    removed <- removed[, c("id", "key_source", "analysis_set", "reason", "field", "current_value")]
    names(removed)[names(removed) == "id"] <- "record_key"
    diagnostics_path <- file.path(output_dir, paste0(prefix, "no_gep_diagnostics.xlsx"))
    write_gep_workbook(
        list(
            Specification_Coding = build_no_gep_specification(results$horizon_months),
            Denominators_Flow = summarize_no_gep_flow(data, results),
            Predictor_Screening = results$predictor_screening,
            Checks = checks,
            Warnings_Limitations = build_no_gep_limitations(results),
            Record_Issues = removed
        ),
        diagnostics_path
    )
    risk_path <- file.path(output_dir, paste0(prefix, "no_gep_group_risk.png"))
    ggplot2::ggsave(risk_path, plot_no_gep_group_risk(results$group_outcomes), width = 14, height = 4.8, dpi = 300, bg = "white")
    auc_path <- file.path(output_dir, paste0(prefix, "no_gep_model_auc.png"))
    ggplot2::ggsave(auc_path, plot_no_gep_model_auc(results$performance), width = 9, height = 8, dpi = 300, bg = "white")
    invisible(list(primary = primary_path, diagnostics = diagnostics_path, figures = c(risk_path, auc_path)))
}

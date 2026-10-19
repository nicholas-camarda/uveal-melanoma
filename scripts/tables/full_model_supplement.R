# Full-model supplementary table (Revision 2, Reviewer 2 comment 1)
#
# Lists every term of each primary adjusted model with its estimate, 95% CI and
# p-value, using the same formatted gtsummary tables the pipeline writes for
# each model, so the supplement cannot disagree with the per-model outputs.

#' Primary adjusted models reported in the full-model supplement
#'
#' Each row names the model, its effect measure, and the path to its result
#' inside the per-cohort objective results returned by `run_my_analysis()`.
#'
#' @return Tibble with `model_label`, `effect_measure`, `objective`, `result`,
#'   `model_field`, and `table_field`.
full_model_supplement_specs <- function() {
    tibble::tribble(
        ~model_label, ~effect_measure, ~objective, ~result, ~model_field, ~table_field,
        "Local recurrence (Cox, full follow-up)", "HR", "objective_1", "recurrence_time_to_event", "cox_model", "cox_table",
        "Metastatic progression (Cox, full follow-up)", "HR", "objective_1", "mets_time_to_event", "cox_model", "cox_table",
        "Overall survival (Cox, 5-year)", "HR", "objective_1", "os_5yr_capped", "model", "table",
        "Progression-free survival (Cox, 5-year)", "HR", "objective_1", "pfs_5yr_capped", "model", "table",
        "VA change, logMAR (linear)", "MD", "objective_2", "vision_changes", "regression_model", "regression_table",
        "Radiation retinopathy (logistic)", "OR", "objective_2", "retinopathy_analysis", "model", "table",
        "Neovascular glaucoma (logistic)", "OR", "objective_2", "nvg_analysis", "model", "table",
        "Serous retinal detachment (logistic)", "OR", "objective_2", "srd_analysis", "model", "table"
    )
}

#' Count fitted observations and events for a supplement model
#'
#' @param model Fitted `coxph`, `glm`, or `lm` object.
#' @return Named list with `fitted_n` and `events` (`NA` for linear models).
full_model_fit_counts <- function(model) {
    fitted_n <- nrow(stats::model.frame(model))
    events <- if (inherits(model, "coxph")) {
        sum(model$y[, ncol(model$y)])
    } else if (inherits(model, "glm")) {
        sum(stats::model.response(stats::model.frame(model)) == 1)
    } else {
        NA_real_
    }
    list(fitted_n = as.integer(fitted_n), events = as.integer(events))
}

#' Format a p-value with the manuscript's precision rule
#'
#' Two decimals at or above 0.01, three decimals from 0.001 to below 0.01, and
#' "<0.001" below that.
#'
#' @param p Numeric vector of p-values.
#' @return Character vector; `NA` stays `NA`.
format_manuscript_pvalue <- function(p) {
    dplyr::case_when(
        is.na(p) ~ NA_character_,
        p < 0.001 ~ "<0.001",
        p < 0.01 ~ formatC(p, format = "f", digits = 3),
        TRUE ~ formatC(p, format = "f", digits = 2)
    )
}

#' Build the full-model supplementary table from objective results
#'
#' A model that was not fitted (for example the restricted-cohort neovascular
#' glaucoma model, which had too few events) appears as one row saying so.
#'
#' @param results_by_cohort Named list (display cohort label -> results of
#'   `run_my_analysis()` for that cohort, containing `objective_1` and
#'   `objective_2`).
#' @return Tibble with cohort, model, effect measure, fitted N, events, term,
#'   estimate, 95% CI, and p-value columns.
build_full_model_supplement <- function(results_by_cohort) {
    specs <- full_model_supplement_specs()
    rows <- list()
    for (cohort_label in names(results_by_cohort)) {
        cohort_results <- results_by_cohort[[cohort_label]]
        for (i in seq_len(nrow(specs))) {
            spec <- specs[i, ]
            result <- cohort_results[[spec$objective]][[spec$result]]
            if (is.null(result)) {
                stop(sprintf("Missing %s$%s for %s.", spec$objective, spec$result, cohort_label), call. = FALSE)
            }
            model <- result[[spec$model_field]]
            if (is.null(model)) {
                rows[[length(rows) + 1]] <- tibble::tibble(
                    cohort = cohort_label, model = spec$model_label, effect_measure = spec$effect_measure,
                    fitted_n = NA_integer_, events = NA_integer_,
                    term = "Model not fitted (insufficient events)", estimate = NA_character_,
                    ci_95 = NA_character_, p_value = NA_character_
                )
                next
            }
            tbl <- result[[spec$table_field]]
            if (!inherits(tbl, "gtsummary")) {
                stop(sprintf("Missing formatted model table for %s, %s.", cohort_label, spec$model_label), call. = FALSE)
            }
            counts <- full_model_fit_counts(model)
            body <- tbl$table_body
            formatted <- gtsummary::as_tibble(tbl, col_labels = FALSE)
            if (nrow(formatted) != nrow(body) || !all(c("estimate", "conf.low", "p.value") %in% names(formatted))) {
                stop(sprintf("Unexpected model table layout for %s, %s.", cohort_label, spec$model_label), call. = FALSE)
            }
            rows[[length(rows) + 1]] <- tibble::tibble(
                cohort = cohort_label,
                model = spec$model_label,
                effect_measure = spec$effect_measure,
                fitted_n = counts$fitted_n,
                events = counts$events,
                term = ifelse(body$row_type == "label", body$label, paste0("    ", body$label)),
                estimate = ifelse(body$reference_row %in% TRUE, "Reference", formatted$estimate),
                ci_95 = formatted$conf.low,
                p_value = format_manuscript_pvalue(body$p.value)
            )
        }
    }
    dplyr::bind_rows(rows)
}

#' Write the full-model supplementary table
#'
#' @param results_by_cohort Named list as in `build_full_model_supplement()`.
#' @param output_dir Directory for `full_model_supplement.xlsx`.
#' @return Path to the written workbook, invisibly.
write_full_model_supplement <- function(results_by_cohort, output_dir = MERGED_TABLES_DIR) {
    supplement <- build_full_model_supplement(results_by_cohort)
    notes <- tibble::tibble(note = c(
        "Each block is one primary adjusted model; terms, estimates, 95% CIs, and p-values are taken from the same formatted table the pipeline writes for that model.",
        "HR, hazard ratio (GKSRS vs PBT and per covariate level); MD, mean difference in logMAR change (positive = less VA loss); OR, odds ratio.",
        "Cox models for local recurrence and metastatic progression use full follow-up; overall and progression-free survival models are administratively censored at 60 months.",
        "Covariates: age at diagnosis (continuous), sex, and tumor location (sparse location categories excluded); the VA model also includes baseline logMAR VA.",
        "fitted_n is the number of patients in the fitted model; events is the number of outcome events in that model (not applicable to the linear model). Both exclude patients in sparse tumor-location categories, so they can be slightly smaller than the cohort counts in Tables 2 and 3.",
        "p-values: two decimals at or above 0.01, three decimals below 0.01, <0.001 below 0.001. A factor's p-value is shown on its level row."
    ))
    path <- file.path(output_dir, "full_model_supplement.xlsx")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write_readable_xlsx(list(Full_models = supplement, Notes = notes), path)
    invisible(path)
}

# Survival estimators shared by the exploratory no-GEP report, the MFS sensitivity
# summary, and the Objective 4 curve figures.

#' Prepare the Shared Incident Post-Treatment MFS KM Dataset
#'
#' Applies the centralized treatment-time eligibility rule and retains the
#' full-follow-up MFS time and event fields used by all Objective 4 KM paths.
#' Raw metastasis columns remain in the returned data for source/audit use but
#' are never used to define the incident KM risk set.
#'
#' @param data Data frame containing the canonical Objective 4 MFS fields.
#' @return Filtered data frame with canonical analysis time and event columns.
prepare_incident_mfs_km_data <- function(data) {
    required_cols <- c(
        "mets_free_at_baseline",
        "tt_mets_months_analysis",
        "mets_event_analysis"
    )
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0L) {
        stop(
            sprintf(
                "Incident MFS KM data is missing required columns: %s",
                paste(missing_cols, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    data %>%
        dplyr::mutate(
            .mfs_time = suppressWarnings(as.numeric(.data$tt_mets_months_analysis)),
            .mfs_event = suppressWarnings(as.integer(.data$mets_event_analysis))
        ) %>%
        dplyr::filter(
            !is.na(.data$mets_free_at_baseline),
            .data$mets_free_at_baseline,
            is.finite(.data$.mfs_time),
            .data$.mfs_time >= 0,
            !is.na(.data$.mfs_event),
            .data$.mfs_event %in% c(0L, 1L)
        ) %>%
        dplyr::mutate(
            tt_mets_months_analysis = .data$.mfs_time,
            mets_event_analysis = .data$.mfs_event
        ) %>%
        dplyr::select(-dplyr::all_of(c(".mfs_time", ".mfs_event")))
}

#' Estimate Kaplan-Meier MFS at a Horizon
#'
#' Calculates a censoring-aware Kaplan-Meier estimate for metastasis-free
#' survival at a requested month horizon and converts that estimate to an
#' observed risk and pseudo-event count on the original sample size scale.
#'
#' @param data Data frame containing the time and event columns.
#' @param timepoint_months Numeric horizon in months.
#' @param time_var Character name of the MFS follow-up column.
#' @param event_var Character name of the metastasis event column.
#'
#' @return Named list with Kaplan-Meier survival, risk, confidence bounds,
#'   pseudo-event counts, and raw event bookkeeping.
estimate_mfs_km_at_horizon <- function(data,
                                       timepoint_months,
                                       time_var = "tt_mets_months",
                                       event_var = "mets_event") {
    valid_data <- data %>%
        dplyr::filter(
            !is.na(.data[[time_var]]),
            !is.na(.data[[event_var]])
        )

    raw_events_by_horizon <- sum(
        valid_data[[event_var]] == 1 & valid_data[[time_var]] <= timepoint_months,
        na.rm = TRUE
    )

    if (nrow(valid_data) == 0) {
        return(list(
            n = 0L,
            survival = NA_real_,
            survival_ci_lower = NA_real_,
            survival_ci_upper = NA_real_,
            risk = NA_real_,
            risk_ci_lower = NA_real_,
            risk_ci_upper = NA_real_,
            observed_events = NA_real_,
            observed_events_ci_lower = NA_real_,
            observed_events_ci_upper = NA_real_,
            raw_events_by_horizon = raw_events_by_horizon
        ))
    }

    km_fit <- survival::survfit(
        survival::Surv(valid_data[[time_var]], valid_data[[event_var]]) ~ 1,
        data = valid_data
    )
    km_summary <- summary(km_fit, times = timepoint_months, extend = TRUE)

    km_survival <- km_summary$surv[1] %||% NA_real_
    km_survival_lower <- km_summary$lower[1] %||% NA_real_
    km_survival_upper <- km_summary$upper[1] %||% NA_real_

    km_risk <- if (!is.na(km_survival)) 1 - km_survival else NA_real_
    km_risk_lower <- if (!is.na(km_survival_upper)) max(0, 1 - km_survival_upper) else NA_real_
    km_risk_upper <- if (!is.na(km_survival_lower)) min(1, 1 - km_survival_lower) else NA_real_

    list(
        n = nrow(valid_data),
        survival = km_survival,
        survival_ci_lower = km_survival_lower,
        survival_ci_upper = km_survival_upper,
        risk = km_risk,
        risk_ci_lower = km_risk_lower,
        risk_ci_upper = km_risk_upper,
        observed_events = if (!is.na(km_risk)) nrow(valid_data) * km_risk else NA_real_,
        observed_events_ci_lower = if (!is.na(km_risk_lower)) nrow(valid_data) * km_risk_lower else NA_real_,
        observed_events_ci_upper = if (!is.na(km_risk_upper)) nrow(valid_data) * km_risk_upper else NA_real_,
        raw_events_by_horizon = raw_events_by_horizon
    )
}

#' Estimate MSS Cumulative Incidence at a Horizon
#'
#' Calculates a censoring-aware Aalen-Johansen cumulative incidence estimate for
#' melanoma-specific death at a requested month horizon and converts that
#' estimate to an observed pseudo-event count on the original sample size scale.
#'
#' @param data Data frame containing the time and event columns.
#' @param timepoint_months Numeric horizon in months.
#' @param time_var Character name of the MSS follow-up column.
#' @param melanoma_event_var Character name of the melanoma-death indicator.
#' @param competing_event_var Character name of the competing-death indicator.
#'
#' @return Named list with cumulative incidence, pseudo-event counts, and raw
#'   event bookkeeping.
estimate_mss_cif_at_horizon <- function(data,
                                        timepoint_months,
                                        time_var = "tt_death_months",
                                        melanoma_event_var = "melanoma_death_event",
                                        competing_event_var = "competing_death_event") {
    valid_data <- data %>%
        dplyr::mutate(
            .time_value = suppressWarnings(as.numeric(.data[[time_var]])),
            .melanoma_value = suppressWarnings(as.integer(.data[[melanoma_event_var]])),
            .competing_value = suppressWarnings(as.integer(.data[[competing_event_var]]))
        ) %>%
        dplyr::filter(
            is.finite(.data$.time_value),
            !is.na(.data$.melanoma_value),
            !is.na(.data$.competing_value)
        )

    raw_events_by_horizon <- sum(
        valid_data$.melanoma_value == 1 & valid_data$.time_value <= timepoint_months,
        na.rm = TRUE
    )

    if (nrow(valid_data) == 0) {
        return(list(
            n = 0L,
            cif = NA_real_,
            observed_events = NA_real_,
            raw_events_by_horizon = raw_events_by_horizon,
            observed_method = "aalen_johansen_cif_at_horizon"
        ))
    }

    event_status <- dplyr::case_when(
        valid_data$.melanoma_value == 1 ~ 1L,
        valid_data$.competing_value == 1 ~ 2L,
        TRUE ~ 0L
    )

    if (!any(event_status == 1L)) {
        return(list(
            n = nrow(valid_data),
            cif = 0,
            observed_events = 0,
            raw_events_by_horizon = raw_events_by_horizon,
            observed_method = "aalen_johansen_cif_at_horizon"
        ))
    }

    cif_fit <- cmprsk::cuminc(
        ftime = valid_data$.time_value,
        fstatus = event_status,
        cencode = 0
    )

    curve_names <- names(cif_fit)
    melanoma_curve_name <- curve_names[grepl("(^1$| 1$)", curve_names)][1]
    cif_curve <- if (!is.na(melanoma_curve_name)) cif_fit[[melanoma_curve_name]] else NULL
    cif_value <- if (!is.null(cif_curve) && length(cif_curve$time) > 0) {
        eligible_idx <- which(cif_curve$time <= timepoint_months)
        if (length(eligible_idx) == 0) {
            0
        } else {
            as.numeric(cif_curve$est[max(eligible_idx)])
        }
    } else {
        NA_real_
    }

    list(
        n = nrow(valid_data),
        cif = cif_value,
        observed_events = if (!is.na(cif_value)) nrow(valid_data) * cif_value else NA_real_,
        raw_events_by_horizon = raw_events_by_horizon,
        observed_method = "aalen_johansen_cif_at_horizon"
    )
}

#' Summarize fixed-horizon out-of-fold performance by target population
#'
#' Uses one keyed OOF prediction table and its unchanged outer-training-fold
#' IPCW weights to report the overall complete-predictor validation population
#' and the prespecified no-GEP target population. The two operational no-GEP
#' groups are counted descriptively; they are not promoted to separate AUC
#' scopes.
#'
#' @param oof_predictions Row-level OOF predictions containing `stable_id`,
#'   `repeat_id`, `prediction`, `horizon_event`, and `ipcw_weight`.
#' @param group_var Column identifying exploratory GEP groups.
#' @return One row per repeat and performance scope with support counts, method
#'   identifiers, and IPCW-weighted AUC, Brier, and calibration results.
summarize_scoped_ipcw_oof_performance <- function(
    oof_predictions,
    group_var = "exploratory_gep_group"
) {
    required <- c(
        "stable_id", "repeat_id", "prediction", "horizon_event",
        "ipcw_weight", group_var
    )
    if (!is.data.frame(oof_predictions) || any(!required %in% names(oof_predictions))) {
        stop(
            "OOF predictions must contain stable keys, repeat IDs, predictions, horizon outcomes, weights, and the declared group.",
            call. = FALSE
        )
    }
    if (anyNA(oof_predictions$stable_id) ||
        anyDuplicated(oof_predictions[c("repeat_id", "stable_id")])) {
        stop("OOF predictions must be uniquely keyed by repeat_id and stable_id.", call. = FALSE)
    }
    if (anyNA(oof_predictions$repeat_id) || anyNA(oof_predictions[[group_var]])) {
        stop("OOF repeat IDs and performance-scope groups must be observed.", call. = FALSE)
    }

    no_gep_levels <- c("GEP Failed/Indeterminate", "GEP Not Tested")
    repeat_values <- sort(unique(oof_predictions$repeat_id))
    scopes <- list(
        Overall = rep(TRUE, nrow(oof_predictions)),
        `No GEP` = as.character(oof_predictions[[group_var]]) %in% no_gep_levels
    )

    purrr::imap_dfr(scopes, function(in_scope, scope_name) {
        scope_data <- oof_predictions[in_scope, , drop = FALSE]
        purrr::map_dfr(repeat_values, function(repeat_value) {
            repeat_data <- scope_data[scope_data$repeat_id == repeat_value, , drop = FALSE]
            auc <- calculate_ipcw_auc(
                repeat_data$horizon_event,
                repeat_data$prediction,
                repeat_data$ipcw_weight
            )
            brier <- calculate_ipcw_brier(
                repeat_data$horizon_event,
                repeat_data$prediction,
                repeat_data$ipcw_weight
            )
            calibration <- summarize_ipcw_calibration(
                repeat_data$horizon_event,
                pmin(pmax(repeat_data$prediction, 1e-6), 1 - 1e-6),
                repeat_data$ipcw_weight
            )
            group_value <- as.character(repeat_data[[group_var]])
            positive_weight <- repeat_data$ipcw_weight > 0

            tibble::tibble(
                repeat_id = as.integer(repeat_value),
                performance_scope = scope_name,
                scope_n = as.integer(length(unique(repeat_data$stable_id))),
                positive_weight_n = as.integer(auc$n_positive_weight),
                case_n = as.integer(sum(
                    positive_weight & repeat_data$horizon_event == 1L,
                    na.rm = TRUE
                )),
                control_n = as.integer(sum(
                    positive_weight & repeat_data$horizon_event == 0L,
                    na.rm = TRUE
                )),
                weighted_cases = as.numeric(auc$weighted_cases),
                weighted_controls = as.numeric(auc$weighted_controls),
                failed_indeterminate_n = as.integer(length(unique(
                    repeat_data$stable_id[group_value == "GEP Failed/Indeterminate"]
                ))),
                not_tested_n = as.integer(length(unique(
                    repeat_data$stable_id[group_value == "GEP Not Tested"]
                ))),
                auc_status = auc$status,
                cv_auc = auc$auc,
                brier_status = brier$status,
                cv_brier = brier$brier,
                calibration_status = calibration$status,
                cv_calibration_intercept = calibration$intercept,
                cv_calibration_slope = calibration$slope,
                evaluation_method = "outer-training-fold IPCW weighted OOF AUC/Brier/calibration",
                censoring_method = "outer-training-fold Kaplan-Meier IPCW"
            )
        })
    })
}

# Definitive-GEP 5-year validation
#
# Patients with a definitive DecisionDx Class 1/2 call and PRAME reported carry
# one of four imported predictions (Class x PRAME). The honest check is therefore
# observed 5-year risk in each of the four groups against the imported value,
# plus tie-aware discrimination of that four-level predictor. No model is fit and
# no value is extrapolated beyond the last observed follow-up.

#' Assign the Class x PRAME validation group
#'
#' @param data Data frame with `gep_group` and `prame_status`.
#' @return `data` with a factor column `validation_group` (levels
#'   `GEP_VALIDATION_GROUP_LEVELS`); rows outside the four groups get `NA`.
assign_gep_validation_group <- function(data) {
    detail <- classify_gep_detail(data$gep_group, data$prame_status)
    data$validation_group <- factor(
        ifelse(detail %in% GEP_VALIDATION_GROUP_LEVELS, detail, NA_character_),
        levels = GEP_VALIDATION_GROUP_LEVELS
    )
    data
}

#' Describe a record's field values as text for Record_Issues
#'
#' @param rows Data frame of records.
#' @param columns Character vector of column names present in `rows`.
#' @return Character vector with one string per row, `column = value; ...`, with
#'   missing values written as `missing`.
format_record_values <- function(rows, columns) {
    vapply(seq_len(nrow(rows)), function(i) {
        parts <- vapply(columns, function(column) {
            value <- rows[[column]][i]
            paste0(column, " = ", if (is.na(value)) "missing" else as.character(value))
        }, character(1))
        paste(parts, collapse = "; ")
    }, character(1))
}

#' Percentile bootstrap interval with a valid-draw requirement
#'
#' One rule for every bootstrap in the GEP analyses: draws that are not finite count as
#' failed, and the 95% percentile interval is reported only when at least
#' `GEP_MIN_VALID_BOOTSTRAP_SHARE` of the draws are valid.
#'
#' @param draws Numeric vector of bootstrap statistics (`NA` for failed draws).
#' @param n_boot Integer number of draws requested; `0` returns no interval.
#' @return List: `ci_lower`, `ci_upper`, `valid_draws`, `failed_draws`.
bootstrap_percentile_interval <- function(draws, n_boot) {
    valid <- sum(is.finite(draws))
    bounds <- c(NA_real_, NA_real_)
    if (n_boot > 0L && valid >= GEP_MIN_VALID_BOOTSTRAP_SHARE * n_boot) {
        bounds <- stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    }
    list(ci_lower = bounds[[1]], ci_upper = bounds[[2]], valid_draws = valid, failed_draws = n_boot - valid)
}

#' Prepare the analysis set for one validation endpoint
#'
#' Keeps the patients flagged eligible for the endpoint, assigns the validation
#' group, and creates uniform `time_months`, `status` (0 censored, 1 event of
#' interest, 2 competing death; MFS never uses 2), and `published_risk` columns.
#'
#' @param data Cohort data frame from the analytic dataset.
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @return Data frame with columns `validation_group`, `time_months`, `status`,
#'   and `published_risk` (1 minus the imported 5-year survival probability).
prepare_gep_validation_endpoint <- function(data, endpoint = c("MFS", "MSS")) {
    endpoint <- match.arg(endpoint)
    spec <- switch(
        endpoint,
        MFS = list(eligible = "mfs_analysis_eligible", time = "tt_mets_months_analysis",
                   status = "mfs_event_type", expected = "biopsy1_gep_mfs"),
        MSS = list(eligible = "mss_analysis_eligible", time = "tt_death_months",
                   status = "mss_event_type", expected = "biopsy1_gep_mss")
    )
    required <- c(spec$eligible, spec$time, spec$status, spec$expected, "gep_group", "prame_status", "gep_test_timing")
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0L) {
        stop(sprintf("%s validation requires columns: %s", endpoint, paste(missing_cols, collapse = ", ")), call. = FALSE)
    }

    prepared <- data %>%
        dplyr::filter(.data[[spec$eligible]] %in% TRUE) %>%
        assign_gep_validation_group() %>%
        dplyr::mutate(
            time_months = as.numeric(.data[[spec$time]]),
            status = as.integer(.data[[spec$status]]),
            published_risk = 1 - as.numeric(.data[[spec$expected]])
        )
    # An eligible patient must be analysable; a silent drop would leave a patient the
    # diagnostics workbook cannot account for.
    unusable <- is.na(prepared$validation_group) | !is.finite(prepared$time_months) | is.na(prepared$status)
    if (any(unusable)) {
        stop(sprintf(
            "%d patient(s) eligible for the %s endpoint have no validation group, follow-up time or status.",
            sum(unusable), endpoint
        ), call. = FALSE)
    }
    if (identical(endpoint, "MSS")) {
        if (!"mss_event_type_unknown_cause_melanoma" %in% names(prepared)) {
            stop("MSS validation requires column: mss_event_type_unknown_cause_melanoma", call. = FALSE)
        }
        prepared$status_unknown_melanoma <- as.integer(prepared$mss_event_type_unknown_cause_melanoma)
    }
    prepared
}

#' Kaplan-Meier risk at a horizon without extrapolation
#'
#' Reports the censoring-aware metastasis risk at `horizon_months`. The estimate
#' is `NA` when the group is not followed to the horizon, rather than carrying the
#' last value forward. The interval is `NA` when no event occurred by the horizon
#' because the Kaplan-Meier interval is then degenerate.
#'
#' @param time_months Numeric follow-up times in months.
#' @param status Integer status, 1 = event, anything else censored.
#' @param horizon_months Numeric horizon in months.
#' @return Named list: `n`, `events`, `censored_before`, `event_free_through`,
#'   `max_followup_months`, `risk`, `ci_lower`, `ci_upper`, `estimable`.
estimate_km_risk_at_horizon <- function(time_months, status, horizon_months) {
    n <- length(time_months)
    event <- status == 1L
    result <- list(
        n = n,
        events = sum(event & time_months <= horizon_months),
        censored_before = sum(!event & time_months < horizon_months),
        competing = 0L,
        event_free_through = NA_integer_,
        max_followup_months = if (n > 0L) max(time_months) else NA_real_,
        risk = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        estimable = FALSE
    )
    result$event_free_through <- n - result$events - result$censored_before
    if (n == 0L || max(time_months) < horizon_months) {
        return(result)
    }
    fit <- survival::survfit(survival::Surv(time_months, event) ~ 1, conf.type = "log-log")
    at_horizon <- summary(fit, times = horizon_months)
    result$risk <- 1 - at_horizon$surv
    result$estimable <- TRUE
    if (result$events > 0L && length(at_horizon$lower) == 1L && is.finite(at_horizon$lower)) {
        result$ci_lower <- 1 - at_horizon$upper
        result$ci_upper <- 1 - at_horizon$lower
    }
    result
}

#' Aalen-Johansen cumulative incidence of the event of interest at a horizon
#'
#' Competing deaths remain in the risk set as competing events. The estimate is
#' `NA` when the group is not followed to the horizon. The interval uses the
#' log-log transform of the cumulative-incidence variance and is `NA` when no
#' event of interest occurred by the horizon.
#'
#' @param time_months Numeric follow-up times in months.
#' @param status Integer status: 0 censored, 1 event of interest, 2 competing.
#' @param horizon_months Numeric horizon in months.
#' @return Same list shape as `estimate_km_risk_at_horizon()`, with `competing`
#'   counting competing deaths by the horizon.
estimate_cif_risk_at_horizon <- function(time_months, status, horizon_months) {
    n <- length(time_months)
    result <- list(
        n = n,
        events = sum(status == 1L & time_months <= horizon_months),
        censored_before = sum(status == 0L & time_months < horizon_months),
        competing = sum(status == 2L & time_months <= horizon_months),
        event_free_through = NA_integer_,
        max_followup_months = if (n > 0L) max(time_months) else NA_real_,
        risk = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        estimable = FALSE
    )
    result$event_free_through <- n - result$events - result$censored_before - result$competing
    if (n == 0L || max(time_months) < horizon_months) {
        return(result)
    }
    result$estimable <- TRUE
    if (result$events == 0L) {
        result$risk <- 0
        return(result)
    }
    fit <- cmprsk::cuminc(ftime = time_months, fstatus = status, cencode = 0)
    point <- cmprsk::timepoints(fit, times = horizon_months)
    cause_row <- grep("1$", rownames(point$est))[1]
    cif <- as.numeric(point$est[cause_row, 1])
    cif_var <- as.numeric(point$var[cause_row, 1])
    result$risk <- cif
    if (is.finite(cif_var) && cif_var > 0 && cif > 0 && cif < 1) {
        log_log <- log(-log(1 - cif))
        se_log_log <- sqrt(cif_var) / ((1 - cif) * (-log(1 - cif)))
        z <- stats::qnorm(0.975)
        result$ci_lower <- 1 - exp(-exp(log_log - z * se_log_log))
        result$ci_upper <- 1 - exp(-exp(log_log + z * se_log_log))
    }
    result
}

#' Estimate horizon risk for one set of patients
#'
#' @param time_months Numeric follow-up times in months.
#' @param status Integer status vector (see the two estimators).
#' @param horizon_months Numeric horizon in months.
#' @param endpoint Character, `"MFS"` (Kaplan-Meier) or `"MSS"` (Aalen-Johansen).
#' @return List from `estimate_km_risk_at_horizon()` or `estimate_cif_risk_at_horizon()`.
estimate_endpoint_risk <- function(time_months, status, horizon_months, endpoint) {
    if (identical(endpoint, "MFS")) {
        estimate_km_risk_at_horizon(time_months, status, horizon_months)
    } else {
        estimate_cif_risk_at_horizon(time_months, status, horizon_months)
    }
}

#' Case-mix-standardised overall observed risk
#'
#' Combines the four group-level observed risks weighted by group size, so the
#' overall observed risk does not depend on a pooled Kaplan-Meier curve that is
#' dominated by whichever group is followed longest (Class 1 patients are
#' censored much earlier than Class 2 patients). No interval is reported: the
#' group-level intervals carry the uncertainty, and resampled groups are often not
#' followed to the horizon.
#'
#' @param group_estimates List of estimates from `estimate_endpoint_risk()`, one per group.
#' @return Numeric overall risk, or `NA` when any group is not followed to the horizon.
standardised_overall_risk <- function(group_estimates) {
    if (!all(vapply(group_estimates, function(e) isTRUE(e$estimable), logical(1)))) {
        return(NA_real_)
    }
    sizes <- vapply(group_estimates, function(e) e$n, numeric(1))
    sum(sizes * vapply(group_estimates, function(e) e$risk, numeric(1))) / sum(sizes)
}

#' Build the observed-versus-imported group table for one endpoint
#'
#' One row per Class x PRAME group plus a case-mix-standardised overall row.
#' `observed_risk` is Kaplan-Meier (MFS) or Aalen-Johansen (MSS) risk at the
#' horizon; `imported_risk` is 1 minus the imported 5-year probability.
#'
#' @param prepared Data frame from `prepare_gep_validation_endpoint()`.
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @param horizon_months Numeric horizon in months.
#' @return Data frame with one row per group and an `Overall` row, with counts at
#'   the horizon, observed and imported risk, their difference in percentage
#'   points, the observed-to-imported ratio, and an `estimable` flag with a note.
build_gep_group_table <- function(prepared, endpoint, horizon_months) {
    group_estimates <- list()
    group_rows <- lapply(levels(prepared$validation_group), function(group_name) {
        g <- prepared[prepared$validation_group == group_name, , drop = FALSE]
        est <- estimate_endpoint_risk(g$time_months, g$status, horizon_months, endpoint)
        group_estimates[[group_name]] <<- est
        imported_values <- unique(round(g$published_risk, 6))
        imported <- if (nrow(g) > 0L) mean(g$published_risk) else NA_real_
        note <- dplyr::case_when(
            nrow(g) == 0L ~ "No patients in this group",
            !est$estimable ~ sprintf("Not followed to %d months (maximum %.0f months)", horizon_months, est$max_followup_months),
            est$events == 0L ~ "No events by the horizon; interval not estimable",
            TRUE ~ ""
        )
        if (length(imported_values) > 1L) {
            note <- trimws(paste(note, "Imported value varies within group"))
        }
        data.frame(
            group = group_name, n = est$n, events_by_horizon = est$events,
            competing_deaths_by_horizon = est$competing,
            censored_before_horizon = est$censored_before,
            event_free_through_horizon = est$event_free_through,
            max_followup_months = est$max_followup_months,
            observed_risk = est$risk, observed_ci_lower = est$ci_lower, observed_ci_upper = est$ci_upper,
            imported_risk = imported, estimable = est$estimable, note = note,
            stringsAsFactors = FALSE
        )
    })
    overall_risk <- standardised_overall_risk(group_estimates)
    groups <- dplyr::bind_rows(group_rows)
    overall_row <- data.frame(
        group = "Overall (case-mix standardised)", n = nrow(prepared),
        events_by_horizon = sum(groups$events_by_horizon),
        competing_deaths_by_horizon = sum(groups$competing_deaths_by_horizon),
        censored_before_horizon = sum(groups$censored_before_horizon),
        event_free_through_horizon = sum(groups$event_free_through_horizon),
        max_followup_months = NA_real_,
        observed_risk = overall_risk, observed_ci_lower = NA_real_, observed_ci_upper = NA_real_,
        imported_risk = mean(prepared$published_risk), estimable = is.finite(overall_risk),
        note = if (is.finite(overall_risk)) {
            "Group-size-weighted mean of group risks; no interval (see group rows)"
        } else {
            "At least one group is not followed to the horizon"
        },
        stringsAsFactors = FALSE
    )
    dplyr::bind_rows(groups, overall_row) %>%
        dplyr::mutate(
            difference_pct_points = 100 * (.data$observed_risk - .data$imported_risk),
            observed_to_imported_ratio = dplyr::if_else(.data$imported_risk > 0, .data$observed_risk / .data$imported_risk, NA_real_)
        )
}

#' Tie-aware Harrell's C at the horizon with a bootstrap interval
#'
#' Follow-up is administratively censored at the horizon. Pairs tied on the
#' imported risk count one half, which matters here because the predictor has
#' only four values. The interval is the percentile interval over patient-level
#' bootstrap draws; draws without an event or with a single risk value are
#' counted as failed.
#'
#' @param time_months Numeric follow-up times in months.
#' @param status Integer status, 1 = event.
#' @param risk Numeric imported risk (higher = worse).
#' @param horizon_months Numeric horizon in months.
#' @param n_boot Integer bootstrap replicates; `0` returns the point estimate only.
#' @param seed Integer seed, applied locally (unused when `n_boot` is 0).
#' @return List: `estimable`, `reason`, `c_index`, `ci_lower`, `ci_upper`,
#'   `events`, `concordant_pairs`, `discordant_pairs`, `tied_risk_pairs`,
#'   `tied_pair_share`, `distinct_risk_values`, `valid_draws`, `failed_draws`.
calculate_gep_harrell_c <- function(time_months, status, risk, horizon_months, n_boot, seed) {
    truncated_time <- pmin(time_months, horizon_months)
    truncated_event <- as.integer(status == 1L & time_months <= horizon_months)
    out <- list(
        estimable = FALSE, reason = NA_character_, c_index = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        events = sum(truncated_event), concordant_pairs = NA_real_, discordant_pairs = NA_real_,
        tied_risk_pairs = NA_real_, tied_pair_share = NA_real_,
        distinct_risk_values = length(unique(risk)), valid_draws = 0L, failed_draws = 0L
    )
    if (length(time_months) < GEP_MIN_SAMPLE_SIZE) {
        out$reason <- sprintf("Fewer than %d patients", GEP_MIN_SAMPLE_SIZE)
        return(out)
    }
    if (out$events < GEP_MIN_DISCRIMINATION_EVENTS) {
        out$reason <- sprintf("Fewer than %d events by %d months", GEP_MIN_DISCRIMINATION_EVENTS, horizon_months)
        return(out)
    }
    if (out$distinct_risk_values < 2L) {
        out$reason <- "Imported risk does not vary"
        return(out)
    }
    concordance_at <- function(idx) {
        if (sum(truncated_event[idx]) == 0L || length(unique(risk[idx])) < 2L) {
            return(NULL)
        }
        survival::concordance(
            survival::Surv(truncated_time[idx], truncated_event[idx]) ~ risk[idx],
            reverse = TRUE
        )
    }
    fit <- concordance_at(seq_along(risk))
    counts <- fit$count
    comparable <- counts[["concordant"]] + counts[["discordant"]] + counts[["tied.x"]]
    out$estimable <- TRUE
    out$c_index <- as.numeric(fit$concordance)
    out$concordant_pairs <- counts[["concordant"]]
    out$discordant_pairs <- counts[["discordant"]]
    out$tied_risk_pairs <- counts[["tied.x"]]
    out$tied_pair_share <- counts[["tied.x"]] / comparable
    draws <- numeric(0)
    if (n_boot > 0L) {
        draws <- withr::with_seed(seed, {
            vapply(seq_len(n_boot), function(i) {
                resampled <- concordance_at(sample.int(length(risk), replace = TRUE))
                if (is.null(resampled)) NA_real_ else as.numeric(resampled$concordance)
            }, numeric(1))
        })
    }
    interval <- bootstrap_percentile_interval(draws, n_boot)
    out[names(interval)] <- interval
    out
}

#' Competing-risk time-dependent AUC for melanoma-specific death at the horizon
#'
#' Cumulative/dynamic AUC with inverse-probability-of-censoring weights, where
#' non-melanoma deaths before the horizon are controls (`timeROC` AUC_2). The
#' interval is the percentile interval over patient-level bootstrap draws.
#'
#' @param time_months Numeric follow-up times in months.
#' @param status Integer status: 0 censored, 1 melanoma death, 2 other death.
#' @param risk Numeric imported risk (higher = worse).
#' @param horizon_months Numeric horizon in months.
#' @param n_boot Integer bootstrap replicates; `0` returns the point estimate only.
#' @param seed Integer seed, applied locally (unused when `n_boot` is 0).
#' @return List: `estimable`, `reason`, `auc`, `ci_lower`, `ci_upper`, `events`,
#'   `valid_draws`, `failed_draws`.
calculate_gep_competing_risk_auc <- function(time_months, status, risk, horizon_months, n_boot, seed) {
    out <- list(
        estimable = FALSE, reason = NA_character_, auc = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        events = sum(status == 1L & time_months <= horizon_months), valid_draws = 0L, failed_draws = 0L
    )
    if (length(time_months) < GEP_MIN_SAMPLE_SIZE) {
        out$reason <- sprintf("Fewer than %d patients", GEP_MIN_SAMPLE_SIZE)
        return(out)
    }
    if (out$events < GEP_MIN_DISCRIMINATION_EVENTS) {
        out$reason <- sprintf("Fewer than %d melanoma-specific deaths by %d months", GEP_MIN_DISCRIMINATION_EVENTS, horizon_months)
        return(out)
    }
    if (length(unique(risk)) < 2L) {
        out$reason <- "Imported risk does not vary"
        return(out)
    }
    auc_at <- function(idx) {
        fit <- timeROC::timeROC(
            T = time_months[idx], delta = status[idx], marker = risk[idx], cause = 1,
            weighting = "marginal", times = horizon_months, iid = FALSE
        )
        # timeROC returns the competing-risk AUCs only when the sample contains a
        # competing death; otherwise `AUC` is the same quantity (no controls to
        # reclassify). It also prepends t = 0 when a patient has zero follow-up,
        # so the requested horizon is always the last element.
        auc_by_time <- if (is.null(fit$AUC_2)) fit$AUC else fit$AUC_2
        auc_by_time[[length(auc_by_time)]]
    }
    # Bootstrap resamples can lack an event or a risk contrast; those draws are
    # counted as failed rather than silently replaced.
    auc_draw <- function(idx) {
        if (sum(status[idx] == 1L & time_months[idx] <= horizon_months) == 0L || length(unique(risk[idx])) < 2L) {
            return(NA_real_)
        }
        value <- tryCatch(auc_at(idx), error = function(e) NA_real_)
        if (length(value) == 1L) as.numeric(value) else NA_real_
    }
    out$auc <- as.numeric(auc_at(seq_along(risk)))
    out$estimable <- is.finite(out$auc)
    if (!out$estimable) {
        out$reason <- "timeROC returned a missing AUC at the horizon"
        return(out)
    }
    draws <- numeric(0)
    if (n_boot > 0L) {
        draws <- withr::with_seed(seed, {
            vapply(seq_len(n_boot), function(i) auc_draw(sample.int(length(risk), replace = TRUE)), numeric(1))
        })
    }
    interval <- bootstrap_percentile_interval(draws, n_boot)
    out[names(interval)] <- interval
    out
}

#' Proportional-hazards check for the imported risk among validation patients
#'
#' Fits a Cox model of the endpoint on the imported risk (one term, so groups
#' without events do not make the fit diverge) and tests scaled Schoenfeld
#' residuals. For MSS, the cause-specific melanoma-death hazard is modelled and
#' competing deaths are censored.
#'
#' @param prepared Data frame from `prepare_gep_validation_endpoint()`.
#' @return One-row data frame: `term`, `chisq`, `df`, `p_value`, `status`
#'   (`no concern`, `concern`, or `not assessed`) and `reason`.
check_gep_proportional_hazards <- function(prepared) {
    if (sum(prepared$status == 1L) < GEP_MIN_DISCRIMINATION_EVENTS || length(unique(prepared$published_risk)) < 2L) {
        return(data.frame(term = "published_risk", chisq = NA_real_, df = NA_real_, p_value = NA_real_,
                          status = "not assessed", reason = "Too few events or no variation in imported risk",
                          stringsAsFactors = FALSE))
    }
    fit <- survival::coxph(survival::Surv(time_months, status == 1L) ~ published_risk, data = prepared)
    zph <- survival::cox.zph(fit)
    row <- as.data.frame(zph$table)["published_risk", ]
    data.frame(
        term = "published_risk", chisq = row$chisq, df = row$df, p_value = row$p,
        status = if (row$p < 0.05) "concern" else "no concern",
        reason = "Schoenfeld residual test of a time-constant effect; p < 0.05 flagged",
        stringsAsFactors = FALSE
    )
}

#' Cause-of-death sensitivity for melanoma-specific risk
#'
#' The primary analysis counts deaths with an unknown recorded cause as competing
#' (non-melanoma) deaths. This recomputes the observed melanoma-specific risk in each
#' group counting them as melanoma deaths, the opposite extreme.
#'
#' @param endpoint_data Data frame with `time_months`, `status`, `status_unknown_melanoma`
#'   and the grouping column.
#' @param group_var Character name of the grouping column.
#' @param horizon_months Numeric horizon in months.
#' @return Data frame per group: patients, deaths reclassified, `primary_risk`,
#'   `unknown_cause_as_melanoma_risk`, and `difference_pct_points`.
compare_unknown_cause_of_death_assumption <- function(endpoint_data, group_var, horizon_months) {
    groups <- levels(endpoint_data[[group_var]])
    rows <- lapply(groups, function(group_name) {
        g <- endpoint_data[endpoint_data[[group_var]] %in% group_name, , drop = FALSE]
        primary <- estimate_cif_risk_at_horizon(g$time_months, g$status, horizon_months)
        alternative <- estimate_cif_risk_at_horizon(g$time_months, g$status_unknown_melanoma, horizon_months)
        data.frame(
            group = group_name, patients = nrow(g),
            deaths_reclassified = sum(g$status != g$status_unknown_melanoma, na.rm = TRUE),
            primary_risk = primary$risk, unknown_cause_as_melanoma_risk = alternative$risk,
            difference_pct_points = 100 * (alternative$risk - primary$risk), stringsAsFactors = FALSE
        )
    })
    dplyr::bind_rows(rows)
}

#' Sensitivity to GEP samples taken long after treatment
#'
#' Recomputes the group tables and discrimination without the patients whose first
#' biopsy was more than `GEP_TEST_TIMING_LATE_DAYS` after treatment (for example at
#' local recurrence), whose sample describes a different tumor state.
#'
#' @param mfs Prepared MFS endpoint data.
#' @param mss Prepared MSS endpoint data.
#' @param horizon_months Numeric horizon in months.
#' @return List: `patients_excluded` (per endpoint) and `table` with `endpoint`,
#'   `quantity`, `primary`, `excluding_late_biopsy`, `difference`, `unit`.
run_gep_late_biopsy_sensitivity <- function(mfs, mss, horizon_months) {
    one_endpoint <- function(prepared, endpoint) {
        keep <- prepared[!(as.character(prepared$gep_test_timing) %in% GEP_TEST_TIMING_LATE), , drop = FALSE]
        primary_table <- build_gep_group_table(prepared, endpoint, horizon_months)
        restricted_table <- build_gep_group_table(keep, endpoint, horizon_months)
        discrimination <- function(d) {
            if (identical(endpoint, "MFS")) {
                r <- calculate_gep_harrell_c(d$time_months, d$status, d$published_risk, horizon_months, n_boot = 0L, seed = 0L)
                list(value = r$c_index, label = "Harrell's C")
            } else {
                r <- calculate_gep_competing_risk_auc(d$time_months, d$status, d$published_risk, horizon_months, n_boot = 0L, seed = 0L)
                list(value = r$auc, label = "Competing-risk AUC")
            }
        }
        primary_disc <- discrimination(prepared)
        restricted_disc <- discrimination(keep)
        risk_rows <- data.frame(
            endpoint = endpoint, quantity = paste(sprintf("Observed %.0f-year risk:", horizon_months / 12), primary_table$group),
            primary = primary_table$observed_risk, excluding_late_biopsy = restricted_table$observed_risk,
            difference = 100 * (restricted_table$observed_risk - primary_table$observed_risk), unit = "percentage points",
            stringsAsFactors = FALSE
        )
        disc_row <- data.frame(
            endpoint = endpoint, quantity = primary_disc$label, primary = primary_disc$value,
            excluding_late_biopsy = restricted_disc$value, difference = restricted_disc$value - primary_disc$value,
            unit = "absolute", stringsAsFactors = FALSE
        )
        list(table = rbind(risk_rows, disc_row), excluded = nrow(prepared) - nrow(keep))
    }
    mfs_result <- one_endpoint(mfs, "MFS")
    mss_result <- one_endpoint(mss, "MSS")
    list(
        patients_excluded = c(MFS = mfs_result$excluded, MSS = mss_result$excluded),
        table = rbind(mfs_result$table, mss_result$table)
    )
}

#' Run the definitive-GEP 5-year validation for one cohort
#'
#' @param data Cohort data frame from the analytic dataset.
#' @param horizon_months Numeric horizon in months (default 60).
#' @param n_boot Integer bootstrap replicates.
#' @param seed Integer seed for all bootstraps in this run.
#' @return List with `horizon_months`, `mfs` and `mss` (each with `prepared`,
#'   `group_table`, `discrimination`, `proportional_hazards`),
#'   `late_biopsy_sensitivity`, and `cause_of_death_sensitivity`.
run_gep_definitive_validation <- function(data,
                                          horizon_months = GEP_VALIDATION_HORIZON_MONTHS,
                                          n_boot = GEP_BOOTSTRAP_ITERATIONS,
                                          seed = GEP_VALIDATION_BOOTSTRAP_SEED) {
    mfs <- prepare_gep_validation_endpoint(data, "MFS")
    mss <- prepare_gep_validation_endpoint(data, "MSS")
    late_biopsy_sensitivity <- run_gep_late_biopsy_sensitivity(mfs, mss, horizon_months)
    list(
        horizon_months = horizon_months,
        late_biopsy_sensitivity = late_biopsy_sensitivity,
        cause_of_death_sensitivity = compare_unknown_cause_of_death_assumption(mss, "validation_group", horizon_months),
        mfs = list(
            prepared = mfs,
            group_table = build_gep_group_table(mfs, "MFS", horizon_months),
            discrimination = calculate_gep_harrell_c(
                mfs$time_months, mfs$status, mfs$published_risk, horizon_months, n_boot, seed
            ),
            proportional_hazards = check_gep_proportional_hazards(mfs)
        ),
        mss = list(
            prepared = mss,
            group_table = build_gep_group_table(mss, "MSS", horizon_months),
            discrimination = calculate_gep_competing_risk_auc(
                mss$time_months, mss$status, mss$published_risk, horizon_months, n_boot, seed
            ),
            proportional_hazards = check_gep_proportional_hazards(mss)
        )
    )
}

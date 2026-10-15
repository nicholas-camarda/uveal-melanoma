# No-GEP analysis
#
# Describes outcomes in every GEP availability group and asks how well baseline
# clinical features rank 5-year risk in patients without a definitive Class 1/2
# call, against two no-fit comparators (T stage alone, tumor diameter alone) and a
# no-information Brier reference. Performance is out-of-fold; one analysis set per
# endpoint carries every table, so denominators agree across outputs.

#' Prepare the cohort for the no-GEP analysis
#'
#' Recodes internal reflectivity "Unknown" to missing (it is recorded as a category
#' for 46% of patients and would otherwise pass the completeness screen).
#'
#' @param data Cohort data frame with `gep_group` and the baseline predictors.
#' @return `data` with `internal_reflectivity` recoded and unused levels dropped.
prepare_no_gep_data <- function(data) {
    required <- c("gep_group", "id", GEP_NO_GEP_CANDIDATE_PREDICTORS)
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0L) {
        stop(sprintf("No-GEP analysis requires columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
    if (anyNA(data$id) || anyDuplicated(data$id)) {
        stop("No-GEP analysis requires a unique, observed patient identifier for deterministic folds.", call. = FALSE)
    }
    reflectivity <- coerce_to_factor_preserving_levels(data$internal_reflectivity)
    reflectivity[as.character(reflectivity) %in% "Unknown"] <- NA
    data$internal_reflectivity <- droplevels(reflectivity)
    if (anyNA(data$gep_group)) {
        stop(sprintf("gep_group is missing for %d patient(s); every patient must be classified.", sum(is.na(data$gep_group))), call. = FALSE)
    }
    data$gep_group <- factor(as.character(data$gep_group), levels = GEP_GROUP_LEVELS)
    data
}

#' Select the analysis set for one endpoint
#'
#' @param data Prepared cohort from `prepare_no_gep_data()`.
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @return Data frame of patients eligible for the endpoint, with `time_months` and
#'   `status` (0 censored, 1 event of interest, 2 competing death).
prepare_no_gep_endpoint <- function(data, endpoint = c("MFS", "MSS")) {
    endpoint <- match.arg(endpoint)
    spec <- switch(
        endpoint,
        MFS = list(time = "tt_mets_months_analysis", status = "mfs_event_type", baseline = "mets_free_at_baseline"),
        MSS = list(time = "tt_death_months", status = "mss_event_type", baseline = NULL)
    )
    required <- c(spec$time, spec$status, spec$baseline)
    missing_cols <- setdiff(required, names(data))
    if (length(missing_cols) > 0L) {
        stop(sprintf("%s analysis requires columns: %s", endpoint, paste(missing_cols, collapse = ", ")), call. = FALSE)
    }
    eligible <- is.finite(data[[spec$time]]) & !is.na(data[[spec$status]])
    if (!is.null(spec$baseline)) {
        eligible <- eligible & (data[[spec$baseline]] %in% TRUE)
    }
    prepared <- data[eligible, , drop = FALSE]
    prepared$time_months <- as.numeric(prepared[[spec$time]])
    prepared$status <- as.integer(prepared[[spec$status]])
    if (identical(endpoint, "MSS")) {
        if (!"mss_event_type_unknown_cause_melanoma" %in% names(prepared)) {
            stop("MSS analysis requires column: mss_event_type_unknown_cause_melanoma", call. = FALSE)
        }
        prepared$status_unknown_melanoma <- as.integer(prepared$mss_event_type_unknown_cause_melanoma)
    }
    prepared
}

#' Screen candidate predictors for completeness and sparse levels
#'
#' A predictor is retained only if, in every analysis set (MFS, MSS, and the
#' definitive Class 1/2 patients used for the surrogate), at least
#' `GEP_NO_GEP_MIN_COMPLETENESS` of values are observed and a factor keeps at least
#' two levels with `GEP_NO_GEP_MIN_LEVEL_COUNT` patients each.
#'
#' @param analysis_sets Named list of data frames (`MFS`, `MSS`, `Surrogate`).
#' @return Data frame with one row per candidate predictor: completeness in each set,
#'   `status` (`retained`/`dropped`), and `reason`.
screen_no_gep_predictors <- function(analysis_sets) {
    rows <- lapply(GEP_NO_GEP_CANDIDATE_PREDICTORS, function(predictor) {
        completeness <- vapply(analysis_sets, function(d) mean(!is.na(d[[predictor]])), numeric(1))
        reason <- NA_character_
        if (any(completeness < GEP_NO_GEP_MIN_COMPLETENESS)) {
            reason <- sprintf(
                "completeness below %.0f%% (%s)", 100 * GEP_NO_GEP_MIN_COMPLETENESS,
                paste(sprintf("%s %.0f%%", names(completeness), 100 * completeness), collapse = ", ")
            )
        } else if (predictor %in% GEP_NO_GEP_FACTOR_PREDICTORS) {
            for (set_name in names(analysis_sets)) {
                sparse <- summarize_sparse_factor_levels(
                    analysis_sets[[set_name]][[predictor]],
                    min_level_count = GEP_NO_GEP_MIN_LEVEL_COUNT,
                    explicit_exclusions = MODELING_LEVEL_EXCLUSIONS[[predictor]] %||% character()
                )
                if (!is.na(sparse$drop_reason)) {
                    reason <- sprintf("%s: %s", set_name, sparse$drop_reason)
                    break
                }
            }
        }
        data.frame(
            predictor = predictor,
            type = if (predictor %in% GEP_NO_GEP_FACTOR_PREDICTORS) "factor" else "numeric",
            completeness_mfs = completeness[["MFS"]], completeness_mss = completeness[["MSS"]],
            completeness_surrogate = completeness[["Surrogate"]],
            status = if (is.na(reason)) "retained" else "dropped",
            reason = if (is.na(reason)) "passes completeness and sparse-level screening" else reason,
            stringsAsFactors = FALSE
        )
    })
    dplyr::bind_rows(rows)
}

#' Build the model data for one analysis set
#'
#' Applies the sparse-level exclusions once and then keeps complete rows, so every
#' downstream table for the set uses the same patients. Removed patients and the
#' reason are returned for the Record_Issues sheet.
#'
#' @param analysis_set Data frame for one endpoint or the surrogate.
#' @param predictors Character vector of retained predictors.
#' @param label Character name of the analysis set.
#' @return List: `data` (model rows), `removed` (data frame with `id`, `reason`,
#'   `field`, `current_value`), and `counts` (input, sparse-level, incomplete, model).
build_no_gep_model_data <- function(analysis_set, predictors, label) {
    factor_predictors <- intersect(predictors, GEP_NO_GEP_FACTOR_PREDICTORS)
    initial_n <- nrow(analysis_set)
    after_sparse <- analysis_set
    if (length(factor_predictors) > 0L) {
        after_sparse <- apply_sparse_level_exclusions(
            analysis_set, variables = factor_predictors, analysis_name = label,
            min_level_count = GEP_NO_GEP_MIN_LEVEL_COUNT, id_col = "id",
            level_exclusions = MODELING_LEVEL_EXCLUSIONS
        )$data
    }
    sparse_removed <- setdiff(analysis_set$id, after_sparse$id)
    complete <- stats::complete.cases(after_sparse[, predictors, drop = FALSE])
    model_data <- after_sparse[complete, , drop = FALSE]
    for (variable in factor_predictors) {
        model_data[[variable]] <- droplevels(coerce_to_factor_preserving_levels(model_data[[variable]]))
    }
    incomplete_removed <- setdiff(after_sparse$id, model_data$id)

    removed_rows <- function(ids, reason, columns) {
        rows <- analysis_set[analysis_set$id %in% ids, , drop = FALSE]
        value <- vapply(seq_len(nrow(rows)), function(i) {
            parts <- vapply(columns, function(column) {
                v <- rows[[column]][i]
                paste0(column, " = ", if (is.na(v)) "missing" else as.character(v))
            }, character(1))
            paste(parts, collapse = "; ")
        }, character(1))
        data.frame(
            id = as.character(rows$id), analysis_set = rep(label, nrow(rows)), reason = rep(reason, nrow(rows)),
            field = rep(paste(columns, collapse = ", "), nrow(rows)), current_value = value, stringsAsFactors = FALSE
        )
    }
    removed <- dplyr::bind_rows(
        removed_rows(sparse_removed, "Factor level with fewer than the minimum patients", factor_predictors),
        removed_rows(incomplete_removed, "Missing baseline predictor", predictors)
    )
    list(
        data = model_data,
        removed = removed,
        counts = data.frame(
            analysis_set = label, eligible = initial_n, removed_sparse_level = length(sparse_removed),
            removed_missing_predictor = length(incomplete_removed), model_n = nrow(model_data),
            stringsAsFactors = FALSE
        )
    )
}

#' Out-of-fold horizon predictions for one endpoint and predictor set
#'
#' Wraps the repeated nested cross-validation (censoring weights estimated inside each
#' outer training fold) and attaches the GEP group and the comparator covariates.
#'
#' @param model_data Model rows from `build_no_gep_model_data()`.
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @param predictors Character vector of predictors for the ridge model.
#' @param horizon_months Numeric horizon in months.
#' @return Data frame of out-of-fold predictions keyed by `stable_id` and `repeat_id`,
#'   with `horizon_event`, `ipcw_weight`, `prediction`, `gep_group`, and the
#'   comparator covariates.
compute_no_gep_oof_predictions <- function(model_data, endpoint, predictors, horizon_months) {
    cv <- cross_validate_horizon_ridge(
        data = model_data, predictors = predictors,
        time_var = "time_months", event_type_var = "status",
        horizon_months = horizon_months, stable_id_var = "id"
    )
    context <- model_data[, c("id", "gep_group", unname(GEP_NO_GEP_COMPARATOR_PREDICTORS)), drop = FALSE]
    context$stable_id <- as.character(context$id)
    context$id <- NULL
    dplyr::left_join(tibble::as_tibble(cv$oof_predictions), context, by = "stable_id")
}

#' Select out-of-fold rows in a performance scope
#'
#' @param oof Out-of-fold predictions with `gep_group`.
#' @param scope Character, `"No GEP"` (every patient without a definitive Class 1/2
#'   call), or `"All patients"`.
#' @return The rows in scope.
filter_no_gep_scope <- function(oof, scope = c("No GEP", "All patients")) {
    scope <- match.arg(scope)
    if (identical(scope, "No GEP")) {
        oof[!(as.character(oof$gep_group) %in% c("Class 1", "Class 2")), , drop = FALSE]
    } else {
        oof
    }
}

#' One performance metric on one repeat's out-of-fold rows
#'
#' @param rows Out-of-fold rows for a single repeat.
#' @param score Numeric score aligned to `rows`.
#' @param metric Character: `"auc"`, `"brier"`, or `"null_brier"` (Brier of the
#'   weighted overall prevalence, the no-information reference).
#' @return Numeric value or `NA_real_` when case-control support is insufficient.
no_gep_repeat_metric <- function(rows, score, metric) {
    if (identical(metric, "auc")) {
        return(calculate_ipcw_auc(rows$horizon_event, score, rows$ipcw_weight)$auc)
    }
    if (identical(metric, "null_brier")) {
        known <- rows$ipcw_weight > 0 & !is.na(rows$horizon_event)
        score <- rep(sum(rows$ipcw_weight[known] * rows$horizon_event[known]) / sum(rows$ipcw_weight[known]), nrow(rows))
    }
    calculate_ipcw_brier(rows$horizon_event, pmin(pmax(score, 0), 1), rows$ipcw_weight)$brier
}

#' Out-of-fold metric with a patient-level bootstrap interval
#'
#' The estimate is the mean over outer-partition repeats. Each bootstrap draw
#' resamples the in-scope patients with replacement and uses one randomly chosen
#' repeat, so the interval reflects both patient sampling and fold assignment. It
#' does not re-run model fitting on bootstrap samples.
#'
#' @param oof Out-of-fold rows already restricted to the scope.
#' @param score_col Character name of the score column.
#' @param metric Character, see `no_gep_repeat_metric()`.
#' @param n_boot Integer bootstrap draws.
#' @param seed Integer seed, applied locally.
#' @param with_interval Logical; compute the bootstrap interval.
#' @return List: `estimate`, `ci_lower`, `ci_upper`, `valid_draws`, `failed_draws`.
summarize_no_gep_oof_metric <- function(oof, score_col, metric, n_boot, seed, with_interval = TRUE) {
    by_repeat <- split(oof, oof$repeat_id)
    values <- vapply(by_repeat, function(rows) {
        suppressWarnings(as.numeric(no_gep_repeat_metric(rows, rows[[score_col]], metric)))
    }, numeric(1))
    out <- list(estimate = mean(values, na.rm = TRUE), ci_lower = NA_real_, ci_upper = NA_real_, valid_draws = 0L, failed_draws = 0L)
    if (!is.finite(out$estimate)) {
        out$estimate <- NA_real_
        return(out)
    }
    if (!isTRUE(with_interval)) {
        return(out)
    }
    draws <- withr::with_seed(seed, {
        vapply(seq_len(n_boot), function(i) {
            rows <- by_repeat[[sample.int(length(by_repeat), 1L)]]
            resampled <- rows[sample.int(nrow(rows), replace = TRUE), , drop = FALSE]
            suppressWarnings(as.numeric(no_gep_repeat_metric(resampled, resampled[[score_col]], metric)))
        }, numeric(1))
    })
    out$valid_draws <- sum(is.finite(draws))
    out$failed_draws <- n_boot - out$valid_draws
    if (out$valid_draws >= 0.9 * n_boot) {
        bounds <- stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
        out$ci_lower <- bounds[[1]]
        out$ci_upper <- bounds[[2]]
    }
    out
}

#' Numeric ordering of the T stage for the no-fit comparator
#'
#' @param stage Factor or character with levels T1 to T4.
#' @return Numeric 1 to 4.
t_stage_to_score <- function(stage) {
    as.numeric(sub("^T", "", as.character(stage)))
}

#' Performance of the ridge models and no-fit comparators for one endpoint
#'
#' @param models Named list: `Full` and `Parsimonious` out-of-fold tables from
#'   `compute_no_gep_oof_predictions()`; the first table also supplies the comparator
#'   covariates (identical patients and weights across models).
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @param n_boot Integer bootstrap draws.
#' @param seed Integer seed.
#' @return Data frame with one row per model and scope: patients, cases, AUC with
#'   interval, Brier, no-information Brier, and Brier skill.
summarize_no_gep_model_performance <- function(models, endpoint, n_boot, seed) {
    rows <- list()
    for (scope in c("No GEP", "All patients")) {
        reference <- filter_no_gep_scope(models[[1]], scope)
        reference_rows <- reference[reference$repeat_id == min(reference$repeat_id), , drop = FALSE]
        scored <- list()
        for (model_name in names(models)) {
            scoped <- filter_no_gep_scope(models[[model_name]], scope)
            scored[[model_name]] <- list(oof = scoped, score = "prediction", brier = TRUE)
        }
        for (label in names(GEP_NO_GEP_COMPARATOR_PREDICTORS)) {
            variable <- GEP_NO_GEP_COMPARATOR_PREDICTORS[[label]]
            comparator <- filter_no_gep_scope(models[[1]], scope)
            comparator$comparator_score <- if (identical(variable, "initial_t_stage_simple")) t_stage_to_score(comparator[[variable]]) else as.numeric(comparator[[variable]])
            scored[[label]] <- list(oof = comparator, score = "comparator_score", brier = FALSE)
        }
        for (model_name in names(scored)) {
            item <- scored[[model_name]]
            auc <- summarize_no_gep_oof_metric(item$oof, item$score, "auc", n_boot, seed)
            brier <- if (isTRUE(item$brier)) summarize_no_gep_oof_metric(item$oof, item$score, "brier", n_boot, seed, with_interval = FALSE) else list(estimate = NA_real_)
            null_brier <- if (isTRUE(item$brier)) summarize_no_gep_oof_metric(item$oof, item$score, "null_brier", n_boot, seed, with_interval = FALSE) else list(estimate = NA_real_)
            first_repeat <- item$oof[item$oof$repeat_id == min(item$oof$repeat_id), , drop = FALSE]
            known <- first_repeat$ipcw_weight > 0 & !is.na(first_repeat$horizon_event)
            rows[[length(rows) + 1L]] <- data.frame(
                endpoint = endpoint, performance_scope = scope, model = model_name,
                patients = nrow(first_repeat), cases_by_horizon = sum(first_repeat$horizon_event[known] == 1L),
                auc = auc$estimate, auc_ci_lower = auc$ci_lower, auc_ci_upper = auc$ci_upper,
                bootstrap_valid_draws = auc$valid_draws, bootstrap_failed_draws = auc$failed_draws,
                brier = brier$estimate, no_information_brier = null_brier$estimate,
                brier_skill = if (is.finite(brier$estimate) && is.finite(null_brier$estimate)) 1 - brier$estimate / null_brier$estimate else NA_real_,
                stringsAsFactors = FALSE
            )
        }
    }
    dplyr::bind_rows(rows)
}

#' Binary AUC with ties counted as one half
#'
#' @param outcome Binary outcome vector coded 0/1.
#' @param score Numeric score; larger values rank higher risk.
#' @return Numeric AUC, or `NA_real_` when the outcome lacks both classes.
calculate_binary_auc <- function(outcome, score) {
    valid <- !is.na(outcome) & !is.na(score)
    outcome <- outcome[valid]
    score <- score[valid]
    n_positive <- sum(outcome == 1)
    n_negative <- sum(outcome == 0)
    if (n_positive == 0 || n_negative == 0) {
        return(NA_real_)
    }
    ranks <- rank(score, ties.method = "average")
    (sum(ranks[outcome == 1]) - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
}

#' Align assessment design columns with the training design
#'
#' Zero-fills a level that is absent from the assessment rows. A level present only in
#' the assessment rows would be scored as the reference level, so it stops the run
#' instead.
#'
#' @param data Rows to encode.
#' @param predictors Character predictor names.
#' @param reference_columns Optional training design columns.
#' @return Numeric design matrix without an intercept.
build_surrogate_design_matrix <- function(data, predictors, reference_columns = NULL) {
    design <- stats::model.matrix(stats::reformulate(predictors), data = data)
    design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
    if (!is.null(reference_columns)) {
        unseen <- setdiff(colnames(design), reference_columns)
        if (length(unseen) > 0L) {
            stop(sprintf("Assessment rows contain levels absent from training: %s", paste(unseen, collapse = ", ")), call. = FALSE)
        }
        absent <- setdiff(reference_columns, colnames(design))
        if (length(absent) > 0L) {
            design <- cbind(design, matrix(0, nrow(design), length(absent), dimnames = list(NULL, absent)))
        }
        design <- design[, reference_columns, drop = FALSE]
    }
    design
}

#' Repeated out-of-fold predictions for the Class 2 surrogate
#'
#' Ridge logistic regression of `Class 2` versus `Class 1` on the baseline predictors
#' among patients with a definitive call, using deterministic outer folds and
#' deterministic inner tuning folds, repeated over partitions.
#'
#' @param model_data Surrogate model rows with `id`, `gep_group`, and predictors.
#' @param predictors Character predictor names.
#' @return Data frame with `stable_id`, `repeat_id`, `prediction`, `outcome`.
compute_surrogate_oof_predictions <- function(model_data, predictors) {
    outcome <- as.integer(model_data$gep_group == "Class 2")
    stable_id <- as.character(model_data$id)
    purrr::map_dfr(seq_len(GEP_EXPLORATORY_CV_REPEATS), function(repeat_id) {
        seed <- as.integer(GEP_EXPLORATORY_CV_SEED + repeat_id - 1L)
        fold <- create_deterministic_fold_ids(strata = outcome, folds = GEP_EXPLORATORY_OUTER_FOLDS, seed = seed, stable_id = stable_id)
        prediction <- rep(NA_real_, nrow(model_data))
        for (k in seq_len(GEP_EXPLORATORY_OUTER_FOLDS)) {
            prediction[fold == k] <- tryCatch(
                {
                    training <- model_data[fold != k, , drop = FALSE]
                    assessment <- model_data[fold == k, , drop = FALSE]
                    x_training <- build_surrogate_design_matrix(training, predictors)
                    x_assessment <- build_surrogate_design_matrix(assessment, predictors, colnames(x_training))
                    y_training <- outcome[fold != k]
                    inner <- create_deterministic_fold_ids(
                        strata = y_training, folds = GEP_EXPLORATORY_INNER_FOLDS,
                        seed = as.integer(seed + 1000L + k), stable_id = as.character(training$id)
                    )
                    fit <- suppressWarnings(glmnet::cv.glmnet(
                        x = x_training, y = y_training, family = "binomial", alpha = 0, foldid = inner,
                        standardize = TRUE, type.measure = "deviance"
                    ))
                    as.numeric(stats::predict(fit, newx = x_assessment, s = "lambda.min", type = "response"))
                },
                error = function(error) {
                    stop(sprintf(
                        "Surrogate nested validation failed at repeat %d, outer fold %d: %s",
                        repeat_id, k, conditionMessage(error)
                    ), call. = FALSE)
                }
            )
        }
        data.frame(stable_id = stable_id, repeat_id = repeat_id, prediction = prediction, outcome = outcome, stringsAsFactors = FALSE)
    })
}

#' Surrogate AUC with a bootstrap interval
#'
#' Same interval construction as `summarize_no_gep_oof_metric()`: patient resampling
#' with one random partition repeat per draw.
#'
#' @param oof Output of `compute_surrogate_oof_predictions()`.
#' @param n_boot Integer bootstrap draws.
#' @param seed Integer seed.
#' @return Data frame row with patients, Class 2 patients, AUC and interval.
summarize_surrogate_performance <- function(oof, n_boot, seed) {
    by_repeat <- split(oof, oof$repeat_id)
    auc_of <- function(rows) calculate_binary_auc(rows$outcome, rows$prediction)
    estimate <- mean(vapply(by_repeat, auc_of, numeric(1)), na.rm = TRUE)
    draws <- withr::with_seed(seed, {
        vapply(seq_len(n_boot), function(i) {
            rows <- by_repeat[[sample.int(length(by_repeat), 1L)]]
            auc_of(rows[sample.int(nrow(rows), replace = TRUE), , drop = FALSE])
        }, numeric(1))
    })
    valid <- sum(is.finite(draws))
    bounds <- if (valid >= 0.9 * n_boot) stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE) else c(NA_real_, NA_real_)
    first <- by_repeat[[1]]
    data.frame(
        endpoint = "Class 2 vs Class 1", performance_scope = "Definitive Class 1/2 patients",
        model = "Class 2-like score (ridge)", patients = nrow(first), cases_by_horizon = sum(first$outcome == 1L),
        auc = estimate, auc_ci_lower = bounds[[1]], auc_ci_upper = bounds[[2]],
        bootstrap_valid_draws = valid, bootstrap_failed_draws = n_boot - valid,
        brier = NA_real_, no_information_brier = NA_real_, brier_skill = NA_real_,
        stringsAsFactors = FALSE
    )
}

#' Observed 5-year risk in every GEP availability group
#'
#' Kaplan-Meier (MFS) or Aalen-Johansen (MSS) risk at the horizon with counts at the
#' horizon, reverse-Kaplan-Meier median follow-up from treatment, and an `estimable`
#' flag. Nothing is carried past a group's last follow-up.
#'
#' @param endpoint_data Data frame from `prepare_no_gep_endpoint()`.
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @param horizon_months Numeric horizon in months.
#' @return Data frame with one row per group in `GEP_GROUP_LEVELS`.
build_no_gep_group_outcomes <- function(endpoint_data, endpoint, horizon_months) {
    rows <- lapply(GEP_GROUP_LEVELS, function(group_name) {
        g <- endpoint_data[endpoint_data$gep_group %in% group_name, , drop = FALSE]
        est <- estimate_endpoint_risk(g$time_months, g$status, horizon_months, endpoint)
        followup <- NA_real_
        if (nrow(g) > 0L) {
            # Reverse Kaplan-Meier: censoring is the "event", so the median is the
            # typical follow-up of patients still under observation.
            reverse <- survival::survfit(survival::Surv(g$time_months, g$status == 0L) ~ 1)
            followup <- unname(summary(reverse)$table["median"])
        }
        note <- dplyr::case_when(
            nrow(g) == 0L ~ "No patients in this group",
            !est$estimable ~ sprintf("Not followed to %d months (maximum %.0f months)", horizon_months, est$max_followup_months),
            est$events == 0L ~ "No events by the horizon; interval not estimable",
            TRUE ~ ""
        )
        data.frame(
            endpoint = endpoint, group = group_name, n = est$n, events_by_horizon = est$events,
            competing_deaths_by_horizon = est$competing, censored_before_horizon = est$censored_before,
            event_free_through_horizon = est$event_free_through, median_followup_months = followup,
            max_followup_months = est$max_followup_months, observed_risk = est$risk,
            observed_ci_lower = est$ci_lower, observed_ci_upper = est$ci_upper,
            estimable = est$estimable, note = note, stringsAsFactors = FALSE
        )
    })
    dplyr::bind_rows(rows)
}

#' Baseline characteristics by GEP availability group
#'
#' Descriptive only: medians with quartiles and counts with percentages, in the same
#' columns for every group and for all patients.
#'
#' @param data Prepared cohort with `gep_group`.
#' @return Long data frame: `characteristic`, `statistic`, and one column per group
#'   plus `All patients`.
build_no_gep_baseline_table <- function(data) {
    groups <- c(list(`All patients` = rep(TRUE, nrow(data))),
                stats::setNames(lapply(GEP_GROUP_LEVELS, function(g) data$gep_group %in% g), GEP_GROUP_LEVELS))
    median_iqr <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0L) return("not available")
        q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
        sprintf("%.1f (%.1f-%.1f)", q[2], q[1], q[3])
    }
    count_pct <- function(flag) {
        flag <- flag[!is.na(flag)]
        if (length(flag) == 0L) return("not available")
        sprintf("%d (%.0f%%)", sum(flag), 100 * mean(flag))
    }
    specs <- list(
        list("Patients", "n", function(d) as.character(nrow(d))),
        list("Age at diagnosis, years", "median (quartiles)", function(d) median_iqr(d$age_at_diagnosis)),
        list("Tumor height, mm", "median (quartiles)", function(d) median_iqr(d$initial_tumor_height)),
        list("Tumor diameter, mm", "median (quartiles)", function(d) median_iqr(d$initial_tumor_diameter)),
        list("T stage 3 or 4", "n (%)", function(d) count_pct(d$initial_t_stage_simple %in% c("T3", "T4"))),
        list("Male", "n (%)", function(d) count_pct(d$sex == "Male")),
        list("Choroidal location", "n (%)", function(d) count_pct(d$location == "Choroidal")),
        list("Subretinal fluid", "n (%)", function(d) count_pct(d$srf == "Yes")),
        list("Optic nerve involvement", "n (%)", function(d) count_pct(d$optic_nerve_involvement == 1)),
        list("Treated with GKSRS", "n (%)", function(d) count_pct(d$treatment_group == "GKSRS")),
        list("Treated in 2013 or later", "n (%)", function(d) count_pct(d$treatment_year >= 2013))
    )
    rows <- lapply(specs, function(spec) {
        values <- vapply(groups, function(in_group) spec[[3]](data[in_group, , drop = FALSE]), character(1))
        data.frame(characteristic = spec[[1]], statistic = spec[[2]], as.list(values), check.names = FALSE, stringsAsFactors = FALSE)
    })
    dplyr::bind_rows(rows)
}

#' Run the no-GEP analysis for the full cohort
#'
#' @param data Cohort data frame (the full cohort, with `gep_group`).
#' @param horizon_months Numeric horizon in months (default 60).
#' @param n_boot Integer bootstrap draws.
#' @param seed Integer seed for all bootstraps.
#' @return List: `group_outcomes`, `cause_of_death_sensitivity`, `baseline`, `predictor_screening`, `predictors`,
#'   `model_flow`, `removed_records`, `performance`, and `design_columns`.
run_no_gep_analysis <- function(data,
                                horizon_months = GEP_VALIDATION_HORIZON_MONTHS,
                                n_boot = GEP_BOOTSTRAP_ITERATIONS,
                                seed = GEP_NO_GEP_BOOTSTRAP_SEED) {
    prepared <- prepare_no_gep_data(data)
    sets <- list(
        MFS = prepare_no_gep_endpoint(prepared, "MFS"),
        MSS = prepare_no_gep_endpoint(prepared, "MSS"),
        Surrogate = prepared[prepared$gep_group %in% c("Class 1", "Class 2"), , drop = FALSE]
    )
    screening <- screen_no_gep_predictors(sets)
    predictors <- screening$predictor[screening$status == "retained"]
    if (length(predictors) == 0L) {
        stop("No baseline predictor passed the completeness and sparse-level screen.", call. = FALSE)
    }
    parsimonious <- GEP_NO_GEP_PARSIMONIOUS_PREDICTORS
    if (!all(parsimonious %in% predictors)) {
        stop(sprintf(
            "Parsimonious predictors were screened out: %s",
            paste(setdiff(parsimonious, predictors), collapse = ", ")
        ), call. = FALSE)
    }

    model_sets <- lapply(names(sets), function(label) build_no_gep_model_data(sets[[label]], predictors, label))
    names(model_sets) <- names(sets)
    # Patients who cannot enter an endpoint (metastatic at baseline, missing follow-up)
    # are listed so the flow reconciles with the cohort.
    not_eligible <- lapply(c(MFS = "MFS", MSS = "MSS"), function(endpoint) setdiff(prepared$id, sets[[endpoint]]$id))
    for (endpoint in names(not_eligible)) {
        ids <- not_eligible[[endpoint]]
        model_sets[[endpoint]]$counts$candidates <- nrow(prepared)
        model_sets[[endpoint]]$counts$not_eligible_for_endpoint <- length(ids)
        model_sets[[endpoint]]$removed <- dplyr::bind_rows(
            data.frame(
                id = as.character(ids), analysis_set = rep(endpoint, length(ids)),
                reason = rep("Not eligible for the endpoint (metastatic at baseline or no follow-up time)", length(ids)),
                field = rep(if (identical(endpoint, "MFS")) "mets_free_at_baseline, tt_mets_months_analysis" else "tt_death_months", length(ids)),
                current_value = rep("see cohort record", length(ids)), stringsAsFactors = FALSE
            ),
            model_sets[[endpoint]]$removed
        )
    }
    model_sets$Surrogate$counts$candidates <- sum(prepared$gep_group %in% c("Class 1", "Class 2"))
    model_sets$Surrogate$counts$not_eligible_for_endpoint <- 0L
    performance <- list()
    design_columns <- list()
    for (endpoint in c("MFS", "MSS")) {
        model_data <- model_sets[[endpoint]]$data
        oof <- list(
            Full = compute_no_gep_oof_predictions(model_data, endpoint, predictors, horizon_months),
            Parsimonious = compute_no_gep_oof_predictions(model_data, endpoint, parsimonious, horizon_months)
        )
        names(oof) <- c("Penalized, all retained predictors", "Penalized, parsimonious predictors")
        performance[[endpoint]] <- summarize_no_gep_model_performance(oof, endpoint, n_boot, seed)
        design_columns[[endpoint]] <- c(
            full = ncol(build_horizon_ridge_design_matrix(model_data, predictors)),
            parsimonious = ncol(build_horizon_ridge_design_matrix(model_data, parsimonious))
        )
    }
    surrogate_oof <- compute_surrogate_oof_predictions(model_sets$Surrogate$data, predictors)
    performance$Surrogate <- summarize_surrogate_performance(surrogate_oof, n_boot, seed)

    list(
        horizon_months = horizon_months,
        group_outcomes = dplyr::bind_rows(
            build_no_gep_group_outcomes(sets$MFS, "MFS", horizon_months),
            build_no_gep_group_outcomes(sets$MSS, "MSS", horizon_months)
        ),
        cause_of_death_sensitivity = compare_unknown_cause_of_death_assumption(
            transform(sets$MSS, gep_group = factor(as.character(gep_group), levels = GEP_GROUP_LEVELS)),
            "gep_group", horizon_months
        ),
        baseline = build_no_gep_baseline_table(prepared),
        predictor_screening = screening,
        predictors = predictors,
        parsimonious_predictors = parsimonious,
        model_flow = dplyr::bind_rows(lapply(model_sets, function(x) x$counts)),
        removed_records = dplyr::bind_rows(lapply(model_sets, function(x) x$removed)),
        performance = dplyr::bind_rows(performance),
        design_columns = design_columns
    )
}

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
    required <- c("gep_group", "id", "treatment_year", GEP_NO_GEP_CANDIDATE_PREDICTORS)
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
    if (anyNA(data$treatment_year)) {
        stop(sprintf("treatment_year is missing for %d patient(s).", sum(is.na(data$treatment_year))), call. = FALSE)
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
    after_sparse <- analysis_set
    if (length(factor_predictors) > 0L) {
        # Removing patients with one sparse level can leave another level below the
        # minimum, so the exclusion repeats until no patient is removed.
        repeat {
            previous_n <- nrow(after_sparse)
            after_sparse <- apply_sparse_level_exclusions(
                after_sparse, variables = factor_predictors, analysis_name = label,
                min_level_count = GEP_NO_GEP_MIN_LEVEL_COUNT, id_col = "id",
                level_exclusions = MODELING_LEVEL_EXCLUSIONS
            )$data
            if (nrow(after_sparse) == previous_n) break
        }
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
        value <- format_record_values(rows, columns)
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
            analysis_set = label, removed_sparse_level = length(sparse_removed),
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
        oof[!(as.character(oof$gep_group) %in% GEP_DEFINITIVE_GROUP_LEVELS), , drop = FALSE]
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
    interval <- bootstrap_percentile_interval(draws, n_boot)
    out[names(interval)] <- interval
    out
}

#' Numeric ordering of the T stage for the no-fit comparator
#'
#' @param stage Factor or character with levels T1 to T4.
#' @return Numeric 1 to 4.
t_stage_to_score <- function(stage) {
    as.numeric(sub("^T", "", as.character(stage)))
}

#' Paired AUC difference between a model and a no-fit comparator
#'
#' Each bootstrap draw resamples the in-scope patients once and one partition repeat,
#' and computes both AUCs on that same resample, so the interval is for the
#' difference itself and not for two separate AUCs.
#'
#' @param oof Out-of-fold rows already restricted to the scope, with `prediction`
#'   and the comparator covariate.
#' @param comparator_variable Character name of the comparator column
#'   (`initial_t_stage_simple` or `initial_tumor_diameter`).
#' @param n_boot Integer bootstrap draws.
#' @param seed Integer seed, applied locally.
#' @return List: `estimate` (model minus comparator), `ci_lower`, `ci_upper`,
#'   `valid_draws`, `failed_draws`.
summarize_no_gep_auc_difference <- function(oof, comparator_variable, n_boot, seed) {
    oof$comparator_score <- if (identical(comparator_variable, "initial_t_stage_simple")) {
        t_stage_to_score(oof[[comparator_variable]])
    } else {
        as.numeric(oof[[comparator_variable]])
    }
    by_repeat <- split(oof, oof$repeat_id)
    difference <- function(rows) {
        model_auc <- calculate_ipcw_auc(rows$horizon_event, rows$prediction, rows$ipcw_weight)$auc
        comparator_auc <- calculate_ipcw_auc(rows$horizon_event, rows$comparator_score, rows$ipcw_weight)$auc
        model_auc - comparator_auc
    }
    estimate <- mean(vapply(by_repeat, function(rows) suppressWarnings(as.numeric(difference(rows))), numeric(1)), na.rm = TRUE)
    draws <- withr::with_seed(seed, {
        vapply(seq_len(n_boot), function(i) {
            rows <- by_repeat[[sample.int(length(by_repeat), 1L)]]
            suppressWarnings(as.numeric(difference(rows[sample.int(nrow(rows), replace = TRUE), , drop = FALSE])))
        }, numeric(1))
    })
    c(list(estimate = estimate), bootstrap_percentile_interval(draws, n_boot))
}

#' AUC differences of the ridge models against each comparator in the no-GEP scope
#'
#' @param models Named list of out-of-fold tables, as in `summarize_no_gep_model_performance()`.
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @param n_boot Integer bootstrap draws.
#' @param seed Integer seed.
#' @return Data frame: endpoint, model, comparator, difference with interval, draws.
summarize_no_gep_auc_differences <- function(models, endpoint, n_boot, seed) {
    rows <- list()
    for (model_name in names(models)) {
        scoped <- filter_no_gep_scope(models[[model_name]], "No GEP")
        for (label in names(GEP_NO_GEP_COMPARATOR_PREDICTORS)) {
            d <- summarize_no_gep_auc_difference(scoped, GEP_NO_GEP_COMPARATOR_PREDICTORS[[label]], n_boot, seed)
            rows[[length(rows) + 1L]] <- data.frame(
                endpoint = endpoint, performance_scope = "No GEP", model = model_name, comparator = label,
                auc_difference = d$estimate, ci_lower = d$ci_lower, ci_upper = d$ci_upper,
                bootstrap_valid_draws = d$valid_draws, bootstrap_failed_draws = d$failed_draws,
                stringsAsFactors = FALSE
            )
        }
    }
    dplyr::bind_rows(rows)
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
                    x_training <- build_horizon_ridge_design_matrix(training, predictors)
                    x_assessment <- build_horizon_ridge_design_matrix(assessment, predictors, colnames(x_training))
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
    interval <- bootstrap_percentile_interval(draws, n_boot)
    first <- by_repeat[[1]]
    data.frame(
        endpoint = "Class 2 vs Class 1", performance_scope = "Definitive Class 1/2 patients",
        model = "Class 2-like score (ridge)", patients = nrow(first), cases_by_horizon = sum(first$outcome == 1L),
        auc = estimate, auc_ci_lower = interval$ci_lower, auc_ci_upper = interval$ci_upper,
        bootstrap_valid_draws = interval$valid_draws, bootstrap_failed_draws = interval$failed_draws,
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

#' Sensitivity of the group risks to treatment era
#'
#' Re-estimates each group's 5-year risk among patients treated in or after
#' `GEP_TESTING_ERA_START_YEAR`, so the tested and untested groups cover the same era
#' (testing appears in this cohort only from that year on). The primary row of each
#' group is the full-group estimate from `build_no_gep_group_outcomes()`.
#'
#' @param endpoint_data Prepared endpoint data with `gep_group`, `treatment_year`,
#'   `time_months` and `status`.
#' @param endpoint Character, `"MFS"` or `"MSS"`.
#' @param horizon_months Numeric horizon in months.
#' @return Data frame, one row per group: patients, patients in and before the era,
#'   primary risk, risk in the era with its interval, the difference in percentage
#'   points, whether the era risk is estimable, and a note.
build_no_gep_treatment_era_sensitivity <- function(endpoint_data, endpoint, horizon_months) {
    rows <- lapply(GEP_GROUP_LEVELS, function(group_name) {
        g <- endpoint_data[endpoint_data$gep_group %in% group_name, , drop = FALSE]
        in_era <- g[g$treatment_year >= GEP_TESTING_ERA_START_YEAR, , drop = FALSE]
        primary <- estimate_endpoint_risk(g$time_months, g$status, horizon_months, endpoint)
        era <- estimate_endpoint_risk(in_era$time_months, in_era$status, horizon_months, endpoint)
        note <- dplyr::case_when(
            nrow(g) == 0L ~ "No patients in this group",
            nrow(in_era) == 0L ~ "No patient treated in the era",
            !era$estimable ~ sprintf("Not followed to %d months in the era (maximum %.0f months)", horizon_months, era$max_followup_months),
            era$events == 0L ~ "No events by the horizon in the era; interval not estimable",
            nrow(g) == nrow(in_era) ~ "Every patient was treated in the era; same as the primary estimate",
            TRUE ~ ""
        )
        data.frame(
            endpoint = endpoint, group = group_name, patients = nrow(g),
            patients_in_era = nrow(in_era), patients_before_era = nrow(g) - nrow(in_era),
            primary_risk = primary$risk, risk_in_era = era$risk,
            ci_lower_in_era = era$ci_lower, ci_upper_in_era = era$ci_upper,
            difference_pct_points = 100 * (era$risk - primary$risk),
            estimable_in_era = era$estimable, note = note, stringsAsFactors = FALSE
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
        list(sprintf("Treated in %d or later", GEP_TESTING_ERA_START_YEAR), "n (%)", function(d) count_pct(d$treatment_year >= GEP_TESTING_ERA_START_YEAR))
    )
    rows <- lapply(specs, function(spec) {
        values <- vapply(groups, function(in_group) {
            # Groups with fewer than GEP_MIN_REPORTED_GROUP_SIZE patients show only their size, so a
            # published table cannot describe one or two individuals.
            if (sum(in_group) < GEP_MIN_REPORTED_GROUP_SIZE && !identical(spec[[1]], "Patients")) {
                return(sprintf("fewer than %d patients", GEP_MIN_REPORTED_GROUP_SIZE))
            }
            spec[[3]](data[in_group, , drop = FALSE])
        }, character(1))
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
#' @return List: `prepared` (the cohort after `prepare_no_gep_data()`), `group_outcomes`, `treatment_era_sensitivity`, `treatment_era_excluded_n`,
#'   `cause_of_death_sensitivity`, `baseline`, `predictor_screening`, `predictors`,
#'   `model_flow`, `removed_records`, `performance`, `auc_differences`, and
#'   `design_columns`.
run_no_gep_analysis <- function(data,
                                horizon_months = GEP_VALIDATION_HORIZON_MONTHS,
                                n_boot = GEP_BOOTSTRAP_ITERATIONS,
                                seed = GEP_NO_GEP_BOOTSTRAP_SEED) {
    prepared <- prepare_no_gep_data(data)
    sets <- list(
        MFS = prepare_no_gep_endpoint(prepared, "MFS"),
        MSS = prepare_no_gep_endpoint(prepared, "MSS"),
        Surrogate = prepared[prepared$gep_group %in% GEP_DEFINITIVE_GROUP_LEVELS, , drop = FALSE]
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
        fields <- intersect(
            if (identical(endpoint, "MFS")) c("mets_free_at_baseline", "tt_mets_months_analysis", "mfs_event_type") else c("tt_death_months", "mss_event_type"),
            names(prepared)
        )
        model_sets[[endpoint]]$removed <- dplyr::bind_rows(
            data.frame(
                id = as.character(ids), analysis_set = rep(endpoint, length(ids)),
                reason = rep("Not eligible for the endpoint (metastatic at baseline or no follow-up time)", length(ids)),
                field = rep(paste(fields, collapse = ", "), length(ids)),
                current_value = format_record_values(prepared[prepared$id %in% ids, , drop = FALSE], fields),
                stringsAsFactors = FALSE
            ),
            model_sets[[endpoint]]$removed
        )
    }
    model_sets$Surrogate$counts$candidates <- sum(prepared$gep_group %in% GEP_DEFINITIVE_GROUP_LEVELS)
    model_sets$Surrogate$counts$not_eligible_for_endpoint <- 0L
    performance <- list()
    auc_differences <- list()
    design_columns <- list()
    for (endpoint in c("MFS", "MSS")) {
        model_data <- model_sets[[endpoint]]$data
        oof <- list(
            compute_no_gep_oof_predictions(model_data, endpoint, predictors, horizon_months),
            compute_no_gep_oof_predictions(model_data, endpoint, parsimonious, horizon_months)
        )
        names(oof) <- unname(GEP_NO_GEP_MODEL_LABELS[c("full", "parsimonious")])
        performance[[endpoint]] <- summarize_no_gep_model_performance(oof, endpoint, n_boot, seed)
        auc_differences[[endpoint]] <- summarize_no_gep_auc_differences(oof, endpoint, n_boot, seed)
        design_columns[[endpoint]] <- c(
            full = ncol(build_horizon_ridge_design_matrix(model_data, predictors)),
            parsimonious = ncol(build_horizon_ridge_design_matrix(model_data, parsimonious))
        )
    }
    surrogate_oof <- compute_surrogate_oof_predictions(model_sets$Surrogate$data, predictors)
    performance$Surrogate <- summarize_surrogate_performance(surrogate_oof, n_boot, seed)

    # Patients treated before the testing era stay in every primary table; they are left out
    # of the treatment-era sensitivity only, and listed so each one can be traced.
    before_era <- prepared$treatment_year < GEP_TESTING_ERA_START_YEAR
    era_removed <- data.frame(
        id = as.character(prepared$id[before_era]),
        analysis_set = rep("Treatment-era sensitivity", sum(before_era)),
        reason = rep(sprintf("Treated before %d: left out of the treatment-era sensitivity only", GEP_TESTING_ERA_START_YEAR), sum(before_era)),
        field = rep("treatment_year", sum(before_era)),
        current_value = as.character(prepared$treatment_year[before_era]),
        stringsAsFactors = FALSE
    )

    list(
        horizon_months = horizon_months,
        prepared = prepared,
        group_outcomes = dplyr::bind_rows(
            build_no_gep_group_outcomes(sets$MFS, "MFS", horizon_months),
            build_no_gep_group_outcomes(sets$MSS, "MSS", horizon_months)
        ),
        treatment_era_sensitivity = dplyr::bind_rows(
            build_no_gep_treatment_era_sensitivity(sets$MFS, "MFS", horizon_months),
            build_no_gep_treatment_era_sensitivity(sets$MSS, "MSS", horizon_months)
        ),
        treatment_era_excluded_n = sum(before_era),
        cause_of_death_sensitivity = compare_unknown_cause_of_death_assumption(
            transform(sets$MSS, gep_group = factor(as.character(gep_group), levels = GEP_GROUP_LEVELS)),
            "gep_group", horizon_months
        ),
        baseline = build_no_gep_baseline_table(prepared),
        predictor_screening = screening,
        predictors = predictors,
        parsimonious_predictors = parsimonious,
        model_flow = dplyr::bind_rows(lapply(model_sets, function(x) x$counts)),
        removed_records = dplyr::bind_rows(c(lapply(model_sets, function(x) x$removed), list(era_removed))),
        performance = dplyr::bind_rows(performance),
        auc_differences = dplyr::bind_rows(auc_differences),
        design_columns = design_columns
    )
}

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

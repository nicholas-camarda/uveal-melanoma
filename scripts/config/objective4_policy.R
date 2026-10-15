# =============================================================================
# OBJECTIVE 4 GEP VALIDATION POLICY
# =============================================================================
# These settings control the GEP validation analysis (Objective 4):
# - Definitive-GEP validation: one 5-year horizon, bootstrap settings, minimum counts
# - Exploratory no-GEP analysis: IPCW and nested cross-validation settings
# - GEP availability groups and the raw laboratory labels that define them

GEP_BOOTSTRAP_ITERATIONS <- 1000

# Definitive-GEP validation (Class 1/2 with PRAME reported): one horizon, one seed.
# Seven- and ten-year results are not reported because follow-up leaves 0-1
# patients at risk at 10 years and the longer horizons were extrapolated.
GEP_VALIDATION_HORIZON_YEARS <- 5
GEP_VALIDATION_HORIZON_MONTHS <- GEP_VALIDATION_HORIZON_YEARS * 12
GEP_VALIDATION_BOOTSTRAP_SEED <- 20261013L
GEP_MIN_DISCRIMINATION_EVENTS <- 3 # Minimum events by the horizon to estimate discrimination
GEP_MIN_AT_RISK_AT_HORIZON <- 5 # Fewer patients known event-free at the horizon is flagged as sparse
GEP_VALIDATION_GROUP_LEVELS <- c(
    "Class 1, PRAME negative",
    "Class 1, PRAME positive",
    "Class 2, PRAME negative",
    "Class 2, PRAME positive"
)

# No-GEP analysis: baseline predictors, comparators, screening, and bootstrap settings.
# Internal reflectivity is excluded because "Unknown" (46% of patients) is recorded as a
# category; it is recoded to missing and fails the completeness screen.
GEP_NO_GEP_CANDIDATE_PREDICTORS <- c(
    "age_at_diagnosis", "sex", "initial_tumor_height", "initial_tumor_diameter", "location",
    "initial_t_stage_simple", "internal_reflectivity", "srf", "initial_vision", "optic_nerve_involvement"
)
GEP_NO_GEP_FACTOR_PREDICTORS <- c("sex", "location", "initial_t_stage_simple", "internal_reflectivity", "srf")
GEP_NO_GEP_PARSIMONIOUS_PREDICTORS <- c("age_at_diagnosis", "initial_tumor_diameter", "location", "initial_t_stage_simple")
GEP_NO_GEP_COMPARATOR_PREDICTORS <- c(
    "T stage alone" = "initial_t_stage_simple",
    "Tumor diameter alone" = "initial_tumor_diameter"
)
GEP_NO_GEP_MIN_COMPLETENESS <- 0.9 # Minimum non-missing share for a candidate predictor
GEP_NO_GEP_MIN_LEVEL_COUNT <- 5 # Minimum patients per retained factor level
GEP_NO_GEP_MIN_EVENTS_PER_PARAMETER <- 10 # Fewer events per design column is flagged
GEP_NO_GEP_MIN_BRIER_SKILL <- 0.05 # Brier skill below this is flagged as little improvement on prevalence
GEP_NO_GEP_BOOTSTRAP_SEED <- 20261015L

# No-GEP exploratory model settings
GEP_MIN_CALIBRATION_EVENTS <- 5 # Minimum events and non-events for stable recalibration fits
GEP_IPCW_WEIGHT_CAP_PROB <- 0.99 # Cap extreme inverse-censoring weights at the 99th percentile
GEP_MAX_CALIBRATION_COEF_ABS <- 10 # Treat larger recalibration coefficients as unstable/quasi-separated fits
GEP_MAX_CALIBRATION_COEF_SE <- 5 # Treat larger recalibration coefficient standard errors as unstable fits

# Deterministic nested cross-validation for the Objective 4 no-GEP analysis
GEP_EXPLORATORY_CV_SEED <- 20260820L
GEP_EXPLORATORY_CV_REPEATS <- 20L
GEP_EXPLORATORY_OUTER_FOLDS <- 5L
GEP_EXPLORATORY_INNER_FOLDS <- 5L

# Sample size requirement for any analysis
GEP_MIN_SAMPLE_SIZE <- 20

GEP_DEFINITIVE_SIMPLE_LEVELS <- c("Class 1", "Class 2")
GEP_CLASS_1_DEFINITIVE_RAW_LEVELS <- c(
    "Class_1A_PRAME_negative",
    "Class_1A_PRAME_positive",
    "Class_1B_PRAME_negative",
    "Class_1B_PRAME_positive"
)
GEP_CLASS_2_DEFINITIVE_RAW_LEVELS <- c(
    "Class_2_PRAME_negative",
    "Class_2_PRAME_positive"
)
GEP_DEFINITIVE_RAW_LEVELS <- c(
    GEP_CLASS_1_DEFINITIVE_RAW_LEVELS,
    GEP_CLASS_2_DEFINITIVE_RAW_LEVELS
)
# Source-string levels that are not definitive Class 1/2 calls. The raw strings
# keep the laboratory's historical 1A/1B wording; reader-facing labels use the
# current Class 1 / Class 2 terminology only.
GEP_CLASS_KNOWN_PRAME_MISSING_RAW_LEVELS <- c(
    "Class_1A_PRAME_not_reported",
    "Class_2_PRAME_not_reported",
    "Class_2_PRAME_Unknown"
)
GEP_DISCORDANT_RAW_LEVELS <- "Class_1A_PRAME_discordant"
GEP_TEST_FAILED_RAW_LEVELS <- c("Failed", "Unknown")
GEP_FAILED_OR_INDETERMINATE_RAW_LEVELS <- c(
    GEP_TEST_FAILED_RAW_LEVELS,
    GEP_CLASS_KNOWN_PRAME_MISSING_RAW_LEVELS,
    GEP_DISCORDANT_RAW_LEVELS
)
GEP_NOT_TESTED_RAW_LEVELS <- c("No", "N/A")

# Reader-facing GEP availability groups. Class 1 / Class 2 are definitive calls
# with PRAME reported (the primary validation set). Patients whose class is known
# but PRAME is missing keep no imputed prediction and are reported separately
# from assay failures.
GEP_GROUP_CLASS_PRAME_MISSING <- "Class known, PRAME not reported"
GEP_GROUP_DISCORDANT <- "Discordant result"
GEP_GROUP_TEST_FAILED <- "Failed or unknown result"
GEP_GROUP_NOT_TESTED <- "Not tested"
GEP_GROUP_LEVELS <- c(
    "Class 1",
    "Class 2",
    GEP_GROUP_CLASS_PRAME_MISSING,
    GEP_GROUP_DISCORDANT,
    GEP_GROUP_TEST_FAILED,
    GEP_GROUP_NOT_TESTED
)

# Recorded causes of death that do not say whether the death was melanoma-specific.
# The primary MSS analysis treats them as competing (non-melanoma) deaths; the
# sensitivity analysis treats them as melanoma-specific deaths.
GEP_UNKNOWN_CAUSE_OF_DEATH_LEVELS <- c("Unknown", "Unknown_Cancer")

# Timing of the first biopsy (the GEP sample) relative to the primary treatment date.
# A sample taken long after treatment (for example at local recurrence) describes a
# different tumor state and starts the GEP-to-outcome clock after follow-up began.
GEP_TEST_TIMING_WINDOW_DAYS <- 31
GEP_TEST_TIMING_LATE_DAYS <- 180
GEP_TEST_TIMING_LEVELS <- c(
    "More than 31 days before treatment",
    "Within 31 days of treatment",
    "31 to 180 days after treatment",
    "More than 180 days after treatment",
    "No biopsy date"
)
GEP_TEST_TIMING_LATE <- GEP_TEST_TIMING_LEVELS[[4]]

#' Classify the biopsy date relative to treatment
#'
#' @param days_from_treatment Numeric days from treatment to biopsy (negative = before).
#' @return Character vector of timing categories (`GEP_TEST_TIMING_LEVELS`).
classify_gep_test_timing <- function(days_from_treatment) {
    dplyr::case_when(
        is.na(days_from_treatment) ~ GEP_TEST_TIMING_LEVELS[[5]],
        days_from_treatment < -GEP_TEST_TIMING_WINDOW_DAYS ~ GEP_TEST_TIMING_LEVELS[[1]],
        days_from_treatment <= GEP_TEST_TIMING_WINDOW_DAYS ~ GEP_TEST_TIMING_LEVELS[[2]],
        days_from_treatment <= GEP_TEST_TIMING_LATE_DAYS ~ GEP_TEST_TIMING_LEVELS[[3]],
        TRUE ~ GEP_TEST_TIMING_LEVELS[[4]]
    )
}

# The laboratory's discordant-result string as typed in the source spreadsheet.
GEP_DISCORDANT_SOURCE_LABEL <- "DISCORDANT CASTLE RESULTS: Class 1A, PRAME not reported"

#' Normalize the source GEP label to the raw level names
#'
#' The only source string that is not already a raw level is the discordant-result
#' text; it maps to `GEP_DISCORDANT_RAW_LEVELS`. Every classification of the raw
#' label goes through this function so the mapping lives in one place.
#'
#' @param label Character or factor vector of source GEP labels.
#' @return Character vector with the discordant source string replaced.
normalize_gep_source_label <- function(label) {
    label <- as.character(label)
    label[label %in% GEP_DISCORDANT_SOURCE_LABEL] <- GEP_DISCORDANT_RAW_LEVELS
    label
}

#' Classify raw GEP source labels into reader-facing availability groups
#'
#' Single mapping from the laboratory source string to `GEP_GROUP_LEVELS`.
#' Every raw level belongs to exactly one group; missing values mean no test.
#'
#' @param raw_label Character or factor vector of raw GEP source labels
#'   (for example `Class_1A_PRAME_negative`).
#' @return Character vector of group labels. Unrecognised non-missing labels
#'   return `NA` so they surface in validation instead of being absorbed.
classify_gep_group <- function(raw_label) {
    raw_label <- as.character(raw_label)
    dplyr::case_when(
        raw_label %in% GEP_CLASS_1_DEFINITIVE_RAW_LEVELS ~ "Class 1",
        raw_label %in% GEP_CLASS_2_DEFINITIVE_RAW_LEVELS ~ "Class 2",
        raw_label %in% GEP_CLASS_KNOWN_PRAME_MISSING_RAW_LEVELS ~ GEP_GROUP_CLASS_PRAME_MISSING,
        raw_label %in% GEP_DISCORDANT_RAW_LEVELS ~ GEP_GROUP_DISCORDANT,
        raw_label %in% GEP_TEST_FAILED_RAW_LEVELS ~ GEP_GROUP_TEST_FAILED,
        raw_label %in% GEP_NOT_TESTED_RAW_LEVELS | is.na(raw_label) ~ GEP_GROUP_NOT_TESTED,
        TRUE ~ NA_character_
    )
}
GEP_INVALID_ANALYSIS_LABELS <- c(
    GEP_FAILED_OR_INDETERMINATE_RAW_LEVELS,
    GEP_NOT_TESTED_RAW_LEVELS,
    "GEP Failed/Indeterminate",
    "GEP Not Tested"
)

# GEP-specific derived variables
GEP_DERIVED_VARIABLES <- c(
    "gep_class_simple", "prame_status", "expected_mfs_5yr", "expected_mss_5yr"
)

GEP_DISPLAY_VARIABLES <- c(
    "biopsy1_gep", "gep_class_simple", "prame_status", "gep12_prame_status"
)

# Objective 4 grouping specs and per-context defaults.
# Change these mappings when reader-facing or technical GEP groupings need to move together.
GEP_GROUPING_SPECS <- list(
    biopsy1_gep = list(
        key = "biopsy1_gep",
        var = "biopsy1_gep",
        label = "GEP Class",
        legend_title = "GEP Class",
        allowed_levels = NULL,
        reader_facing = FALSE
    ),
    biopsy1_gep_model = list(
        key = "biopsy1_gep_model",
        var = "biopsy1_gep_model",
        label = "GEP Class",
        legend_title = "GEP Class",
        allowed_levels = NULL,
        reader_facing = FALSE
    ),
    gep_class_simple = list(
        key = "gep_class_simple",
        var = "gep_class_simple",
        label = "Simplified GEP Class",
        legend_title = "GEP Class (Simple)",
        allowed_levels = c("Class 1", "Class 2"),
        reader_facing = TRUE
    )
)

OBJECTIVE4_GEP_GROUPING <- list(
    mfs = list(
        observed_expected = "biopsy1_gep",
        visuals = "biopsy1_gep",
        visuals_model = "biopsy1_gep_model"
    ),
    mss = list(
        observed_expected = "biopsy1_gep",
        competing_risk = "biopsy1_gep",
        reporting = "biopsy1_gep",
        visuals = "gep_class_simple"
    )
)

get_gep_grouping_spec <- function(grouping_key) {
    if (!grouping_key %in% names(GEP_GROUPING_SPECS)) {
        stop(sprintf("Unknown Objective 4 grouping key: %s", grouping_key))
    }

    GEP_GROUPING_SPECS[[grouping_key]]
}

get_gep_grouping_for_context <- function(outcome, context) {
    if (!outcome %in% names(OBJECTIVE4_GEP_GROUPING)) {
        stop(sprintf("Unknown Objective 4 outcome grouping scope: %s", outcome))
    }

    outcome_grouping <- OBJECTIVE4_GEP_GROUPING[[outcome]]
    if (!context %in% names(outcome_grouping)) {
        stop(sprintf("Unknown Objective 4 grouping context '%s' for outcome '%s'", context, outcome))
    }

    get_gep_grouping_spec(outcome_grouping[[context]])
}

#' Collapse GEP availability groups to the simple class label used by eligibility flags
#'
#' @param gep_group Character vector from `classify_gep_group()`.
#' @return Character vector: `Class 1`, `Class 2`, `GEP Not Tested`, or
#'   `GEP Failed/Indeterminate` (every other group); `NA` stays `NA`.
derive_gep_class_simple <- function(gep_group) {
    dplyr::case_when(
        gep_group %in% c("Class 1", "Class 2") ~ as.character(gep_group),
        gep_group %in% GEP_GROUP_NOT_TESTED ~ "GEP Not Tested",
        !is.na(gep_group) ~ "GEP Failed/Indeterminate",
        TRUE ~ NA_character_
    )
}

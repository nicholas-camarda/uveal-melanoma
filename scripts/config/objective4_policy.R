# =============================================================================
# OBJECTIVE 4 GEP VALIDATION POLICY
# =============================================================================
# These settings control the GEP validation analysis (Objective 4):
# - Definitive-GEP validation: one 5-year horizon, bootstrap settings, minimum counts
# - Exploratory no-GEP analysis: IPCW and nested cross-validation settings
# - GEP availability groups and the raw laboratory labels that define them

GEP_BOOTSTRAP_ITERATIONS <- 1000
GEP_MIN_VALID_BOOTSTRAP_SHARE <- 0.9 # Share of bootstrap draws that must be valid before an interval is reported

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
# Labels of the two penalized models; the analysis names its results with them and the checks look them up.
GEP_NO_GEP_MODEL_LABELS <- c(
    full = "Penalized, all retained predictors",
    parsimonious = "Penalized, parsimonious predictors"
)
GEP_NO_GEP_COMPARATOR_PREDICTORS <- c(
    "T stage alone" = "initial_t_stage_simple",
    "Tumor diameter alone" = "initial_tumor_diameter"
)
GEP_NO_GEP_MIN_COMPLETENESS <- 0.9 # Minimum non-missing share for a candidate predictor
GEP_NO_GEP_MIN_LEVEL_COUNT <- 5 # Minimum patients per retained factor level
GEP_NO_GEP_MIN_EVENTS_PER_PARAMETER <- 10 # Fewer events per design column is flagged
# The no-GEP analysis runs on the full cohort: the restricted and GKSRS-only cohorts are subsets of
# it, so a second run would repeat the same patients.
GEP_NO_GEP_DATASET <- "uveal_melanoma_full_cohort"

# GEP testing appears in this cohort from this treatment year on. The treatment-era sensitivity
# re-estimates every group's risk among patients treated in or after this year, so tested and
# untested groups cover the same era; the baseline table reports the share treated in or after it.
GEP_TESTING_ERA_START_YEAR <- 2013
GEP_MIN_REPORTED_GROUP_SIZE <- 3 # Groups smaller than this show only their size in published no-GEP tables
GEP_NO_GEP_MIN_BRIER_SKILL <- 0.05 # Brier skill below this is flagged as little improvement on prevalence
GEP_NO_GEP_BOOTSTRAP_SEED <- 20261015L

# No-GEP exploratory model settings
GEP_IPCW_WEIGHT_CAP_PROB <- 0.99 # Cap extreme inverse-censoring weights at the 99th percentile

# Deterministic nested cross-validation for the Objective 4 no-GEP analysis
GEP_EXPLORATORY_CV_SEED <- 20260820L
GEP_EXPLORATORY_CV_REPEATS <- 20L
GEP_EXPLORATORY_OUTER_FOLDS <- 5L
GEP_EXPLORATORY_INNER_FOLDS <- 5L

# Sample size requirement for any analysis
GEP_MIN_SAMPLE_SIZE <- 20

GEP_DEFINITIVE_GROUP_LEVELS <- c("Class 1", "Class 2")
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
# Groups with a laboratory result that is not a definitive Class 1/2 call with PRAME
# reported; they are drawn de-emphasised in survival figures.
GEP_GROUPS_WITHOUT_DEFINITIVE_CALL <- c(
    GEP_GROUP_CLASS_PRAME_MISSING,
    GEP_GROUP_DISCORDANT,
    GEP_GROUP_TEST_FAILED
)
# Class x PRAME detail used for the `biopsy1_gep` display factor: the four validation
# groups, then every other availability group under its own label.
GEP_DETAIL_LEVELS <- c(
    GEP_VALIDATION_GROUP_LEVELS,
    setdiff(GEP_GROUP_LEVELS, GEP_DEFINITIVE_GROUP_LEVELS)
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
GEP_TEST_TIMING_BEFORE <- sprintf("More than %d days before treatment", GEP_TEST_TIMING_WINDOW_DAYS)
GEP_TEST_TIMING_AT_TREATMENT <- sprintf("Within %d days of treatment", GEP_TEST_TIMING_WINDOW_DAYS)
GEP_TEST_TIMING_AFTER <- sprintf("%d to %d days after treatment", GEP_TEST_TIMING_WINDOW_DAYS, GEP_TEST_TIMING_LATE_DAYS)
GEP_TEST_TIMING_LATE <- sprintf("More than %d days after treatment", GEP_TEST_TIMING_LATE_DAYS)
GEP_TEST_TIMING_NO_DATE <- "No biopsy date"
GEP_TEST_TIMING_LEVELS <- c(
    GEP_TEST_TIMING_BEFORE,
    GEP_TEST_TIMING_AT_TREATMENT,
    GEP_TEST_TIMING_AFTER,
    GEP_TEST_TIMING_LATE,
    GEP_TEST_TIMING_NO_DATE
)
# A biopsy this many days before treatment (or earlier) may be a date-entry error or a
# long-standing earlier sample; the patient stays in the analysis and is listed for review.
GEP_TEST_REVIEW_BEFORE_DAYS <- 365
GEP_TEST_REVIEW_BEFORE_REASON <- sprintf("More than %d days before treatment: confirm the biopsy date", GEP_TEST_REVIEW_BEFORE_DAYS)

#' Classify the biopsy date relative to treatment
#'
#' @param days_from_treatment Numeric days from treatment to biopsy (negative = before).
#' @return Character vector of timing categories (`GEP_TEST_TIMING_LEVELS`).
classify_gep_test_timing <- function(days_from_treatment) {
    dplyr::case_when(
        is.na(days_from_treatment) ~ GEP_TEST_TIMING_NO_DATE,
        days_from_treatment < -GEP_TEST_TIMING_WINDOW_DAYS ~ GEP_TEST_TIMING_BEFORE,
        days_from_treatment <= GEP_TEST_TIMING_WINDOW_DAYS ~ GEP_TEST_TIMING_AT_TREATMENT,
        days_from_treatment <= GEP_TEST_TIMING_LATE_DAYS ~ GEP_TEST_TIMING_AFTER,
        TRUE ~ GEP_TEST_TIMING_LATE
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

#' Classify availability group and PRAME status into the Class x PRAME detail label
#'
#' Single mapping behind the `biopsy1_gep` display factor and the four validation
#' groups. Class 1 and Class 2 patients carry PRAME by definition of the group, so
#' they split into `GEP_VALIDATION_GROUP_LEVELS`; every other group keeps its own
#' label.
#'
#' @param gep_group Character or factor vector from `classify_gep_group()`.
#' @param prame_status Character or factor vector with `Negative`/`Positive` for
#'   Class 1/2 patients.
#' @return Character vector of `GEP_DETAIL_LEVELS`. A Class 1/2 patient without a
#'   `Negative`/`Positive` PRAME status returns `NA` so it surfaces in validation.
classify_gep_detail <- function(gep_group, prame_status) {
    gep_group <- as.character(gep_group)
    prame_status <- as.character(prame_status)
    dplyr::case_when(
        gep_group == "Class 1" & prame_status %in% "Negative" ~ GEP_VALIDATION_GROUP_LEVELS[[1]],
        gep_group == "Class 1" & prame_status %in% "Positive" ~ GEP_VALIDATION_GROUP_LEVELS[[2]],
        gep_group == "Class 2" & prame_status %in% "Negative" ~ GEP_VALIDATION_GROUP_LEVELS[[3]],
        gep_group == "Class 2" & prame_status %in% "Positive" ~ GEP_VALIDATION_GROUP_LEVELS[[4]],
        gep_group %in% setdiff(GEP_GROUP_LEVELS, GEP_DEFINITIVE_GROUP_LEVELS) ~ gep_group,
        TRUE ~ NA_character_
    )
}

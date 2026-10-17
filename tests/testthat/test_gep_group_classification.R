test_that("every raw GEP source level maps to exactly one availability group", {
    raw_levels <- c(
        GEP_CLASS_1_DEFINITIVE_RAW_LEVELS,
        GEP_CLASS_2_DEFINITIVE_RAW_LEVELS,
        GEP_CLASS_KNOWN_PRAME_MISSING_RAW_LEVELS,
        GEP_DISCORDANT_RAW_LEVELS,
        GEP_TEST_FAILED_RAW_LEVELS,
        GEP_NOT_TESTED_RAW_LEVELS
    )

    # The named vectors partition the raw levels: no level is listed twice.
    expect_false(anyDuplicated(raw_levels) > 0)

    groups <- classify_gep_group(raw_levels)
    expect_false(anyNA(groups))
    expect_true(all(groups %in% GEP_GROUP_LEVELS))
    expect_setequal(unique(groups), GEP_GROUP_LEVELS)
})

test_that("GEP group classification separates class-known, discordant, and failed results", {
    expect_identical(
        classify_gep_group(c(
            "Class_1A_PRAME_negative", "Class_1B_PRAME_positive",
            "Class_2_PRAME_negative", "Class_2_PRAME_positive",
            "Class_1A_PRAME_not_reported", "Class_2_PRAME_not_reported", "Class_2_PRAME_Unknown",
            "Class_1A_PRAME_discordant",
            "Failed", "Unknown",
            "No", NA_character_
        )),
        c(
            "Class 1", "Class 1",
            "Class 2", "Class 2",
            GEP_GROUP_CLASS_PRAME_MISSING, GEP_GROUP_CLASS_PRAME_MISSING, GEP_GROUP_CLASS_PRAME_MISSING,
            GEP_GROUP_DISCORDANT,
            GEP_GROUP_TEST_FAILED, GEP_GROUP_TEST_FAILED,
            GEP_GROUP_NOT_TESTED, GEP_GROUP_NOT_TESTED
        )
    )
    # Factor input (as stored in the cohort) classifies the same way.
    expect_identical(
        classify_gep_group(factor("Class_2_PRAME_Unknown")),
        GEP_GROUP_CLASS_PRAME_MISSING
    )
    # An unrecognised source string is surfaced as NA rather than absorbed.
    expect_true(is.na(classify_gep_group("Class_3_PRAME_negative")))
})

test_that("reader-facing GEP group labels use current Class 1/Class 2 terminology only", {
    expect_false(any(grepl("1A|1B", GEP_GROUP_LEVELS)))
    expect_false(any(grepl("1A|1B", unlist(get_objective0_contract_allowed_values("gep_group")))))
})

test_that("definitive GEP groups are exactly the Class 1/Class 2 calls with PRAME reported", {
    raw <- c(
        GEP_CLASS_1_DEFINITIVE_RAW_LEVELS,
        GEP_CLASS_2_DEFINITIVE_RAW_LEVELS,
        GEP_CLASS_KNOWN_PRAME_MISSING_RAW_LEVELS,
        GEP_DISCORDANT_RAW_LEVELS,
        GEP_TEST_FAILED_RAW_LEVELS,
        GEP_NOT_TESTED_RAW_LEVELS
    )
    group <- classify_gep_group(raw)
    expect_identical(group %in% GEP_DEFINITIVE_GROUP_LEVELS, raw %in% GEP_DEFINITIVE_RAW_LEVELS)
})

test_that("the discordant source string is normalized once and classified as its own group", {
    expect_identical(normalize_gep_source_label(GEP_DISCORDANT_SOURCE_LABEL), GEP_DISCORDANT_RAW_LEVELS)
    expect_identical(normalize_gep_source_label(factor("Class_2_PRAME_negative")), "Class_2_PRAME_negative")
    expect_true(is.na(classify_gep_group(GEP_DISCORDANT_SOURCE_LABEL)))
    expect_identical(
        classify_gep_group(normalize_gep_source_label(GEP_DISCORDANT_SOURCE_LABEL)),
        GEP_GROUP_DISCORDANT
    )

    # The derivation and cohort creation both go through the one normalizer.
    derivation_text <- paste(readLines(here::here("scripts", "data_helper", "data_derivation.R"), warn = FALSE), collapse = "\n")
    creation_text <- paste(readLines(here::here("scripts", "data_helper", "cohort_creation.R"), warn = FALSE), collapse = "\n")
    expect_true(grepl("classify_gep_group(normalize_gep_source_label(biopsy1_gep))", derivation_text, fixed = TRUE))
    expect_true(grepl("normalize_gep_source_label(biopsy1_gep)", creation_text, fixed = TRUE))
    # The source string is defined once, in the policy, not retyped in the derivation code.
    expect_false(grepl("DISCORDANT CASTLE RESULTS", paste(derivation_text, creation_text), fixed = TRUE))
})

test_that("the Class x PRAME detail label splits only the definitive groups", {
    expect_identical(
        classify_gep_detail(
            c("Class 1", "Class 1", "Class 2", "Class 2", GEP_GROUP_CLASS_PRAME_MISSING, GEP_GROUP_DISCORDANT, GEP_GROUP_TEST_FAILED, GEP_GROUP_NOT_TESTED),
            c("Negative", "Positive", "Negative", "Positive", "Unknown", "Not Available", "Not Available", "Not Available")
        ),
        GEP_DETAIL_LEVELS
    )
    # A definitive class without a reported PRAME status surfaces as NA, never as a group.
    expect_true(is.na(classify_gep_detail("Class 1", "Unknown")))
    expect_true(is.na(classify_gep_detail(NA_character_, "Negative")))

    # The detail levels are the validation groups plus every non-definitive availability group.
    expect_identical(GEP_DETAIL_LEVELS[seq_along(GEP_VALIDATION_GROUP_LEVELS)], GEP_VALIDATION_GROUP_LEVELS)
    expect_setequal(
        GEP_DETAIL_LEVELS[-seq_along(GEP_VALIDATION_GROUP_LEVELS)],
        setdiff(GEP_GROUP_LEVELS, GEP_DEFINITIVE_GROUP_LEVELS)
    )
    expect_identical(
        GEP_GROUPS_WITHOUT_DEFINITIVE_CALL,
        setdiff(GEP_GROUP_LEVELS, c(GEP_DEFINITIVE_GROUP_LEVELS, GEP_GROUP_NOT_TESTED))
    )
})

test_that("live code carries one GEP label set: no pooled failed/indeterminate label or retired variables", {
    script_files <- list.files(here::here("scripts"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
    script_text <- vapply(script_files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1))
    retired <- c(
        "GEP Failed/Indeterminate", "GEP Not Tested", "gep_class_simple", "derive_gep_class_simple",
        "GEP_FAILED_OR_INDETERMINATE_RAW_LEVELS", "GEP_INVALID_ANALYSIS_LABELS", "GEP_DEFINITIVE_SIMPLE_LEVELS",
        "GEP_GROUPING_SPECS", "OBJECTIVE4_GEP_GROUPING", "refresh_gep_analysis_flags",
        "Class 1 PRAME Negative", "Class 2 PRAME Positive"
    )
    for (term in retired) {
        offenders <- basename(script_files)[vapply(script_text, grepl, logical(1), pattern = term, fixed = TRUE)]
        expect_identical(offenders, character(0), info = term)
    }
})

test_that("biopsy timing categories use the documented windows and sign", {
    days <- c(-2197, -32, -31, 0, 31, 32, 180, 181, 2627, NA)
    expect_identical(
        classify_gep_test_timing(days),
        c(
            GEP_TEST_TIMING_LEVELS[[1]], GEP_TEST_TIMING_LEVELS[[1]],
            GEP_TEST_TIMING_LEVELS[[2]], GEP_TEST_TIMING_LEVELS[[2]], GEP_TEST_TIMING_LEVELS[[2]],
            GEP_TEST_TIMING_LEVELS[[3]], GEP_TEST_TIMING_LEVELS[[3]],
            GEP_TEST_TIMING_LEVELS[[4]], GEP_TEST_TIMING_LEVELS[[4]],
            GEP_TEST_TIMING_LEVELS[[5]]
        )
    )
    # Days are biopsy minus treatment: a biopsy long before treatment is negative and
    # must not be classified as late.
    expect_false(identical(classify_gep_test_timing(-400), GEP_TEST_TIMING_LATE))
    expect_identical(classify_gep_test_timing(400), GEP_TEST_TIMING_LATE)
    expect_identical(GEP_UNKNOWN_CAUSE_OF_DEATH_LEVELS, c("Unknown", "Unknown_Cancer"))
})

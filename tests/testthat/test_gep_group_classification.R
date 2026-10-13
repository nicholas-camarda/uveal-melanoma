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
    expect_setequal(
        GEP_FAILED_OR_INDETERMINATE_RAW_LEVELS,
        c(
            GEP_TEST_FAILED_RAW_LEVELS,
            GEP_CLASS_KNOWN_PRAME_MISSING_RAW_LEVELS,
            GEP_DISCORDANT_RAW_LEVELS
        )
    )

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

test_that("definitive GEP groups reproduce the primary validation set definition", {
    raw <- c(
        GEP_CLASS_1_DEFINITIVE_RAW_LEVELS,
        GEP_CLASS_2_DEFINITIVE_RAW_LEVELS,
        GEP_FAILED_OR_INDETERMINATE_RAW_LEVELS,
        GEP_NOT_TESTED_RAW_LEVELS
    )
    simple <- dplyr::case_when(
        raw %in% GEP_CLASS_1_DEFINITIVE_RAW_LEVELS ~ "Class 1",
        raw %in% GEP_CLASS_2_DEFINITIVE_RAW_LEVELS ~ "Class 2",
        TRUE ~ NA_character_
    )
    group <- classify_gep_group(raw)
    # Class 1/Class 2 in gep_group are exactly the definitive calls with PRAME reported.
    expect_identical(group %in% GEP_DEFINITIVE_SIMPLE_LEVELS, !is.na(simple))
})

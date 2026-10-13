test_that("Objective 3 PFS-2 endpoint contract is explicit and isolated from Objective 1 PFS", {
    source_text <- paste(
        readLines(here::here("scripts", "data_helper", "data_derivation.R"), warn = FALSE),
        collapse = "\n"
    )

    expect_true(grepl("pfs2_second_recurrence_observed", source_text, fixed = TRUE))
    expect_true(grepl("recurrence2", source_text, fixed = TRUE))
    expect_true(grepl("pfs2_event", source_text, fixed = TRUE))
    expect_true(grepl("tt_pfs2_months", source_text, fixed = TRUE))
    expect_false(grepl("pfs2_event\\s*=\\s*if_else\\([^\\n]*(mets_event|death_event)", source_text))
})

test_that("Objective 4 GEP MFS/MSS endpoints remain separate from Objective 1 PFS", {
    validation_text <- paste(
        readLines(here::here("scripts", "gep", "validation", "gep_definitive_validation.R"), warn = FALSE),
        collapse = "\n"
    )

    # The validation reads the canonical MFS/MSS follow-up and event-type fields.
    expect_true(grepl('time = "tt_mets_months_analysis"', validation_text, fixed = TRUE))
    expect_true(grepl('status = "mfs_event_type"', validation_text, fixed = TRUE))
    expect_true(grepl('time = "tt_death_months"', validation_text, fixed = TRUE))
    expect_true(grepl('status = "mss_event_type"', validation_text, fixed = TRUE))
    expect_false(grepl("pfs_event|tt_pfs_months", validation_text))
})

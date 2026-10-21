test_that("project R sources keep non-ASCII text out of code", {
    files <- list.files(
        c(here::here("scripts"), here::here("tests")),
        pattern = "[.]R$", recursive = TRUE, full.names = TRUE
    )
    offenders <- character()
    for (file in files) {
        parsed <- utils::getParseData(parse(file, keep.source = TRUE))
        code <- parsed$text[parsed$terminal & parsed$token != "COMMENT"]
        if (any(grepl("[\x80-\xff]", code, perl = TRUE, useBytes = TRUE))) {
            offenders <- c(offenders, file)
        }
    }

    # A literal such as the "greater than or equal" sign is read in the session's
    # encoding; written as ≥ it is the same string in every locale.
    expect_identical(offenders, character())
})

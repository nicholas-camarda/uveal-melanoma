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

    # A literal such as the "greater than or equal" sign is read in the
    # session's encoding; written as ≥ it is the same string in every locale.
    expect_identical(offenders, character())
})

test_that("loading the project sets a UTF-8 session locale", {
    # Label lookups on text such as "≥ 80 years" depend on a UTF-8 session.
    expect_true(isTRUE(l10n_info()[["UTF-8"]]))
})

test_that("loading the project fixes one locale whatever locale R started in", {
    # testthat resets LC_COLLATE to C inside every test, so the contract is
    # checked in fresh processes started in locales that differ from the
    # project locale. Under C collation capitals sort before lower case
    # (A, B, a, b); the published outputs were sorted in en_US.UTF-8 order.
    script <- withr::local_tempfile(fileext = ".R")
    writeLines(
        c(
            sprintf("source(%s)", dQuote(here::here("scripts", "load_all.R"), FALSE)),
            "cat(Sys.getlocale('LC_CTYPE'), Sys.getlocale('LC_COLLATE'),",
            "    paste(sort(c('b', 'A', 'a', 'B')), collapse = ''), '\\n')"
        ),
        script
    )
    for (start_locale in c("C", "C.UTF-8")) {
        output <- suppressWarnings(system2(
            file.path(R.home("bin"), "Rscript"),
            script,
            stdout = TRUE,
            stderr = FALSE,
            env = c(
                paste0("LC_ALL=", start_locale),
                "OCULAR_REQUIRE_TEST_MANIFEST=false"
            )
        ))
        expect_identical(
            trimws(tail(output, 1L)),
            "en_US.UTF-8 en_US.UTF-8 aAbB",
            info = paste("started in", start_locale)
        )
    }
})

test_that("a missing project locale stops loading with instructions", {
    # Run the locale block of load_all.R with a locale name no computer has,
    # in a fresh process, and check the message tells the reader what to do.
    load_all_lines <- readLines(here::here("scripts", "load_all.R"), warn = FALSE)
    block_end <- grep("ANALYSIS SETTINGS", load_all_lines, fixed = TRUE)[[1]] - 2L
    block <- load_all_lines[seq_len(block_end)]
    block <- sub("^PROJECT_LOCALE <- .*$", 'PROJECT_LOCALE <- "xx_XX.UTF-8"', block)
    script <- withr::local_tempfile(fileext = ".R")
    writeLines(block, script)

    output <- suppressWarnings(system2(
        file.path(R.home("bin"), "Rscript"),
        script,
        stdout = TRUE,
        stderr = TRUE
    ))
    message_text <- paste(output, collapse = "\n")

    expect_gt(attr(output, "status") %||% 0L, 0L)
    expect_match(message_text, "LOCALE ERROR", fixed = TRUE)
    expect_match(message_text, "xx_XX.UTF-8", fixed = TRUE)
    expect_match(message_text, "locale -a", fixed = TRUE)
    expect_match(message_text, "locale-gen", fixed = TRUE)
    expect_match(message_text, "scripts/load_all.R", fixed = TRUE)
})

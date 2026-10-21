source(here::here("scripts", "ci", "run_integration_suite.R"))

#' Create a stand-in repository whose test runner records what it receives
#'
#' @param exit_status Exit status the stub `run_testthat.R` returns.
#' @return A list with the stand-in repository root and the record file path.
make_stub_runner_repo <- function(exit_status = 0L) {
    repo_root <- withr::local_tempdir(.local_envir = parent.frame())
    record_path <- file.path(repo_root, "record.txt")
    dir.create(file.path(repo_root, "scripts", "ci"), recursive = TRUE)
    writeLines(
        c(
            sprintf("record <- %s", dQuote(record_path, FALSE)),
            "writeLines(c(",
            "    commandArgs(trailingOnly = TRUE),",
            "    Sys.getenv('OCULAR_INTEGRATION_RAW_DATA_DIR'),",
            "    Sys.getenv('OCULAR_INTEGRATION_PROCESSED_DATA_DIR')",
            "), record)",
            sprintf("quit(status = %d)", exit_status)
        ),
        file.path(repo_root, "scripts", "ci", "run_testthat.R")
    )
    list(repo_root = repo_root, record_path = record_path)
}

test_that("integration runner names both inputs for the child process only", {
    withr::local_envvar(
        OCULAR_INTEGRATION_RAW_DATA_DIR = NA,
        OCULAR_INTEGRATION_PROCESSED_DATA_DIR = NA
    )
    child <- build_integration_suite_command("/raw", "/processed", "/repo")

    expect_identical(
        child$env,
        c(
            "OCULAR_INTEGRATION_RAW_DATA_DIR='/raw'",
            "OCULAR_INTEGRATION_PROCESSED_DATA_DIR='/processed'"
        )
    )
    expect_identical(
        child$args,
        shQuote(c("/repo/scripts/ci/run_testthat.R", "/repo/tests/integration"))
    )
    expect_identical(Sys.getenv("OCULAR_INTEGRATION_RAW_DATA_DIR"), "")
    expect_identical(Sys.getenv("OCULAR_INTEGRATION_PROCESSED_DATA_DIR"), "")
})

test_that("integration runner passes inputs and the test directory to run_testthat", {
    stub <- make_stub_runner_repo(0L)
    raw_dir <- withr::local_tempdir()
    processed_dir <- withr::local_tempdir()

    expect_true(run_integration_suite(raw_dir, processed_dir, stub$repo_root))

    received <- readLines(stub$record_path)
    expect_identical(received[[1]], file.path(stub$repo_root, "tests", "integration"))
    expect_identical(received[[2]], raw_dir)
    expect_identical(received[[3]], processed_dir)
})

test_that("integration runner handles data directories containing spaces", {
    stub <- make_stub_runner_repo(0L)
    parent <- withr::local_tempdir()
    raw_dir <- file.path(parent, "Project Vault", "Original Files")
    processed_dir <- file.path(parent, "Analytic Dataset")
    dir.create(raw_dir, recursive = TRUE)
    dir.create(processed_dir)

    expect_true(run_integration_suite(raw_dir, processed_dir, stub$repo_root))

    received <- readLines(stub$record_path)
    expect_identical(received[[2]], raw_dir)
    expect_identical(received[[3]], processed_dir)
})

test_that("integration runner fails when the integration tests fail", {
    stub <- make_stub_runner_repo(1L)

    expect_error(
        run_integration_suite(
            withr::local_tempdir(),
            withr::local_tempdir(),
            stub$repo_root
        ),
        "integration suite failed"
    )
})

test_that("integration runner fails before testing when the configured raw input is missing", {
    withr::local_envvar(RAW_DATA_DIR = withr::local_tempdir())

    expect_error(
        run_integration_suite(repo_root = here::here()),
        "Required raw input path checks failed"
    )
})

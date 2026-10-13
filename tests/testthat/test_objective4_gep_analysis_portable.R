build_objective4_output_dirs <- function(test_output_dir) {
    list(
        obj4_mfs = file.path(test_output_dir, "04_GEP_Validation", "a_mfs_validation"),
        obj4_mss = file.path(test_output_dir, "04_GEP_Validation", "b_mss_validation"),
        obj4_ph_diagnostics = file.path(test_output_dir, "04_GEP_Validation", "c_proportional_hazards_diagnostics")
    )
}

test_that("Two-cohort MFS poster KM stack writer saves combined panel", {
    output_dir <- file.path(TEST_OUTPUT_DIR, "objective4_simple_binary_poster_km_stack")
    withr::defer(unlink(output_dir, recursive = TRUE), envir = parent.frame())

    paths <- write_mfs_simple_binary_poster_km_stack(
        full_data = create_test_dataset(),
        gksrs_data = create_test_dataset(),
        output_dir = output_dir,
        filename = "test_stack.png"
    )

    expect_true(file.exists(paths$png))
    expect_equal(basename(paths$png), "test_stack.png")
})

test_that("Objective run state fails closed for fatal issues", {
    expect_identical(
        determine_run_state(fatal_issues = "mss_validation: forced failure"),
        "failed"
    )
})

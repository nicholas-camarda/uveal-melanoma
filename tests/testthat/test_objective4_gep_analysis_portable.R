build_objective4_output_dirs <- function(test_output_dir) {
    list(
        obj4_gep = file.path(test_output_dir, "04_GEP_Validation"),
        obj4_no_gep = file.path(test_output_dir, "04_GEP_Validation", "d_no_gep")
    )
}

test_that("Objective run state fails closed for fatal issues", {
    expect_identical(
        determine_run_state(fatal_issues = "mss_validation: forced failure"),
        "failed"
    )
})

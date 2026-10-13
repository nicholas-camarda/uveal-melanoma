#' Run Objective 4: GEP Validation
#'
#' Validates the imported DecisionDx class-level 5-year predictions in patients with
#' a definitive Class 1/2 call and PRAME reported, then runs the exploratory no-GEP
#' report (full cohort only) and the MFS sensitivity summary. A failing component
#' stops the objective; the caller records the failure.
#'
#' @param data Data frame with analytic dataset
#' @param dataset_name Character string for dataset name
#' @param output_dirs List of output directories
#' @param prefix Character string for file prefix
#' @param confounders Unused by Objective 4; kept so every objective runner shares one interface.
#' @return List with `validation` (results of `run_gep_definitive_validation()`),
#'   `output_files`, `exploratory_no_gep_results`, and `mfs_sensitivity_results`.
run_objective_4 <- function(data, dataset_name, output_dirs, prefix, confounders = NULL) {
    step4_start_time <- Sys.time()
    display_name <- tools::toTitleCase(gsub("_", " ", gsub("uveal_melanoma_|_cohort", "", dataset_name)))
    log_phase(paste("STEP 4: GEP PREDICTIVE ACCURACY VALIDATION", display_name, sep = " - "))

    logger::log_info(formatted("Executing run_gep_definitive_validation: 5-year observed vs imported risk", indent = 1))
    validation <- run_gep_definitive_validation(data)
    output_files <- write_gep_validation_outputs(validation, data, dataset_name, output_dirs, prefix)
    logger::log_info(formatted("Definitive-GEP validation completed", indent = 1))

    exploratory_no_gep_results <- NULL
    if (identical(dataset_name, "uveal_melanoma_full_cohort")) {
        logger::log_info(formatted("Executing exploratory no-GEP summary integration for the full cohort", indent = 1))
        no_gep_output_dir <- file.path(dirname(output_dirs$obj4_mfs), "d_exploratory_no_gep")
        collected_no_gep_results <- collect_exploratory_no_gep_analysis(
            data = data,
            dataset_name = dataset_name,
            verify_km_fix = FALSE
        )
        exploratory_no_gep_results <- run_exploratory_no_gep_report(
            dataset_name = dataset_name,
            output_dir = no_gep_output_dir,
            verify_km_fix = FALSE,
            data = data,
            collected_results = collected_no_gep_results
        )
    }

    logger::log_info(formatted("Executing Objective 4 MFS sensitivity summary", indent = 1))
    mfs_sensitivity_results <- run_objective4_mfs_sensitivity_summary(
        data = data,
        dataset_name = dataset_name,
        output_dir = file.path(dirname(output_dirs$obj4_mfs), "unified_summary"),
        prefix = prefix
    )

    logger::log_info(sprintf(">>> COMPLETED %s (Duration: %.1f seconds)",
        "STEP 4: GEP PREDICTIVE ACCURACY VALIDATION",
        as.numeric(difftime(Sys.time(), step4_start_time, units = "secs"))
    ))

    list(
        validation = validation,
        output_files = output_files,
        mfs_sensitivity_results = mfs_sensitivity_results,
        exploratory_no_gep_results = exploratory_no_gep_results
    )
}

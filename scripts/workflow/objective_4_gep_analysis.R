#' Run Objective 4: GEP Validation
#'
#' Validates the imported DecisionDx class-level 5-year predictions in patients with
#' a definitive Class 1/2 call and PRAME reported, then (full cohort only) describes
#' outcomes in every GEP availability group and evaluates baseline clinical features
#' in patients without a definitive call. A failing component stops the objective;
#' the caller records the failure.
#'
#' @param data Data frame with analytic dataset
#' @param dataset_name Character string for dataset name
#' @param output_dirs List of output directories
#' @param prefix Character string for file prefix
#' @param confounders Unused by Objective 4; kept so every objective runner shares one interface.
#' @return List with `validation` (results of `run_gep_definitive_validation()`),
#'   `output_files`, `no_gep_results`, and `no_gep_files` (the last two are `NULL`
#'   outside the full cohort).
run_objective_4 <- function(data, dataset_name, output_dirs, prefix, confounders = NULL) {
    step4_start_time <- Sys.time()
    display_name <- tools::toTitleCase(gsub("_", " ", gsub("uveal_melanoma_|_cohort", "", dataset_name)))
    log_phase(paste("STEP 4: GEP PREDICTIVE ACCURACY VALIDATION", display_name, sep = " - "))

    logger::log_info(formatted("Executing run_gep_definitive_validation: 5-year observed vs imported risk", indent = 1))
    validation <- run_gep_definitive_validation(data)
    output_files <- write_gep_validation_outputs(validation, data, dataset_name, output_dirs, prefix)
    logger::log_info(formatted("Definitive-GEP validation completed", indent = 1))

    no_gep_results <- NULL
    no_gep_files <- NULL
    if (identical(dataset_name, GEP_NO_GEP_DATASET)) {
        logger::log_info(formatted("Executing run_no_gep_analysis for the full cohort", indent = 1))
        no_gep_results <- run_no_gep_analysis(data)
        no_gep_files <- write_no_gep_outputs(
            no_gep_results, no_gep_results$prepared,
            output_dir = get_output_dir(output_dirs, "obj4_no_gep"), prefix = prefix
        )
        logger::log_info(formatted("No-GEP analysis completed", indent = 1))
    }

    logger::log_info(sprintf(">>> COMPLETED %s (Duration: %.1f seconds)",
        "STEP 4: GEP PREDICTIVE ACCURACY VALIDATION",
        as.numeric(difftime(Sys.time(), step4_start_time, units = "secs"))
    ))

    list(
        validation = validation,
        output_files = output_files,
        no_gep_results = no_gep_results,
        no_gep_files = no_gep_files
    )
}

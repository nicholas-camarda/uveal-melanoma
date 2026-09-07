test_that("ridge reporting preserves original-scale coefficients and their signs", {
    set.seed(42)
    data <- data.frame(
        age = seq(20, 80, length.out = 120),
        sex = factor(rep(c("Male", "Female"), 60), levels = c("Male", "Female"))
    )
    x <- model.matrix(~ age + sex, data)[, -1]
    y <- rbinom(120, 1, plogis(-3 + 0.035 * x[, "age"] - x[, "sexFemale"]))
    fit <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = 0,
        foldid = rep(1:5, 24), standardize = TRUE)
    coefficients <- extract_binary_model_coefficients(fit, c("age", "sex"), data)
    reported <- summarize_predictor_contributions(coefficients, "Fixture")

    # A rescaling or sign loss at the reporting boundary would no longer
    # reconstruct the fitted model on the original input units.
    expect_true(all(c("coefficient", "abs_coefficient") %in% names(reported)))
    expect_equal(reported$coefficient, coefficients$estimate[match(reported$dominant_term, coefficients$term)])
    expect_equal(reported$abs_coefficient, abs(reported$coefficient))
    expect_identical(unname(reported$reference_level[reported$predictor == "sex"]), "Male")
    expect_identical(reported$dominant_term[reported$predictor == "sex"], "sexFemale")
    beta <- coefficients$estimate[match(colnames(x), coefficients$term)]
    intercept <- coefficients$estimate[coefficients$term == "(Intercept)"]
    expect_equal(as.numeric(plogis(intercept + x %*% beta)),
        as.numeric(predict(fit, newx = x, s = "lambda.min", type = "response")))
})

test_that("coefficient plots retain all predictors and the actual coded factor levels", {
    model <- "Direct 5-Year MFS Risk"
    coefficients <- tibble::tibble(
        term = c("sexFemale", "srfNo", "initial_t_stage_simpleT3", "initial_tumor_diameter", "age_at_diagnosis", "optic_nerve_involvement"),
        predictor = c("sex", "srf", "initial_t_stage_simple", "initial_tumor_diameter", "age_at_diagnosis", "optic_nerve_involvement"),
        estimate = c(-0.6, -0.4, 0.5, 0.1, 0.01, -0.2),
        reference_level = c("Male", "Yes", "T1", NA, NA, NA)
    )
    contribution <- summarize_predictor_contributions(coefficients, model) %>%
        dplyr::mutate(section = "model_contribution")
    path <- tempfile(fileext = ".png")
    withr::defer(unlink(path))
    result <- create_exploratory_no_gep_direct_model_contributions_plot(
        contribution, model, path, return_plot = TRUE)
    expect_setequal(result$plot_data$dominant_term, coefficients$term)
    expect_equal(result$plot_data$coefficient,
        coefficients$estimate[match(result$plot_data$dominant_term, coefficients$term)])
    labels <- as.character(result$plot_data$display_label)
    expect_true(any(grepl("Female", labels)))
    expect_true(any(grepl("Female vs Male", labels)))
    expect_true(any(grepl("No", labels)))
    expect_true(any(grepl("T3", labels)))
    expect_true(any(grepl("mm", labels)))
    expect_true(any(grepl("year", labels)))
    expect_match(result$plot$labels$y, "original", ignore.case = TRUE)
    expect_match(result$plot$labels$caption, "not.*importance", ignore.case = TRUE)
    expect_true(file.exists(path))
})

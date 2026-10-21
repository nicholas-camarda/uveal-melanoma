#' Build a deterministic table whose simulated Fisher p-value depends on the RNG
#'
#' `fisher.test(simulate.p.value = TRUE)` only simulates for tables larger than
#' 2 x 2, so the fixture crosses two treatment groups with a sparse four-level
#' and a sparse three-level categorical variable. No random numbers are used to
#' build the data.
#'
#' @return Tibble with `treatment_group`, `grade` and `site`.
create_simulated_test_fixture <- function() {
    tibble::tibble(
        treatment_group = factor(
            rep(c("PBT", "GKSRS"), times = c(40L, 30L)),
            levels = c("PBT", "GKSRS")
        ),
        grade = factor(
            c(
                rep(c("A", "B", "C", "D"), times = c(22L, 12L, 5L, 1L)),
                rep(c("A", "B", "C", "D"), times = c(10L, 9L, 7L, 4L))
            ),
            levels = c("A", "B", "C", "D")
        ),
        site = factor(
            c(
                rep(c("x", "y", "z"), times = c(25L, 13L, 2L)),
                rep(c("x", "y", "z"), times = c(14L, 11L, 5L))
            ),
            levels = c("x", "y", "z")
        )
    )
}

#' Extract the non-missing p-values of a gtsummary table in table order
#'
#' @param tbl A gtsummary table with a p-value column.
#' @return Numeric vector of p-values.
extract_table_p_values <- function(tbl) {
    p_values <- tbl$table_body$p.value
    p_values[!is.na(p_values)]
}

#' Run an expression after two different amounts of unrelated RNG use
#'
#' Evaluates `fn` once after a short stream and once after a long stream, under
#' different global seeds, as unrelated earlier code would leave it.
#'
#' @param fn Function of no arguments returning the value to compare.
#' @return List with the results after each RNG history.
evaluate_after_unrelated_draws <- function(fn) {
    set.seed(1)
    invisible(stats::runif(3))
    first <- fn()

    set.seed(2024)
    invisible(stats::rnorm(5000))
    second <- fn()

    list(first = first, second = second)
}

test_that("the fixture is sensitive to the random-number stream", {
    fixture <- create_simulated_test_fixture()
    cross_table <- table(fixture$treatment_group, fixture$grade)

    set.seed(1)
    p_stream_one <- stats::fisher.test(cross_table, simulate.p.value = TRUE)$p.value
    set.seed(2)
    p_stream_two <- stats::fisher.test(cross_table, simulate.p.value = TRUE)$p.value

    # If this fails the tests below could pass without the seed doing anything.
    expect_false(identical(p_stream_one, p_stream_two))
})

test_that("with_local_seed gives the same draws and leaves the caller's stream untouched", {
    draws <- evaluate_after_unrelated_draws(function() {
        with_local_seed(SIMULATED_FISHER_SEED, stats::runif(5))
    })
    expect_identical(draws$first, draws$second)

    set.seed(99)
    stream_before <- .Random.seed
    invisible(with_local_seed(SIMULATED_FISHER_SEED, stats::runif(5)))
    expect_identical(.Random.seed, stream_before)

    rm(".Random.seed", envir = globalenv())
    invisible(with_local_seed(SIMULATED_FISHER_SEED, stats::runif(5)))
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("with_local_seed does not depend on the session's generator kinds", {
    default_draws <- with_local_seed(SIMULATED_FISHER_SEED, stats::runif(5))

    previous_kinds <- RNGkind()
    withr::defer(do.call(RNGkind, as.list(previous_kinds)))
    suppressWarnings(RNGkind("Wichmann-Hill", "Box-Muller", "Rounding"))
    other_kind_draws <- suppressWarnings(
        with_local_seed(SIMULATED_FISHER_SEED, stats::runif(5))
    )

    expect_identical(other_kind_draws, default_draws)
    expect_identical(RNGkind(), c("Wichmann-Hill", "Box-Muller", "Rounding"))
})

test_that("Objective 2 simulated Fisher p-values ignore earlier RNG use", {
    fixture <- create_simulated_test_fixture()
    compute_p <- function() {
        calculate_objective2_fisher_p_value(fixture, "treatment_group", "grade")
    }

    p_values <- evaluate_after_unrelated_draws(compute_p)
    expect_identical(p_values$first, p_values$second)
    expect_equal(p_values$first, 0.112443778110945, tolerance = 1e-12)

    set.seed(99)
    stream_before <- .Random.seed
    invisible(compute_p())
    expect_identical(.Random.seed, stream_before)
})

test_that("gtsummary simulated Fisher p-values ignore earlier RNG use", {
    fixture <- create_simulated_test_fixture()
    build_table <- function() {
        fixture %>%
            gtsummary::tbl_summary(by = treatment_group) %>%
            add_p_seeded(
                test = list(gtsummary::all_categorical() ~ "fisher.test"),
                test.args = list(
                    gtsummary::all_categorical() ~ list(simulate.p.value = TRUE)
                )
            )
    }

    p_values <- evaluate_after_unrelated_draws(function() {
        extract_table_p_values(build_table())
    })
    expect_length(p_values$first, 2L)
    expect_identical(p_values$first, p_values$second)
    expect_equal(
        p_values$first,
        c(0.112443778110945, 0.216391804097951),
        tolerance = 1e-12
    )

    set.seed(99)
    stream_before <- .Random.seed
    invisible(build_table())
    expect_identical(.Random.seed, stream_before)
})

test_that("safe_add_p_to_summary applies the simulated-test seed", {
    fixture <- create_simulated_test_fixture()
    build_table <- function() {
        fixture %>%
            dplyr::select(treatment_group, grade) %>%
            gtsummary::tbl_summary(by = treatment_group) %>%
            safe_add_p_to_summary(
                context_label = "RNG independence fixture",
                test = list(gtsummary::all_categorical() ~ "fisher.test"),
                test.args = list(
                    gtsummary::all_categorical() ~ list(simulate.p.value = TRUE)
                )
            )
    }

    p_values <- evaluate_after_unrelated_draws(function() {
        extract_table_p_values(build_table())
    })
    expect_identical(p_values$first, p_values$second)
    expect_equal(p_values$first, 0.112443778110945, tolerance = 1e-12)
})

test_that("tumor-size jitter positions ignore earlier RNG use", {
    tumor_data <- tibble::tibble(
        treatment_group = factor(
            rep(c("PBT", "GKSRS"), each = 12L),
            levels = c("PBT", "GKSRS")
        ),
        initial_tumor_height = c(seq(2, 8.5, length.out = 12L), seq(3, 9.5, length.out = 12L))
    )
    jitter_positions <- function() {
        summary <- summarize_tumor_size_by_treatment(tumor_data)
        built <- ggplot2::ggplot_build(summary$plot)
        # Layer 2 is the point layer; layer 1 is the boxplot.
        built$data[[2L]]$x
    }

    positions <- evaluate_after_unrelated_draws(jitter_positions)
    expect_identical(positions$first, positions$second)
    # The points are jittered at all (not drawn at the integer group positions).
    expect_true(any(positions$first != round(positions$first)))

    set.seed(99)
    stream_before <- .Random.seed
    invisible(jitter_positions())
    expect_identical(.Random.seed, stream_before)
})

#' Walk a parsed expression and visit every call with its enclosing call names
#'
#' @param expr Parsed expression (call, symbol or constant).
#' @param visit Function of `(call, name, ancestors)` invoked for each call node,
#'   where `name` is the called function and `ancestors` are the names of the
#'   enclosing calls.
#' @param ancestors Character vector of enclosing call names.
#' @return `NULL`, invisibly.
walk_calls <- function(expr, visit, ancestors = character()) {
    if (!is.call(expr)) {
        return(invisible(NULL))
    }
    head <- expr[[1L]]
    name <- if (is.symbol(head)) {
        as.character(head)
    } else if (is.call(head) && as.character(head[[1L]]) %in% c("::", ":::")) {
        as.character(head[[3L]])
    } else {
        ""
    }
    visit(expr, name, ancestors)
    for (argument in as.list(expr)[-1L]) {
        if (!missing(argument)) {
            walk_calls(argument, visit, c(ancestors, name))
        }
    }
    invisible(NULL)
}

test_that("every random draw in scripts/ runs under a local seed", {
    seeded_wrappers <- c("with_local_seed", "with_seed", "add_p_seeded", "safe_add_p_to_summary")
    random_draws <- c("sample", "sample.int", "runif", "rnorm", "rbinom", "rexp", "rpois", "rmultinom")
    script_files <- list.files(
        here::here("scripts"),
        pattern = "\\.R$",
        recursive = TRUE,
        full.names = TRUE
    )
    # Definitions of the wrappers themselves are the seeded call sites.
    script_files <- script_files[!grepl("rng_utilities\\.R$", script_files)]

    unseeded <- character()
    scanned_calls <- c(simulated = 0L, jitter = 0L, draws = 0L)

    for (path in script_files) {
        label <- sub(paste0("^", here::here(), "/"), "", path)
        visit <- function(call, name, ancestors) {
            argument_names <- names(call)
            is_seeded <- any(seeded_wrappers %in% c(name, ancestors))
            if (!is.null(argument_names) && "simulate.p.value" %in% argument_names) {
                scanned_calls[["simulated"]] <<- scanned_calls[["simulated"]] + 1L
                if (!is_seeded) {
                    unseeded <<- c(unseeded, paste(label, "simulate.p.value outside a seeded wrapper"))
                }
            }
            if (name == "geom_jitter") {
                scanned_calls[["jitter"]] <<- scanned_calls[["jitter"]] + 1L
                if (!"position" %in% argument_names || any(c("width", "height") %in% argument_names)) {
                    unseeded <<- c(unseeded, paste(label, "geom_jitter without a seeded position_jitter"))
                }
            }
            if (name == "position_jitter") {
                scanned_calls[["jitter"]] <<- scanned_calls[["jitter"]] + 1L
                seed_argument <- if ("seed" %in% argument_names) call[["seed"]] else NULL
                if (is.null(seed_argument) || (is.logical(seed_argument) && is.na(seed_argument))) {
                    unseeded <<- c(unseeded, paste(label, "position_jitter without a fixed seed"))
                }
            }
            if (name %in% random_draws) {
                scanned_calls[["draws"]] <<- scanned_calls[["draws"]] + 1L
                if (!is_seeded) {
                    unseeded <<- c(unseeded, paste(label, name, "outside with_local_seed or withr::with_seed"))
                }
            }
            if (name == "set.seed") {
                unseeded <<- c(unseeded, paste(label, "set.seed on the global stream; use with_local_seed"))
            }
            if (name == "cv.glmnet" && !"foldid" %in% argument_names) {
                unseeded <<- c(unseeded, paste(label, "cv.glmnet without foldid draws its own folds"))
            }
        }
        for (expression in parse(path, keep.source = FALSE)) {
            walk_calls(expression, visit)
        }
    }

    # The scan must actually see the call sites, or it would pass vacuously.
    expect_gte(scanned_calls[["simulated"]], 10L)
    expect_gte(scanned_calls[["jitter"]], 2L)
    expect_gte(scanned_calls[["draws"]], 4L)
    expect_identical(unseeded, character())
})

test_that("cross-validation fold assignment ignores earlier RNG use", {
    strata <- rep(c("event", "no_event"), times = c(12L, 24L))
    stable_id <- sprintf("patient_%02d", seq_along(strata))
    compute_folds <- function() {
        create_deterministic_fold_ids(strata, folds = 3L, seed = GEP_EXPLORATORY_CV_SEED, stable_id = stable_id)
    }

    folds <- evaluate_after_unrelated_draws(compute_folds)
    expect_identical(folds$first, folds$second)

    set.seed(99)
    stream_before <- .Random.seed
    invisible(compute_folds())
    expect_identical(.Random.seed, stream_before)
})

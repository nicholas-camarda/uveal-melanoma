# Random-number utilities
#
# Helpers that scope a fixed seed to one calculation so reported values do not
# depend on the global random-number stream. Seed values are defined in
# scripts/config/rng_policy.R.

#' Evaluate an expression under a fixed, locally scoped random seed
#'
#' Sets `seed` with the default Mersenne-Twister, Inversion and Rejection
#' generators, evaluates `expr`, and then restores the caller's generator kinds
#' and `.Random.seed` (or its absence). The result therefore depends only on
#' `seed` and `expr`, not on earlier random draws or on generator settings
#' changed elsewhere in the session, and the calculation does not advance the
#' caller's stream.
#'
#' @param seed Integer scalar random seed.
#' @param expr Expression to evaluate under the seed.
#' @return The value produced by `expr`.
with_local_seed <- function(seed, expr) {
    withr::with_seed(
        seed,
        expr,
        .rng_kind = "Mersenne-Twister",
        .rng_normal_kind = "Inversion",
        .rng_sample_kind = "Rejection"
    )
}

#' Add p-values to a gtsummary table under the central simulated-test seed
#'
#' Calls `gtsummary::add_p()` with `seed` set for the whole call, so Monte Carlo
#' tests (`fisher.test(simulate.p.value = TRUE)`) return the same p-values for
#' the same table specification and data. All tests in one call share one
#' stream, in table order.
#'
#' @param x A gtsummary table.
#' @param ... Arguments passed to `gtsummary::add_p()`.
#' @param seed Integer scalar random seed; defaults to `SIMULATED_FISHER_SEED`.
#' @return The gtsummary table with a p-value column.
add_p_seeded <- function(x, ..., seed = SIMULATED_FISHER_SEED) {
    with_local_seed(seed, gtsummary::add_p(x, ...))
}

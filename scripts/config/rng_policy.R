# =============================================================================
# RANDOM-NUMBER SEEDS
# =============================================================================
# Every output that depends on random numbers draws them under one of the fixed
# seeds below, scoped to the calculation (see `with_local_seed()` in
# scripts/utils/rng_utilities.R). The seeds make the displayed values a function
# of the data and the code alone: they do not depend on which analysis ran
# earlier, on the run mode, or on the session's global random-number stream.
# Changing a seed changes the affected published numbers once.

# Monte Carlo (simulate.p.value = TRUE) Fisher exact tests for tables larger
# than 2 x 2, in every baseline, outcome and Objective 2 descriptive table.
SIMULATED_FISHER_SEED <- 20260422L

# Jitter of the points drawn over the tumor-size boxplots.
FIGURE_JITTER_SEED <- 20261021L

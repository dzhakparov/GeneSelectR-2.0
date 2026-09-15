# Outcome shuffling provides a gene-specific reference under joint outcome
# independence. The resulting ratios are ranking measurements.


compute_null_selection_frequencies <- function(
    X,
    y,
    B = 20L,
    n_permutations = 20L,
    alpha = 0.5,
    k_folds = 5L,
    random_seed = 123L,
    workers = 1L,
    max_failed_fraction = 0.1,
    verbose = TRUE
) {
    n_genes <- ncol(X)
    null_frequencies <- matrix(
        0,
        nrow = n_permutations, ncol = n_genes,
        dimnames = list(NULL, colnames(X))
    )
    null_contribution <- null_frequencies
    diagnostics <- vector("list", n_permutations)

    for (permutation_index in seq_len(n_permutations)) {
        permutation_seed <- random_seed + 10000L * permutation_index
        shuffled_y <- withr::with_seed(permutation_seed, sample(y))
        splits <- create_subsamples(
            shuffled_y,
            B = B,
            random_seed = permutation_seed,
            k_folds = k_folds
        )
        fits <- .run_resamples(
            X, shuffled_y, splits,
            alpha = alpha, workers = workers,
            base_seed = permutation_seed, verbose = FALSE
        )
        diagnostics[[permutation_index]] <- .check_resample_results(
            fits, max_failed_fraction
        )
        aggregated <- .aggregate_resamples(fits, n_genes)
        contribution <- .compute_predictive_contribution(
            X, splits, aggregated,
            verbose = FALSE
        )
        null_frequencies[permutation_index, ] <- aggregated$recurrence
        null_contribution[permutation_index, ] <- contribution$raw

        if (isTRUE(verbose) &&
            (permutation_index %% 5L == 0L ||
                permutation_index == n_permutations)) {
            message(sprintf(
                "Completed %d of %d shuffled outcomes.",
                permutation_index, n_permutations
            ))
        }
    }

    list(
        null_frequencies = null_frequencies,
        null_contribution = null_contribution,
        diagnostics = .summarise_null_diagnostics(diagnostics)
    )
}


.summarise_null_diagnostics <- function(diagnostics) {
    n_fits <- sum(vapply(
        diagnostics, function(x) x$n_resamples, integer(1)
    ))
    n_fit_failures <- sum(vapply(
        diagnostics, function(x) x$n_fit_failures, integer(1)
    ))
    n_auc_failures <- sum(vapply(
        diagnostics, function(x) x$n_auc_failures, integer(1)
    ))
    list(
        n_permutations = length(diagnostics),
        n_fits = n_fits,
        n_fit_failures = n_fit_failures,
        fit_failure_fraction = n_fit_failures / n_fits,
        n_auc_failures = n_auc_failures,
        auc_failure_fraction = n_auc_failures / n_fits
    )
}

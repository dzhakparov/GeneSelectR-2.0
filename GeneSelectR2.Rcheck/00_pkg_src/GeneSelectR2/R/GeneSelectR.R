#' Fit the GeneSelectR ranking model
#'
#' Repeated five-fold elastic-net models provide two gene-level measurements.
#' Recurrence is the fraction of models with a non-zero coefficient. Predictive
#' contribution combines exact linear SHAP measurements on excluded samples
#' with mutual information. Both measurements are compared with values from
#' models fitted after shuffling the outcome. Their adjusted values receive
#' equal weight in a geometric mean.
#'
#' Biological annotation is performed after ranking with the package's
#' biological interpretation functions.
#'
#' @param X Numeric matrix with samples in rows and named genes in columns.
#' @param y Two-level factor. The second level is the modelled class.
#' @param gene_names Optional gene names. Column names of `X` are used by
#'   default.
#' @param B Number of repeated model fits.
#' @param alpha Elastic-net mixing parameter between zero and one.
#' @param k_folds Number of stratified folds. The benchmark used five.
#' @param permutations Number of shuffled-outcome references.
#' @param null_B Number of repeated fits for each shuffled outcome.
#' @param recurrence_epsilon Smoothing constant for the recurrence ratio.
#' @param contribution_epsilon Optional smoothing constant for the predictive
#'   contribution ratio. The 95th percentile of positive observed values is
#'   used by default.
#' @param ratio_limit Maximum absolute log2 ratio. The benchmark used four.
#' @param workers Number of parallel fitting workers.
#' @param random_seed Random seed, restored when the function returns.
#' @param max_failed_fraction Maximum permitted fraction of failed fits or
#'   excluded-sample evaluations.
#' @param verbose Logical value controlling progress messages.
#'
#' @return A `geneselectr_result` list. `gene_scores` contains the complete
#'   ranking and the separate gene-level measurements. Ratios are ranking
#'   quantities; they are not probabilities, p-values, or false-discovery-rate
#'   estimates.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 12), nrow = 40)
#' colnames(X) <- paste0("gene", seq_len(ncol(X)))
#' y <- factor(
#'     rep(c("control", "case"), each = 20),
#'     levels = c("control", "case")
#' )
#' signal_columns <- seq_len(2L)
#' X[y == "case", signal_columns] <- X[y == "case", signal_columns] + 1
#' fit <- geneselectr2_fit(
#'     X, y,
#'     B = 2, permutations = 1, null_B = 2, verbose = FALSE
#' )
#' head(fit$gene_scores)
#' @export
geneselectr2_fit <- function(
    X,
    y,
    gene_names = NULL,
    B = 50L,
    alpha = 0.5,
    k_folds = 5L,
    permutations = 20L,
    null_B = 20L,
    recurrence_epsilon = 0.01,
    contribution_epsilon = NULL,
    ratio_limit = 4,
    workers = 1L,
    random_seed = 123L,
    max_failed_fraction = 0.1,
    verbose = TRUE
) {
    start_time <- Sys.time()
    withr::local_seed(random_seed)

    inputs <- .prepare_fit_inputs(X, y, gene_names)
    X <- inputs$X
    y <- inputs$y
    gene_names <- inputs$gene_names
    .validate_fit_settings(
        B, alpha, k_folds, permutations, null_B, recurrence_epsilon,
        contribution_epsilon, ratio_limit, workers, random_seed,
        max_failed_fraction
    )

    .fit_message(
        verbose,
        "Fitting %d elastic-net models for %d samples and %d genes.",
        B, nrow(X), ncol(X)
    )
    resampling_start <- Sys.time()
    subsamples <- create_subsamples(
        y,
        B = B, random_seed = random_seed, k_folds = k_folds
    )
    fits <- .run_resamples(
        X, y, subsamples,
        alpha = alpha, workers = workers,
        base_seed = random_seed, verbose = verbose
    )
    diagnostics <- .check_resample_results(fits, max_failed_fraction)
    aggregated <- .aggregate_resamples(fits, ncol(X))
    resampling_seconds <- .elapsed_seconds(resampling_start)

    contribution <- .compute_predictive_contribution(
        X, subsamples, aggregated,
        verbose = verbose
    )

    .fit_message(
        verbose,
        "Fitting %d shuffled-outcome references with %d models each.",
        permutations, min(B, null_B)
    )
    adjustment_start <- Sys.time()
    null <- compute_null_selection_frequencies(
        X, y,
        B = min(B, null_B), n_permutations = permutations,
        alpha = alpha, k_folds = k_folds, random_seed = random_seed,
        workers = workers, max_failed_fraction = max_failed_fraction,
        verbose = verbose
    )
    adjusted <- .adjust_ranking_components(
        aggregated$recurrence,
        contribution$raw,
        null,
        recurrence_epsilon,
        contribution_epsilon,
        ratio_limit
    )
    adjustment_seconds <- .elapsed_seconds(adjustment_start)

    combined_score <- exp(
        0.5 * log(adjusted$recurrence + 1e-10) +
            0.5 * log(adjusted$contribution + 1e-10)
    )
    final_score <- percentile01(combined_score)
    gene_scores <- .build_gene_score_table(
        gene_names, final_score, combined_score, aggregated$recurrence,
        contribution, adjusted
    )
    stability <- compute_nogueira_stability(aggregated$selection_matrix)

    result <- list(
        gene_scores = gene_scores,
        stability = c(
            stability,
            list(selection_matrix = aggregated$selection_matrix)
        ),
        instance_importance = contribution$shap_matrix,
        cv_results = .summarise_auc(fits),
        diagnostics = list(
            resampling = diagnostics,
            permutation = null$diagnostics
        ),
        parameters = list(
            B = as.integer(B),
            k_folds = as.integer(k_folds),
            alpha = alpha,
            permutations = as.integer(permutations),
            null_B = as.integer(min(B, null_B)),
            recurrence_epsilon = recurrence_epsilon,
            contribution_epsilon = adjusted$contribution_epsilon,
            ratio_limit = ratio_limit,
            workers = as.integer(workers),
            random_seed = as.integer(random_seed),
            positive_class = levels(y)[2]
        ),
        timing = list(
            total_seconds = .elapsed_seconds(start_time),
            resampling_seconds = resampling_seconds,
            adjustment_seconds = adjustment_seconds
        )
    )
    class(result) <- "geneselectr_result"
    result
}


.prepare_fit_inputs <- function(X, y, gene_names) {
    if (!is.null(gene_names)) {
        if (!is.character(gene_names) || length(gene_names) != ncol(X)) {
            stop("gene_names must contain one name for each column of X")
        }
        colnames(X) <- gene_names
    }
    validate_inputs(X, y, colnames(X), min_per_class = 5L)
    list(X = X, y = y, gene_names = colnames(X))
}


.validate_fit_settings <- function(
    B, alpha, k_folds, permutations, null_B, recurrence_epsilon,
    contribution_epsilon, ratio_limit, workers, random_seed,
    max_failed_fraction
) {
    .validate_positive_integer(B, "B")
    .validate_positive_integer(k_folds, "k_folds", minimum = 2L)
    .validate_positive_integer(permutations, "permutations")
    .validate_positive_integer(null_B, "null_B")
    .validate_positive_integer(workers, "workers")
    .validate_integer(random_seed, "random_seed")
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
        alpha < 0 || alpha > 1) {
        stop("alpha must be one finite value between zero and one")
    }
    .validate_positive_number(recurrence_epsilon, "recurrence_epsilon")
    if (!is.null(contribution_epsilon)) {
        .validate_positive_number(
            contribution_epsilon, "contribution_epsilon"
        )
    }
    .validate_positive_number(ratio_limit, "ratio_limit")
    if (!is.numeric(max_failed_fraction) ||
        length(max_failed_fraction) != 1L ||
        !is.finite(max_failed_fraction) || max_failed_fraction < 0 ||
        max_failed_fraction >= 1) {
        stop("max_failed_fraction must be one value in [0, 1)")
    }
    invisible(TRUE)
}


.validate_integer <- function(value, name, minimum = NULL) {
    valid <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value == as.integer(value)
    if (!is.null(minimum)) {
        valid <- valid && value >= minimum
    }
    if (!valid) {
        stop(name, " must be an integer", call. = FALSE)
    }
    invisible(TRUE)
}


.validate_positive_integer <- function(value, name, minimum = 1L) {
    .validate_integer(value, name, minimum)
}


.validate_positive_number <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
        stop(name, " must be one positive finite number", call. = FALSE)
    }
    invisible(TRUE)
}


.fit_one_resample <- function(X, y, split, alpha) {
    training <- split$train
    excluded <- split$oob
    fit <- tryCatch(
        fit_regularized_model(X[training, , drop = FALSE], y[training], alpha),
        error = identity
    )
    if (inherits(fit, "error")) {
        return(.failed_resample(ncol(X), conditionMessage(fit)))
    }

    prediction <- tryCatch(
        fit$predict_fn(X[excluded, , drop = FALSE]),
        error = identity
    )
    if (inherits(prediction, "error")) {
        auc <- NA_real_
        auc_error <- conditionMessage(prediction)
    } else {
        auc <- tryCatch(compute_auc(y[excluded], prediction), error = identity)
        auc_error <- if (inherits(auc, "error")) {
            conditionMessage(auc)
        } else {
            NA_character_
        }
        if (inherits(auc, "error")) {
            auc <- NA_real_
        }
    }

    list(
        selected = fit$selected,
        coefficients = fit$full_coef_vector,
        mutual_information = compute_mi_vectorized(
            X[training, , drop = FALSE], y[training]
        ),
        auc = auc,
        fit_ok = TRUE,
        fit_error = NA_character_,
        auc_ok = is.finite(auc),
        auc_error = auc_error
    )
}


.failed_resample <- function(n_genes, message) {
    list(
        selected = integer(),
        coefficients = numeric(n_genes),
        mutual_information = numeric(n_genes),
        auc = NA_real_,
        fit_ok = FALSE,
        fit_error = message,
        auc_ok = FALSE,
        auc_error = "Evaluation was skipped because model fitting failed"
    )
}


.run_resamples <- function(
    X, y, subsamples, alpha, workers, base_seed, verbose
) {
    run_one <- function(split) {
        .fit_one_resample(X, y, split, alpha)
    }
    if (workers == 1L) {
        return(lapply(subsamples, run_one))
    }

    cluster <- tryCatch(
        parallel::makeCluster(workers, type = "PSOCK"),
        error = identity
    )
    if (inherits(cluster, "error")) {
        warning(
            "Parallel workers could not be created. Serial fitting was used.",
            call. = FALSE
        )
        return(lapply(subsamples, run_one))
    }
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    .fit_message(verbose, "Using %d parallel workers.", workers)
    parallel::clusterSetRNGStream(cluster, iseed = base_seed)
    parallel::parLapply(cluster, subsamples, run_one)
}


.check_resample_results <- function(results, max_failed_fraction) {
    fit_ok <- vapply(results, function(x) isTRUE(x$fit_ok), logical(1))
    auc_ok <- vapply(results, function(x) isTRUE(x$auc_ok), logical(1))
    fit_fraction <- mean(!fit_ok)
    auc_fraction <- mean(!auc_ok)
    fit_errors <- unique(vapply(
        results[!fit_ok], function(x) x$fit_error, character(1)
    ))
    auc_errors <- unique(vapply(
        results[!auc_ok], function(x) x$auc_error, character(1)
    ))
    if (fit_fraction > max_failed_fraction) {
        failure <- .failure_message("Model fitting", fit_ok, fit_errors)
        stop(failure, call. = FALSE)
    }
    if (auc_fraction > max_failed_fraction) {
        stop(.failure_message("Excluded-sample evaluation", auc_ok, auc_errors),
            call. = FALSE
        )
    }
    if (any(!fit_ok)) {
        warning(.failure_message("Model fitting", fit_ok, fit_errors),
            call. = FALSE
        )
    }
    if (any(!auc_ok)) {
        warning(.failure_message(
            "Excluded-sample evaluation", auc_ok, auc_errors
        ), call. = FALSE)
    }
    list(
        n_resamples = length(results),
        n_fit_failures = sum(!fit_ok),
        fit_failure_fraction = fit_fraction,
        fit_errors = fit_errors,
        n_auc_failures = sum(!auc_ok),
        auc_failure_fraction = auc_fraction,
        auc_errors = auc_errors,
        max_failed_fraction = max_failed_fraction
    )
}


.failure_message <- function(stage, ok, messages) {
    details <- paste(
        utils::head(messages[!is.na(messages) & nzchar(messages)], 3L),
        collapse = " | "
    )
    sprintf(
        "%s failed in %d of %d resamples: %s",
        stage, sum(!ok), length(ok), details
    )
}


.aggregate_resamples <- function(results, n_genes) {
    n_resamples <- length(results)
    selection_matrix <- matrix(FALSE, nrow = n_genes, ncol = n_resamples)
    coefficient_matrix <- matrix(0, nrow = n_genes, ncol = n_resamples)
    mi_matrix <- matrix(0, nrow = n_genes, ncol = n_resamples)
    for (index in seq_len(n_resamples)) {
        selected <- results[[index]]$selected
        if (length(selected)) {
            selection_matrix[selected, index] <- TRUE
        }
        coefficient_matrix[, index] <- results[[index]]$coefficients
        mi_matrix[, index] <- results[[index]]$mutual_information
    }
    list(
        recurrence = rowMeans(selection_matrix),
        mutual_information = rowMeans(mi_matrix),
        selection_matrix = selection_matrix,
        coefficient_matrix = coefficient_matrix
    )
}


.compute_predictive_contribution <- function(
    X, subsamples, aggregated, verbose
) {
    n_genes <- ncol(X)
    candidates <- which(aggregated$recurrence > 0)
    shap_frequency <- numeric(n_genes)
    shap_matrix <- matrix(
        numeric(),
        nrow = nrow(X), ncol = 0L,
        dimnames = list(rownames(X), character())
    )
    if (length(candidates)) {
        .fit_message(
            verbose,
            "Calculating excluded-sample contributions for %d genes.",
            length(candidates)
        )
        shap <- compute_shap_utility(
            X, subsamples, candidates, aggregated$coefficient_matrix
        )
        shap_frequency[candidates] <- shap$instance_frequency
        shap_matrix <- shap$shap_matrix
    }
    raw <- sqrt(
        shap_frequency * aggregated$mutual_information + 1e-10
    )
    list(
        shap_frequency = shap_frequency,
        mutual_information = aggregated$mutual_information,
        raw = raw,
        shap_matrix = shap_matrix
    )
}


.adjust_ranking_components <- function(
    recurrence, contribution, null, recurrence_epsilon,
    contribution_epsilon, ratio_limit
) {
    if (is.null(contribution_epsilon)) {
        positive <- contribution[is.finite(contribution) & contribution > 0]
        contribution_epsilon <- if (length(positive)) {
            as.numeric(stats::quantile(positive, 0.95, names = FALSE))
        } else {
            1e-6
        }
    }
    list(
        recurrence = calibrate_by_null(
            recurrence, null$null_frequencies,
            epsilon = recurrence_epsilon, winsorize_at = ratio_limit
        ),
        contribution = calibrate_by_null(
            contribution, null$null_contribution,
            epsilon = contribution_epsilon, winsorize_at = ratio_limit
        ),
        contribution_epsilon = contribution_epsilon
    )
}


.build_gene_score_table <- function(
    genes, final_score, combined_score, recurrence, contribution, adjusted
) {
    scores <- data.frame(
        gene = genes,
        final_score = final_score,
        combined_score = combined_score,
        recurrence = recurrence,
        shap_frequency = contribution$shap_frequency,
        mutual_information = contribution$mutual_information,
        predictive_contribution = contribution$raw,
        adjusted_recurrence = adjusted$recurrence,
        adjusted_contribution = adjusted$contribution,
        stringsAsFactors = FALSE
    )
    score_order <- order(scores$final_score, decreasing = TRUE)
    scores <- scores[score_order, , drop = FALSE]
    rownames(scores) <- NULL
    scores
}


.summarise_auc <- function(results) {
    auc <- vapply(results, function(x) x$auc, numeric(1))
    list(
        auc_scores = auc,
        mean_auc = mean(auc, na.rm = TRUE),
        sd_auc = stats::sd(auc, na.rm = TRUE),
        median_auc = stats::median(auc, na.rm = TRUE),
        n_resamples = length(results)
    )
}


.fit_message <- function(verbose, format, ...) {
    if (isTRUE(verbose)) {
        message(sprintf(format, ...))
    }
    invisible(NULL)
}


.elapsed_seconds <- function(start) {
    as.numeric(difftime(Sys.time(), start, units = "secs"))
}


#' Print a GeneSelectR result
#'
#' @param x A `geneselectr_result` object.
#' @param ... Additional arguments, currently ignored.
#' @return `x`, invisibly.
#' @export
print.geneselectr_result <- function(x, ...) {
    cat("GeneSelectR result\n")
    cat(sprintf("  Genes ranked: %d\n", nrow(x$gene_scores)))
    cat(sprintf("  Mean excluded-sample AUC: %.4f\n", x$cv_results$mean_auc))
    cat(sprintf(
        "  Nogueira stability: %.4f\n",
        x$stability$nogueira_index
    ))
    cat(sprintf("  Elastic-net alpha: %.2f\n", x$parameters$alpha))
    invisible(x)
}

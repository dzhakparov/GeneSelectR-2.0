#' Validate GeneSelectR inputs
#'
#' @param X Numeric expression matrix with samples in rows.
#' @param y Two-level factor.
#' @param gene_names Gene names.
#' @param min_per_class Minimum observations per outcome level.
#' @return `TRUE`, invisibly.
#' @keywords internal
validate_inputs <- function(X, y, gene_names = NULL, min_per_class = 5L) {
    if (!is.matrix(X) || !is.numeric(X) || !nrow(X) || !ncol(X)) {
        stop("X must be a non-empty numeric matrix")
    }
    if (!is.factor(y) || nlevels(y) != 2L || anyNA(y) ||
        any(table(y) == 0L)) {
        stop("y must be a two-level factor without missing values")
    }
    if (nrow(X) != length(y)) {
        stop("The number of rows in X must equal the length of y")
    }
    if (is.null(colnames(X))) {
        stop("X must have column names containing gene names")
    }
    if (anyNA(colnames(X)) || any(!nzchar(colnames(X))) ||
        anyDuplicated(colnames(X))) {
        stop("X column names must be unique, non-missing gene names")
    }
    if (!is.null(gene_names) && length(gene_names) != ncol(X)) {
        stop("gene_names must contain one name for each column of X")
    }
    if (any(table(y) < min_per_class)) {
        stop(sprintf(
            "Each outcome level must contain at least %d observations",
            min_per_class
        ))
    }
    if (any(!is.finite(X))) {
        stop("X must not contain missing or non-finite values")
    }
    zero_variance <- apply(X, 2L, stats::sd) == 0
    if (any(zero_variance)) {
        warning(
            sprintf("%d genes have zero variance", sum(zero_variance)),
            call. = FALSE
        )
    }
    invisible(TRUE)
}


#' Convert positive values to percentile scores
#'
#' Zero values remain zero. Positive values are assigned their fractional rank.
#'
#' @param x Numeric vector.
#' @return Numeric vector between zero and one.
#' @keywords internal
percentile01 <- function(x) {
    if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0)) {
        stop("x must contain finite non-negative values")
    }
    positive <- x > 0
    if (!any(positive) || length(unique(x)) == 1L) {
        return(numeric(length(x)))
    }
    result <- numeric(length(x))
    result[positive] <- rank(
        x[positive],
        ties.method = "average"
    ) / sum(positive)
    result
}


adaptive_bin_count <- function(n, maximum = 5L) {
    max(3L, min(as.integer(maximum), floor(sqrt(n))))
}


discretize_continuous <- function(x, n_bins = 5L) {
    if (length(unique(x)) < 2L) {
        return(integer(length(x)))
    }
    n_bins <- adaptive_bin_count(length(x), n_bins)
    breaks <- unique(stats::quantile(
        x,
        probs = seq(0, 1, length.out = n_bins + 1L), names = FALSE
    ))
    if (length(breaks) < 3L) {
        return(integer(length(x)))
    }
    as.integer(cut(
        x,
        breaks = breaks, labels = FALSE, include.lowest = TRUE
    ))
}


#' Compute mutual information with a binary outcome
#'
#' Expression values are divided into quantile-based intervals. The returned
#' value is measured in nats.
#'
#' @param x Numeric expression vector.
#' @param y Two-level factor.
#' @param n_bins Maximum number of intervals.
#' @return Mutual information in nats.
#' @examples
#' x <- c(0, 0, 1, 1, 2, 2)
#' y <- factor(c("a", "a", "a", "b", "b", "b"))
#' compute_mutual_information(x, y)
#' @export
compute_mutual_information <- function(x, y, n_bins = 5L) {
    if (!is.numeric(x) || any(!is.finite(x)) || length(x) != length(y)) {
        stop("x must be a finite numeric vector aligned with y")
    }
    if (!is.factor(y) || nlevels(y) != 2L || anyNA(y)) {
        stop("y must be a two-level factor without missing values")
    }
    binned <- discretize_continuous(x, n_bins)
    if (length(unique(binned)) < 2L) {
        return(0)
    }
    joint <- table(binned, y) / length(y)
    marginal_x <- rowSums(joint)
    marginal_y <- colSums(joint)
    information <- 0
    for (row in seq_along(marginal_x)) {
        for (column in seq_along(marginal_y)) {
            if (joint[row, column] > 0) {
                information <- information + joint[row, column] * log(
                    joint[row, column] /
                        (marginal_x[row] * marginal_y[column])
                )
            }
        }
    }
    information
}


compute_mi_vectorized <- function(X, y, n_bins = 5L) {
    vapply(
        seq_len(ncol(X)),
        function(column) compute_mutual_information(X[, column], y, n_bins),
        numeric(1)
    )
}


#' Create repeated stratified cross-validation splits
#'
#' Repeated K-fold training sets are used by the implemented GeneSelectR
#' workflow. The training sets overlap, so recurrence is a descriptive ranking
#' measurement.
#'
#' @param y Two-level factor.
#' @param B Number of training and excluded-sample splits.
#' @param random_seed Random seed, restored on return.
#' @param k_folds Number of folds.
#' @return List of training and excluded-sample index vectors.
#' @examples
#' y <- factor(
#'     rep(c("control", "case"), each = 10),
#'     levels = c("control", "case")
#' )
#' create_subsamples(y, B = 5, k_folds = 5)
#' @export
create_subsamples <- function(
    y, B = 50L, random_seed = 123L, k_folds = 5L
) {
    if (!is.factor(y) || nlevels(y) != 2L || anyNA(y)) {
        stop("y must be a two-level factor without missing values")
    }
    .validate_positive_integer(B, "B")
    .validate_positive_integer(k_folds, "k_folds", minimum = 2L)
    if (k_folds > min(table(y))) {
        stop("k_folds cannot exceed the size of the smaller outcome level")
    }
    withr::local_seed(random_seed)
    class_indices <- split(seq_along(y), y)
    repeats <- ceiling(B / k_folds)
    splits <- vector("list", 0L)
    for (repeat_index in seq_len(repeats)) {
        fold_id <- integer(length(y))
        for (level in names(class_indices)) {
            indices <- class_indices[[level]]
            fold_id[indices] <- sample(rep(
                seq_len(k_folds),
                length.out = length(indices)
            ))
        }
        for (fold in seq_len(k_folds)) {
            if (length(splits) >= B) {
                break
            }
            excluded <- which(fold_id == fold)
            splits[[length(splits) + 1L]] <- list(
                train = setdiff(seq_along(y), excluded),
                oob = excluded,
                subsample_id = length(splits) + 1L,
                repeat_id = repeat_index,
                fold_id = fold
            )
        }
    }
    splits
}


fit_regularized_model <- function(X_train, y_train, alpha = 0.5) {
    class_counts <- table(y_train)
    if (length(class_counts) != 2L || min(class_counts) < 3L) {
        stop(
            "Model fitting requires three observations per outcome level",
            call. = FALSE
        )
    }
    cv_folds <- min(5L, min(class_counts))
    fold_id <- NULL
    if (cv_folds < 5L) {
        fold_id <- integer(length(y_train))
        for (level in levels(y_train)) {
            indices <- which(y_train == level)
            fold_id[indices] <- sample(rep(
                seq_len(cv_folds),
                length.out = length(indices)
            ))
        }
    }
    fit <- glmnet::cv.glmnet(
        X_train,
        as.integer(y_train) - 1L,
        family = "binomial",
        alpha = alpha,
        nfolds = cv_folds,
        foldid = fold_id,
        type.measure = "auc"
    )
    coefficients <- as.vector(
        glmnet::coef.glmnet(fit, s = "lambda.min")
    )
    intercept <- coefficients[1L]
    coefficients <- coefficients[-1L]
    selected <- which(coefficients != 0)
    predict_fn <- function(X_new) {
        stats::plogis(as.numeric(X_new %*% coefficients) + intercept)
    }
    list(
        selected = selected,
        full_coef_vector = coefficients,
        intercept = intercept,
        lambda = fit$lambda.min,
        predict_fn = predict_fn
    )
}


compute_auc <- function(y_true, predictions) {
    if (!is.factor(y_true) || nlevels(y_true) != 2L) {
        stop("y_true must be a two-level factor")
    }
    if (!is.numeric(predictions) || length(predictions) != length(y_true) ||
        any(!is.finite(predictions))) {
        stop("predictions must be finite and aligned with y_true")
    }
    response <- as.integer(y_true) - 1L
    as.numeric(pROC::auc(pROC::roc(
        response, predictions,
        quiet = TRUE, direction = "<"
    )))
}


#' Compute the Nogueira feature-selection stability index
#'
#' The index is calculated from repeated selected sets. Repeated K-fold
#' training sets overlap, so the value is descriptive for GeneSelectR results.
#'
#' @param selection_matrix Logical matrix with genes in rows and fits in
#'   columns.
#' @return List containing `nogueira_index`.
#' @examples
#' selections <- matrix(
#'     c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE),
#'     nrow = 4
#' )
#' compute_nogueira_stability(selections)
#' @export
compute_nogueira_stability <- function(selection_matrix) {
    if (!is.matrix(selection_matrix) || ncol(selection_matrix) < 2L ||
        !is.logical(selection_matrix) || anyNA(selection_matrix)) {
        stop(
            "selection_matrix must be logical with at least two columns",
            call. = FALSE
        )
    }
    n_fits <- ncol(selection_matrix)
    n_genes <- nrow(selection_matrix)
    mean_size <- mean(colSums(selection_matrix))
    denominator <- mean_size * (1 - mean_size / n_genes)
    if (denominator <= 0) {
        return(list(nogueira_index = NA_real_))
    }
    frequency <- rowMeans(selection_matrix)
    variance <- frequency * (1 - frequency) * n_fits / (n_fits - 1)
    list(nogueira_index = 1 - sum(variance) / denominator)
}


compute_shap_utility <- function(
    X, subsamples, candidate_index, coefficient_matrix,
    threshold_quantile = 0.75
) {
    n_samples <- nrow(X)
    n_candidates <- length(candidate_index)
    shap_sum <- matrix(0, nrow = n_samples, ncol = n_candidates)
    observation_count <- integer(n_samples)
    for (fit_index in seq_along(subsamples)) {
        excluded <- subsamples[[fit_index]]$oob
        training <- subsamples[[fit_index]]$train
        coefficients <- coefficient_matrix[candidate_index, fit_index]
        observation_count[excluded] <- observation_count[excluded] + 1L
        if (all(coefficients == 0)) {
            next
        }
        training_mean <- colMeans(
            X[training, candidate_index, drop = FALSE]
        )
        centred <- sweep(
            X[excluded, candidate_index, drop = FALSE],
            2L, training_mean, "-"
        )
        shap_sum[excluded, ] <- shap_sum[excluded, ] + abs(sweep(
            centred, 2L, coefficients, "*"
        ))
    }
    observed <- observation_count > 0L
    shap_mean <- shap_sum
    shap_mean[observed, ] <- shap_sum[observed, , drop = FALSE] /
        observation_count[observed]
    colnames(shap_mean) <- colnames(X)[candidate_index]

    substantial <- matrix(FALSE, nrow = n_samples, ncol = n_candidates)
    for (sample_index in which(observed)) {
        values <- shap_mean[sample_index, ]
        positive <- values[values > 0]
        if (!length(positive)) {
            next
        }
        threshold <- stats::quantile(
            positive, threshold_quantile,
            names = FALSE
        )
        substantial[sample_index, ] <- values >= threshold
    }
    instance_frequency <- if (any(observed)) {
        colMeans(substantial[observed, , drop = FALSE])
    } else {
        numeric(n_candidates)
    }
    names(instance_frequency) <- colnames(X)[candidate_index]
    list(
        instance_frequency = instance_frequency,
        shap_matrix = shap_mean
    )
}

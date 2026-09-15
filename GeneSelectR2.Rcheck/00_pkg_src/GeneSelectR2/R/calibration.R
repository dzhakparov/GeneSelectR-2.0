#' Compare observed gene scores with shuffled-outcome references
#'
#' Each observed value is divided by the mean shuffled-outcome value for the
#' same gene after adding a smoothing constant to both. Ratios are bounded on
#' the log2 scale. They are ranking measurements and provide no probability,
#' p-value, or false-discovery-rate interpretation.
#'
#' @param observed Numeric vector of observed scores.
#' @param null_matrix Matrix with shuffled outcomes in rows and genes in
#'   columns.
#' @param epsilon Positive smoothing constant.
#' @param winsorize_at Maximum absolute log2 ratio.
#' @return Numeric vector of adjusted scores.
#' @keywords internal
calibrate_by_null <- function(
    observed, null_matrix, epsilon, winsorize_at = 4
) {
    if (!is.numeric(observed) || any(!is.finite(observed)) ||
        any(observed < 0)) {
        stop("observed must contain finite non-negative values")
    }
    if (!is.matrix(null_matrix) || !is.numeric(null_matrix) ||
        ncol(null_matrix) != length(observed) ||
        any(!is.finite(null_matrix)) || any(null_matrix < 0)) {
        stop("null_matrix must contain finite non-negative gene scores")
    }
    .validate_positive_number(epsilon, "epsilon")
    .validate_positive_number(winsorize_at, "winsorize_at")

    reference <- colMeans(null_matrix)
    ratio <- (observed + epsilon) / (reference + epsilon)
    lower <- 2^(-winsorize_at)
    upper <- 2^winsorize_at
    pmax(lower, pmin(upper, ratio))
}

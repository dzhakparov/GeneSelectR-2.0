#' Plot gene-level GeneSelectR measurements
#'
#' Display the recurrence, permutation-adjusted predictive contribution and
#' final ranking score for the leading genes in a fitted GeneSelectR result.
#'
#' @param x A GeneSelectR result.
#' @param n Number of leading genes to display.
#' @param colours Three colours used for recurrence, contribution and final
#'   score.
#'
#' @return The displayed rows of `x$gene_scores`, invisibly.
#' @examples
#' data(asthma_example)
#' fit <- select_genes(asthma_example, "severity",
#'     n_genes = 10,
#'     alpha = 0.5, B = 2, permutations = 1, null_B = 2
#' )
#' plot_gene_ranking(fit, n = 8)
#' @export
plot_gene_ranking <- function(x, n = 15L,
                                colours = c("#087E82", "#7A5BB5", "#D27342")) {
    scores <- .validate_gene_ranking_plot(x, n)
    if (!is.character(colours) || length(colours) != 3L || anyNA(colours)) {
        stop("colours must contain three valid colour specifications")
    }
    try(grDevices::col2rgb(colours), silent = FALSE)

    scores <- scores[rev(seq_len(nrow(scores))), , drop = FALSE]
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::par(
        mfrow = c(1, 3), mar = c(4.2, 1.2, 3.1, 0.8),
        oma = c(0, 6.5, 2.2, 0), xaxs = "i"
    )

    .gene_dot_plot(scores$recurrence, scores$gene,
        "Within-training recurrence", colours[[1]],
        limits = c(0, 1), show_labels = TRUE
    )
    .gene_dot_plot(scores$adjusted_contribution, scores$gene,
        "Adjusted contribution", colours[[2]],
        reference = 1
    )
    .gene_dot_plot(scores$final_score, scores$gene,
        "Final ranking score", colours[[3]],
        reference = 1
    )
    graphics::mtext("Gene-level ranking measurements",
        outer = TRUE,
        side = 3, line = 0.4, font = 2, cex = 1.05
    )
    invisible(scores)
}

.validate_gene_ranking_plot <- function(x, n) {
    if (!is.list(x) || !is.data.frame(x$gene_scores)) {
        stop("x must contain a gene_scores data frame")
    }
    required <- c(
        "gene", "recurrence", "adjusted_contribution", "final_score"
    )
    missing <- setdiff(required, names(x$gene_scores))
    if (length(missing)) {
        stop(
            sprintf(
                "gene_scores is missing: %s",
                paste(missing, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
        n < 1 || n != as.integer(n)) {
        stop("n must be a positive integer")
    }
    n <- min(as.integer(n), nrow(x$gene_scores))
    scores <- utils::head(x$gene_scores, n)
    for (column in required[-1L]) {
        valid_column <- is.numeric(scores[[column]]) &&
            all(is.finite(scores[[column]]))
        if (!valid_column) {
            stop(sprintf(
                "%s must contain finite numeric values", column
            ), call. = FALSE)
        }
    }
    scores
}

.gene_dot_plot <- function(values, labels, title, colour, limits = NULL,
                            reference = NULL, show_labels = FALSE) {
    if (is.null(limits)) {
        upper <- max(c(values, reference), na.rm = TRUE)
        limits <- c(0, upper * 1.08)
    }
    y <- seq_along(values)
    graphics::plot(values, y,
        type = "n", xlim = limits,
        ylim = c(0.5, length(values) + 0.5), yaxt = "n",
        ylab = "", xlab = title, bty = "n"
    )
    graphics::abline(v = pretty(limits), col = "#E6E8EB", lwd = 0.7)
    if (!is.null(reference)) {
        graphics::abline(v = reference, col = "#60656D", lty = 2, lwd = 1)
    }
    graphics::segments(0, y, values, y, col = "#C7CCD1", lwd = 1.4)
    graphics::points(values, y, pch = 19, cex = 1.15, col = colour)
    if (show_labels) {
        graphics::axis(2,
            at = y, labels = labels, las = 1,
            tick = FALSE, hadj = 1, cex.axis = 0.82
        )
    }
    graphics::box(bty = "l", col = "#7A7F86")
}

#' Plot gene-level predictive and biological evidence
#'
#' Display how often genes entered selected sets, their within-training
#' recurrence, permutation-adjusted predictive contribution and a continuous
#' disease-association score. Missing association scores are labelled.
#'
#' @param x Gene-level data containing selection counts, recurrence,
#'   contribution and disease association.
#' @param comparison_label Label for the optional comparison-count series.
#' @param colours Colours for GeneSelectR, the comparison series, recurrence,
#'   contribution and association.
#'
#' @return The displayed data, invisibly.
#' @examples
#' data(asthma_case_study)
#' plot_gene_evidence(asthma_case_study, comparison_label = "DGE")
#' @export
plot_gene_evidence <- function(
    x, comparison_label = "Comparison",
    colours = c("#087E82", "#9AA1A8", "#087E82", "#7A5BB5", "#D27342")
) {
    required <- c(
        "gene", "selection_count", "recurrence", "contribution",
        "association"
    )
    if (!is.data.frame(x) || length(setdiff(required, names(x)))) {
        input_message <- paste(
            "x must contain gene, selection_count, recurrence,",
            "contribution and association"
        )
        stop(
            input_message,
            call. = FALSE
        )
    }
    if (!nrow(x) || anyDuplicated(x$gene) || anyNA(x$gene)) {
        stop("x must contain unique gene names")
    }
    numeric_columns <- setdiff(required, "gene")
    for (column in numeric_columns) {
        if (!is.numeric(x[[column]]) ||
            any(!is.finite(x[[column]][!is.na(x[[column]])]))) {
            stop(sprintf("%s must contain numeric values", column),
                call. = FALSE
            )
        }
    }
    if (length(colours) != 5L) {
        stop("colours must contain five colour specifications")
    }
    try(grDevices::col2rgb(colours), silent = FALSE)

    x <- x[order(-x$selection_count, x$gene), , drop = FALSE]
    x <- x[rev(seq_len(nrow(x))), , drop = FALSE]
    y <- seq_len(nrow(x))
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::par(
        mfrow = c(1, 4), mar = c(4.3, 1.1, 3.1, 0.7),
        oma = c(0, 6.5, 2.2, 0), xaxs = "i"
    )

    count_max <- max(c(x$selection_count, x$comparison_count), na.rm = TRUE)
    graphics::plot(x$selection_count, y,
        type = "n", xlim = c(0, count_max * 1.08),
        ylim = c(0.5, nrow(x) + 0.5), yaxt = "n", ylab = "",
        xlab = "Outer selection count", bty = "n"
    )
    graphics::abline(v = pretty(c(0, count_max)), col = "#E6E8EB", lwd = 0.7)
    if ("comparison_count" %in% names(x)) {
        graphics::segments(x$comparison_count, y - 0.12, x$selection_count,
            y + 0.12,
            col = "#C7CCD1", lwd = 1
        )
        graphics::points(x$comparison_count, y - 0.12,
            pch = 1, cex = 1,
            col = colours[[2]]
        )
    }
    graphics::points(x$selection_count, y + 0.12,
        pch = 19, cex = 1,
        col = colours[[1]]
    )
    graphics::axis(2,
        at = y, labels = x$gene, las = 1, tick = FALSE,
        hadj = 1, cex.axis = 0.82
    )
    graphics::legend("bottomright",
        legend = c("GeneSelectR", comparison_label),
        pch = c(19, 1), col = colours[seq_len(2L)],
        bty = "n", cex = 0.72
    )
    graphics::box(bty = "l", col = "#7A7F86")

    .gene_dot_plot(x$recurrence, x$gene, "Within-training recurrence",
        colours[[3]],
        limits = c(0, 1)
    )
    .gene_dot_plot(x$contribution, x$gene, "Adjusted contribution",
        colours[[4]],
        reference = 1
    )

    association_limit <- max(x$association, na.rm = TRUE)
    association_limit <- if (
        is.finite(association_limit) && association_limit > 0
    ) {
        association_limit * 1.18
    } else {
        1
    }
    graphics::plot(x$association, y,
        type = "n", xlim = c(0, association_limit),
        ylim = c(0.5, nrow(x) + 0.5), yaxt = "n", ylab = "",
        xlab = "Disease-association score", bty = "n"
    )
    graphics::abline(
        v = pretty(c(0, association_limit)), col = "#E6E8EB",
        lwd = 0.7
    )
    available <- !is.na(x$association)
    graphics::segments(0, y[available], x$association[available], y[available],
        col = "#C7CCD1", lwd = 1.4
    )
    graphics::points(x$association[available], y[available],
        pch = 19,
        cex = 1.1, col = colours[[5]]
    )
    if (any(!available)) {
        graphics::text(0, y[!available],
            labels = "NA", pos = 4,
            cex = 0.76, col = "#555A60"
        )
    }
    graphics::box(bty = "l", col = "#7A7F86")
    graphics::mtext("Asthma gene-level output",
        outer = TRUE, side = 3,
        line = 0.4, font = 2, cex = 1.05
    )
    invisible(x)
}

#' Plot biological assessment against matched random gene sets
#'
#' Draw a heatmap of selected-to-random ratios for one feature-selection
#' method. A ratio of one denotes the mean of same-size random gene sets.
#' Biological measurements remain separate and are not averaged.
#'
#' @param x Biological ratios by dataset and method.
#' @param method Method to display.
#' @param measures Named character vector mapping labels to columns of `x`.
#'   `NULL` uses all biological measures in `benchmark_biology`.
#' @param dataset_order Optional order for dataset labels.
#' @param colour_limit Maximum absolute log2 ratio used for the colour scale.
#'
#' @return The displayed method-specific data, invisibly.
#' @examples
#' data(benchmark_biology)
#' plot_biology_comparison(benchmark_biology, method = "GeneSelectR")
#' @export
plot_biology_comparison <- function(
    x, method = "GeneSelectR", measures = NULL,
    dataset_order = NULL, colour_limit = NULL
) {
    if (is.null(measures)) {
        measures <- c(
            "GO similarity" = "go_ratio",
            "Hallmark sharing" = "hallmark_ratio",
            "Open Targets >=0.05" = "open_targets_05_ratio",
            "Open Targets >=0.10" = "open_targets_10_ratio"
        )
    }
    required <- c("dataset", "method", unname(measures))
    missing <- setdiff(required, names(x))
    if (!is.data.frame(x) || length(missing)) {
        stop(
            sprintf("x is missing: %s", paste(missing, collapse = ", ")),
            call. = FALSE
        )
    }
    if (!is.character(method) || length(method) != 1L || is.na(method)) {
        stop("method must be one method name")
    }
    shown <- x[x$method == method, required, drop = FALSE]
    if (!nrow(shown)) {
        stop("method is absent from x")
    }
    if (anyDuplicated(shown$dataset)) {
        stop("x must contain one row per dataset and method")
    }
    values <- as.matrix(shown[, unname(measures), drop = FALSE])
    storage.mode(values) <- "double"
    if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
        stop("biological ratios must be finite and non-negative")
    }
    rownames(values) <- shown$dataset
    colnames(values) <- names(measures)
    if (!is.null(dataset_order)) {
        if (!setequal(dataset_order, rownames(values))) {
            stop("dataset_order must contain each displayed dataset once")
        }
        values <- values[rev(dataset_order), , drop = FALSE]
    }

    log_values <- log2(pmax(values, 1 / 16))
    if (is.null(colour_limit)) {
        colour_limit <- max(1, ceiling(max(abs(log_values))))
    }
    if (length(colour_limit) != 1L || !is.finite(colour_limit) ||
        colour_limit <= 0) {
        stop("colour_limit must be one positive number")
    }
    clipped <- matrix(
        pmax(-colour_limit, pmin(colour_limit, log_values)),
        nrow = nrow(log_values), dimnames = dimnames(log_values)
    )
    palette <- grDevices::colorRampPalette(
        c("#D98978", "#F6F5F2", "#2F8F88")
    )(201)
    colour_index <- round(
        (clipped + colour_limit) / (2 * colour_limit) * 200
    ) + 1L

    nr <- nrow(values)
    nc <- ncol(values)
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::par(mar = c(6.2, 8.5, 3.8, 1.1), xaxs = "i", yaxs = "i")
    graphics::plot(c(0.5, nc + 0.5), c(0.5, nr + 0.5),
        type = "n",
        axes = FALSE, xlab = "", ylab = "", bty = "n"
    )
    for (row in seq_len(nr)) {
        for (column in seq_len(nc)) {
            graphics::rect(column - 0.5, row - 0.5, column + 0.5, row + 0.5,
                col = palette[colour_index[row, column]],
                border = "white", lwd = 1.4
            )
            label_colour <- if (
                abs(clipped[row, column]) > colour_limit * 0.55
            ) {
                "white"
            } else {
                "#17324D"
            }
            graphics::text(column, row, sprintf("%.2fx", values[row, column]),
                col = label_colour, cex = 0.82, font = 2
            )
        }
    }
    source_labels <- vapply(
        colnames(values),
        function(label) paste(strwrap(label, width = 17), collapse = "\n"),
        character(1)
    )
    graphics::axis(1,
        at = seq_len(nc), labels = source_labels,
        las = 1, tick = FALSE, cex.axis = 0.78, line = -0.2
    )
    graphics::axis(2,
        at = seq_len(nr), labels = rownames(values),
        las = 1, tick = FALSE, cex.axis = 0.86
    )
    graphics::title(main = paste(method, "biological assessment"))
    graphics::mtext("Selected gene set / matched-random mean",
        side = 1,
        line = 4.8, cex = 0.82
    )
    invisible(shown)
}

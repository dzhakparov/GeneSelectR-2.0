test_that("gene ranking plot uses fitted score columns", {
    score_table <- data.frame(
        gene = LETTERS[1:4], recurrence = c(0.9, 0.7, 0.5, 0.2),
        adjusted_contribution = c(3, 2, 1.2, 0.8),
        final_score = c(1, 0.75, 0.5, 0.25)
    )
    fit <- list(gene_scores = score_table)
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file)
    shown <- plot_gene_ranking(fit, n = 3)
    grDevices::dev.off()
    expect_equal(nrow(shown), 3)
    expect_true(file.info(file)$size > 0)
    expect_error(plot_gene_ranking(list(), n = 3), "gene_scores")
})

test_that("gene evidence plot preserves missing annotations", {
    data("asthma_case_study", package = "GeneSelectR2")
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file, width = 10, height = 5)
    shown <- plot_gene_evidence(asthma_case_study, comparison_label = "DGE")
    grDevices::dev.off()
    expect_equal(sum(is.na(shown$association)), 1)
    expect_true(file.info(file)$size > 0)
})

test_that("biological assessment plot uses separate ratios", {
    data("benchmark_biology", package = "GeneSelectR2")
    file <- tempfile(fileext = ".pdf")
    grDevices::pdf(file, width = 8, height = 6)
    shown <- plot_biology_comparison(benchmark_biology, method = "GeneSelectR")
    grDevices::dev.off()
    expect_equal(nrow(shown), 7)
    expect_true(file.info(file)$size > 0)
    expect_error(
        plot_biology_comparison(benchmark_biology, method = "missing"),
        "absent"
    )
})

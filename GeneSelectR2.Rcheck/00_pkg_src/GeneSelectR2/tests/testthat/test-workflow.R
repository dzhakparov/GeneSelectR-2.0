test_that("the asthma example has aligned, documented inputs", {
    data("asthma_example", package = "GeneSelectR2")
    expect_s4_class(asthma_example, "SummarizedExperiment")
    expect_equal(dim(asthma_example), c(300, 96))
    expect_equal(as.integer(table(asthma_example$severity)), c(48, 48))
    expect_identical(colnames(asthma_example),
                     rownames(SummarizedExperiment::colData(asthma_example)))
})

test_that("the workflow validates sample alignment and set size", {
    data("asthma_example", package = "GeneSelectR2")
    expect_error(select_genes(asthma_example, "absent"), "colData")
    expect_error(select_genes(asthma_example, "severity", n_genes = 301),
                 "n_genes")
    expect_error(select_genes(matrix(1, 4, 4), factor(c(0, 0, 1, 1))),
                 "gene names")
})
test_that("the high-level workflow matches the configured core fit", {
    data("asthma_example", package = "GeneSelectR2")
    example <- asthma_example[seq_len(30), ]
    result <- select_genes(example, "severity", n_genes = 10,
                           alpha = 0.5, B = 5, permutations = 2, null_B = 5)
    reference <- geneselectr2_fit(
        t(as.matrix(SummarizedExperiment::assay(example))), example$severity,
        B = 5, alpha = 0.5, permutations = 2, null_B = 5,
        workers = 1, random_seed = 123, verbose = FALSE
    )
    expect_equal(result$gene_scores, reference$gene_scores)
    expect_identical(result$selected_genes,
                     head(reference$gene_scores$gene, 10))
})

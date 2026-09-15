test_that("the core fit returns the implemented ranking measurements", {
    set.seed(3)
    X <- matrix(rnorm(80 * 16), nrow = 80)
    colnames(X) <- paste0("gene", seq_len(ncol(X)))
    y <- factor(rep(c("control", "case"), each = 40),
                levels = c("control", "case"))
    X[y == "case", 1:3] <- X[y == "case", 1:3] + 1

    fit <- geneselectr2_fit(
        X, y, B = 5, permutations = 2, null_B = 5, verbose = FALSE
    )

    expected <- c(
        "gene", "final_score", "combined_score", "recurrence",
        "shap_frequency", "mutual_information", "predictive_contribution",
        "adjusted_recurrence", "adjusted_contribution"
    )
    expect_s3_class(fit, "geneselectr_result")
    expect_identical(names(fit$gene_scores), expected)
    expect_equal(nrow(fit$gene_scores), ncol(X))
    expect_true(all(diff(fit$gene_scores$final_score) <= 0))
    expect_true(all(fit$gene_scores$recurrence >= 0))
    expect_true(all(fit$gene_scores$recurrence <= 1))
    expect_equal(
        fit$gene_scores$combined_score,
        exp(
            0.5 * log(fit$gene_scores$adjusted_recurrence + 1e-10) +
                0.5 * log(fit$gene_scores$adjusted_contribution + 1e-10)
        )
    )
})


test_that("random_seed controls the complete core fit", {
    set.seed(4)
    X <- matrix(rnorm(80 * 12), nrow = 80)
    colnames(X) <- paste0("gene", seq_len(ncol(X)))
    y <- factor(rep(c("control", "case"), each = 40),
                levels = c("control", "case"))
    X[y == "case", 1:3] <- X[y == "case", 1:3] + 1
    arguments <- list(
        X = X, y = y, B = 5, permutations = 2, null_B = 5,
        random_seed = 17, verbose = FALSE
    )

    first <- do.call(geneselectr2_fit, arguments)
    runif(20)
    second <- do.call(geneselectr2_fit, arguments)

    expect_identical(first$gene_scores$gene, second$gene_scores$gene)
    expect_equal(first$gene_scores$final_score,
                 second$gene_scores$final_score, tolerance = 0)
})


test_that("discarded experimental settings are absent from the interface", {
    removed <- c(
        "gate_method", "pfer", "q_max", "regularization_method",
        "gene_groups", "bio_mode", "components", "score_formula",
        "score_weights", "selection_method"
    )
    expect_false(any(removed %in% names(formals(geneselectr2_fit))))
})

test_that("repeated splits are stratified and use every sample", {
    y <- factor(rep(c("control", "case"), each = 10),
                levels = c("control", "case"))
    splits <- create_subsamples(y, B = 7, k_folds = 5, random_seed = 2)
    expect_length(splits, 7)
    for (split in splits) {
        expect_length(intersect(split$train, split$oob), 0)
        expect_setequal(c(split$train, split$oob), seq_along(y))
        expect_identical(levels(droplevels(y[split$oob])), levels(y))
    }
})


test_that("mutual information is zero for a constant gene", {
    y <- factor(rep(c("control", "case"), each = 10),
                levels = c("control", "case"))
    expect_equal(compute_mutual_information(rep(1, 20), y), 0)
    expect_gt(compute_mutual_information(as.numeric(y), y), 0)
})


test_that("the Nogueira index is one for identical selected sets", {
    selections <- matrix(FALSE, nrow = 5, ncol = 4)
    selections[1:2, ] <- TRUE
    result <- compute_nogueira_stability(selections)
    expect_equal(result$nogueira_index, 1)
})


test_that("outcome adjustment is gene specific and bounded", {
    observed <- c(0.8, 0.8)
    null <- rbind(c(0.2, 0.6), c(0.2, 0.6))
    adjusted <- GeneSelectR:::calibrate_by_null(
        observed, null, epsilon = 0.01, winsorize_at = 2
    )
    expect_gt(adjusted[1], adjusted[2])
    expect_true(all(adjusted >= 0.25 & adjusted <= 4))
})


test_that("AUC preserves the specified class direction", {
    y <- factor(
        c("control", "control", "case", "case"),
        levels = c("control", "case")
    )
    expect_equal(GeneSelectR2:::compute_auc(y, c(0.9, 0.8, 0.2, 0.1)), 0)
    expect_equal(GeneSelectR2:::compute_auc(y, c(0.1, 0.2, 0.8, 0.9)), 1)
})

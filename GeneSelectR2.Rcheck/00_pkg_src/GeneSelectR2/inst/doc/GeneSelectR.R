## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----installation, eval=FALSE-------------------------------------------------
# BiocManager::install("GeneSelectR2")

## ----data---------------------------------------------------------------------
library(GeneSelectR2)

set.seed(10)
n <- 40
p <- 20
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("gene", seq_len(p))
y <- factor(
    rep(c("control", "case"), each = n / 2),
    levels = c("control", "case")
)
X[y == "case", 1:3] <- X[y == "case", 1:3] + 1

## ----fit----------------------------------------------------------------------
fit <- geneselectr2_fit(
    X,
    y,
    B = 5,
    permutations = 2,
    null_B = 5,
    verbose = FALSE
)

## ----results------------------------------------------------------------------
head(fit$gene_scores[, c(
    "gene", "recurrence", "adjusted_contribution", "final_score"
)])
fit$stability["nogueira_index"]
fit$cv_results[c("mean_auc", "sd_auc", "n_resamples")]

## ----provenance---------------------------------------------------------------
fit$parameters
sessionInfo()


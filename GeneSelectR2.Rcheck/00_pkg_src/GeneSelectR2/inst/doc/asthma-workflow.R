## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----installation, eval=FALSE-------------------------------------------------
# BiocManager::install("GeneSelectR2")

## ----input--------------------------------------------------------------------
library(GeneSelectR2)
library(SummarizedExperiment)
data(asthma_example)
dim(asthma_example)
table(asthma_example$severity)
metadata(asthma_example)

## ----fit----------------------------------------------------------------------
fit <- select_genes(
    asthma_example, outcome = "severity", n_genes = 20,
    alpha = 0.5, B = 5, permutations = 2, null_B = 5,
    workers = 1, seed = 123
)
fit$alpha_comparison
fit$selected_genes

## ----full-settings, eval=FALSE------------------------------------------------
# fit <- select_genes(training_experiment, "severity", n_genes = 20)

## ----output-------------------------------------------------------------------
head(fit$gene_scores)
fit$parameters

## ----ranking-plot, fig.width=9, fig.height=5----------------------------------
plot_gene_ranking(fit, n = 12)

## ----biology-heatmap, fig.width=8.5, fig.height=6.3---------------------------
data(benchmark_biology)
plot_biology_comparison(
    benchmark_biology, method = "GeneSelectR",
    dataset_order = c(
        "Atopic dermatitis", "Bladder cancer", "Asthma", "Sepsis",
        "Psoriasis", "Tuberculosis", "Crohn disease"
    )
)

## ----case-study, fig.width=11, fig.height=5.5---------------------------------
data(asthma_case_study)
plot_gene_evidence(asthma_case_study, comparison_label = "DGE")

## ----session------------------------------------------------------------------
sessionInfo()


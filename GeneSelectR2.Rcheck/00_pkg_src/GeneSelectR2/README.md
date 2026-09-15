# GeneSelectR

GeneSelectR is an R workflow for predictive gene selection in binary-outcome
transcriptomic studies. `select_genes()` combines repeated elastic-net
inclusion and excluded-sample predictive contribution after comparison with
shuffled outcomes. The output retains both measurements for inspection.
Biological interpretation follows gene ranking.

## Installation

The package is under development. A local source checkout can be installed with:

```r
install.packages(".", repos = NULL, type = "source")
```

Following Bioconductor acceptance, the release version will be installed with:

```r
BiocManager::install("GeneSelectR")
```

## Asthma example

```r
library(GeneSelectR)

data(asthma_example)
fit <- select_genes(
    asthma_example, "severity", n_genes = 20,
    alpha = 0.5, B = 5, permutations = 2, null_B = 5
)
fit$selected_genes
head(fit$gene_scores)
```

The evaluated asthma vignette uses 96 samples and 300 genes from GSE69683.
It demonstrates the interface with reduced computation and does not reproduce
the full benchmark. The standard workflow uses 50 fits, 20 permutations and
20 fits per permutation, and compares alpha 0.5 and 1. One worker is used by
default. No external database download is needed to run the example.

Input can be a `SummarizedExperiment` or a samples-by-genes numeric matrix.
Output includes the complete gene score table, the selected gene names, the
internal alpha comparison, and model summaries. Internal AUC is a tuning
measurement. Separate test samples are required to estimate predictive
performance.

The ranking table reports recurrence, SHAP contribution frequency, mutual
information, raw predictive contribution, the two shuffled-outcome-adjusted
measurements, and the final ranking score. Gene Ontology, Hallmark gene sets,
Open Targets and STRING functions support downstream interpretation.

## Development checks

Source submissions should pass `R CMD build`, `R CMD check`, and
`BiocCheck::BiocCheck()` under the current Bioconductor development release.

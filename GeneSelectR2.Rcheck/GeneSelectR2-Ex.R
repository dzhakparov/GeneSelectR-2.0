pkgname <- "GeneSelectR2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('GeneSelectR2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("asthma_example")
### * asthma_example

flush(stderr()); flush(stdout())

### Name: asthma_example
### Title: Compact U-BIOPRED asthma expression example
### Aliases: asthma_example
### Keywords: datasets

### ** Examples

data(asthma_example)
dim(asthma_example)
table(SummarizedExperiment::colData(asthma_example)$severity)



cleanEx()
nameEx("biological_scorer")
### * biological_scorer

flush(stderr()); flush(stdout())

### Name: biological_scorer
### Title: Compute Biological Relevance Scores
### Aliases: biological_scorer

### ** Examples

if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    biological_scorer(
        c("TP53", "BRCA1"),
        mode = "supervised",
        target_terms = "GO:0006915"
    )
}




cleanEx()
nameEx("cache_info")
### * cache_info

flush(stderr()); flush(stdout())

### Name: cache_info
### Title: Inspect cached annotation data
### Aliases: cache_info

### ** Examples

cache_info()



cleanEx()
nameEx("clear_cache")
### * clear_cache

flush(stderr()); flush(stdout())

### Name: clear_cache
### Title: Clear All Cached Data
### Aliases: clear_cache

### ** Examples

if (interactive()) {
    clear_cache(confirm = FALSE)
}



cleanEx()
nameEx("compute_mutual_information")
### * compute_mutual_information

flush(stderr()); flush(stdout())

### Name: compute_mutual_information
### Title: Compute mutual information with a binary outcome
### Aliases: compute_mutual_information

### ** Examples

x <- c(0, 0, 1, 1, 2, 2)
y <- factor(c("a", "a", "a", "b", "b", "b"))
compute_mutual_information(x, y)



cleanEx()
nameEx("compute_nogueira_stability")
### * compute_nogueira_stability

flush(stderr()); flush(stdout())

### Name: compute_nogueira_stability
### Title: Compute the Nogueira feature-selection stability index
### Aliases: compute_nogueira_stability

### ** Examples

selections <- matrix(
    c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE),
    nrow = 4
)
compute_nogueira_stability(selections)



cleanEx()
nameEx("create_subsamples")
### * create_subsamples

flush(stderr()); flush(stdout())

### Name: create_subsamples
### Title: Create repeated stratified cross-validation splits
### Aliases: create_subsamples

### ** Examples

y <- factor(
    rep(c("control", "case"), each = 10),
    levels = c("control", "case")
)
create_subsamples(y, B = 5, k_folds = 5)



cleanEx()
nameEx("evaluate_gene_set_coherence")
### * evaluate_gene_set_coherence

flush(stderr()); flush(stdout())

### Name: evaluate_gene_set_coherence
### Title: Evaluate the Semantic Coherence of a Selected Gene Set
### Aliases: evaluate_gene_set_coherence

### ** Examples

if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    evaluate_gene_set_coherence(c("TP53", "BRCA1"), verbose = FALSE)
}



cleanEx()
nameEx("geneselectr2_fit")
### * geneselectr2_fit

flush(stderr()); flush(stdout())

### Name: geneselectr2_fit
### Title: Fit the GeneSelectR ranking model
### Aliases: geneselectr2_fit

### ** Examples

set.seed(1)
X <- matrix(rnorm(40 * 12), nrow = 40)
colnames(X) <- paste0("gene", seq_len(ncol(X)))
y <- factor(
    rep(c("control", "case"), each = 20),
    levels = c("control", "case")
)
signal_columns <- seq_len(2L)
X[y == "case", signal_columns] <- X[y == "case", signal_columns] + 1
fit <- geneselectr2_fit(
    X, y,
    B = 2, permutations = 1, null_B = 2, verbose = FALSE
)
head(fit$gene_scores)



cleanEx()
nameEx("get_disease_seeds_opentargets")
### * get_disease_seeds_opentargets

flush(stderr()); flush(stdout())

### Name: get_disease_seeds_opentargets
### Title: Derive Disease Seed Genes from Open Targets
### Aliases: get_disease_seeds_opentargets

### ** Examples

if (interactive()) {
    get_disease_seeds_opentargets(
        "EFO_0000274",
        max_seeds = 5, verbose = FALSE
    )
}



cleanEx()
nameEx("load_go_cache")
### * load_go_cache

flush(stderr()); flush(stdout())

### Name: load_go_cache
### Title: Load GO Annotations with Caching
### Aliases: load_go_cache

### ** Examples

if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    go_cache <- load_go_cache("human")
    go_cache[["TP53"]]
}




cleanEx()
nameEx("multilayer_bio_scorer")
### * multilayer_bio_scorer

flush(stderr()); flush(stdout())

### Name: multilayer_bio_scorer
### Title: Multi-Layer Biological Scorer
### Aliases: multilayer_bio_scorer

### ** Examples

if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    multilayer_bio_scorer(
        c("TP53", "BRCA1"),
        layers = "go",
        target_terms = "GO:0006915", verbose = FALSE
    )
}




cleanEx()
nameEx("plot_biology_comparison")
### * plot_biology_comparison

flush(stderr()); flush(stdout())

### Name: plot_biology_comparison
### Title: Plot biological assessment against matched random gene sets
### Aliases: plot_biology_comparison

### ** Examples

data(benchmark_biology)
plot_biology_comparison(benchmark_biology, method = "GeneSelectR")



cleanEx()
nameEx("plot_gene_evidence")
### * plot_gene_evidence

flush(stderr()); flush(stdout())

### Name: plot_gene_evidence
### Title: Plot gene-level predictive and biological evidence
### Aliases: plot_gene_evidence

### ** Examples

data(asthma_case_study)
plot_gene_evidence(asthma_case_study, comparison_label = "DGE")



cleanEx()
nameEx("plot_gene_ranking")
### * plot_gene_ranking

flush(stderr()); flush(stdout())

### Name: plot_gene_ranking
### Title: Plot gene-level GeneSelectR measurements
### Aliases: plot_gene_ranking

### ** Examples

data(asthma_example)
fit <- select_genes(asthma_example, "severity",
    n_genes = 10,
    alpha = 0.5, B = 2, permutations = 1, null_B = 2
)
plot_gene_ranking(fit, n = 8)



cleanEx()
nameEx("score_network_layer")
### * score_network_layer

flush(stderr()); flush(stdout())

### Name: score_network_layer
### Title: Score Genes by Network Proximity to Disease Seeds
### Aliases: score_network_layer

### ** Examples

if (interactive()) {
    score_network_layer(
        c("IL6", "STAT3"), "EFO_0000274",
        verbose = FALSE
    )
}



cleanEx()
nameEx("score_semantic_layer")
### * score_semantic_layer

flush(stderr()); flush(stdout())

### Name: score_semantic_layer
### Title: Score Genes by Semantic Similarity to Disease Biology
### Aliases: score_semantic_layer

### ** Examples

if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    score_semantic_layer(
        c("TP53", "BRCA1"), "GO:0006915",
        verbose = FALSE
    )
}



cleanEx()
nameEx("select_genes")
### * select_genes

flush(stderr()); flush(stdout())

### Name: select_genes
### Title: Run the predictive GeneSelectR workflow
### Aliases: select_genes

### ** Examples

data(asthma_example)
fit <- select_genes(asthma_example, "severity",
    n_genes = 10,
    alpha = 0.5, B = 2, permutations = 1, null_B = 2
)
head(fit$selected_genes)



cleanEx()
nameEx("set_cache_options")
### * set_cache_options

flush(stderr()); flush(stdout())

### Name: set_cache_options
### Title: Set Cache Options
### Aliases: set_cache_options

### ** Examples

set_cache_options(verbose = FALSE)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

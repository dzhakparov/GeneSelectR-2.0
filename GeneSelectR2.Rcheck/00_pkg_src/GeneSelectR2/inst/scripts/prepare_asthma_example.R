# Reconstruct the compact example from the benchmark's prepared GEO files.
# Run from the benchmark repository root. No download or model fitting occurs.
config <- list(
    expression = "data/GSE69683/expression_prepared.csv",
    metadata = "data/GSE69683/metadata_prepared.csv",
    destination = "package/GeneSelectR/data/asthma_example.rda"
)
expression <- as.matrix(read.csv(config$expression, row.names = 1,
                                check.names = FALSE))
metadata <- read.csv(config$metadata, stringsAsFactors = FALSE)
stopifnot(identical(colnames(expression), metadata$sample_id))
# Sampling uses identifiers and gene order only, with no fitted importance or
# outcome-association criterion. The subset is for usage examples only.
samples <- unlist(lapply(c("moderate", "severe"), function(level) {
    head(sort(metadata$sample_id[metadata$outcome == level]), 48)
}), use.names = FALSE)
genes <- sort(rownames(expression))
genes <- genes[unique(round(seq(1, length(genes), length.out = 300)))]
meta <- metadata[match(samples, metadata$sample_id), ]
asthma_example <- SummarizedExperiment::SummarizedExperiment(
    assays = list(log_expression = expression[genes, samples]),
    colData = S4Vectors::DataFrame(
        severity = factor(meta$outcome, levels = c("moderate", "severe")),
        row.names = samples
    ),
    metadata = list(
        accession = "GSE69683", platform = "GPL13158",
        source = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69683",
        scale = "Submitter-normalized log-scale expression",
        selection = "48 sample IDs per class; 300 evenly spaced sorted symbols",
        purpose = "Executable usage example; not the full benchmark"
    )
)
dir.create(dirname(config$destination), recursive = TRUE, showWarnings = FALSE)
save(asthma_example, file = config$destination, compress = "xz")

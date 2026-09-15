CONFIG <- list(
    case_study = "output/biostec_2027/latex/verification/data/asthma_eight_completed_OT.csv",
    biology = "output/biostec_2027/latex/display_options/biology_summary_source.csv",
    open_targets = "output/manuscript/2026-09-07/source_data/dev7_dense_OT_by_dataset.csv",
    output = "package/GeneSelectR/data"
)

case_source <- read.csv(CONFIG$case_study, check.names = FALSE)
asthma_case_study <- data.frame(
    gene = case_source$gene,
    selection_count = case_source$n_outer_top20,
    comparison_count = case_source$DGE_top20_count,
    recurrence = case_source$mean_geneselectr_internal_recurrence,
    adjusted_recurrence = case_source$mean_geneselectr_calibrated_recurrence,
    contribution = case_source$mean_geneselectr_calibrated_utility,
    association = case_source$OT_association_score,
    stringsAsFactors = FALSE
)

biology_source <- read.csv(CONFIG$biology, check.names = FALSE)
open_targets_source <- read.csv(CONFIG$open_targets, check.names = FALSE)
open_targets_wide <- reshape(
    open_targets_source[, c("dataset", "method", "cutoff",
                            "open_targets_enrichment")],
    idvar = c("dataset", "method"), timevar = "cutoff", direction = "wide"
)
names(open_targets_wide)[names(open_targets_wide) ==
    "open_targets_enrichment.0.05"] <- "open_targets_05_ratio"
names(open_targets_wide)[names(open_targets_wide) ==
    "open_targets_enrichment.0.1"] <- "open_targets_10_ratio"

benchmark_biology <- merge(
    biology_source[, c("dataset", "method", "nogueira_stability",
                       "GO_semantic_enrichment", "hallmark_enrichment")],
    open_targets_wide, by = c("dataset", "method"), all = FALSE, sort = FALSE
)
names(benchmark_biology) <- c(
    "dataset", "method", "stability", "go_ratio", "hallmark_ratio",
    "open_targets_05_ratio", "open_targets_10_ratio"
)
method_names <- c(
    GS_full_ungrouped = "GeneSelectR", RF_importance = "Random forest",
    ElasticNet = "Elastic net", Boruta = "Boruta", DGE = "DGE",
    LASSO = "LASSO", mRMR = "mRMR"
)
dataset_names <- c(
    GSE101794 = "Crohn disease", GSE107994 = "Tuberculosis",
    GSE13355 = "Psoriasis", GSE65682 = "Sepsis", GSE69683 = "Asthma",
    imvigor210 = "Bladder cancer", sosall = "Atopic dermatitis"
)
benchmark_biology$method <- unname(method_names[benchmark_biology$method])
benchmark_biology$dataset <- unname(dataset_names[benchmark_biology$dataset])
stopifnot(!anyNA(benchmark_biology), nrow(benchmark_biology) == 49L)

dir.create(CONFIG$output, recursive = TRUE, showWarnings = FALSE)
save(asthma_case_study, file = file.path(CONFIG$output, "asthma_case_study.rda"),
     compress = "xz")
save(benchmark_biology, file = file.path(CONFIG$output, "benchmark_biology.rda"),
     compress = "xz")

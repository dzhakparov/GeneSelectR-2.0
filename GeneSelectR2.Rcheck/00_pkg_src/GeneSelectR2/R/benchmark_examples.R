#' Verified asthma gene-level example
#'
#' Gene-level measurements for the eight genes most often present in the
#' GeneSelectR top-20 sets in the GSE69683 outer validation analysis. The data
#' contain selection counts over 15 outer divisions, mean within-training
#' recurrence, mean permutation-adjusted predictive contribution and continuous
#' Open Targets asthma association scores. `NA` records an unavailable Open
#' Targets score. `comparison_count` gives the corresponding DGE top-20 count.
#'
#' @format A data frame with 8 rows and 7 columns.
#' @source GeneSelectR benchmark result extract; Open Targets release 26.06,
#'   retrieved 10 September 2026.
"asthma_case_study"

#' Verified biological assessment example
#'
#' Biological and stability summaries for seven transcriptomic datasets and
#' seven feature-selection methods. GO, Hallmark and Open Targets values are
#' observed-to-random ratios calculated with same-size gene sets from the
#' training-specific candidate genes. Stability is the mean Nogueira measure.
#' Biological information was applied after predictive ranking.
#'
#' @format A data frame with 49 rows and 7 columns: dataset, method, stability,
#'   go_ratio, hallmark_ratio, open_targets_05_ratio and
#'   open_targets_10_ratio.
#' @source GeneSelectR benchmark result tables verified 10 September 2026.
"benchmark_biology"

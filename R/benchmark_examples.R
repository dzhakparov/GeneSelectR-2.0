#' Asthma gene-level result example
#'
#' Gene-level measurements for the eight genes most often present in the
#' GeneSelectR 2.0 top-20 sets in the GSE69683 outer validation analysis.
#'
#' @format A data frame with 8 rows and 7 columns:
#' \describe{
#'   \item{gene}{Gene symbol.}
#'   \item{selection_count}{Number of GeneSelectR 2.0 top-20 sets containing the
#'     gene across 15 outer validation divisions.}
#'   \item{comparison_count}{Number of differential-expression top-20 sets
#'     containing the gene across 15 outer validation divisions.}
#'   \item{recurrence}{Mean within-training selection recurrence.}
#'   \item{adjusted_recurrence}{Mean recurrence divided by its shuffled-outcome
#'     reference.}
#'   \item{contribution}{Mean predictive contribution divided by its
#'     shuffled-outcome reference.}
#'   \item{association}{Open Targets asthma association score from release
#'     26.06. Missing annotation is recorded as `NA`.}
#' }
#' @source GeneSelectR 2.0 benchmark result extract; Open Targets release 26.06,
#'   retrieved 10 September 2026.
"asthma_case_study"

#' Biological assessment example
#'
#' Biological and stability summaries for seven transcriptomic datasets and
#' seven feature-selection methods. GO, Hallmark and Open Targets values are
#' observed-to-random ratios calculated with same-size gene sets from the
#' training-specific candidate genes. Stability is the mean Nogueira measure.
#' Biological information was applied after predictive ranking.
#'
#' @format A data frame with 49 rows and 7 columns:
#' \describe{
#'   \item{dataset}{Study phenotype or disease label.}
#'   \item{method}{Feature-selection method.}
#'   \item{stability}{Mean Nogueira feature-selection stability.}
#'   \item{go_ratio}{GO semantic coherence divided by the matched random-set
#'     mean.}
#'   \item{hallmark_ratio}{Hallmark sharing divided by the matched random-set
#'     mean.}
#'   \item{open_targets_05_ratio}{Open Targets gene count at score 0.05 divided
#'     by the matched random-set mean.}
#'   \item{open_targets_10_ratio}{Open Targets gene count at score 0.10 divided
#'     by the matched random-set mean.}
#' }
#' @source GeneSelectR 2.0 benchmark result tables extracted 10 September 2026.
"benchmark_biology"

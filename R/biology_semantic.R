## GO semantic similarity against externally specified disease terms.
## Rectangular matrices support disease scoring. Square matrices support
## downstream gene-set coherence calculations.

#' Build a candidate-term by target-term similarity matrix
#'
#' Computes the exact rectangular block used by disease-target semantic
#' scoring. Similarities among candidate terms are omitted because they do not
#' enter the best-match-average calculation against external target terms.
#'
#' @param candidate_terms GO terms annotating the candidate genes
#' @param target_terms Externally specified disease GO terms
#' @param ic_scores Named numeric vector of information content
#' @param ancestor_map Named list of term to ancestor terms
#' @param method Similarity metric (resnik/lin/jiang/rel)
#' @param ontology Ontology branch, used in the cache key
#' @param use_cache Read and write the disk cache
#' @param verbose Print progress
#' @return Numeric matrix with candidate terms in rows and target terms in
#'   columns
#' @keywords internal
build_target_similarity_matrix <- function(candidate_terms,
                                            target_terms,
                                            ic_scores, ancestor_map,
                                            method = "resnik",
                                            ontology = "BP",
                                            use_cache = TRUE,
                                            verbose = TRUE) {
    row_terms <- sort(unique(candidate_terms[
        candidate_terms %in% names(ic_scores)
    ]))
    column_terms <- sort(unique(target_terms[
        target_terms %in% names(ic_scores)
    ]))

    if (length(row_terms) == 0L || length(column_terms) == 0L) {
        return(matrix(0,
            nrow = length(row_terms), ncol = length(column_terms),
            dimnames = list(row_terms, column_terms)
        ))
    }

    cache_file <- NULL
    if (use_cache) {
        if (!requireNamespace("digest", quietly = TRUE)) {
            warning(
                "digest is unavailable; target matrices will not be cached",
                call. = FALSE
            )
            use_cache <- FALSE
        } else {
            relevant_terms <- sort(unique(c(row_terms, column_terms)))
            key_hash <- digest::digest(
                list(
                    method = method,
                    ontology = ontology,
                    row_terms = row_terms,
                    column_terms = column_terms,
                    ic_scores = ic_scores[relevant_terms],
                    ancestors = ancestor_map[relevant_terms]
                ),
                algo = "md5"
            )
            cache_file <- file.path(
                get_cache_dir(), sprintf("term_sim_target_v1_%s.rds", key_hash)
            )
        }

        if (!is.null(cache_file)) {
            cached <- .read_cache_file(cache_file)
            if (!is.null(cached) &&
                identical(rownames(cached), row_terms) &&
                identical(colnames(cached), column_terms)) {
                if (verbose) {
                    message("    Loading the cached target-term matrix\n")
                }
                return(cached)
            }
            if (verbose) message("    Target-term cache mismatch; rebuilding\n")
        }
    }

    if (verbose) {
        message(sprintf(
            "    Building %d x %d candidate-target term matrix (%d pairs)\n",
            length(row_terms), length(column_terms),
            length(row_terms) * length(column_terms)
        ))
    }

    similarity_matrix <- matrix(
        0,
        nrow = length(row_terms), ncol = length(column_terms),
        dimnames = list(row_terms, column_terms)
    )
    for (i in seq_along(row_terms)) {
        for (j in seq_along(column_terms)) {
            similarity_matrix[i, j] <- compute_semantic_similarity(
                row_terms[i], column_terms[j], ic_scores, ancestor_map, method
            )
        }
    }

    if (use_cache && !is.null(cache_file)) {
        .write_cache_file(similarity_matrix, cache_file)
        if (verbose) {
            message("    Cached target-term similarity matrix to disk\n")
        }
    }

    similarity_matrix
}

#' Build a square term-by-term similarity matrix
#'
#' Computes pairwise semantic similarity across a term universe. Downstream
#' comparisons use matrix subsets for the required term sets. The cache key
#' includes the term universe, ontology, similarity method, information content,
#' and ancestor relationships.
#'
#' @param term_universe Character vector of all GO term IDs that any comparison
#'   might touch (candidate-gene terms, plus any disease target terms).
#' @param ic_scores Named numeric vector of information content
#' @param ancestor_map Named list of term -> ancestors
#' @param method Similarity metric (resnik/lin/jiang/rel)
#' @param ontology Ontology branch, used only in the cache key
#' @param use_cache Read/write the disk cache
#' @param verbose Print progress
#' @return Symmetric numeric matrix with dimnames = list(terms, terms)
#' @keywords internal
build_term_similarity_matrix <- function(term_universe,
                                            ic_scores, ancestor_map,
                                            method = "resnik",
                                            ontology = "BP",
                                            use_cache = TRUE,
                                            verbose = TRUE) {
    terms <- sort(unique(term_universe[term_universe %in% names(ic_scores)]))
    n_terms <- length(terms)

    if (n_terms == 0) {
        return(matrix(0, 0, 0))
    }

    cache_file <- NULL
    if (use_cache) {
        if (!requireNamespace("digest", quietly = TRUE)) {
            warning(
                "digest is unavailable; term matrices will not be cached",
                call. = FALSE
            )
            use_cache <- FALSE
        } else {
            key_hash <- digest::digest(
                list(
                    method = method,
                    ontology = ontology,
                    terms = terms,
                    ic_scores = ic_scores[terms],
                    ancestors = ancestor_map[terms]
                ),
                algo = "md5"
            )

            cache_file <- file.path(
                get_cache_dir(),
                sprintf("term_sim_square_v2_%s.rds", key_hash)
            )
        }

        if (!is.null(cache_file)) {
            if (verbose) message("    Loading cached term similarity matrix\n")
            cached <- .read_cache_file(cache_file)
            if (!is.null(cached) && nrow(cached) == n_terms &&
                identical(rownames(cached), terms)) {
                return(cached)
            }
            if (verbose) message("    Cache mismatch, rebuilding\n")
        }
    }

    n_pairs <- n_terms * (n_terms - 1) / 2
    if (verbose) {
        message(sprintf(
            "    Building %d x %d square term matrix (%.0f unique pairs)\n",
            n_terms, n_terms, n_pairs
        ))
        message("    The computed matrix is stored in the package cache.\n")
    }

    sim_matrix <- matrix(0, n_terms, n_terms, dimnames = list(terms, terms))
    diag(sim_matrix) <- 1

    ## Similarity is symmetric, so only the upper triangle is computed.
    for (i in seq_len(n_terms - 1)) {
        for (j in (i + 1):n_terms) {
            pair_similarity <- compute_semantic_similarity(
                terms[i], terms[j],
                ic_scores, ancestor_map, method
            )
            sim_matrix[i, j] <- pair_similarity
            sim_matrix[j, i] <- pair_similarity
        }
        if (verbose && i %% 250 == 0) {
            message(sprintf(
                "      %d/%d terms (%.1f%%)\n",
                i, n_terms, 100 * i / n_terms
            ))
        }
    }

    if (use_cache && !is.null(cache_file)) {
        .write_cache_file(sim_matrix, cache_file)
        if (verbose) message("    Cached term similarity matrix to disk\n")
    }

    sim_matrix
}


#' Best-Match Average From a Precomputed Term Matrix
#'
#' Calculates the mean of row maxima and column maxima from the matrix subset
#' defined by two term sets.
#'
#' @param terms_a Terms indexing the matrix rows
#' @param terms_b Terms indexing the matrix columns
#' @param term_sim_matrix Precomputed similarity matrix
#' @return Numeric BMA similarity
#' @keywords internal
bma_from_matrix <- function(terms_a, terms_b, term_sim_matrix) {
    terms_a <- terms_a[terms_a %in% rownames(term_sim_matrix)]
    terms_b <- terms_b[terms_b %in% colnames(term_sim_matrix)]

    if (length(terms_a) == 0 || length(terms_b) == 0) {
        return(0)
    }

    sub_matrix <- term_sim_matrix[terms_a, terms_b, drop = FALSE]

    best_a_to_b <- apply(sub_matrix, 1, max)
    best_b_to_a <- apply(sub_matrix, 2, max)

    mean(c(best_a_to_b, best_b_to_a))
}

.load_semantic_resources <- function(organism, ontology, use_cache) {
    go_cache <- load_go_cache(organism = organism)
    go_cache <- filter_go_cache_by_ontology(go_cache, ontology = ontology)
    list(
        go_cache = go_cache,
        ic_scores = load_ic_cache(go_cache),
        ancestor_map = load_ancestor_map(
            organism = organism,
            use_cache = use_cache
        )
    )
}

.annotated_go_genes <- function(genes, go_cache) {
    genes <- intersect(genes, names(go_cache))
    genes[vapply(
        go_cache[genes],
        function(terms) !is.null(terms) && length(terms) > 0L,
        logical(1)
    )]
}

.score_gene_terms <- function(genes, gene_terms, target_terms, term_matrix) {
    scores <- stats::setNames(rep(0, length(genes)), genes)
    for (gene in names(gene_terms)) {
        scores[gene] <- bma_from_matrix(
            gene_terms[[gene]],
            target_terms,
            term_matrix
        )
    }
    scores
}

.report_semantic_scores <- function(scores, n_annotated, n_targets) {
    message(sprintf(
        "  Semantic: %d annotated genes vs %d target terms\n",
        n_annotated,
        n_targets
    ))
    message(sprintf(
        "    %d/%d genes received similarity scores\n",
        sum(scores > 0),
        length(scores)
    ))
}


#' Score Genes by Semantic Similarity to Disease Biology
#'
#' Scores each candidate gene by the best-match-average GO semantic similarity
#' between its own GO terms and a user-supplied set of disease target terms.
#'
#' Targets are specified externally. Similarity within a selected gene set is
#' evaluated separately with \code{\link{evaluate_gene_set_coherence}}.
#'
#' @param genes Character vector of candidate gene symbols to score
#' @param target_terms Character vector of externally specified disease GO term
#'   identifiers.
#' @param ontology GO ontology branch (default "BP")
#' @param sim_method Term-pair similarity metric (default "resnik")
#' @param organism Organism for the GO cache (default "human")
#' @param use_cache Logical - cache the term matrix to disk
#' @param verbose Logical
#' @return Named numeric vector of percentile-normalised similarity scores
#' @examples
#' if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'     score_semantic_layer(
#'         c("TP53", "BRCA1"), "GO:0006915",
#'         verbose = FALSE
#'     )
#' }
#' @export
score_semantic_layer <- function(genes,
                                    target_terms,
                                    ontology = "BP",
                                    sim_method = "resnik",
                                    organism = "human",
                                    use_cache = TRUE,
                                    verbose = TRUE) {
    n_genes <- length(genes)
    zero_scored <- stats::setNames(rep(0, n_genes), genes)

    if (is.null(target_terms) || length(target_terms) == 0) {
        stop(
            "score_semantic_layer requires disease GO term IDs in target_terms",
            call. = FALSE
        )
    }

    resources <- .load_semantic_resources(organism, ontology, use_cache)
    annotated_genes <- .annotated_go_genes(genes, resources$go_cache)
    if (length(annotated_genes) == 0L) {
        warning("No candidate genes have GO annotations; returning zeros.")
        return(zero_scored)
    }

    gene_terms <- resources$go_cache[annotated_genes]
    candidate_terms <- unique(unlist(gene_terms, use.names = FALSE))
    term_sim <- build_target_similarity_matrix(
        candidate_terms = candidate_terms,
        target_terms = target_terms,
        ic_scores = resources$ic_scores,
        ancestor_map = resources$ancestor_map,
        method = sim_method, ontology = ontology,
        use_cache = use_cache, verbose = verbose
    )
    if (nrow(term_sim) == 0) {
        return(zero_scored)
    }

    raw_scores <- .score_gene_terms(
        genes, gene_terms, target_terms, term_sim
    )

    if (verbose) {
        .report_semantic_scores(
            raw_scores,
            length(annotated_genes),
            length(target_terms)
        )
    }

    percentile01(raw_scores)
}

.gene_set_coherence <- function(gene_set, go_cache, term_matrix) {
    n_genes <- length(gene_set)
    if (n_genes < 2L) {
        return(NA_real_)
    }
    term_sets <- go_cache[gene_set]
    pair_values <- numeric(n_genes * (n_genes - 1L) / 2L)
    pair_index <- 1L
    for (i in seq_len(n_genes - 1L)) {
        for (j in (i + 1L):n_genes) {
            pair_values[pair_index] <- bma_from_matrix(
                term_sets[[i]], term_sets[[j]], term_matrix
            )
            pair_index <- pair_index + 1L
        }
    }
    mean(pair_values, na.rm = TRUE)
}

.coherence_null <- function(
    background_genes,
    selected_genes,
    go_cache,
    term_matrix,
    n_permutations
) {
    background <- setdiff(
        .annotated_go_genes(background_genes, go_cache),
        selected_genes
    )
    if (length(background) < length(selected_genes)) {
        return(list(expected = NA_real_, empirical_p = NA_real_))
    }
    null_values <- replicate(
        n_permutations,
        .gene_set_coherence(
            sample(background, length(selected_genes)),
            go_cache,
            term_matrix
        )
    )
    list(
        expected = mean(null_values, na.rm = TRUE),
        null_values = null_values
    )
}

.prepare_coherence_context <- function(
    gene_set,
    background_genes,
    resources,
    sim_method,
    ontology,
    use_cache,
    verbose
) {
    selected <- .annotated_go_genes(gene_set, resources$go_cache)
    if (length(selected) < 2L) {
        warning(
            "Fewer than two selected genes have GO annotations.",
            call. = FALSE
        )
        return(NULL)
    }
    universe_genes <- if (is.null(background_genes)) {
        selected
    } else {
        unique(c(
            selected,
            .annotated_go_genes(background_genes, resources$go_cache)
        ))
    }
    term_universe <- unique(unlist(
        resources$go_cache[universe_genes],
        use.names = FALSE
    ))
    term_matrix <- build_term_similarity_matrix(
        term_universe = term_universe,
        ic_scores = resources$ic_scores,
        ancestor_map = resources$ancestor_map,
        method = sim_method,
        ontology = ontology,
        use_cache = use_cache,
        verbose = verbose
    )
    if (nrow(term_matrix) == 0L) {
        return(NULL)
    }
    list(selected = selected, term_matrix = term_matrix)
}


#' Evaluate the Semantic Coherence of a Selected Gene Set
#'
#' Computes mean pairwise best-match-average GO semantic similarity among
#' selected genes. This downstream assessment does not change the ranking.
#'
#' To contextualise the value, an empirical null is estimated by drawing random
#' gene sets of the same size from the background and computing their coherence.
#'
#' @param gene_set Character vector containing the selected genes.
#' @param background_genes Character vector to draw the null from (typically all
#'   candidate genes). If NULL, no null is computed.
#' @param n_permutations Number of random gene sets used for comparison.
#' @param ontology GO ontology branch (default "BP")
#' @param sim_method Term-pair similarity metric (default "resnik")
#' @param organism Organism for the GO cache (default "human")
#' @param use_cache Logical
#' @param random_seed Integer seed for the random-gene-set comparison.
#' @param verbose Logical
#' @return List with: mean_coherence, n_annotated, expected_coherence,
#'   empirical_p (NA if no background supplied)
#' @examples
#' if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'     evaluate_gene_set_coherence(c("TP53", "BRCA1"), verbose = FALSE)
#' }
#' @export
evaluate_gene_set_coherence <- function(
    gene_set,
    background_genes = NULL,
    n_permutations = 50,
    ontology = "BP",
    sim_method = "resnik",
    organism = "human",
    use_cache = TRUE,
    random_seed = 123,
    verbose = TRUE
) {
    withr::local_seed(random_seed)

    empty_result <- list(
        mean_coherence = NA_real_, n_annotated = 0L,
        expected_coherence = NA_real_, empirical_p = NA_real_
    )

    resources <- .load_semantic_resources(organism, ontology, use_cache)
    context <- .prepare_coherence_context(
        gene_set, background_genes, resources, sim_method,
        ontology, use_cache, verbose
    )
    if (is.null(context)) {
        return(empty_result)
    }
    selected_annotated <- context$selected
    term_sim <- context$term_matrix

    if (verbose) {
        message(sprintf(
            "  Gene-set coherence: %d annotated genes\n",
            length(selected_annotated)
        ))
    }
    observed <- .gene_set_coherence(
        selected_annotated,
        resources$go_cache,
        term_sim
    )

    expected <- NA_real_
    empirical_p <- NA_real_

    if (!is.null(background_genes)) {
        if (verbose) {
            message(sprintf(
                "  Comparing with %d random gene sets...\n",
                n_permutations
            ))
        }
        null_result <- .coherence_null(
            background_genes,
            selected_annotated,
            resources$go_cache,
            term_sim,
            n_permutations
        )
        expected <- null_result$expected
        if (!is.null(null_result$null_values)) {
            empirical_p <- (
                sum(null_result$null_values >= observed, na.rm = TRUE) + 1
            ) / (n_permutations + 1)
        }
    }

    list(
        mean_coherence = observed,
        n_annotated = length(selected_annotated),
        expected_coherence = expected,
        empirical_p        = empirical_p
    )
}

## Gene Ontology semantic similarity, enrichment, and MSigDB membership.

.validate_biological_inputs <- function(
    genes,
    mode,
    ontology,
    sim_method,
    enrich_fdr,
    min_term_freq,
    max_enriched_terms,
    n_top_sims,
    ic_quantile
) {
    if (!is.character(genes) || length(genes) == 0L || anyNA(genes)) {
        stop(
            "genes must be a non-empty character vector without missing values"
        )
    }
    mode <- match.arg(mode, c("supervised", "data_driven"))
    ontology <- match.arg(ontology, c("BP", "MF", "CC"), several.ok = TRUE)
    sim_method <- match.arg(sim_method, c("resnik", "lin", "jiang", "rel"))
    .validate_range(enrich_fdr, "enrich_fdr", 0, 1, lower_open = TRUE)
    .validate_range(ic_quantile, "ic_quantile", 0, 1, upper_open = TRUE)
    .validate_positive(n_top_sims, "n_top_sims")
    .validate_positive(max_enriched_terms, "max_enriched_terms")
    if (!is.null(min_term_freq)) {
        .validate_positive(min_term_freq, "min_term_freq")
    }
    list(mode = mode, ontology = ontology, sim_method = sim_method)
}

.validate_positive <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 1) {
        stop(name, " must be one positive number", call. = FALSE)
    }
}

.validate_range <- function(
    value,
    name,
    lower,
    upper,
    lower_open = FALSE,
    upper_open = FALSE
) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        stop(name, " must be one finite number", call. = FALSE)
    }
    lower_invalid <- if (lower_open) value <= lower else value < lower
    upper_invalid <- if (upper_open) value >= upper else value > upper
    if (lower_invalid || upper_invalid) {
        interval <- sprintf(
            "%s%s, %s%s",
            if (lower_open) "(" else "[",
            lower,
            upper,
            if (upper_open) ")" else "]"
        )
        stop(name, " must be one finite number in ", interval, call. = FALSE)
    }
}

.prepare_biological_resources <- function(
    genes,
    go_cache,
    ontology,
    organism,
    use_cache
) {
    if (is.null(go_cache)) {
        go_cache <- if (use_cache) {
            load_go_cache(organism = organism)
        } else {
            download_go_annotations(organism)
        }
    }
    annotated <- intersect(genes, names(go_cache))
    if (length(annotated) == 0L) {
        warning(
            "No GO annotations were found; returning zero evidence scores",
            call. = FALSE
        )
        return(NULL)
    }
    if (length(annotated) < length(genes) * 0.5) {
        warning(sprintf(
            "Only %.1f%% of genes have GO annotations (%d/%d)",
            100 * length(annotated) / length(genes),
            length(annotated),
            length(genes)
        ))
    }
    go_cache <- filter_go_cache_by_ontology(go_cache, ontology)
    if (!any(genes %in% names(go_cache))) {
        warning(sprintf(
            "No genes have GO annotations in ontology [%s]",
            paste(ontology, collapse = ", ")
        ))
        return(NULL)
    }
    list(
        go_cache = go_cache,
        ic_scores = compute_information_content(go_cache),
        ancestor_map = load_ancestor_map(organism, use_cache),
        similarity_cache = create_similarity_cache()
    )
}

.compute_biological_scores <- function(
    genes,
    target_terms,
    enrichment_genes,
    mode,
    resources,
    settings
) {
    if (mode == "supervised") {
        return(compute_supervised_scores(
            genes = genes,
            target_terms = target_terms,
            go_cache = resources$go_cache,
            ic_scores = resources$ic_scores,
            similarity_cache = resources$similarity_cache,
            ancestor_map = resources$ancestor_map,
            sim_method = settings$sim_method,
            n_top_sims = settings$n_top_sims
        ))
    }
    compute_data_driven_scores(
        genes = genes,
        enrichment_genes = enrichment_genes,
        go_cache = resources$go_cache,
        ic_scores = resources$ic_scores,
        similarity_cache = resources$similarity_cache,
        ancestor_map = resources$ancestor_map,
        sim_method = settings$sim_method,
        enrich_fdr = settings$enrich_fdr,
        min_term_freq = settings$min_term_freq,
        max_enriched_terms = settings$max_enriched_terms,
        n_top_sims = settings$n_top_sims,
        ic_quantile = settings$ic_quantile
    )
}


#' Compute Biological Relevance Scores
#'
#' Scores genes by GO semantic similarity to specified target terms or to terms
#' identified by an enrichment analysis.
#'
#' @param genes Character vector of gene symbols to score.
#' @param mode "supervised" or "data_driven"
#' @param target_terms GO term IDs for supervised mode (e.g., "GO:0006955")
#' @param ontology Character vector of GO ontologies to use. Any subset of
#'   c("BP", "MF", "CC"). Default: "BP".
#' @param sim_method Semantic similarity metric. One of: "resnik", "lin",
#'   "jiang", "rel".
#' @param enrich_fdr FDR threshold for Fisher enrichment in data-driven mode.
#' @param min_term_freq Minimum annotation frequency for candidate terms.
#' @param max_enriched_terms Maximum number of enriched terms used as targets.
#' @param n_top_sims Top-k similarities to average per gene.
#' @param ic_quantile Specificity filter quantile for IC.
#' @param enrichment_genes Character vector of gene symbols used to identify
#'   enriched GO terms in data-driven mode. The default is `genes`.
#' @param go_cache Pre-loaded GO cache (gene -> GO terms). If NULL, loaded.
#' @param use_cache Use disk/memory caching.
#' @param organism Organism name.
#'
#' @return Numeric vector of biological scores (one per gene, same order as
#'   input), percentile-normalized. A vector of zeros indicates that the
#'   requested genes have no annotation evidence.
#' @examples
#' if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'     biological_scorer(
#'         c("TP53", "BRCA1"),
#'         mode = "supervised",
#'         target_terms = "GO:0006915"
#'     )
#' }
#'
#' @export
biological_scorer <- function(
    genes,
    mode = "supervised",
    target_terms = NULL,
    ontology = "BP",
    sim_method = "resnik",
    enrich_fdr = 0.05,
    min_term_freq = NULL,
    max_enriched_terms = 100,
    n_top_sims = 5,
    ic_quantile = 0.5,
    enrichment_genes = NULL,
    go_cache = NULL,
    use_cache = TRUE,
    organism = "human"
) {
    validated <- .validate_biological_inputs(
        genes, mode, ontology, sim_method, enrich_fdr,
        min_term_freq, max_enriched_terms, n_top_sims, ic_quantile
    )
    mode <- validated$mode
    ontology <- validated$ontology
    sim_method <- validated$sim_method
    if (mode == "supervised" &&
        (is.null(target_terms) || length(target_terms) == 0L)) {
        stop("target_terms required for supervised mode")
    }
    if (is.null(enrichment_genes)) {
        enrichment_genes <- genes
    }
    resources <- .prepare_biological_resources(
        genes, go_cache, ontology, organism, use_cache
    )
    if (is.null(resources)) {
        return(rep(0, length(genes)))
    }

    settings <- list(
        sim_method = sim_method,
        enrich_fdr = enrich_fdr,
        min_term_freq = min_term_freq,
        max_enriched_terms = max_enriched_terms,
        n_top_sims = n_top_sims,
        ic_quantile = ic_quantile
    )
    scores <- .compute_biological_scores(
        genes, target_terms, enrichment_genes, mode, resources, settings
    )
    percentile01(scores)
}


#' Filter GO Cache by Ontology
#'
#' Removes GO terms that don't belong to the specified ontologies.
#' Requires GO.db and terminates if the requested ontology cannot be resolved.
#'
#' @param go_cache Named list of gene -> GO terms
#' @param ontology Character vector of ontologies to keep ("BP", "MF", "CC")
#' @return Filtered go_cache (genes with zero remaining terms are dropped)
#' @keywords internal
filter_go_cache_by_ontology <- function(go_cache, ontology = "BP") {
    if (length(ontology) == 3 &&
        all(c("BP", "MF", "CC") %in% ontology)) {
        return(go_cache)
    }

    if (!requireNamespace("GO.db", quietly = TRUE)) {
        stop(
            "GO.db is required to filter GO annotations by ontology. ",
            "Install with: BiocManager::install('GO.db')"
        )
    }

    all_terms <- unique(unlist(go_cache))
    if (length(all_terms) == 0) {
        return(go_cache)
    }

    tryCatch(
        {
            term_info <- AnnotationDbi::select(
                GO.db::GO.db,
                keys = all_terms,
                columns = c("GOID", "ONTOLOGY"),
                keytype = "GOID"
            )
            keep_terms <- term_info$GOID[term_info$ONTOLOGY %in% ontology]

            go_cache_out <- lapply(go_cache, function(terms) {
                intersect(terms, keep_terms)
            })

            go_cache_out <- go_cache_out[lengths(go_cache_out) > 0]

            if (getOption("geneselectr2.verbose", FALSE)) {
                message(sprintf(
                    paste(
                        "  Ontology filter [%s]: %d -> %d terms,",
                        "%d -> %d annotated genes"
                    ),
                    paste(ontology, collapse = "+"),
                    length(all_terms), length(keep_terms),
                    length(go_cache), length(go_cache_out)
                ))
            }

            return(go_cache_out)
        },
        error = function(error) {
            stop(
                "GO.db ontology lookup failed: ", conditionMessage(error),
                call. = FALSE
            )
        }
    )
}


#' Compute Supervised Scores
#'
#' Scores each gene by its mean-of-top-k semantic similarity to a set of
#' user-specified target GO terms.
#'
#' @param genes Gene symbols
#' @param target_terms Target GO term IDs
#' @param go_cache GO annotations (ontology-filtered)
#' @param ic_scores Information content scores
#' @param similarity_cache Cache environment for similarities
#' @param ancestor_map Named list of GO term -> ancestors
#' @param sim_method Similarity metric name
#' @param n_top_sims Number of top similarities to average
#' @return Numeric vector of similarity scores (one per gene)
#' @keywords internal
compute_supervised_scores <- function(genes, target_terms, go_cache, ic_scores,
                                        similarity_cache, ancestor_map = NULL,
                                        sim_method = "resnik",
                                        n_top_sims = 5) {
    n_genes <- length(genes)
    scores <- numeric(n_genes)

    for (i in seq_len(n_genes)) {
        gene <- genes[i]
        gene_terms <- go_cache[[gene]]

        if (is.null(gene_terms) || length(gene_terms) == 0) {
            next
        }
        scores[i] <- .mean_top_similarity(
            gene_terms, target_terms, ic_scores, similarity_cache,
            ancestor_map, sim_method, n_top_sims
        )
    }

    return(scores)
}

.mean_top_similarity <- function(
    source_terms,
    target_terms,
    ic_scores,
    similarity_cache,
    ancestor_map,
    sim_method,
    n_top_sims
) {
    if (length(source_terms) == 0L || length(target_terms) == 0L) {
        return(0)
    }
    similarities <- numeric(length(source_terms) * length(target_terms))
    similarity_index <- 1L
    for (source_term in source_terms) {
        for (target_term in target_terms) {
            similarities[similarity_index] <- get_or_compute_similarity(
                source_term, target_term, ic_scores, similarity_cache,
                ancestor_map, sim_method
            )
            similarity_index <- similarity_index + 1L
        }
    }
    similarities <- sort(similarities, decreasing = TRUE)
    mean(utils::head(similarities, n_top_sims))
}

.specific_go_annotations <- function(
    enrichment_genes,
    go_cache,
    ic_scores,
    ic_quantile
) {
    enrichment_genes <- intersect(enrichment_genes, names(go_cache))
    if (length(enrichment_genes) < 10L) {
        warning(
            "Fewer than 10 enrichment genes are annotated; returning zeros",
            call. = FALSE
        )
        return(NULL)
    }
    all_terms <- unique(unlist(go_cache[enrichment_genes]))
    if (length(all_terms) == 0L) {
        warning(
            "No enrichment genes remain after the ontology filter",
            call. = FALSE
        )
        return(NULL)
    }
    available_ic <- ic_scores[names(ic_scores) %in% all_terms]
    if (length(available_ic) > 0L && ic_quantile > 0) {
        threshold <- stats::quantile(
            available_ic, probs = ic_quantile, na.rm = TRUE
        )
        specific_terms <- names(available_ic[available_ic >= threshold])
        go_cache <- lapply(
            go_cache,
            function(terms) intersect(terms, specific_terms)
        )
        go_cache <- go_cache[lengths(go_cache) > 0L]
    }
    annotated_genes <- intersect(enrichment_genes, names(go_cache))
    if (length(annotated_genes) < 10L) {
        warning(
            "Fewer than 10 enrichment genes pass the specificity filter",
            call. = FALSE
        )
        return(NULL)
    }
    list(go_cache = go_cache, annotated_genes = annotated_genes)
}

.select_data_driven_terms <- function(
    annotated_genes,
    go_cache,
    min_term_freq,
    enrich_fdr,
    max_enriched_terms
) {
    term_counts <- table(unlist(go_cache[annotated_genes]))
    minimum_frequency <- if (is.null(min_term_freq)) {
        max(5, ceiling(0.01 * length(annotated_genes)))
    } else {
        min_term_freq
    }
    target_terms <- names(term_counts[term_counts >= minimum_frequency])
    if (length(target_terms) == 0L) {
        keep <- seq_len(min(20L, length(term_counts)))
        target_terms <- names(sort(term_counts, decreasing = TRUE))[keep]
    }

    background <- setdiff(names(go_cache), annotated_genes)
    if (length(annotated_genes) >= 30L && length(background) >= 50L) {
        enrichment <- test_go_enrichment(
            selected_genes = annotated_genes,
            background_genes = background,
            go_cache = go_cache
        )
        significant <- enrichment$term[enrichment$p_adj < enrich_fdr]
        if (length(significant) >= 5L) {
            target_terms <- significant
        }
    }
    if (length(target_terms) > max_enriched_terms) {
        frequencies <- term_counts[target_terms]
        target_terms <- names(sort(
            frequencies,
            decreasing = TRUE
        ))[seq_len(max_enriched_terms)]
    }
    target_terms
}

.score_data_driven_genes <- function(
    genes,
    target_terms,
    go_cache,
    ic_scores,
    similarity_cache,
    ancestor_map,
    sim_method,
    n_top_sims
) {
    scores <- numeric(length(genes))
    for (i in seq_along(genes)) {
        gene_terms <- go_cache[[genes[i]]]
        if (is.null(gene_terms) || length(gene_terms) == 0L) {
            next
        }
        comparison_terms <- setdiff(target_terms, gene_terms)
        if (length(comparison_terms) == 0L) {
            scores[i] <- min(length(gene_terms) / length(target_terms), 0.5)
            next
        }
        scores[i] <- .mean_top_similarity(
            gene_terms, comparison_terms, ic_scores, similarity_cache,
            ancestor_map, sim_method, n_top_sims
        )
    }
    scores
}

#' Compute data-driven biological relevance scores from GO annotations
#'
#' @param genes Character vector of gene identifiers to score.
#' @param enrichment_genes Genes used to select enriched terms.
#' @param go_cache Named list mapping gene -> GO terms.
#' @param ic_scores Named numeric vector of GO term information content.
#' @param similarity_cache Cache used for semantic similarity computation.
#' @param ancestor_map Optional GO ancestor mapping.
#' @param sim_method Semantic similarity method (e.g. "resnik").
#' @param enrich_fdr Adjusted p-value threshold for enrichment.
#' @param min_term_freq Minimum term frequency threshold.
#' @param max_enriched_terms Maximum number of target terms retained.
#' @param n_top_sims Number of top similarities averaged per gene.
#' @param ic_quantile IC quantile threshold for specificity filtering.
#'
#' @return Numeric vector of biological relevance scores.
#'
#' @keywords internal
#'
compute_data_driven_scores <- function(genes,
                                        enrichment_genes,
                                        go_cache,
                                        ic_scores,
                                        similarity_cache,
                                        ancestor_map = NULL,
                                        sim_method = "resnik",
                                        enrich_fdr = 0.05,
                                        min_term_freq = NULL,
                                        max_enriched_terms = 100,
                                        n_top_sims = 5,
                                        ic_quantile = 0.5) {
    if (is.null(enrichment_genes)) {
        enrichment_genes <- genes
    }
    specific <- .specific_go_annotations(
        enrichment_genes, go_cache, ic_scores, ic_quantile
    )
    if (is.null(specific)) {
        return(rep(0, length(genes)))
    }
    target_terms <- .select_data_driven_terms(
        specific$annotated_genes, specific$go_cache,
        min_term_freq, enrich_fdr, max_enriched_terms
    )
    if (length(target_terms) == 0L) {
        return(rep(0, length(genes)))
    }
    .score_data_driven_genes(
        genes, target_terms, specific$go_cache, ic_scores,
        similarity_cache, ancestor_map, sim_method, n_top_sims
    )
}


#' Compute Semantic Similarity Between Two GO Terms
#'
#' Dispatcher that calls the appropriate metric function.
#' All metrics are based on the Most Informative Common Ancestor (MICA):
#'   the shared ancestor of the two terms that has the highest IC.
#'
#' @param term1 First GO term ID
#' @param term2 Second GO term ID
#' @param ic_scores Named vector of Information Content scores
#' @param ancestor_map Named list of GO term -> ancestor vectors
#' @param method One of "resnik", "lin", "jiang", "rel"
#' @return Numeric similarity score between 0 and 1
#' @keywords internal
compute_semantic_similarity <- function(term1, term2, ic_scores,
                                        ancestor_map = NULL,
                                        method = "resnik") {
    if (term1 == term2) {
        return(1.0)
    }

    ancestors1 <- get_go_ancestors(term1, ancestor_map)
    ancestors2 <- get_go_ancestors(term2, ancestor_map)

    common_ancestors <- intersect(ancestors1, ancestors2)
    if (length(common_ancestors) == 0) {
        return(0.0)
    }

    common_with_ic <- common_ancestors[common_ancestors %in% names(ic_scores)]
    if (length(common_with_ic) == 0) {
        return(0.0)
    }

    mica_ic <- max(ic_scores[common_with_ic], na.rm = TRUE)

    ic1 <- ic_scores[term1]
    ic2 <- ic_scores[term2]
    if (is.na(ic1) || is.na(ic2)) {
        return(0.0)
    }
    if (ic1 == 0 && ic2 == 0) {
        return(0.0)
    }

    sim <- switch(method,
        resnik = sim_resnik(mica_ic, ic1, ic2),
        lin    = sim_lin(mica_ic, ic1, ic2),
        jiang  = sim_jiang(mica_ic, ic1, ic2),
        rel    = sim_rel(mica_ic, ic1, ic2),
        stop(sprintf("Unknown similarity method: '%s'", method))
    )

    return(max(0, min(sim, 1.0)))
}


#' Resnik Similarity (max-normalized)
#'
#' sim = IC(MICA) / max(IC(t1), IC(t2))
#' Range: (0, 1)
#'
#' @keywords internal
sim_resnik <- function(mica_ic, ic1, ic2) {
    max_ic <- max(ic1, ic2)
    if (is.infinite(max_ic) || max_ic == 0) {
        return(0.0)
    }
    mica_ic / max_ic
}

#' Lin Similarity
#'
#' sim = 2 * IC(MICA) / (IC(t1) + IC(t2))
#' Range: (0, 1)
#'
#' @keywords internal
sim_lin <- function(mica_ic, ic1, ic2) {
    denom <- ic1 + ic2
    if (is.infinite(denom) || denom == 0) {
        return(0.0)
    }
    2 * mica_ic / denom
}

#' Jiang-Conrath Similarity
#'
#' distance = IC(t1) + IC(t2) - 2 * IC(MICA)
#' sim = 1 / (1 + distance)
#' Range: (0, 1)
#'
#' @keywords internal
sim_jiang <- function(mica_ic, ic1, ic2) {
    distance <- ic1 + ic2 - 2 * mica_ic
    distance <- max(0, distance)
    1 / (1 + distance)
}

#' Relevance Similarity (Schlicker et al.)
#'
#' sim = Lin(t1, t2) * (1 - p(MICA))
#' where p(MICA) = exp(-IC(MICA)) is the annotation probability.
#' Penalizes broad common ancestors.
#' Range: (0, 1).
#'
#' @keywords internal
sim_rel <- function(mica_ic, ic1, ic2) {
    lin_val <- sim_lin(mica_ic, ic1, ic2)
    p_mica <- exp(-mica_ic)
    lin_val * (1 - p_mica)
}


#' Get GO Ancestors
#'
#' Retrieves all ancestor terms of a GO term by traversing the GO DAG.
#' Uses a pre-built ancestor map if available, otherwise falls back to GO.db.
#'
#' @param term GO term ID (e.g., "GO:0006955")
#' @param ancestor_map Named list mapping GO terms to their ancestor vectors.
#' @return Character vector of ancestor terms (always includes the term itself)
#' @keywords internal
get_go_ancestors <- function(term, ancestor_map = NULL) {
    if (!is.null(ancestor_map)) {
        return(unique(c(term, ancestor_map[[term]])))
    }
    if (!requireNamespace("GO.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
        return(term)
    }
    ontology <- tryCatch(
        AnnotationDbi::select(
            GO.db::GO.db, keys = term,
            columns = "ONTOLOGY", keytype = "GOID"
        )$ONTOLOGY[[1L]],
        error = function(error) NA_character_
    )
    if (is.na(ontology)) {
        return(term)
    }
    ancestors <- tryCatch(
        AnnotationDbi::as.list(get(
            paste0("GO", ontology, "ANCESTOR"),
            envir = asNamespace("GO.db")
        ))[[term]],
        error = function(error) NULL
    )
    unique(c(term, ancestors[ancestors != "all"]))
}


#' Build Ancestor Map from GO.db
#'
#' Pre-computes ancestor relationships for all GO terms to avoid repeated
#' lookups. Results are cached to disk.
#'
#' @param organism Character, organism name
#' @param use_cache Logical, use disk caching
#' @return Named list mapping each GO term to its ancestor terms
#' @keywords internal
load_ancestor_map <- function(organism = "human", use_cache = TRUE) {
    if (!requireNamespace("GO.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
        stop("GO semantic scoring requires GO.db and AnnotationDbi")
    }
    go_version <- as.character(utils::packageVersion("GO.db"))
    cache_key <- paste0("go_ancestor_map_v2_", organism, "_", go_version)

    if (use_cache && exists(cache_key, envir = .geneselectr2_cache)) {
        return(get(cache_key, envir = .geneselectr2_cache))
    }

    if (use_cache) {
        cache_file <- file.path(get_cache_dir(), paste0(cache_key, ".rds"))
        ancestor_map <- .read_cache_file(cache_file)
        if (!is.null(ancestor_map)) {
            assign(cache_key, ancestor_map, envir = .geneselectr2_cache)
            return(ancestor_map)
        }
    }

    if (getOption("geneselectr2.verbose", TRUE)) {
        message("Building the GO ancestor map from GO.db...")
    }
    ancestor_map <- list()
    for (ontology in c("BP", "MF", "CC")) {
        object_name <- paste0("GO", ontology, "ANCESTOR")
        ontology_map <- tryCatch(
            AnnotationDbi::as.list(get(
                object_name, envir = asNamespace("GO.db")
            )),
            error = function(error) {
                stop(
                    "Could not load ", object_name, ": ",
                    conditionMessage(error), call. = FALSE
                )
            }
        )
        ontology_map <- lapply(ontology_map, setdiff, "all")
        ancestor_map[names(ontology_map)] <- ontology_map
    }
    if (getOption("geneselectr2.verbose", TRUE)) {
        message(sprintf(
            "  Loaded ancestors for %d GO terms",
            length(ancestor_map)
        ))
    }

    if (length(ancestor_map) == 0L) {
        stop("GO.db returned an empty ancestor map")
    }

    if (use_cache && length(ancestor_map) > 0) {
        cache_file <- file.path(get_cache_dir(), paste0(cache_key, ".rds"))
        .write_cache_file(ancestor_map, cache_file)
        assign(cache_key, ancestor_map, envir = .geneselectr2_cache)
    }
    ancestor_map
}


#' Compute Information Content
#'
#' IC(term) = -log( (count(term) + 1) / (n_genes + 1) )
#' Laplace-smoothed to avoid log(0).
#'
#' @param go_cache GO annotations (gene -> terms), already ontology-filtered
#' @return Named numeric vector of IC scores
#' @keywords internal
compute_information_content <- function(go_cache) {
    all_terms <- unique(unlist(go_cache))

    if (length(all_terms) == 0) {
        return(numeric(0))
    }

    term_counts <- table(unlist(go_cache))
    n_genes <- length(go_cache)

    ic_scores <- -log((term_counts + 1) / (n_genes + 1))

    return(ic_scores)
}


#' Test GO Enrichment via Fisher's Exact Test
#'
#' For each GO term annotated to any selected gene, tests whether the term
#' is over-represented in the selected set vs the background.
#'
#' @param selected_genes Character vector of selected gene names
#' @param background_genes Character vector of background gene names
#' @param go_cache GO annotations
#' @return Data frame with columns: term, p_value, odds_ratio, n_selected,
#'   n_background, p_adj (BH-adjusted)
#' @keywords internal
test_go_enrichment <- function(selected_genes, background_genes, go_cache) {
    selected_terms <- unique(unlist(go_cache[selected_genes]))

    if (length(selected_terms) == 0) {
        return(data.frame(
            term = character(0),
            p_value = numeric(0),
            odds_ratio = numeric(0),
            n_selected = integer(0),
            n_background = integer(0),
            p_adj = numeric(0),
            stringsAsFactors = FALSE
        ))
    }

    n_sel <- length(selected_genes)
    n_bg <- length(background_genes)

    results <- lapply(selected_terms, function(term) {
        n_sel_with <- sum(vapply(
            go_cache[selected_genes],
            function(x) term %in% x, logical(1)
        ))
        n_bg_with <- sum(vapply(
            go_cache[background_genes],
            function(x) term %in% x, logical(1)
        ))

        contingency <- matrix(c(
            n_sel_with, n_sel - n_sel_with,
            n_bg_with,  n_bg - n_bg_with
        ), nrow = 2)

        test_result <- stats::fisher.test(
            contingency, alternative = "greater"
        )

        data.frame(
            term = term,
            p_value = test_result$p.value,
            odds_ratio = as.numeric(test_result$estimate),
            n_selected = n_sel_with,
            n_background = n_bg_with,
            stringsAsFactors = FALSE
        )
    })

    results_df <- do.call(rbind, results)
    results_df$p_adj <- stats::p.adjust(results_df$p_value, method = "BH")
    results_df <- results_df[order(results_df$p_value), ]

    return(results_df)
}


#' Multi-Layer Biological Scorer
#'
#' Combines Gene Ontology semantic similarity and MSigDB gene-set membership.
#' Each active source contributes equally through a geometric mean.
#'
#' @param genes Character vector of gene symbols to score
#' @param layers Character vector, subset of \code{c("go", "msigdb")}.
#'   Controls which evidence layers are active.
#' @param disease_term Character string describing the disease (e.g.,
#'   "breast cancer", "ulcerative colitis"). Used by the MSigDB source.
#' @param target_terms Character vector of GO term IDs (for the GO layer)
#' @param go_mode Character, "supervised" or "data_driven" (for the GO layer)
#' @param msigdb_categories Character vector of MSigDB collections to search.
#'   Default: c("C2", "C7", "H") = curated, immunologic and Hallmark gene
#'   sets.
#' @param organism Character, species name for msigdbr. Default: "Homo sapiens".
#' @param use_cache Logical, cache results to disk.
#' @param verbose Logical, print progress.
#' @param ... Additional arguments passed to the GO layer
#'   (\code{biological_scorer}).
#'
#' @return Numeric vector of combined biological scores (one per gene),
#'   percentile-normalized between 0 and 1.
#' @examples
#' if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'     multilayer_bio_scorer(
#'         c("TP53", "BRCA1"),
#'         layers = "go",
#'         target_terms = "GO:0006915", verbose = FALSE
#'     )
#' }
#'
#' @export
multilayer_bio_scorer <- function(
    genes,
    layers = c("go", "msigdb"),
    disease_term = NULL,
    target_terms = NULL,
    go_mode = "supervised",
    msigdb_categories = c("C2", "C7", "H"),
    organism = "Homo sapiens",
    use_cache = TRUE,
    verbose = TRUE,
    ...
) {
    layers <- match.arg(layers, c("go", "msigdb"),
        several.ok = TRUE
    )
    n_genes <- length(genes)

    if ("msigdb" %in% layers && is.null(disease_term)) {
        stop("disease_term is required for the MSigDB source")
    }

    if ("go" %in% layers && go_mode == "supervised" &&
        (is.null(target_terms) || length(target_terms) == 0)) {
        stop("target_terms required for GO layer in supervised mode")
    }

    layer_scores <- list()

    if ("go" %in% layers) {
        if (verbose) {
            message("  GO layer: computing semantic similarity...\n")
        }

        go_raw <- biological_scorer(
            genes = genes,
            mode = go_mode,
            target_terms = target_terms,
            use_cache = use_cache,
            ...
        )
        layer_scores[["go"]] <- go_raw

        if (verbose) {
            message(sprintf(
                "    %d/%d genes scored (mean = %.3f)\n",
                sum(go_raw > 0), n_genes, mean(go_raw, na.rm = TRUE)
            ))
        }
    }

    if ("msigdb" %in% layers) {
        if (verbose) {
            message(sprintf(
                "  Bio layer [MSigDB]: searching '%s' in [%s]...\n",
                disease_term, paste(msigdb_categories, collapse = ", ")
            ))
        }

        msigdb_raw <- score_msigdb_layer(
            genes = genes,
            disease_keyword = disease_term,
            categories = msigdb_categories,
            organism = organism,
            verbose = verbose
        )
        layer_scores[["msigdb"]] <- percentile01(msigdb_raw)

        if (verbose) {
            message(sprintf(
                "    %d/%d genes found in disease-relevant sets\n",
                sum(msigdb_raw > 0), n_genes
            ))
        }
    }

    if (length(layer_scores) == 0) {
        stop("No requested biological layer produced a score vector")
    }

    if (length(layer_scores) == 1) {
        combined <- layer_scores[[1]]
    } else {
        eps <- 1e-10
        n_layers <- length(layer_scores)
        log_sum <- numeric(n_genes)
        for (layer_name in names(layer_scores)) {
            log_sum <- log_sum + log(layer_scores[[layer_name]] + eps)
        }
        combined <- exp(log_sum / n_layers)

        ## The offset keeps log() finite. Genes without evidence remain zero.
        layer_mat <- do.call(cbind, layer_scores)
        no_evidence <- rowSums(layer_mat > 0) == 0
        combined[no_evidence] <- 0

        if (verbose && any(no_evidence)) {
            message(sprintf(
                "    %d/%d genes have no evidence in any layer (b = 0)\n",
                sum(no_evidence), n_genes
            ))
        }
    }

    combined <- percentile01(combined)

    if (verbose) {
        message(sprintf(
            "  Bio combined [%s]: mean = %.3f\n",
            paste(names(layer_scores), collapse = " + "),
            mean(combined, na.rm = TRUE)
        ))
    }

    attr(combined, "layer_scores") <- layer_scores
    attr(combined, "layers_used") <- names(layer_scores)

    return(combined)
}


#' Score Genes by MSigDB Gene Set Membership
#'
#' For a given disease keyword, finds all matching gene sets in MSigDB's
#' curated collections and scores each gene by the number of matching sets that
#' contain it.
#'
#' @param genes Character vector of gene symbols
#' @param disease_keyword Character, disease name to search in gene set names
#'   (e.g., "breast cancer", "alzheimer", "ulcerative colitis").
#'   Matching is case-insensitive and uses partial matching on gene set names.
#' @param categories Character vector of MSigDB collection codes.
#'   Default: c("C2", "C7", "H"). C2 = curated gene sets (includes CGP:
#'   chemical and genetic perturbations). C7 = immunologic gene sets. H =
#'   Hallmark gene sets (50 biological states).
#' @param organism Character, species name for msigdbr. Default: "Homo sapiens".
#' @param verbose Logical, print progress.
#' @return Numeric vector (one per gene). The number of matching gene sets is
#'   transformed as `log2(count + 1)`.
#'
#' @keywords internal
score_msigdb_layer <- function(genes,
                                disease_keyword,
                                categories = c("C2", "C7", "H"),
                                organism = "Homo sapiens",
                                verbose = TRUE) {
    n_genes <- length(genes)

    if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop(
            "msigdbr is required for the MSigDB layer. ",
            "Install with: install.packages('msigdbr')"
        )
    }

    all_sets <- do.call(rbind, lapply(categories, function(category) {
        tryCatch(
            {
                ## msigdbr 10.0.0 renamed `category` to `collection`.
                msigdb_args <- list(species = organism)
                if ("collection" %in% names(formals(msigdbr::msigdbr))) {
                    msigdb_args$collection <- category
                } else {
                    msigdb_args$category <- category
                }
                do.call(msigdbr::msigdbr, msigdb_args)
            },
            error = function(e) {
                stop(
                    sprintf(
                        "MSigDB category '%s' failed: %s",
                        category, conditionMessage(e)
                    ),
                    call. = FALSE
                )
            }
        )
    }))

    if (is.null(all_sets) || nrow(all_sets) == 0) {
        stop(
            "No MSigDB gene sets were retrieved for these categories",
            call. = FALSE
        )
    }

    keyword_pattern <- gsub("\\s+", ".*", tolower(disease_keyword))
    name_match <- grepl(keyword_pattern, tolower(all_sets$gs_name))

    if ("gs_description" %in% colnames(all_sets)) {
        desc_match <- grepl(keyword_pattern, tolower(all_sets$gs_description))
        disease_rows <- name_match | desc_match
    } else {
        disease_rows <- name_match
    }

    disease_sets <- all_sets[disease_rows, ]
    n_matching_sets <- length(unique(disease_sets$gs_name))

    if (n_matching_sets == 0) {
        warning(sprintf(
            "No MSigDB gene sets match keyword '%s' in categories [%s]. ",
            disease_keyword, paste(categories, collapse = ", ")
        ))
        return(rep(0, n_genes))
    }

    if (verbose) {
        message(sprintf(
            "    Found %d gene sets matching '%s'\n",
            n_matching_sets, disease_keyword
        ))
    }

    gene_set_counts <- table(disease_sets$gene_symbol)

    scores <- numeric(n_genes)
    names(scores) <- genes
    matched <- intersect(genes, names(gene_set_counts))
    scores[matched] <- as.numeric(gene_set_counts[matched])

    scores <- log2(scores + 1)

    return(scores)
}

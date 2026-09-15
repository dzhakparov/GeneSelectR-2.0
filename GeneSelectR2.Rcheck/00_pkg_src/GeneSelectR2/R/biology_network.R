# =============================================================================
# GeneSelectR 2.0 - Network-Propagation Biology Scorer
# =============================================================================
#
# Open Targets supplies disease-associated genes. STRING random walk with
# restart assigns network-proximity scores to candidate genes. Retrieved seed
# sets are cached with their query settings.


# -----------------------------------------------------------------------------
#  Open Targets seed derivation
# -----------------------------------------------------------------------------

#' Resolve a Disease Name to an EFO ID via Open Targets Search
#'
#' @param disease_term Character, free-text disease name, such as
#'   "atopic dermatitis".
#' @param verbose Logical, print resolution progress
#' @return Character EFO ID (e.g. "EFO_0000274"), or NULL if not resolved
#' @keywords internal
resolve_efo_id <- function(disease_term, verbose = FALSE) {
    if (!requireNamespace("httr", quietly = TRUE) ||
        !requireNamespace("jsonlite", quietly = TRUE)) {
        dependency_message <- paste(
            "Network biology requires 'httr' and 'jsonlite'.",
            "Install both packages before running this function."
        )
        stop(
            dependency_message,
            call. = FALSE
        )
    }

    # If the term already looks like an EFO/MONDO/ontology ID, pass it through.
    if (grepl("^(EFO|MONDO|HP|Orphanet|DOID)[:_]", disease_term)) {
        if (verbose) {
            message("    The term is already an ontology ID; using it as-is.\n")
        }
        return(gsub(":", "_", disease_term))
    }

    endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"

    query <- '
    query resolve($q: String!) {
        search(
        queryString: $q,
        entityNames: ["disease"],
        page: {index: 0, size: 1}
        ) {
        hits { id name entity }
        }
    }'

    response <- tryCatch(
        httr::POST(endpoint,
            body = list(
                query = query,
                variables = list(q = disease_term)
            ),
            encode = "json"
        ),
        error = function(e) NULL
    )
    if (is.null(response)) stop("Open Targets disease search request failed")
    if (httr::status_code(response) != 200) {
        stop(
            sprintf(
                "Open Targets disease search returned HTTP %d",
                httr::status_code(response)
            ),
            call. = FALSE
        )
    }

    parsed <- tryCatch(
        jsonlite::fromJSON(httr::content(response,
            as = "text",
            encoding = "UTF-8"
        )),
        error = function(e) NULL
    )
    if (is.null(parsed)) {
        stop("Open Targets disease search returned invalid JSON")
    }
    if (!is.null(parsed$errors)) {
        stop("Open Targets rejected the disease search query", call. = FALSE)
    }
    hits <- parsed$data$search$hits
    if (is.null(hits) || length(hits) == 0 || nrow(hits) == 0) {
        stop(
            sprintf(
                "Open Targets found no disease matching '%s'",
                disease_term
            ),
            call. = FALSE
        )
    }

    if (verbose) {
        message(sprintf(
            "    Search matched '%s' to %s.\n",
            hits$name[1], hits$id[1]
        ))
    }
    # Attach the resolved disease name so callers can record (not just print)
    # what the query actually resolved to. A top-1 search hit can resolve to a
    # neighbouring disease, and silent mis-resolution is a wrong-prior fault.
    id <- hits$id[1]
    attr(id, "resolved_name") <- hits$name[1]
    id
}


#' Derive Disease Seed Genes from Open Targets
#'
#' Queries the Open Targets Platform for the top targets associated with a
#' disease, ranked by overall association score. Results are frozen to the
#' package cache keyed by (EFO ID, max_seeds, min_score) so the API is called
#' at most once per parameter combination - subsequent runs load from disk and
#' are fully reproducible offline.
#'
#' @param disease_term Disease name or EFO ID
#' @param max_seeds Integer, take at most this many top-scoring targets
#'   (default: 100)
#' @param min_score Numeric, drop associations below this overall score
#'   (default: 0.1). Applied after the top-N cut, so a disease with weak
#'   associations yields fewer seeds rather than padding with noise.
#' @param use_cache Logical, read/write the frozen seed file (default: TRUE)
#' @param force_refresh Logical, ignore any cached file and submit a new
#'   request.
#' @param verbose Logical
#' @return Data frame with columns: ensembl_id, symbol, score. An empty data
#'   frame indicates that the resolved disease has no qualifying associations.
#' @examples
#' if (interactive()) {
#'     get_disease_seeds_opentargets(
#'         "EFO_0000274",
#'         max_seeds = 5, verbose = FALSE
#'     )
#' }
#' @export
get_disease_seeds_opentargets <- function(disease_term,
                                            max_seeds = 100,
                                            min_score = 0.1,
                                            use_cache = TRUE,
                                            force_refresh = FALSE,
                                            verbose = TRUE) {
    if (!requireNamespace("httr", quietly = TRUE) ||
        !requireNamespace("jsonlite", quietly = TRUE)) {
        dependency_message <- paste(
            "Network biology requires 'httr' and 'jsonlite'.",
            "Install both packages before running this function."
        )
        stop(
            dependency_message,
            call. = FALSE
        )
    }

    # A frozen cache records the submitted disease term and resolved
    # ontology identifier. Search those files before an online term-resolution
    # request. This makes a previously fetched seed set reproducible offline.
    if (use_cache && !force_refresh &&
        !grepl("^(EFO|MONDO|HP|Orphanet|DOID)[:_]", disease_term)) {
        seed_pattern <- sprintf(
            "^ot_seeds_.*_n%d_s%s\\.rds$",
            max_seeds, format(min_score, trim = TRUE)
        )
        cached_files <- list.files(get_cache_dir(),
            pattern = seed_pattern,
            full.names = TRUE
        )
        for (candidate_file in cached_files) {
            cached_seeds <- tryCatch(readRDS(candidate_file),
                error = function(e) NULL
            )
            if (!is.null(cached_seeds) &&
                identical(attr(cached_seeds, "disease"), disease_term) &&
                all(c("ensembl_id", "symbol", "score") %in%
                    colnames(cached_seeds))) {
                if (verbose) {
                    cache_format <- paste(
                        "  [Open Targets] loaded %d frozen seeds for",
                        "'%s' from cache (%s)\n"
                    )
                    message(sprintf(
                        cache_format,
                        nrow(cached_seeds), disease_term,
                        basename(candidate_file)
                    ))
                }
                return(cached_seeds)
            }
        }
    }

    # --- Resolve disease to EFO ID ---
    if (verbose) {
        message(sprintf(
            "  [Open Targets] resolving disease term: '%s'\n",
            disease_term
        ))
    }
    efo_id <- resolve_efo_id(disease_term, verbose = verbose)
    if (verbose) {
        message(sprintf("  [Open Targets] resolved to EFO ID: %s\n", efo_id))
    }
    resolved_name <- attr(efo_id, "resolved_name")

    # --- Cache key: frozen seed file ---
    cache_file <- file.path(
        get_cache_dir(),
        sprintf(
            "ot_seeds_%s_n%d_s%g.rds",
            gsub("[^A-Za-z0-9]", "", efo_id), max_seeds, min_score
        )
    )

    if (use_cache && !force_refresh && file.exists(cache_file)) {
        cached_seeds <- readRDS(cache_file)
        if (verbose) {
            message(sprintf(
                "  [Open Targets] loaded %d frozen seeds from cache (%s)\n",
                nrow(cached_seeds), basename(cache_file)
            ))
        }
        return(cached_seeds)
    }

    # --- Query Open Targets for associated targets ---
    if (verbose) {
        message(sprintf(
            "  [Open Targets] querying associated targets (top %d)...\n",
            max_seeds
        ))
    }

    endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
    query <- "
    query assoc($efoId: String!, $size: Int!) {
        disease(efoId: $efoId) {
        id
        name
        associatedTargets(page: {index: 0, size: $size}) {
            count
            rows {
            target { id approvedSymbol }
            score
            }
        }
        }
    }"

    response <- tryCatch(
        httr::POST(endpoint,
            body = list(
                query = query,
                variables = list(
                    efoId = efo_id,
                    size = max_seeds
                )
            ),
            encode = "json"
        ),
        error = function(e) NULL
    )
    if (is.null(response)) {
        stop("Open Targets association request failed")
    }

    status <- httr::status_code(response)
    if (verbose) message(sprintf("  [Open Targets] HTTP status: %d\n", status))
    if (status != 200) {
        stop(
            sprintf(
                "Open Targets association request returned HTTP %d",
                status
            ),
            call. = FALSE
        )
    }

    parsed <- tryCatch(
        jsonlite::fromJSON(httr::content(response,
            as = "text",
            encoding = "UTF-8"
        )),
        error = function(e) NULL
    )
    if (is.null(parsed)) {
        stop("Open Targets association response contained invalid JSON")
    }

    # Report GraphQL errors explicitly - the API returns HTTP 200 with an
    # `errors` field when the query itself is malformed or the EFO is unknown.
    if (!is.null(parsed$errors)) {
        error_messages <- tryCatch(
            paste(parsed$errors$message, collapse = "; "),
            error = function(e) "unparseable error object"
        )
        stop(
            sprintf(
                "Open Targets rejected the association query: %s",
                error_messages
            ),
            call. = FALSE
        )
    }

    disease_name <- tryCatch(parsed$data$disease$name, error = function(e) NA)
    total_assoc <- tryCatch(parsed$data$disease$associatedTargets$count,
        error = function(e) NA
    )
    if (verbose && !is.na(disease_name)) {
        message(sprintf(
            "  [Open Targets] matched disease: '%s' (%s total associations)\n",
            disease_name, format(total_assoc)
        ))
    }

    rows <- parsed$data$disease$associatedTargets$rows
    if (is.null(rows) || length(rows) == 0 ||
        is.null(rows$target) || nrow(rows$target) == 0) {
        warning(
            sprintf("No associated targets were returned for %s", efo_id),
            call. = FALSE
        )
        return(data.frame(
            ensembl_id = character(0), symbol = character(0),
            score = numeric(0), stringsAsFactors = FALSE
        ))
    }

    seeds <- data.frame(
        ensembl_id = rows$target$id,
        symbol = rows$target$approvedSymbol,
        score = rows$score,
        stringsAsFactors = FALSE
    )

    n_before_filter <- nrow(seeds)

    # Apply the soft minimum score (top-N already enforced by the query size).
    keep <- !is.na(seeds$score) & seeds$score >= min_score
    seeds <- seeds[keep, , drop = FALSE]
    seeds <- seeds[order(seeds$score, decreasing = TRUE), , drop = FALSE]
    rownames(seeds) <- NULL

    if (verbose) {
        message(sprintf(
            "  [Open Targets] %d targets returned, %d pass score >= %g\n",
            n_before_filter, nrow(seeds), min_score
        ))
        if (nrow(seeds) > 0) {
            top_preview <- utils::head(seeds, 5)
            message(sprintf(
                "  [Open Targets] top seeds: %s\n",
                paste(sprintf(
                    "%s(%.2f)", top_preview$symbol,
                    top_preview$score
                ), collapse = ", ")
            ))
        }
    }

    # --- Freeze to cache ---
    if (use_cache) {
        attr(seeds, "efo_id") <- efo_id
        attr(seeds, "disease") <- disease_term
        attr(seeds, "resolved_name") <- resolved_name
        attr(seeds, "fetch_date") <- Sys.Date()
        saveRDS(seeds, cache_file)
    }

    seeds
}


# -----------------------------------------------------------------------------
#  Network propagation (random walk with restart)
# -----------------------------------------------------------------------------

#' Score Genes by Network Proximity to Disease Seeds
#'
#' Places the disease seed genes on the STRING functional interaction network
#' and propagates relevance via random walk with restart (RWR). Each candidate
#' gene receives a score equal to its stationary RWR probability - high for
#' genes that are seeds themselves or sit close to many seeds in the network,
#' low for genes far from disease biology.
#'
#' @param genes Character vector of candidate gene symbols to score
#' @param disease_term Disease name or EFO ID (passed to Open Targets)
#' @param string_score_threshold Integer, minimum STRING combined score for an
#'   edge to be included (default: 400, STRING's "medium confidence")
#' @param restart_prob Numeric, RWR restart probability (default: 0.5). Higher
#'   retains more probability near the seeds; lower values allow wider
#'   propagation.
#' @param max_seeds,min_score Passed to \code{get_disease_seeds_opentargets}
#' @param organism STRING species id (default: 9606, human)
#' @param string_version STRING version (default: "12.0")
#' @param use_cache Logical
#' @param verbose Logical
#' @return Numeric vector (length = length(genes)) of percentile-normalised
#'   network relevance scores between 0 and 1, named by gene.
#' @examples
#' if (interactive()) {
#'     score_network_layer(
#'         c("IL6", "STAT3"), "EFO_0000274",
#'         verbose = FALSE
#'     )
#' }
#' @export
score_network_layer <- function(genes,
                                disease_term,
                                string_score_threshold = 400,
                                restart_prob = 0.5,
                                max_seeds = 100,
                                min_score = 0.1,
                                organism = 9606,
                                string_version = "12.0",
                                use_cache = TRUE,
                                verbose = TRUE) {
    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
        stop(
            "Network biology requires 'STRINGdb'. ",
            "Install with: BiocManager::install('STRINGdb')"
        )
    }
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(
            "Network biology requires 'igraph'. ",
            "Install with: install.packages('igraph')"
        )
    }

    n_genes <- length(genes)
    zero_scored <- stats::setNames(rep(0, n_genes), genes)

    # --- 1. Disease seeds from Open Targets ---
    seeds <- get_disease_seeds_opentargets(
        disease_term,
        max_seeds = max_seeds, min_score = min_score,
        use_cache = use_cache, verbose = verbose
    )
    if (nrow(seeds) == 0) {
        warning("No disease seeds available; network layer returns zeros.")
        return(zero_scored)
    }

    # --- 2. Build STRING network over the union of seeds + candidate genes ---
    if (verbose) message("  Building STRING network...\n")

    # STRING reference files are large. Store them in a persistent directory
    # and allow enough time for the initial download.
    string_cache_dir <- getOption(
        "GeneSelectR.string_cache",
        file.path(path.expand("~"), ".cache", "GeneSelectR", "stringdb")
    )
    dir.create(string_cache_dir, recursive = TRUE, showWarnings = FALSE)

    previous_timeout <- getOption("timeout")
    if (is.numeric(previous_timeout) && previous_timeout < 3600) {
        options(timeout = 3600)
        on.exit(options(timeout = previous_timeout), add = TRUE)
    }

    string_db <- tryCatch(
        STRINGdb::STRINGdb$new(
            version = string_version,
            species = organism,
            score_threshold = string_score_threshold,
            input_directory = string_cache_dir
        ),
        error = function(e) NULL
    )
    all_symbols <- unique(c(genes, seeds$symbol))

    # STRINGdb validates its version through an API request even when every
    # reference file is present locally. A frozen, score-filtered edge cache and
    # the official protein-info table provide an equivalent offline path when
    # that request is unavailable.
    using_offline_string <- is.null(string_db)
    if (using_offline_string) {
        offline_key <- paste0(
            "string_offline_", organism, "_v", string_version,
            "_s", string_score_threshold
        )
        offline_resources <- if (
            exists(offline_key, envir = .geneselectr_cache)
        ) {
            get(offline_key, envir = .geneselectr_cache)
        } else {
            info_file <- file.path(
                string_cache_dir,
                sprintf("%s.protein.info.v%s.txt.gz", organism, string_version)
            )
            edge_file <- file.path(
                string_cache_dir,
                sprintf(
                    "%s.protein.links.score%d.v%s.rds", organism,
                    string_score_threshold, string_version
                )
            )
            if (!file.exists(info_file) || !file.exists(edge_file)) {
                cache_message <- paste(
                    "STRING initialization failed and the offline cache",
                    "is incomplete. Required files: %s and %s"
                )
                stop(
                    sprintf(
                        cache_message,
                        basename(info_file), basename(edge_file)
                    ),
                    call. = FALSE
                )
            }
            protein_info <- utils::read.delim(
                gzfile(info_file),
                skip = 1L, header = FALSE, quote = "",
                stringsAsFactors = FALSE,
                col.names = c(
                    "STRING_id", "preferred_name", "protein_size",
                    "annotation"
                )
            )
            edge_table <- readRDS(edge_file)
            required_edge_columns <- c("protein1", "protein2", "combined_score")
            if (!all(required_edge_columns %in% colnames(edge_table)) ||
                any(!is.finite(edge_table$combined_score)) ||
                any(edge_table$combined_score < string_score_threshold)) {
                stop("The offline STRING edge cache failed validation.")
            }
            offline_graph <- igraph::graph_from_data_frame(
                edge_table[, required_edge_columns],
                directed = FALSE
            )
            resources <- list(
                symbol_to_id = stats::setNames(
                    protein_info$STRING_id,
                    protein_info$preferred_name
                ),
                graph = offline_graph
            )
            assign(offline_key, resources, envir = .geneselectr_cache)
            resources
        }
        mapped_ids <- unname(offline_resources$symbol_to_id[all_symbols])
        mapping <- data.frame(
            gene = all_symbols[!is.na(mapped_ids)],
            STRING_id = mapped_ids[!is.na(mapped_ids)],
            stringsAsFactors = FALSE
        )
        g <- offline_resources$graph
        if (verbose) {
            message("  STRING API unavailable; using frozen local v12 files\n")
        }
    } else {
        mapping <- tryCatch(
            string_db$map(
                data.frame(gene = all_symbols, stringsAsFactors = FALSE),
                "gene",
                removeUnmappedRows = TRUE
            ),
            error = function(e) NULL
        )
    }
    if (is.null(mapping) || nrow(mapping) < 5) {
        stop("Fewer than five candidate or seed genes mapped to STRING")
    }

    if (verbose) {
        n_seed_mapped <- sum(mapping$gene %in% seeds$symbol)
        message(sprintf(
            "  STRING mapping: %d/%d symbols mapped, %d seeds mapped\n",
            nrow(mapping), length(all_symbols), n_seed_mapped
        ))
    }

    # Load the complete STRING network. PageRank on an induced candidate/seed
    # subgraph excludes intermediate proteins and makes proximity depend on the
    # candidate pool. The full graph retains all paths allowed by the configured
    # STRING score threshold.
    pull_graph <- function() {
        tryCatch(string_db$get_graph(), error = function(e) NULL)
    }
    if (!using_offline_string) g <- pull_graph()

    # A truncated links file in the cache produces NA edges here. That is a
    # corrupt download, not an absence of interactions, so delete the fragment
    # and try once more before giving up.
    if (is.null(g) && !using_offline_string) {
        links <- list.files(string_cache_dir,
            pattern = "protein\\.links",
            full.names = TRUE
        )
        if (length(links) > 0) {
            warning(
                sprintf(
                    "Removing incomplete STRING files and retrying: %s",
                    paste(basename(links), collapse = ", ")
                ),
                call. = FALSE
            )
            unlink(links)
            string_db <- tryCatch(
                STRINGdb::STRINGdb$new(
                    version = string_version, species = organism,
                    score_threshold = string_score_threshold,
                    input_directory = string_cache_dir
                ),
                error = function(e) NULL
            )
            if (!is.null(string_db)) g <- pull_graph()
        }
    }
    if (is.null(g) || igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
        stop("The complete STRING graph could not be loaded")
    }

    # --- 3. Run random walk with restart on the complete graph ---
    g <- igraph::simplify(g, edge.attr.comb = "max")

    graph_nodes <- igraph::V(g)$name

    # Give each seed probability in proportion to its disease-association score.
    # to its Open Targets association score.
    #
    # Map symbols to STRING identifiers while keeping scores aligned.
    seed_score_by_symbol <- stats::setNames(seeds$score, seeds$symbol)

    # The mapping contains one row per symbol and STRING identifier.
    mapping_scores <- seed_score_by_symbol[mapping$gene]
    is_seed_row <- !is.na(mapping_scores)

    seed_scores_by_id <- stats::setNames(
        mapping_scores[is_seed_row],
        mapping$STRING_id[is_seed_row]
    )

    # Collapse duplicate STRING ids (several symbols can map to one id) by max.
    if (any(duplicated(names(seed_scores_by_id)))) {
        seed_scores_by_id <- tapply(
            seed_scores_by_id,
            names(seed_scores_by_id), max
        )
    }

    restart_vector <- stats::setNames(rep(0, length(graph_nodes)), graph_nodes)
    present_seeds <- intersect(names(seed_scores_by_id), graph_nodes)
    if (length(present_seeds) == 0) {
        warning("No seed genes present in the STRING graph; returning zeros.")
        return(zero_scored)
    }

    # Seed scores can contain NA (a mapped STRING id whose symbol didn't match
    # back to the seed table) or zeros. Either would make the normalised restart
    # vector contain NA/NaN, which igraph's PageRank rejects outright. Keep only
    # finite, positive seed weights.
    seed_weights <- seed_scores_by_id[present_seeds]
    valid <- is.finite(seed_weights) & seed_weights > 0
    present_seeds <- present_seeds[valid]
    seed_weights <- seed_weights[valid]

    if (length(present_seeds) == 0 || sum(seed_weights) <= 0) {
        seed_message <- paste(
            "No finite, positive seed weights map into the graph;",
            "returning zeros."
        )
        warning(
            seed_message,
            call. = FALSE
        )
        return(zero_scored)
    }

    restart_vector[present_seeds] <- seed_weights
    restart_vector <- restart_vector / sum(restart_vector)

    # Final guard: personalization must be all-finite and sum to 1.
    if (any(!is.finite(restart_vector))) {
        restart_message <- paste(
            "The restart vector contains non-finite values after",
            "normalisation; returning zeros."
        )
        warning(
            restart_message,
            call. = FALSE
        )
        return(zero_scored)
    }

    if (verbose) {
        graph_format <- paste(
            "  RWR on complete STRING graph: %d nodes, %d edges,",
            "%d seeds (restart=%.2f)\n"
        )
        message(sprintf(
            graph_format,
            length(graph_nodes), igraph::ecount(g),
            length(present_seeds), restart_prob
        ))
    }

    # Random walk with restart via igraph's personalized PageRank, which is
    # exactly RWR with the personalization vector as the restart distribution.
    edge_attributes <- igraph::edge_attr_names(g)
    weight_name <- intersect(
        c("combined_score", "score", "weight"),
        edge_attributes
    )[1]
    graph_weights <- if (is.na(weight_name)) {
        NULL
    } else {
        igraph::edge_attr(g, weight_name)
    }
    if (!is.null(graph_weights) && any(!is.finite(graph_weights))) {
        stop("STRING graph contains non-finite edge weights")
    }

    rwr <- igraph::page_rank(
        g,
        damping = 1 - restart_prob,
        personalized = restart_vector[graph_nodes],
        weights = graph_weights
    )$vector

    # --- 4. Map RWR scores back to candidate gene symbols ---
    # Each candidate symbol -> its STRING id -> its RWR score.
    symbol_to_id <- stats::setNames(mapping$STRING_id, mapping$gene)
    candidate_ids <- symbol_to_id[genes]

    raw_scores <- rep(0, n_genes)
    in_graph <- !is.na(candidate_ids) & candidate_ids %in% names(rwr)
    raw_scores[in_graph] <- rwr[candidate_ids[in_graph]]
    names(raw_scores) <- genes

    if (verbose) {
        nonzero <- raw_scores[raw_scores > 0]
        score_format <- paste(
            "  [network] %d/%d candidate genes scored;",
            "RWR score range [%.2e, %.2e]\n"
        )
        message(sprintf(
            score_format,
            sum(in_graph), n_genes,
            if (length(nonzero)) min(nonzero) else 0,
            if (length(nonzero)) max(nonzero) else 0
        ))
        if (length(unique(round(nonzero, 10))) <= 1 && length(nonzero) > 1) {
            message(
                "  [network] All scores are identical; RWR may be degenerate\n"
            )
        }
    }

    # Percentile-normalise to [0,1] for combination with the MPO framework.
    percentile01(raw_scores)
}

## Open Targets supplies disease-associated genes. A random walk with restart
## on STRING assigns network-proximity scores to candidate genes.

.is_ontology_id <- function(disease_term) {
    grepl("^(EFO|MONDO|HP|Orphanet|DOID)[:_]", disease_term)
}

.require_packages <- function(packages, context) {
    missing_packages <- packages[!vapply(
        packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )]
    if (length(missing_packages) > 0L) {
        stop(
            sprintf(
                "%s requires: %s",
                context,
                paste(sprintf("'%s'", missing_packages), collapse = ", ")
            ),
            call. = FALSE
        )
    }
}

.open_targets_request <- function(query, variables, request_name) {
    endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
    response <- tryCatch(
        httr::POST(
            endpoint,
            body = list(query = query, variables = variables),
            encode = "json"
        ),
        error = identity
    )
    if (inherits(response, "error")) {
        stop(
            sprintf("%s failed: %s", request_name, conditionMessage(response)),
            call. = FALSE
        )
    }

    status <- httr::status_code(response)
    if (status != 200L) {
        stop(
            sprintf("%s returned HTTP %d", request_name, status),
            call. = FALSE
        )
    }

    parsed <- tryCatch(
        jsonlite::fromJSON(
            httr::content(response, as = "text", encoding = "UTF-8")
        ),
        error = identity
    )
    if (inherits(parsed, "error")) {
        stop(
            sprintf("%s returned invalid JSON", request_name),
            call. = FALSE
        )
    }
    if (!is.null(parsed$errors)) {
        messages <- tryCatch(
            paste(parsed$errors$message, collapse = "; "),
            error = function(e) "unparseable GraphQL error"
        )
        stop(
            sprintf("%s was rejected: %s", request_name, messages),
            call. = FALSE
        )
    }
    parsed
}

#' Resolve a Disease Name to an EFO ID via Open Targets Search
#'
#' @param disease_term Character, free-text disease name, such as
#'   "atopic dermatitis".
#' @param verbose Logical, print resolution progress
#' @return Character EFO ID (e.g. "EFO_0000274"), or NULL if not resolved
#' @keywords internal
resolve_efo_id <- function(disease_term, verbose = FALSE) {
    if (.is_ontology_id(disease_term)) {
        if (verbose) {
            message("    The submitted term is an ontology identifier.\n")
        }
        return(gsub(":", "_", disease_term))
    }
    .require_packages(c("httr", "jsonlite"), "Open Targets disease search")

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
    parsed <- .open_targets_request(
        query,
        list(q = disease_term),
        "Open Targets disease search"
    )
    hits <- parsed$data$search$hits
    if (is.null(hits) || length(hits) == 0L || nrow(hits) == 0L) {
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

    efo_id <- hits$id[1]
    attr(efo_id, "resolved_name") <- hits$name[1]
    efo_id
}

.empty_seed_table <- function() {
    data.frame(
        ensembl_id = character(),
        symbol = character(),
        score = numeric(),
        stringsAsFactors = FALSE
    )
}

.seed_cache_file <- function(efo_id, max_seeds, min_score) {
    file.path(
        get_cache_dir(),
        sprintf(
            "ot_seeds_%s_n%d_s%g.rds",
            gsub("[^A-Za-z0-9]", "", efo_id),
            max_seeds,
            min_score
        )
    )
}

.find_seed_cache <- function(disease_term, max_seeds, min_score, verbose) {
    seed_pattern <- sprintf(
        "^ot_seeds_.*_n%d_s%s\\.rds$",
        max_seeds,
        format(min_score, trim = TRUE)
    )
    cached_files <- list.files(
        get_cache_dir(),
        pattern = seed_pattern,
        full.names = TRUE
    )
    for (candidate_file in cached_files) {
        seeds <- .read_cache_file(candidate_file)
        valid_columns <- !is.null(seeds) &&
            all(c("ensembl_id", "symbol", "score") %in% colnames(seeds))
        if (valid_columns &&
            identical(attr(seeds, "disease"), disease_term)) {
            if (verbose) {
                message(sprintf(
                    "  [Open Targets] loaded %d seeds for '%s' from %s\n",
                    nrow(seeds),
                    disease_term,
                    basename(candidate_file)
                ))
            }
            return(seeds)
        }
    }
    NULL
}

.query_open_targets_seeds <- function(efo_id, max_seeds, verbose) {
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
    parsed <- .open_targets_request(
        query,
        list(efoId = efo_id, size = max_seeds),
        "Open Targets association request"
    )
    disease <- parsed$data$disease
    if (is.null(disease)) {
        return(.empty_seed_table())
    }
    if (verbose) {
        message(sprintf(
            "  [Open Targets] matched '%s' (%s total associations)\n",
            disease$name,
            format(disease$associatedTargets$count)
        ))
    }

    rows <- disease$associatedTargets$rows
    if (is.null(rows) || length(rows) == 0L ||
        is.null(rows$target) || nrow(rows$target) == 0L) {
        return(.empty_seed_table())
    }
    seeds <- data.frame(
        ensembl_id = rows$target$id,
        symbol = rows$target$approvedSymbol,
        score = rows$score,
        stringsAsFactors = FALSE
    )
    attr(seeds, "resolved_name") <- disease$name
    seeds
}

.filter_and_report_seeds <- function(seeds, min_score, verbose) {
    n_returned <- nrow(seeds)
    keep <- !is.na(seeds$score) & seeds$score >= min_score
    seeds <- seeds[keep, , drop = FALSE]
    seeds <- seeds[order(seeds$score, decreasing = TRUE), , drop = FALSE]
    rownames(seeds) <- NULL

    if (verbose) {
        message(sprintf(
            "  [Open Targets] %d targets returned, %d pass score >= %g\n",
            n_returned,
            nrow(seeds),
            min_score
        ))
        if (nrow(seeds) > 0L) {
            preview <- utils::head(seeds, 5L)
            labels <- sprintf("%s(%.2f)", preview$symbol, preview$score)
            message(sprintf(
                "  [Open Targets] top seeds: %s\n",
                paste(labels, collapse = ", ")
            ))
        }
    }
    seeds
}

.validate_seed_request <- function(disease_term, max_seeds, min_score) {
    if (!is.character(disease_term) || length(disease_term) != 1L ||
        is.na(disease_term) || !nzchar(disease_term)) {
        stop("disease_term must be one non-empty character value")
    }
    if (!is.numeric(max_seeds) || length(max_seeds) != 1L ||
        !is.finite(max_seeds) || max_seeds < 1 ||
        max_seeds != as.integer(max_seeds)) {
        stop("max_seeds must be one positive integer")
    }
    if (!is.numeric(min_score) || length(min_score) != 1L ||
        !is.finite(min_score) || min_score < 0 || min_score > 1) {
        stop("min_score must be one finite number in [0, 1]")
    }
    invisible(NULL)
}

#' Derive Disease Seed Genes from Open Targets
#'
#' Queries the Open Targets Platform for the top targets associated with a
#' disease, ranked by overall association score. Results are stored in the
#' package cache with the ontology identifier and query parameters. Cached
#' results support repeated analyses without a new request.
#'
#' @param disease_term Disease name or EFO ID
#' @param max_seeds Integer, take at most this many top-scoring targets
#'   (default: 100)
#' @param min_score Numeric, drop associations below this overall score
#'   (default: 0.1). Applied after the top-N query.
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
get_disease_seeds_opentargets <- function(
    disease_term,
    max_seeds = 100,
    min_score = 0.1,
    use_cache = TRUE,
    force_refresh = FALSE,
    verbose = TRUE
) {
    .require_packages(c("httr", "jsonlite"), "Open Targets seed retrieval")
    .validate_seed_request(disease_term, max_seeds, min_score)

    if (use_cache && !force_refresh && !.is_ontology_id(disease_term)) {
        cached <- .find_seed_cache(
            disease_term, max_seeds, min_score, verbose
        )
        if (!is.null(cached)) {
            return(cached)
        }
    }

    efo_id <- resolve_efo_id(disease_term, verbose = verbose)
    cache_file <- .seed_cache_file(efo_id, max_seeds, min_score)
    if (use_cache && !force_refresh) {
        cached <- .read_cache_file(cache_file)
        if (!is.null(cached)) {
            return(cached)
        }
    }

    seeds <- .query_open_targets_seeds(efo_id, max_seeds, verbose)
    query_name <- attr(seeds, "resolved_name")
    seeds <- .filter_and_report_seeds(seeds, min_score, verbose)
    attr(seeds, "efo_id") <- unname(efo_id)
    attr(seeds, "disease") <- disease_term
    resolved_name <- attr(efo_id, "resolved_name")
    if (is.null(resolved_name)) {
        resolved_name <- query_name
    }
    attr(seeds, "resolved_name") <- resolved_name
    attr(seeds, "fetch_date") <- Sys.Date()
    if (use_cache) {
        .write_cache_file(seeds, cache_file)
    }
    seeds
}

.initialize_string_db <- function(
    string_version,
    organism,
    string_score_threshold,
    string_cache_dir
) {
    tryCatch(
        STRINGdb::STRINGdb$new(
            version = string_version,
            species = organism,
            score_threshold = string_score_threshold,
            input_directory = string_cache_dir
        ),
        error = function(e) NULL
    )
}

.read_offline_string <- function(
    info_file,
    edge_file,
    string_score_threshold
) {
    protein_info <- utils::read.delim(
        gzfile(info_file),
        skip = 1L,
        header = FALSE,
        quote = "",
        stringsAsFactors = FALSE,
        col.names = c(
            "STRING_id", "preferred_name", "protein_size", "annotation"
        )
    )
    edge_table <- .read_cache_file(edge_file)
    required_columns <- c("protein1", "protein2", "combined_score")
    valid_edges <- !is.null(edge_table) &&
        all(required_columns %in% colnames(edge_table)) &&
        all(is.finite(edge_table$combined_score)) &&
        all(edge_table$combined_score >= string_score_threshold)
    if (!valid_edges) {
        stop("The offline STRING edge cache failed validation")
    }
    list(
        symbol_to_id = stats::setNames(
            protein_info$STRING_id,
            protein_info$preferred_name
        ),
        graph = igraph::graph_from_data_frame(
            edge_table[, required_columns],
            directed = FALSE
        )
    )
}

.load_offline_string <- function(
    string_cache_dir,
    organism,
    string_version,
    string_score_threshold
) {
    cache_key <- paste0(
        "string_offline_", organism, "_v", string_version,
        "_s", string_score_threshold
    )
    if (exists(cache_key, envir = .geneselectr2_cache)) {
        return(get(cache_key, envir = .geneselectr2_cache))
    }

    info_file <- file.path(
        string_cache_dir,
        sprintf("%s.protein.info.v%s.txt.gz", organism, string_version)
    )
    edge_file <- file.path(
        string_cache_dir,
        sprintf(
            "%s.protein.links.score%d.v%s.rds",
            organism,
            string_score_threshold,
            string_version
        )
    )
    if (!file.exists(info_file) || !file.exists(edge_file)) {
        cache_error <- paste0(
            "STRING initialization failed; the offline cache requires ",
            "%s and %s"
        )
        stop(
            sprintf(
                cache_error,
                basename(info_file),
                basename(edge_file)
            ),
            call. = FALSE
        )
    }

    resources <- .read_offline_string(
        info_file,
        edge_file,
        string_score_threshold
    )
    assign(cache_key, resources, envir = .geneselectr2_cache)
    resources
}

.prepare_string_resources <- function(
    string_db,
    all_symbols,
    string_cache_dir,
    organism,
    string_version,
    string_score_threshold,
    verbose
) {
    if (is.null(string_db)) {
        offline <- .load_offline_string(
            string_cache_dir,
            organism,
            string_version,
            string_score_threshold
        )
        mapped_ids <- unname(offline$symbol_to_id[all_symbols])
        mapping <- data.frame(
            gene = all_symbols[!is.na(mapped_ids)],
            STRING_id = mapped_ids[!is.na(mapped_ids)],
            stringsAsFactors = FALSE
        )
        if (verbose) {
            message("  STRING API unavailable; using cached STRING files\n")
        }
        return(list(
            mapping = mapping,
            graph = offline$graph,
            offline = TRUE
        ))
    }

    mapping <- tryCatch(
        string_db$map(
            data.frame(gene = all_symbols, stringsAsFactors = FALSE),
            "gene",
            removeUnmappedRows = TRUE
        ),
        error = function(e) NULL
    )
    list(mapping = mapping, graph = NULL, offline = FALSE)
}

.load_string_graph <- function(
    resources,
    string_db,
    string_cache_dir,
    string_version,
    organism,
    string_score_threshold
) {
    if (resources$offline) {
        return(resources$graph)
    }
    graph <- tryCatch(string_db$get_graph(), error = function(e) NULL)
    if (!is.null(graph)) {
        return(graph)
    }

    links <- list.files(
        string_cache_dir,
        pattern = "protein\\.links",
        full.names = TRUE
    )
    if (length(links) > 0L) {
        warning(
            sprintf(
                "Removing incomplete STRING files and retrying: %s",
                paste(basename(links), collapse = ", ")
            ),
            call. = FALSE
        )
        unlink(links)
    }
    string_db <- .initialize_string_db(
        string_version,
        organism,
        string_score_threshold,
        string_cache_dir
    )
    if (is.null(string_db)) {
        return(NULL)
    }
    tryCatch(string_db$get_graph(), error = function(e) NULL)
}

.prepare_restart_vector <- function(seeds, mapping, graph_nodes) {
    score_by_symbol <- stats::setNames(seeds$score, seeds$symbol)
    mapping_scores <- score_by_symbol[mapping$gene]
    seed_scores <- stats::setNames(
        mapping_scores[!is.na(mapping_scores)],
        mapping$STRING_id[!is.na(mapping_scores)]
    )
    if (anyDuplicated(names(seed_scores))) {
        seed_scores <- tapply(seed_scores, names(seed_scores), max)
    }

    present_ids <- intersect(names(seed_scores), graph_nodes)
    weights <- seed_scores[present_ids]
    valid <- is.finite(weights) & weights > 0
    present_ids <- present_ids[valid]
    weights <- weights[valid]
    if (length(present_ids) == 0L || sum(weights) <= 0) {
        return(NULL)
    }

    restart <- stats::setNames(rep(0, length(graph_nodes)), graph_nodes)
    restart[present_ids] <- weights
    restart <- restart / sum(restart)
    if (any(!is.finite(restart))) {
        return(NULL)
    }
    restart
}

.run_string_pagerank <- function(graph, restart_vector, restart_prob) {
    edge_attributes <- igraph::edge_attr_names(graph)
    weight_name <- intersect(
        c("combined_score", "score", "weight"),
        edge_attributes
    )
    graph_weights <- if (length(weight_name) == 0L) {
        NULL
    } else {
        igraph::edge_attr(graph, weight_name[1])
    }
    if (!is.null(graph_weights) && any(!is.finite(graph_weights))) {
        stop("STRING graph contains non-finite edge weights")
    }
    igraph::page_rank(
        graph,
        damping = 1 - restart_prob,
        personalized = restart_vector[igraph::V(graph)$name],
        weights = graph_weights
    )$vector
}

.report_string_run <- function(
    mapping,
    all_symbols,
    seeds,
    graph,
    restart_vector,
    restart_prob
) {
    n_seed_mapped <- sum(mapping$gene %in% seeds$symbol)
    message(sprintf(
        "  STRING mapping: %d/%d symbols mapped, %d seeds mapped\n",
        nrow(mapping),
        length(all_symbols),
        n_seed_mapped
    ))
    message(sprintf(
        "  RWR on STRING: %d nodes, %d edges, %d seeds (restart=%.2f)\n",
        igraph::vcount(graph),
        igraph::ecount(graph),
        sum(restart_vector > 0),
        restart_prob
    ))
}

.map_candidate_network_scores <- function(
    genes,
    mapping,
    network_scores,
    verbose
) {
    symbol_to_id <- stats::setNames(mapping$STRING_id, mapping$gene)
    candidate_ids <- symbol_to_id[genes]
    scores <- stats::setNames(rep(0, length(genes)), genes)
    in_graph <- !is.na(candidate_ids) &
        candidate_ids %in% names(network_scores)
    scores[in_graph] <- network_scores[candidate_ids[in_graph]]
    if (verbose) {
        positive <- scores[scores > 0]
        message(sprintf(
            "  [network] %d/%d candidate genes scored; range [%.2e, %.2e]\n",
            sum(in_graph),
            length(genes),
            if (length(positive)) min(positive) else 0,
            if (length(positive)) max(positive) else 0
        ))
    }
    scores
}

.prepare_string_context <- function(
    genes,
    seeds,
    string_score_threshold,
    organism,
    string_version,
    verbose
) {
    cache_dir <- getOption(
        "GeneSelectR2.string_cache",
        file.path(get_cache_dir(), "stringdb")
    )
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    string_db <- .initialize_string_db(
        string_version, organism, string_score_threshold, cache_dir
    )
    all_symbols <- unique(c(genes, seeds$symbol))
    resources <- .prepare_string_resources(
        string_db, all_symbols, cache_dir, organism,
        string_version, string_score_threshold, verbose
    )
    if (is.null(resources$mapping) || nrow(resources$mapping) < 5L) {
        stop("Fewer than five candidate or seed genes mapped to STRING")
    }
    graph <- .load_string_graph(
        resources, string_db, cache_dir, string_version,
        organism, string_score_threshold
    )
    if (is.null(graph) || igraph::vcount(graph) == 0L ||
        igraph::ecount(graph) == 0L) {
        stop("The complete STRING graph could not be loaded")
    }
    list(
        mapping = resources$mapping,
        all_symbols = all_symbols,
        graph = igraph::simplify(graph, edge.attr.comb = "max")
    )
}

.score_string_network <- function(
    genes,
    seeds,
    string_score_threshold,
    restart_prob,
    organism,
    string_version,
    verbose
) {
    context <- .prepare_string_context(
        genes, seeds, string_score_threshold,
        organism, string_version, verbose
    )
    graph <- context$graph
    mapping <- context$mapping
    restart <- .prepare_restart_vector(seeds, mapping, igraph::V(graph)$name)
    if (is.null(restart)) {
        warning(
            "No finite, positive seed weights map into the STRING graph",
            call. = FALSE
        )
        return(stats::setNames(rep(0, length(genes)), genes))
    }
    if (verbose) {
        .report_string_run(
            mapping, context$all_symbols, seeds, graph, restart, restart_prob
        )
    }

    network_scores <- .run_string_pagerank(graph, restart, restart_prob)
    raw_scores <- .map_candidate_network_scores(
        genes, mapping, network_scores, verbose
    )
    percentile01(raw_scores)
}

#' Score Genes by Network Proximity to Disease Seeds
#'
#' Places disease seed genes on the STRING functional interaction network and
#' propagates their Open Targets association scores with a random walk with
#' restart. Each candidate gene receives its stationary probability from the
#' network propagation.
#'
#' @param genes Character vector of candidate gene symbols to score
#' @param disease_term Disease name or EFO ID (passed to Open Targets)
#' @param string_score_threshold Integer, minimum STRING combined score for an
#'   edge to be included (default: 400, STRING's "medium confidence")
#' @param restart_prob Numeric, RWR restart probability (default: 0.5). Higher
#'   values retain more probability near the seeds.
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
score_network_layer <- function(
    genes,
    disease_term,
    string_score_threshold = 400,
    restart_prob = 0.5,
    max_seeds = 100,
    min_score = 0.1,
    organism = 9606,
    string_version = "12.0",
    use_cache = TRUE,
    verbose = TRUE
) {
    .require_packages(c("STRINGdb", "igraph"), "Network biology")
    if (!is.character(genes) || length(genes) == 0L || anyNA(genes)) {
        stop("genes must be a non-empty character vector")
    }
    if (!is.numeric(restart_prob) || length(restart_prob) != 1L ||
        !is.finite(restart_prob) || restart_prob <= 0 || restart_prob >= 1) {
        stop("restart_prob must be one finite number in (0, 1)")
    }

    seeds <- get_disease_seeds_opentargets(
        disease_term,
        max_seeds = max_seeds,
        min_score = min_score,
        use_cache = use_cache,
        verbose = verbose
    )
    if (nrow(seeds) == 0L) {
        warning("No disease seeds are available; returning zero scores")
        return(stats::setNames(rep(0, length(genes)), genes))
    }

    previous_timeout <- getOption("timeout")
    if (is.numeric(previous_timeout) && previous_timeout < 3600) {
        options(timeout = 3600)
        on.exit(options(timeout = previous_timeout), add = TRUE)
    }
    .score_string_network(
        genes, seeds, string_score_threshold, restart_prob,
        organism, string_version, verbose
    )
}

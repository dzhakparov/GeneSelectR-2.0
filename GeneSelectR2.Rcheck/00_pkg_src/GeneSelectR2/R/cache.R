#' Get Cache Directory (CRAN-compliant)
#'
#' Uses tools::R_user_dir() for proper cache location
#'
#' @return Character path to cache directory
#' @keywords internal
get_cache_dir <- function() {
    # CRAN-compliant cache directory
    cache_dir <- tools::R_user_dir("GeneSelectR", "cache")

    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }

    return(cache_dir)
}

# Global cache environment (package-level)
.geneselectr_cache <- new.env(parent = emptyenv())


#' Load GO Annotations with Caching
#'
#' Loads GO annotations (gene symbol -> GO term IDs) from memory cache,
#' disk cache, or org.Hs.eg.db. Results are cached for fast reuse.
#'
#' IMPORTANT: This function guards against empty caches. If a previous
#' download failed and saved an empty list to disk, the stale file is
#' removed and the download is retried.
#'
#' @param organism Character, organism name (default: "human")
#' @param force_reload Logical, ignore cache and reload (default: FALSE)
#' @return Named list mapping gene symbols to GO term ID vectors
#'
#' @examples
#' if (interactive() && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'     go_cache <- load_go_cache("human")
#'     go_cache[["TP53"]]
#' }
#'
#' @export
load_go_cache <- function(organism = "human", force_reload = FALSE) {
    org_pkg <- switch(organism,
        human = "org.Hs.eg.db",
        mouse = "org.Mm.eg.db",
        stop(sprintf("Unsupported organism: %s", organism), call. = FALSE)
    )
    annotation_version <- if (requireNamespace(org_pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(org_pkg))
    } else {
        "missing"
    }
    cache_key <- paste0("go_v2_", organism, "_", annotation_version)
    cache_file <- file.path(get_cache_dir(), paste0(cache_key, ".rds"))

    # Check memory cache first (fastest)
    if (!force_reload && exists(cache_key, envir = .geneselectr_cache)) {
        cached <- get(cache_key, envir = .geneselectr_cache)
        # Guard: don't return empty cache from a previous failed download
        if (length(cached) > 0) {
            if (getOption("geneselectr.verbose", TRUE)) {
                message(sprintf(
                    "Using GO annotations from memory cache (%d genes)",
                    length(cached)
                ))
            }
            return(cached)
        }
    }

    # Check disk cache (fast)
    if (!force_reload && file.exists(cache_file)) {
        go_data <- readRDS(cache_file)

        # Guard: don't use empty disk cache
        if (length(go_data) > 0) {
            if (getOption("geneselectr.verbose", TRUE)) {
                message(sprintf(
                    "Loading GO annotations from disk cache (%d genes)",
                    length(go_data)
                ))
            }
            assign(cache_key, go_data, envir = .geneselectr_cache)
            return(go_data)
        } else {
            # Remove stale empty cache so we retry the download
            if (getOption("geneselectr.verbose", TRUE)) {
                message(
                    "Found an empty GO cache on disk; rebuilding the cache."
                )
            }
            file.remove(cache_file)
        }
    }

    # Build from database
    if (getOption("geneselectr.verbose", TRUE)) {
        message("Building GO annotation cache (one-time operation)...")
    }

    go_data <- download_go_annotations(organism)

    # Only cache if we got real data
    if (length(go_data) > 0) {
        if (getOption("geneselectr.verbose", TRUE)) {
            message(sprintf(
                "Caching GO annotations for %d genes to disk",
                length(go_data)
            ))
        }
        saveRDS(go_data, cache_file)
        assign(cache_key, go_data, envir = .geneselectr_cache)
    } else {
        annotation_message <- paste(
            "GO annotation loading returned no data for %s.",
            "Install and verify %s."
        )
        stop(
            sprintf(
                annotation_message,
                organism, org_pkg
            ),
            call. = FALSE
        )
    }

    return(go_data)
}


#' Download GO Annotations from Bioconductor Annotation Database
#'
#' Extracts gene symbol -> GO term mappings from org.Hs.eg.db (human) or
#' org.Mm.eg.db (mouse). Returns a named list where each element is a
#' character vector of GO term IDs for that gene.
#'
#' @param organism Character, "human" or "mouse"
#' @importFrom stats median
#' @return Named list: gene symbol -> character vector of GO term IDs.
#' @details Missing annotation packages and extraction failures terminate with
#'   an error so that an operational failure cannot enter a biological score.
#' @keywords internal
download_go_annotations <- function(organism = "human") {
    # --- Determine which OrgDb to use ---
    org_pkg <- switch(organism,
        human = "org.Hs.eg.db",
        mouse = "org.Mm.eg.db",
        stop(
            "Unsupported organism: '", organism,
            "'. For another organism, pass a pre-built go_cache ",
            "to biological_scorer()."
        )
    )

    if (!requireNamespace(org_pkg, quietly = TRUE)) {
        stop(
            org_pkg, " is required for GO annotation scoring. Install with: ",
            "BiocManager::install('", org_pkg, "')"
        )
    }

    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
        stop(
            "AnnotationDbi is required for GO annotation scoring. ",
            "Install with: BiocManager::install('AnnotationDbi')"
        )
    }

    # --- Extract all gene-to-GO mappings ---
    orgdb <- getExportedValue(org_pkg, org_pkg)

    if (getOption("geneselectr.verbose", TRUE)) {
        message(sprintf("  Extracting GO annotations from %s...", org_pkg))
    }

    tryCatch(
        {
            # Get all gene symbols with GO annotations via select()
            # AnnotationDbi returns SYMBOL, GO, EVIDENCE and ONTOLOGY columns.
            all_symbols <- AnnotationDbi::keys(orgdb, keytype = "SYMBOL")

            if (getOption("geneselectr.verbose", TRUE)) {
                message(sprintf(
                    "  Found %d gene symbols in %s",
                    length(all_symbols), org_pkg
                ))
            }

            # Query in chunks to avoid memory issues with very large databases
            chunk_size <- 5000
            n_chunks <- ceiling(length(all_symbols) / chunk_size)
            go_table_list <- vector("list", n_chunks)

            for (i in seq_len(n_chunks)) {
                start_idx <- (i - 1) * chunk_size + 1
                end_idx <- min(i * chunk_size, length(all_symbols))
                chunk_symbols <- all_symbols[start_idx:end_idx]

                chunk_result <- AnnotationDbi::select(
                    orgdb,
                    keys = chunk_symbols,
                    columns = c("SYMBOL", "GO", "EVIDENCE"),
                    keytype = "SYMBOL"
                )

                go_table_list[[i]] <- chunk_result
            }

            go_table <- do.call(rbind, go_table_list)

            # Remove rows with NA GO terms
            go_table <- go_table[!is.na(go_table$GO), ]

            # The ND code denotes missing biological annotation.
            if ("EVIDENCE" %in% colnames(go_table)) {
                go_table <- go_table[go_table$EVIDENCE != "ND", ]
            }

            if (nrow(go_table) == 0) {
                stop(
                    sprintf("No GO annotations were found in %s", org_pkg),
                    call. = FALSE
                )
            }

            # Split into named list: gene symbol -> unique GO term IDs
            go_cache <- split(go_table$GO, go_table$SYMBOL)
            go_cache <- lapply(go_cache, unique)

            # Remove genes with no terms (safety check)
            go_cache <- go_cache[lengths(go_cache) > 0]

            if (getOption("geneselectr.verbose", TRUE)) {
                n_genes <- length(go_cache)
                n_terms <- length(unique(unlist(go_cache)))
                median_terms <- median(lengths(go_cache))
                cache_format <- paste(
                    "  Loaded %d genes with GO annotations",
                    "(%d unique terms; median %d terms per gene)"
                )
                message(sprintf(
                    cache_format,
                    n_genes, n_terms, median_terms
                ))
            }

            return(go_cache)
        },
        error = function(e) {
            stop(
                sprintf(
                    "Failed to extract GO annotations from %s: %s",
                    org_pkg, e$message
                ),
                call. = FALSE
            )
        }
    )
}


#' Load Information Content Scores with Caching
#'
#' @param go_cache GO annotation cache from load_go_cache()
#' @param force_reload Logical, ignore cache and recompute
#' @return Named vector of IC scores
#'
#' @keywords internal
load_ic_cache <- function(go_cache, force_reload = FALSE) {
    # IC depends on the annotation corpus it is computed from, and the caller
    # passes the ONTOLOGY-FILTERED cache. A constant key ("ic_scores") therefore
    # returned BP information content to a caller asking for MF whenever both
    # ontologies were used in one session. Key the cache on the corpus itself:
    # Include the actual term assignments. Gene names and term counts can remain
    # unchanged across annotation releases while the assigned terms change.
    corpus_hash <- if (requireNamespace("digest", quietly = TRUE)) {
        ordered_cache <- go_cache[sort(names(go_cache))]
        ordered_cache <- lapply(ordered_cache, sort)
        digest::digest(ordered_cache, algo = "md5")
    } else {
        # Fallback without digest: never share across corpora.
        paste0("nocache_", as.numeric(Sys.time()))
    }
    cache_key <- paste0("ic_scores_", corpus_hash)
    cache_file <- file.path(get_cache_dir(), paste0(cache_key, ".rds"))

    # Check memory cache
    if (!force_reload && exists(cache_key, envir = .geneselectr_cache)) {
        return(get(cache_key, envir = .geneselectr_cache))
    }

    # Check disk cache
    if (!force_reload && file.exists(cache_file)) {
        ic_scores <- readRDS(cache_file)
        assign(cache_key, ic_scores, envir = .geneselectr_cache)
        return(ic_scores)
    }

    # Compute IC scores
    if (getOption("geneselectr.verbose", TRUE)) {
        message("Computing information content scores...")
    }

    ic_scores <- compute_information_content(go_cache)

    # Save to cache
    saveRDS(ic_scores, cache_file)
    assign(cache_key, ic_scores, envir = .geneselectr_cache)

    return(ic_scores)
}

#' Create Similarity Cache Environment
#'
#' @return Environment for storing GO similarity scores
#' @keywords internal
create_similarity_cache <- function() {
    cache <- new.env(parent = emptyenv())
    return(cache)
}

#' Get or Compute GO Similarity with Caching
#'
#' Checks the cache for a previously computed similarity between two GO terms.
#' If not found, computes it using the specified method and stores it.
#'
#' The cache key includes the similarity method, so switching methods
#' (e.g., from "resnik" to "lin") does not return stale cached values.
#'
#' @param term1 Character, first GO term
#' @param term2 Character, second GO term
#' @param ic_scores Named vector of IC scores
#' @param cache Environment for caching similarities
#' @param ancestor_map Named list mapping GO terms to their ancestors
#' @param sim_method Character, similarity metric: "resnik", "lin", "jiang"
#'   or "rel".
#' @return Numeric similarity score
#'
#' @keywords internal
get_or_compute_similarity <- function(term1, term2, ic_scores, cache,
                                        ancestor_map = NULL,
                                        sim_method = "resnik") {
    # Cache key includes the method so different metrics don't collide
    cache_key <- paste(sim_method, paste(sort(c(term1, term2)), collapse = "_"),
        sep = ":"
    )

    # Check cache
    if (exists(cache_key, envir = cache)) {
        return(get(cache_key, envir = cache))
    }

    # Compute using the dispatcher that supports all four metrics
    similarity <- compute_semantic_similarity(
        term1, term2, ic_scores, ancestor_map,
        method = sim_method
    )

    # Store in cache
    assign(cache_key, similarity, envir = cache)

    return(similarity)
}

#' Clear All Cached Data
#'
#' Removes both memory and disk caches.
#' Use this if you want to force fresh downloads.
#'
#' @param confirm Logical, require confirmation (default: TRUE)
#'
#' @examples
#' if (interactive()) {
#'     clear_cache(confirm = FALSE)
#' }
#' @return `NULL`, invisibly.
#'
#' @export
clear_cache <- function(confirm = TRUE) {
    if (confirm) {
        response <- readline("Clear all cached data? (yes/no): ")
        if (tolower(response) != "yes") {
            message("Cache clearing cancelled")
            return(invisible(NULL))
        }
    }

    # Clear disk cache
    cache_dir <- get_cache_dir()
    if (dir.exists(cache_dir)) {
        unlink(cache_dir, recursive = TRUE)
        message("Disk cache cleared")
    }

    # Clear memory cache
    rm(list = ls(envir = .geneselectr_cache), envir = .geneselectr_cache)
    message("Memory cache cleared")

    message("All caches cleared successfully")
    invisible(NULL)
}

#' Inspect cached annotation data
#'
#' @return A list containing the cache directory, disk files, total disk size,
#'   and memory-cache keys.
#' @examples
#' cache_info()
#' @export
cache_info <- function() {
    cache_dir <- get_cache_dir()
    files <- list.files(cache_dir, full.names = TRUE)
    file_table <- data.frame(
        file = basename(files),
        size_bytes = as.numeric(file.size(files)),
        modified = as.POSIXct(file.mtime(files)),
        stringsAsFactors = FALSE
    )
    list(
        directory = cache_dir,
        files = file_table,
        total_size_bytes = sum(file_table$size_bytes),
        memory_keys = ls(envir = .geneselectr_cache)
    )
}

#' Set Cache Options
#'
#' Configure caching behavior
#'
#' @param verbose Logical, print cache messages
#' @return `NULL`, invisibly.
#' @examples
#' set_cache_options(verbose = FALSE)
#'
#' @export
set_cache_options <- function(verbose = TRUE) {
    options(geneselectr.verbose = verbose)

    invisible(NULL)
}

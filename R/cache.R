#' Locate the package cache directory
#'
#' Uses `tools::R_user_dir()` to obtain a platform-specific cache location.
#'
#' @return Character path to the cache directory.
#' @keywords internal
get_cache_dir <- function() {
    cache_dir <- tools::R_user_dir("GeneSelectR2", "cache")
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    cache_dir
}

## In-memory objects are restricted to the package namespace.
.geneselectr2_cache <- new.env(parent = emptyenv())

.read_cache_file <- function(path) {
    if (!file.exists(path)) {
        return(NULL)
    }
    tryCatch(
        readRDS(path),
        error = function(error) {
            unlink(path)
            NULL
        }
    )
}

.write_cache_file <- function(value, path) {
    temporary_file <- tempfile(
        pattern = paste0(basename(path), "."),
        tmpdir = dirname(path)
    )
    on.exit(unlink(temporary_file), add = TRUE)
    saveRDS(value, temporary_file)
    installed <- file.rename(temporary_file, path)
    if (!installed) {
        installed <- file.copy(temporary_file, path, overwrite = TRUE)
    }
    if (!installed) {
        warning("The cache file could not be written: ", path, call. = FALSE)
    }
    invisible(installed)
}


#' Load GO annotations
#'
#' Loads gene-symbol-to-GO-term mappings from memory, the package cache, or an
#' organism annotation package. Empty cache entries are removed and rebuilt.
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

    if (!force_reload && exists(cache_key, envir = .geneselectr2_cache)) {
        cached <- get(cache_key, envir = .geneselectr2_cache)
        if (length(cached) > 0) {
            if (getOption("geneselectr2.verbose", TRUE)) {
                message(sprintf(
                    "Using GO annotations from memory cache (%d genes)",
                    length(cached)
                ))
            }
            return(cached)
        }
    }

    if (!force_reload) {
        go_data <- .read_cache_file(cache_file)
        if (!is.null(go_data) && length(go_data) > 0L) {
            if (getOption("geneselectr2.verbose", TRUE)) {
                message(sprintf(
                    "Loading GO annotations from disk cache (%d genes)",
                    length(go_data)
                ))
            }
            assign(cache_key, go_data, envir = .geneselectr2_cache)
            return(go_data)
        }
        if (!is.null(go_data)) {
            if (getOption("geneselectr2.verbose", TRUE)) {
                message(
                    "Found an empty GO cache on disk; rebuilding the cache."
                )
            }
            file.remove(cache_file)
        }
    }

    if (getOption("geneselectr2.verbose", TRUE)) {
        message("Building the GO annotation cache...")
    }

    go_data <- download_go_annotations(organism)

    if (length(go_data) > 0) {
        if (getOption("geneselectr2.verbose", TRUE)) {
            message(sprintf(
                "Caching GO annotations for %d genes to disk",
                length(go_data)
            ))
        }
        .write_cache_file(go_data, cache_file)
        assign(cache_key, go_data, envir = .geneselectr2_cache)
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

    go_data
}


#' Extract GO annotations from a Bioconductor annotation package
#'
#' Extracts gene symbol -> GO term mappings from org.Hs.eg.db (human) or
#' org.Mm.eg.db (mouse). Returns a named list where each element is a
#' character vector of GO term IDs for that gene.
#'
#' @param organism Character, "human" or "mouse"
#' @return Named list: gene symbol -> character vector of GO term IDs.
#' @details Missing annotation packages and extraction failures terminate with
#'   an error so that an operational failure cannot enter a biological score.
#' @keywords internal
download_go_annotations <- function(organism = "human") {
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

    orgdb <- getExportedValue(org_pkg, org_pkg)

    if (getOption("geneselectr2.verbose", TRUE)) {
        message(sprintf("  Extracting GO annotations from %s...", org_pkg))
    }

    tryCatch(
        {
            all_symbols <- AnnotationDbi::keys(orgdb, keytype = "SYMBOL")

            if (getOption("geneselectr2.verbose", TRUE)) {
                message(sprintf(
                    "  Found %d gene symbols in %s",
                    length(all_symbols), org_pkg
                ))
            }

            ## Chunked queries limit memory use for large annotation databases.
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

            go_table <- go_table[!is.na(go_table$GO), ]

            ## The ND evidence code denotes an absence of biological data.
            if ("EVIDENCE" %in% colnames(go_table)) {
                go_table <- go_table[go_table$EVIDENCE != "ND", ]
            }

            if (nrow(go_table) == 0) {
                stop(
                    sprintf("No GO annotations were found in %s", org_pkg),
                    call. = FALSE
                )
            }

            go_cache <- split(go_table$GO, go_table$SYMBOL)
            go_cache <- lapply(go_cache, unique)
            go_cache <- go_cache[lengths(go_cache) > 0]

            if (getOption("geneselectr2.verbose", TRUE)) {
                n_genes <- length(go_cache)
                n_terms <- length(unique(unlist(go_cache)))
                median_terms <- stats::median(lengths(go_cache))
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
    ## The corpus hash prevents reuse across ontology branches and annotation
    ## releases with different gene-to-term assignments.
    corpus_hash <- if (requireNamespace("digest", quietly = TRUE)) {
        ordered_cache <- go_cache[sort(names(go_cache))]
        ordered_cache <- lapply(ordered_cache, sort)
        digest::digest(ordered_cache, algo = "md5")
    } else {
        ## Without digest, each call receives a distinct cache key.
        paste0("nocache_", as.numeric(Sys.time()))
    }
    cache_key <- paste0("ic_scores_", corpus_hash)
    cache_file <- file.path(get_cache_dir(), paste0(cache_key, ".rds"))

    if (!force_reload && exists(cache_key, envir = .geneselectr2_cache)) {
        return(get(cache_key, envir = .geneselectr2_cache))
    }

    if (!force_reload) {
        ic_scores <- .read_cache_file(cache_file)
    } else {
        ic_scores <- NULL
    }
    if (!is.null(ic_scores)) {
        assign(cache_key, ic_scores, envir = .geneselectr2_cache)
        return(ic_scores)
    }

    if (getOption("geneselectr2.verbose", TRUE)) {
        message("Computing information content scores...")
    }

    ic_scores <- compute_information_content(go_cache)

    .write_cache_file(ic_scores, cache_file)
    assign(cache_key, ic_scores, envir = .geneselectr2_cache)
    ic_scores
}

#' Create Similarity Cache Environment
#'
#' @return Environment for storing GO similarity scores
#' @keywords internal
create_similarity_cache <- function() {
    new.env(parent = emptyenv())
}

#' Get or Compute GO Similarity with Caching
#'
#' Retrieves a previously computed GO-term similarity or computes and stores a
#' new value. The cache key includes the similarity method.
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
    cache_key <- paste(sim_method, paste(sort(c(term1, term2)), collapse = "_"),
        sep = ":"
    )

    if (exists(cache_key, envir = cache)) {
        return(get(cache_key, envir = cache))
    }

    similarity <- compute_semantic_similarity(
        term1, term2, ic_scores, ancestor_map,
        method = sim_method
    )

    assign(cache_key, similarity, envir = cache)

    similarity
}

#' Clear All Cached Data
#'
#' Removes the package's memory and disk caches.
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

    cache_dir <- get_cache_dir()
    if (dir.exists(cache_dir)) {
        unlink(cache_dir, recursive = TRUE)
        message("Disk cache cleared")
    }

    rm(list = ls(envir = .geneselectr2_cache), envir = .geneselectr2_cache)
    message("Memory cache cleared")
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
        memory_keys = ls(envir = .geneselectr2_cache)
    )
}

#' Set cache options
#'
#' Configures diagnostic messages from cache operations.
#'
#' @param verbose Logical, print cache messages
#' @return `NULL`, invisibly.
#' @examples
#' set_cache_options(verbose = FALSE)
#'
#' @export
set_cache_options <- function(verbose = TRUE) {
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("verbose must be TRUE or FALSE", call. = FALSE)
    }
    options(geneselectr2.verbose = verbose)
    invisible(NULL)
}

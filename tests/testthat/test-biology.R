test_that("information content reflects annotation frequency", {
    go_annotations <- list(
        gene_a = c("GO:ROOT", "GO:RARE"),
        gene_b = "GO:ROOT"
    )

    information_content <- GeneSelectR2:::compute_information_content(
        go_annotations
    )

    expect_equal(unname(information_content["GO:ROOT"]), 0)
    expect_gt(information_content["GO:RARE"], information_content["GO:ROOT"])
    expect_true(all(is.finite(information_content)))
})

test_that("semantic similarity uses the most informative common ancestor", {
    information_content <- c(
        "GO:ROOT" = 0,
        "GO:PARENT" = 2,
        "GO:CHILD" = 3,
        "GO:OTHER" = 4
    )
    ancestors <- list(
        "GO:ROOT" = character(),
        "GO:PARENT" = "GO:ROOT",
        "GO:CHILD" = c("GO:PARENT", "GO:ROOT"),
        "GO:OTHER" = "GO:ROOT"
    )

    parent_child <- GeneSelectR2:::compute_semantic_similarity(
        "GO:PARENT", "GO:CHILD",
        information_content, ancestors,
        method = "resnik"
    )
    unrelated <- GeneSelectR2:::compute_semantic_similarity(
        "GO:CHILD", "GO:OTHER",
        information_content, ancestors,
        method = "resnik"
    )

    expect_equal(parent_child, 2 / 3)
    expect_equal(unrelated, 0)
    expect_equal(GeneSelectR2:::sim_lin(2, 2, 3), 0.8)
    expect_equal(GeneSelectR2:::sim_jiang(2, 2, 3), 0.5)
})

test_that("term matrices preserve dimensions and symmetry", {
    information_content <- c(
        "GO:ROOT" = 0,
        "GO:A" = 2,
        "GO:B" = 3
    )
    ancestors <- list(
        "GO:ROOT" = character(),
        "GO:A" = "GO:ROOT",
        "GO:B" = c("GO:A", "GO:ROOT")
    )

    target_matrix <- GeneSelectR2:::build_target_similarity_matrix(
        candidate_terms = c("GO:A", "GO:B"),
        target_terms = "GO:B",
        ic_scores = information_content,
        ancestor_map = ancestors,
        use_cache = FALSE,
        verbose = FALSE
    )
    square_matrix <- GeneSelectR2:::build_term_similarity_matrix(
        term_universe = c("GO:ROOT", "GO:A", "GO:B"),
        ic_scores = information_content,
        ancestor_map = ancestors,
        use_cache = FALSE,
        verbose = FALSE
    )

    expect_identical(dim(target_matrix), c(2L, 1L))
    expect_equal(target_matrix["GO:B", "GO:B"], 1)
    expect_equal(square_matrix, t(square_matrix))
    expect_equal(unname(diag(square_matrix)), rep(1, 3))
    expect_equal(
        GeneSelectR2:::bma_from_matrix("GO:A", "GO:B", target_matrix),
        target_matrix["GO:A", "GO:B"]
    )
})

test_that("GO enrichment returns adjusted probabilities", {
    go_annotations <- list(
        selected_a = c("GO:1", "GO:2"),
        selected_b = "GO:1",
        background_a = "GO:2",
        background_b = "GO:3"
    )

    enrichment <- GeneSelectR2:::test_go_enrichment(
        selected_genes = c("selected_a", "selected_b"),
        background_genes = c("background_a", "background_b"),
        go_cache = go_annotations
    )

    expect_named(enrichment, c(
        "term", "p_value", "odds_ratio", "n_selected",
        "n_background", "p_adj"
    ))
    expect_true(all(enrichment$p_adj >= 0 & enrichment$p_adj <= 1))
    expect_equal(enrichment$n_selected[enrichment$term == "GO:1"], 2)
})

test_that("data-driven scoring retains annotated genes", {
    gene_names <- paste0("gene", seq_len(12))
    go_annotations <- stats::setNames(
        rep(list(c("GO:PARENT", "GO:CHILD")), length(gene_names)),
        gene_names
    )
    information_content <- GeneSelectR2:::compute_information_content(
        go_annotations
    )
    ancestors <- list(
        "GO:PARENT" = character(),
        "GO:CHILD" = "GO:PARENT"
    )

    scores <- GeneSelectR2:::compute_data_driven_scores(
        genes = gene_names,
        enrichment_genes = gene_names,
        go_cache = go_annotations,
        ic_scores = information_content,
        similarity_cache = GeneSelectR2:::create_similarity_cache(),
        ancestor_map = ancestors,
        min_term_freq = 1,
        max_enriched_terms = 5,
        n_top_sims = 1,
        ic_quantile = 0
    )

    expect_length(scores, length(gene_names))
    expect_true(all(is.finite(scores)))
    expect_true(all(scores > 0))
})

test_that("biological scorer rejects invalid control parameters", {
    expect_error(
        biological_scorer(character(), enrich_fdr = 0.05),
        "genes must be"
    )
    expect_error(
        biological_scorer("TP53", mode = "unsupported"),
        "arg"
    )
    expect_error(
        biological_scorer("TP53", enrich_fdr = 2),
        "enrich_fdr"
    )
    expect_error(
        biological_scorer("TP53", n_top_sims = 0),
        "n_top_sims"
    )
})

test_that("ontology identifiers bypass free-text resolution", {
    efo_id <- GeneSelectR2:::resolve_efo_id("EFO:0000274")

    expect_identical(efo_id, "EFO_0000274")
})

test_that("GO.db fallback resolves ontology ancestors", {
    skip_if_not_installed("GO.db")

    ancestors <- suppressMessages(
        GeneSelectR2:::get_go_ancestors("GO:0006915")
    )

    expect_contains(ancestors, "GO:0006915")
    expect_gt(length(ancestors), 1L)
})

test_that("STRING restart weights retain mapped association scores", {
    seeds <- data.frame(
        symbol = c("A", "B"),
        score = c(0.2, 0.8)
    )
    mapping <- data.frame(
        gene = c("A", "B", "C"),
        STRING_id = c("id_a", "id_b", "id_c")
    )

    restart <- GeneSelectR2:::.prepare_restart_vector(
        seeds,
        mapping,
        graph_nodes = c("id_a", "id_b", "id_c")
    )

    expect_equal(sum(restart), 1)
    expect_equal(unname(restart[c("id_a", "id_b", "id_c")]), c(0.2, 0.8, 0))

    seeds$score <- c(NA, 0)
    expect_null(GeneSelectR2:::.prepare_restart_vector(
        seeds,
        mapping,
        graph_nodes = c("id_a", "id_b", "id_c")
    ))
})

test_that("disk cache helpers recover from invalid files", {
    cache_file <- tempfile(fileext = ".rds")
    on.exit(unlink(cache_file), add = TRUE)
    cached_value <- list(gene = "TP53", score = 1)

    expect_true(GeneSelectR2:::.write_cache_file(cached_value, cache_file))
    expect_identical(
        GeneSelectR2:::.read_cache_file(cache_file),
        cached_value
    )

    writeLines("invalid RDS data", cache_file)
    expect_null(GeneSelectR2:::.read_cache_file(cache_file))
    expect_false(file.exists(cache_file))
})

test_that("cache verbosity requires one logical value", {
    expect_error(set_cache_options(1), "TRUE or FALSE")
    expect_error(set_cache_options(NA), "TRUE or FALSE")
})

# GeneSelectR 2.0 (version 0.99.4)

- Extended both vignettes with study-design, Bioconductor interoperability,
  external-validation, application-scope, and provenance guidance.
- Documented all fields in the packaged benchmark result examples.
- Moved the default STRING cache under the platform-specific R package cache.
- Consolidated disk-cache reads and writes and added recovery from invalid
  cache files.
- Simplified Gene Ontology ancestor retrieval and removed redundant namespace
  checks.
- Added deterministic tests for semantic similarity, enrichment, and cache
  behavior.
- Reduced implementation comments to scientific assumptions and operational
  constraints, and clarified internal variable names.
- Removed committed package-check and source-archive output and corrected
  submission metadata.

# GeneSelectR 2.0 (version 0.99.3)

- Restricted the package implementation to the workflow evaluated in the
  manuscript: repeated elastic net, excluded-sample contribution,
  shuffled-outcome adjustment, and equal-weight score combination.
- Removed discarded selection gates, grouped-model variants, module fitting,
  alternative score formulas, and Rashomon gene-set generation.
- Removed the unused PubTator literature scoring path.
- Simplified the result table and documentation to use explicit gene-level
  measurement names.

# GeneSelectR 2.0 (version 0.99.2)

- Added package functions for plotting gene-level ranking measurements,
  gene-level predictive and disease evidence, and biological assessment
  against matched random gene sets.
- Added verified asthma case-study and biological benchmark example data.
- Extended the asthma vignette with reproducible biological and gene-level
  figures.

# GeneSelectR 2.0 (version 0.99.1)

- Added a predictive workflow interface with SummarizedExperiment support.
- Added a documented GSE69683 subset and an evaluated asthma vignette.
- Added input validation and workflow tests.
- Replaced speculative implementation comments with technical descriptions.

# GeneSelectR 0.99.3

- Restricted the package implementation to the workflow evaluated in the
  manuscript: repeated elastic net, excluded-sample contribution,
  shuffled-outcome adjustment, and equal-weight score combination.
- Removed discarded selection gates, grouped-model variants, module fitting,
  alternative score formulas, and Rashomon gene-set generation.
- Removed the unused PubTator literature scoring path.
- Simplified the result table and documentation to use explicit gene-level
  measurement names.

# GeneSelectR 0.99.2

- Added package functions for plotting gene-level ranking measurements,
  gene-level predictive and disease evidence, and biological assessment
  against matched random gene sets.
- Added verified asthma case-study and biological benchmark example data.
- Extended the asthma vignette with reproducible biological and gene-level
  figures.

# GeneSelectR 0.99.1

- Added a predictive workflow interface with SummarizedExperiment support.
- Added a documented GSE69683 subset and an evaluated asthma vignette.
- Added input validation and workflow tests.
- Replaced speculative implementation comments with technical descriptions.

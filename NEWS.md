# mttk 0.2.6

## New features

- `readMTTKExperiment()`: construct an `MTTKExperiment` directly from
  delimited files (TSV/CSV). Accepts a gene count matrix, gene annotations,
  sample metadata, genome counts, genome metadata, a Newick tree, and
  optional KO-to-module/pathway link tables. Separator is auto-detected from
  file extension. Link tables (`gene_to_genome`, `gene_to_ko`) are built
  automatically from the annotation file when `genome_id` and `ko_id` columns
  are present.
- `fromPhyloseq()`: convert one or two `phyloseq` objects into an
  `MTTKExperiment`. The primary object supplies gene-level RNA counts,
  sample metadata, and gene annotations (via `tax_table`). An optional
  `genomePhyloseq` argument provides genome-level DNA counts, genome
  metadata, and a genome phylogeny. Requires the `phyloseq` package
  (Bioconductor).

---

# mttk 0.2.5

## New features

- `fitKODispersion()`: fits a KO-level negative-binomial mixed model with a
  variable dispersion term (`dispformula = ~ variable`), testing whether
  expression variability (overdispersion) differs across conditions or along a
  continuous gradient. A positive estimate means higher phi (more
  Poisson-like, less variable) in the contrast group; a negative estimate means
  more overdispersed expression.
- `varianceDecomposition()`: extracts per-feature variance components from
  stored backend model objects in an `MTTKFit` (requires `keepFits = TRUE`).
  Returns genome random-effect variance, optional sample-block variance, the
  NB2 dispersion parameter phi, an approximated latent residual variance
  (`log(1 + 1/phi)`), total variance, and genome/block intraclass correlation
  coefficients (ICC) on the log link scale.

---

# mttk 0.2.4

## Bug fixes

- Fixed `library(mttk)` failing in Binder: `binder/install.R` now resolves the
  repo path via the `REPO_DIR` environment variable set by repo2docker, so
  `remotes::install_local()` finds the package regardless of working directory.
- Fixed plot functions (`plotVolcano`, `plotEffects`, `plotEffectHeatmap`,
  `plotGenomeEffects`, `plotCladeScan`) not being found after `library(mttk)`:
  the five exports were missing from `NAMESPACE` because `devtools::document()`
  had not been run after `R/MTTKPlot-fit.R` was added.
- Added `view`, `columns`, and `exclude` arguments to `fitTable()` for concise
  result tables; default `view = "standard"` returns feature ID, effect label,
  estimate, std error, p-value, q-value, and status only.

---

# mttk 0.2.3

## Bug fixes

- Fixed Binder environment: `binder/install.R` now installs mttk from the
  locally cloned repository files (`remotes::install_local()`) instead of
  re-downloading from GitHub, so `library(mttk)` works correctly in the
  Binder session.
- Removed redundant build-time warning from README (the Binder image is
  cached per repository and shared across all users).

---

# mttk 0.2.2

## Documentation

- Added `mttk_walkthrough.Rmd`: a self-contained, chunk-by-chunk interactive
  tutorial covering gene-, genome-, KO-, module-, and clade-level workflows on
  the built-in showcase dataset. Designed to be opened directly in the Binder
  RStudio session.
- Updated README to direct Binder users to `mttk_walkthrough.Rmd`.

---

# mttk 0.2.1

## Documentation

- Added Binder support: the README badge launches a browser-based RStudio
  session with mttk and all dependencies pre-installed, so users can run
  examples without a local R setup.

---

# mttk 0.2.0

## New features

### Result extraction

- `coefTable()` returns the full fixed-effect coefficient table stored inside an
  `MTTKFit`, with optional filtering by term, feature, and sort order. BH
  multiple-testing correction is applied automatically.
- `termFits()` splits a multi-term fit into a named list of per-term `MTTKFit`
  objects without refitting any models, reusing the stored backend objects and
  coefficient table.
- `koGenomeEffects()`, `moduleGenomeEffects()`, and `pathwayGenomeEffects()`
  return the genome-specific conditional effects stored in random-slope fits.

### Visualisation (`MTTKPlot-fit`)

Five `ggplot2`-based plotting functions are now available for `MTTKFit` objects.
All return plain `ggplot` objects that can be extended with additional layers.

- `plotVolcano()` — volcano plot (effect estimate vs. −log₁₀ q) for any
  `MTTKFit`, with optional top-hit labels.
- `plotEffects()` — horizontal ranked lollipop or bar chart of effect estimates
  with optional 95 % CI error bars.
- `plotEffectHeatmap()` — feature × genome heatmap of genome-specific
  conditional effects from random-slope fits, with optional hierarchical
  clustering of rows and columns.
- `plotGenomeEffects()` — phylogenetic tree overlaid with per-genome effect
  estimates; tips are coloured by direction and scaled by magnitude. Uses
  `ape` for the tree layout (no `ggtree` dependency).
- `plotCladeScan()` — phylogenetic tree where tips belonging to significant
  clades (from `scanGenomeClades()` or `scanKOClades()`) are highlighted by
  response direction.

### Workflow discovery

- `findWorkflow()` has been substantially expanded with a richer interactive
  decision tree, more biological questions, and ready-to-edit code templates
  for each recommended workflow.

### Internal

- Multi-term model fits now store a coefficient-level summary table in
  `metadata(x)$mttk_fit$coefficients`, enabling `coefTable()` and `termFits()`
  without re-fitting.
- `MTTKFit` validity checks extended to cover the `groupEffects` slot.

## Dependencies

- `ggplot2` added to `Imports`.

---

# mttk 0.1.0

- Initial release.

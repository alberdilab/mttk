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

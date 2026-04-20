#' MTTKExperiment Class
#'
#' `MTTKExperiment` stores genome-resolved metatranscriptomic measurements and
#' the metadata needed to analyze them across nested biological and functional
#' hierarchies.
#'
#' The class extends
#' [TreeSummarizedExperiment::TreeSummarizedExperiment-class], so it behaves like
#' a standard Bioconductor assay container while leaving room for tree-structured
#' relationships among features or samples.
#'
#' @section Stored data:
#' The core assay container is inherited from
#' `TreeSummarizedExperiment`:
#'
#' - `geneExperiment(x)` names the gene-level layer of the object. In practice,
#'   the `MTTKExperiment` itself is the gene-level experiment.
#' - `geneAssays(x)` hold gene-level measurements such as RNA counts,
#'   normalized RNA values, offsets, fitted values, or residuals.
#' - `rowData` stores feature-level metadata for genes or transcripts.
#' - `colData` stores sample-level metadata.
#' - `genomeExperiment(x)` stores genome-level assays such as genome DNA counts
#'   or aggregated genome RNA counts in a dedicated companion experiment with
#'   one row per genome.
#' - `genomeAssays(x)` names the assay collection stored in
#'   `genomeExperiment(x)`.
#' - `rowTree` and `rowLinks` can store optional named hierarchies when a
#'   relationship is naturally tree-structured.
#'
#' @section MTTK metadata contract:
#' MTTK-specific components are kept close to existing Bioconductor
#' infrastructure:
#'
#' - `genomeData(x)` is stored as `rowData(genomeExperiment(x))`.
#' - explicit assay names such as `"rna_gene_counts"`, `"rna_genome_counts"`,
#'   and `"dna_genome_counts"` make it clear which data layer an assay belongs
#'   to.
#' - `links`: a named `S4Vectors::SimpleList` of mapping tables, such as
#'   gene-to-genome links, gene-to-annotation links, or parent-child annotation
#'   relationships.
#' - `activeHierarchies`: a character vector naming the hierarchies currently
#'   in use.
#'
#' @section Design notes:
#' Only `MTTKExperiment` is defined at this stage. Dedicated helper classes such
#' as `MTTKHierarchy`, `MTTKFit`, and `MTTKNetwork` are intentionally deferred
#' until their storage contracts are clearer.
#'
#' @name MTTKExperiment-class
#' @rdname MTTKExperiment-class
#' @aliases MTTKExperiment-class
#' @exportClass MTTKExperiment
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
NULL

methods::setClass(
    "MTTKExperiment",
    contains = "TreeSummarizedExperiment"
)

methods::setValidity("MTTKExperiment", .validate_mttk_experiment)

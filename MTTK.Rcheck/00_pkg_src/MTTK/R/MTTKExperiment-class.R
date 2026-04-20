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
#' - `assays` can hold RNA counts, DNA counts, normalized values, offsets,
#'   fitted values, residuals, or other assay matrices.
#' - `rowData` stores feature-level metadata for genes or transcripts.
#' - `colData` stores sample-level metadata.
#' - `rowTree` and `rowLinks` can store optional named hierarchies when a
#'   relationship is naturally tree-structured.
#'
#' @section MTTK metadata contract:
#' The first version keeps MTTK-specific components in `metadata(x)$mttk` rather
#' than adding custom slots. This keeps the class close to existing
#' Bioconductor infrastructure and makes the stored information easy to inspect.
#' The following entries are reserved:
#'
#' - `genomeData`: an `S4Vectors::DataFrame` with one row per genome.
#' - `links`: a named `S4Vectors::SimpleList` of mapping tables, such as
#'   gene-to-genome links, gene-to-annotation links, or parent-child annotation
#'   relationships.
#' - `activeHierarchies`: a character vector naming the hierarchies currently in
#'   use.
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

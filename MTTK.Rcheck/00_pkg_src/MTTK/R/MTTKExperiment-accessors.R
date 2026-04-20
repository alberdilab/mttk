#' Genome-Level Metadata
#'
#' `genomeData()` returns the genome-level metadata stored in an
#' `MTTKExperiment`.
#'
#' @param x An `MTTKExperiment`.
#' @param value For replacement methods, an `S4Vectors::DataFrame`,
#'   `data.frame`, or `NULL`.
#'
#' @return `genomeData()` returns an `S4Vectors::DataFrame`. The replacement
#'   method returns an updated `MTTKExperiment`.
#'
#' @rdname genomeData
#' @export
methods::setMethod("genomeData", "MTTKExperiment", function(x) {
    .mttk_state(x)$genomeData
})

#' @rdname genomeData
#' @export
methods::setReplaceMethod("genomeData", "MTTKExperiment", function(x, value) {
    .replace_mttk_component(x, "genomeData", value)
})

#' Explicit Mapping Tables
#'
#' `links()` returns the named collection of mapping tables stored in an
#' `MTTKExperiment`. These tables can represent relationships such as
#' gene-to-genome, gene-to-annotation, or annotation-to-parent mappings.
#'
#' @param x An `MTTKExperiment`.
#' @param value For replacement methods, a named `S4Vectors::SimpleList`, a
#'   named `list`, or `NULL`.
#'
#' @return `links()` returns an `S4Vectors::SimpleList`. The replacement method
#'   returns an updated `MTTKExperiment`.
#'
#' @rdname links
#' @export
methods::setMethod("links", "MTTKExperiment", function(x) {
    .mttk_state(x)$links
})

#' @rdname links
#' @export
methods::setReplaceMethod("links", "MTTKExperiment", function(x, value) {
    .replace_mttk_component(x, "links", value)
})

#' Active Hierarchies
#'
#' `activeHierarchies()` returns the names of the hierarchies currently in use
#' for an `MTTKExperiment`.
#'
#' @param x An `MTTKExperiment`.
#' @param value For replacement methods, a character vector or `NULL`.
#'
#' @return `activeHierarchies()` returns a character vector. The replacement
#'   method returns an updated `MTTKExperiment`.
#'
#' @rdname activeHierarchies
#' @export
methods::setMethod("activeHierarchies", "MTTKExperiment", function(x) {
    .mttk_state(x)$activeHierarchies
})

#' @rdname activeHierarchies
#' @export
methods::setReplaceMethod(
    "activeHierarchies",
    "MTTKExperiment",
    function(x, value) {
        .replace_mttk_component(x, "activeHierarchies", value)
    }
)

#' RNA Count Assay
#'
#' `rnaCounts()` provides a convenient accessor for the assay named
#' `"rna_counts"`.
#'
#' @param x An `MTTKExperiment`.
#' @param withDimnames Logical; should row and column names be included in the
#'   returned assay.
#' @param value For replacement methods, a matrix-like object or `NULL`.
#'
#' @return `rnaCounts()` returns the `"rna_counts"` assay when present, or
#'   `NULL` otherwise. The replacement method returns an updated
#'   `MTTKExperiment`.
#'
#' @rdname rnaCounts
#' @export
methods::setMethod("rnaCounts", "MTTKExperiment", function(x, withDimnames = TRUE) {
    .named_assay(x, assay_name = "rna_counts", withDimnames = withDimnames)
})

#' @rdname rnaCounts
#' @export
methods::setReplaceMethod("rnaCounts", "MTTKExperiment", function(x, value) {
    .set_named_assay(x, assay_name = "rna_counts", value = value)
})

#' DNA Count Assay
#'
#' `dnaCounts()` provides a convenient accessor for the assay named
#' `"dna_counts"`.
#'
#' @param x An `MTTKExperiment`.
#' @param withDimnames Logical; should row and column names be included in the
#'   returned assay.
#' @param value For replacement methods, a matrix-like object or `NULL`.
#'
#' @return `dnaCounts()` returns the `"dna_counts"` assay when present, or
#'   `NULL` otherwise. The replacement method returns an updated
#'   `MTTKExperiment`.
#'
#' @rdname dnaCounts
#' @export
methods::setMethod("dnaCounts", "MTTKExperiment", function(x, withDimnames = TRUE) {
    .named_assay(x, assay_name = "dna_counts", withDimnames = withDimnames)
})

#' @rdname dnaCounts
#' @export
methods::setReplaceMethod("dnaCounts", "MTTKExperiment", function(x, value) {
    .set_named_assay(x, assay_name = "dna_counts", value = value)
})

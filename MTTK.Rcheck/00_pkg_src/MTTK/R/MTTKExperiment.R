#' Construct an MTTKExperiment
#'
#' `MTTKExperiment()` creates the main container used by MTTK. It wraps
#' [TreeSummarizedExperiment::TreeSummarizedExperiment()] and adds reserved
#' metadata fields for genome-level annotations, explicit mapping tables, and
#' the names of active hierarchies.
#'
#' @param ... Arguments passed to
#'   [TreeSummarizedExperiment::TreeSummarizedExperiment()], typically including
#'   `assays`, `rowData`, `colData`, and `metadata`.
#' @param rowTree,colTree Optional row and column trees passed to
#'   `TreeSummarizedExperiment()`.
#' @param rowNodeLab,colNodeLab Optional node labels for `rowTree` and
#'   `colTree`.
#' @param referenceSeq Optional reference sequence information passed to
#'   `TreeSummarizedExperiment()`.
#' @param genomeData Genome-level metadata as an `S4Vectors::DataFrame` or
#'   `data.frame`. Use `NULL` to create an empty table.
#' @param links A named `S4Vectors::SimpleList` or named `list` of mapping
#'   tables. Each entry must be an `S4Vectors::DataFrame` or `data.frame`.
#' @param activeHierarchies Character vector naming the hierarchies currently in
#'   use.
#' @param metadata Additional experiment-level metadata. Existing
#'   `metadata$mttk` entries are respected when the corresponding dedicated
#'   arguments are left as `NULL`.
#'
#' @return A valid `MTTKExperiment`.
#'
#' @export
MTTKExperiment <- function(
    ...,
    rowTree = NULL,
    colTree = NULL,
    rowNodeLab = NULL,
    colNodeLab = NULL,
    referenceSeq = NULL,
    genomeData = NULL,
    links = NULL,
    activeHierarchies = NULL,
    metadata = list()
) {
    tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
        ...,
        rowTree = rowTree,
        colTree = colTree,
        rowNodeLab = rowNodeLab,
        colNodeLab = colNodeLab,
        referenceSeq = referenceSeq,
        metadata = metadata
    )

    out <- methods::as(tse, "MTTKExperiment")
    out <- .update_mttk_state(
        out,
        genomeData = genomeData,
        links = links,
        activeHierarchies = activeHierarchies
    )

    methods::validObject(out)
    out
}

methods::setMethod("show", "MTTKExperiment", function(object) {
    assay_names <- names(SummarizedExperiment::assays(object, withDimnames = FALSE))
    row_data_names <- names(SummarizedExperiment::rowData(object))
    col_data_names <- names(SummarizedExperiment::colData(object))
    genome_count <- nrow(genomeData(object))
    link_names <- names(links(object))
    hierarchy_names <- activeHierarchies(object)
    row_tree_names <- TreeSummarizedExperiment::rowTreeNames(object)
    col_tree_names <- TreeSummarizedExperiment::colTreeNames(object)

    cat(
        class(object)[1L],
        "with",
        nrow(object),
        "features and",
        ncol(object),
        "samples\n"
    )
    cat("assays(", length(assay_names), "): ", .format_preview(assay_names), "\n", sep = "")
    cat(
        "rowData names(",
        length(row_data_names),
        "): ",
        .format_preview(row_data_names),
        "\n",
        sep = ""
    )
    cat(
        "colData names(",
        length(col_data_names),
        "): ",
        .format_preview(col_data_names),
        "\n",
        sep = ""
    )
    cat("genomeData rows(", genome_count, "): ", genome_count, "\n", sep = "")
    cat("links(", length(link_names), "): ", .format_preview(link_names), "\n", sep = "")
    cat(
        "activeHierarchies(",
        length(hierarchy_names),
        "): ",
        .format_preview(hierarchy_names),
        "\n",
        sep = ""
    )
    cat(
        "rowTreeNames(",
        length(row_tree_names),
        "): ",
        .format_preview(row_tree_names),
        "\n",
        sep = ""
    )
    cat(
        "colTreeNames(",
        length(col_tree_names),
        "): ",
        .format_preview(col_tree_names),
        "\n",
        sep = ""
    )
})

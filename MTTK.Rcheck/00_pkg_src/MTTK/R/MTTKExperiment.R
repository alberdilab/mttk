#' Construct an MTTKExperiment
#'
#' `MTTKExperiment()` creates the main container used by MTTK. It wraps
#' [TreeSummarizedExperiment::TreeSummarizedExperiment()] and adds MTTK-specific
#' support for an explicit gene-level layer, a genome-level companion
#' experiment, explicit mapping tables, and active hierarchy labels.
#'
#' @param ... Arguments passed to
#'   [TreeSummarizedExperiment::TreeSummarizedExperiment()], typically including
#'   gene-level `assays`, `rowData`, `colData`, and `metadata`.
#' @param rowTree,colTree Optional row and column trees passed to
#'   `TreeSummarizedExperiment()`.
#' @param rowNodeLab,colNodeLab Optional node labels for `rowTree` and
#'   `colTree`.
#' @param referenceSeq Optional reference sequence information passed to
#'   `TreeSummarizedExperiment()`.
#' @param genomeExperiment Optional genome-level `SummarizedExperiment` with one
#'   row per genome and one column per sample.
#' @param genomeAssays Optional named list of genome-level assays. This is a
#'   convenient way to provide genome DNA abundances or aggregated genome RNA
#'   assays without building `genomeExperiment` by hand.
#' @param genomeData Genome-level metadata as an `S4Vectors::DataFrame` or
#'   `data.frame`. When genome-level assays are present, this becomes
#'   `rowData(genomeExperiment(x))`.
#' @param links A named `S4Vectors::SimpleList` or named `list` of mapping
#'   tables. Each entry must be an `S4Vectors::DataFrame` or `data.frame`.
#' @param activeHierarchies Character vector naming the hierarchies currently in
#'   use.
#' @param metadata Additional experiment-level metadata. Existing
#'   `metadata$mttk` entries for `links` and `activeHierarchies` are respected
#'   when the corresponding dedicated arguments are left as `NULL`.
#'
#' Gene-level and genome-level assays should use explicit names such as
#' `"rna_gene_counts"`, `"rna_genome_counts"`, and `"dna_genome_counts"`.
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
    genomeExperiment = NULL,
    genomeAssays = NULL,
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

    assay_names <- names(SummarizedExperiment::assays(tse, withDimnames = FALSE))
    invalid_gene_assays <- intersect(
        assay_names,
        c("rna_counts", "dna_counts", "rna_genome_counts", "dna_genome_counts")
    )

    if (length(invalid_gene_assays) > 0L) {
        stop(
            "Gene-level assays must use explicit gene-level names. Invalid assay name(s): ",
            paste(invalid_gene_assays, collapse = ", "),
            ". Use names such as 'rna_gene_counts'.",
            call. = FALSE
        )
    }

    out <- .as_mttk_experiment(tse)
    genome_experiment <- .build_genome_experiment(
        col_data = SummarizedExperiment::colData(out),
        genomeExperiment = genomeExperiment,
        genomeAssays = genomeAssays,
        genomeData = genomeData
    )

    if (!is.null(genome_experiment)) {
        out <- .set_genome_experiment(out, genome_experiment)
    }

    out <- .update_mttk_state(
        out,
        links = links,
        activeHierarchies = activeHierarchies
    )

    methods::validObject(out)
    out
}

methods::setMethod("show", "MTTKExperiment", function(object) {
    gene_assay_names <- names(SummarizedExperiment::assays(object, withDimnames = FALSE))
    row_data_names <- names(SummarizedExperiment::rowData(object))
    col_data_names <- names(SummarizedExperiment::colData(object))
    genome_experiment <- genomeExperiment(object)
    genome_count <- nrow(genomeData(object))
    genome_assay_names <- if (is.null(genome_experiment)) {
        character()
    } else {
        names(SummarizedExperiment::assays(genome_experiment, withDimnames = FALSE))
    }
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
    cat(
        "geneAssays(",
        length(gene_assay_names),
        "): ",
        .format_preview(gene_assay_names),
        "\n",
        sep = ""
    )
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
    cat(
        "genomeAssays(",
        length(genome_assay_names),
        "): ",
        .format_preview(genome_assay_names),
        "\n",
        sep = ""
    )
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

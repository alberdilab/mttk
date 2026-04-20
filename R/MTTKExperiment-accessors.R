#' Gene-Level Experiment
#'
#' `geneExperiment()` returns the gene-level
#' `TreeSummarizedExperiment` view of an `MTTKExperiment`. In the current
#' design, the `MTTKExperiment` itself is the gene-level experiment, and this
#' accessor provides an explicit name for that layer.
#'
#' @param x An `MTTKExperiment`.
#' @param value For replacement methods, a `TreeSummarizedExperiment`,
#'   `SummarizedExperiment`, or `MTTKExperiment`.
#'
#' @return `geneExperiment()` returns a `TreeSummarizedExperiment`. The
#'   replacement method returns an updated `MTTKExperiment`.
#'
#' @rdname geneExperiment
#' @export
methods::setMethod("geneExperiment", "MTTKExperiment", function(x) {
    .as_tree_summarized_experiment(x)
})

#' @rdname geneExperiment
#' @export
methods::setReplaceMethod("geneExperiment", "MTTKExperiment", function(x, value) {
    .set_gene_experiment(x, value)
})

#' Genome-Level Experiment
#'
#' `genomeExperiment()` returns the genome-level
#' `SummarizedExperiment` stored inside an `MTTKExperiment`. This companion
#' experiment contains one row per genome and stores genome-level assays such as
#' genome DNA abundances or aggregated genome RNA counts.
#'
#' @param x An `MTTKExperiment`.
#' @param value For replacement methods, a `SummarizedExperiment` or `NULL`.
#'
#' @return `genomeExperiment()` returns a `SummarizedExperiment` or `NULL`. The
#'   replacement method returns an updated `MTTKExperiment`.
#'
#' @rdname genomeExperiment
#' @export
methods::setMethod("genomeExperiment", "MTTKExperiment", function(x) {
    .get_genome_experiment(x)
})

#' @rdname genomeExperiment
#' @export
methods::setReplaceMethod("genomeExperiment", "MTTKExperiment", function(x, value) {
    .set_genome_experiment(x, value)
})

#' Gene-Level Assays
#'
#' `geneAssays()` returns the assay collection stored in the gene-level layer of
#' an `MTTKExperiment`.
#'
#' @param x An `MTTKExperiment`.
#' @param withDimnames Logical; should row and column names be included in the
#'   returned assays.
#' @param value For replacement methods, a named list of assays, an
#'   `S4Vectors::SimpleList`, or `NULL`.
#'
#' @return `geneAssays()` returns the gene-level assay collection. The
#'   replacement method returns an updated `MTTKExperiment`.
#'
#' @rdname geneAssays
#' @export
methods::setMethod("geneAssays", "MTTKExperiment", function(x, withDimnames = TRUE) {
    SummarizedExperiment::assays(x, withDimnames = withDimnames)
})

#' @rdname geneAssays
#' @export
methods::setReplaceMethod("geneAssays", "MTTKExperiment", function(x, value) {
    value <- .normalize_assay_collection(value, "geneAssays")
    SummarizedExperiment::assays(x, withDimnames = FALSE) <- value
    methods::validObject(x)
    x
})

#' Genome-Level Assays
#'
#' `genomeAssays()` returns the assay collection stored in `genomeExperiment(x)`.
#'
#' @param x An `MTTKExperiment`.
#' @param withDimnames Logical; should row and column names be included in the
#'   returned assays.
#' @param value For replacement methods, a named list of assays, an
#'   `S4Vectors::SimpleList`, or `NULL`.
#'
#' @return `genomeAssays()` returns the genome-level assay collection. The
#'   replacement method returns an updated `MTTKExperiment`.
#'
#' @rdname genomeAssays
#' @export
methods::setMethod("genomeAssays", "MTTKExperiment", function(x, withDimnames = TRUE) {
    genome_experiment <- genomeExperiment(x)

    if (is.null(genome_experiment)) {
        return(.empty_assay_list())
    }

    SummarizedExperiment::assays(genome_experiment, withDimnames = withDimnames)
})

#' @rdname genomeAssays
#' @export
methods::setReplaceMethod("genomeAssays", "MTTKExperiment", function(x, value) {
    value <- .normalize_assay_collection(value, "genomeAssays")
    genome_experiment <- genomeExperiment(x)

    if (is.null(genome_experiment)) {
        genome_experiment <- .build_genome_experiment(
            col_data = SummarizedExperiment::colData(x),
            genomeAssays = value,
            genomeData = genomeData(x)
        )

        if (is.null(genome_experiment)) {
            methods::validObject(x)
            return(x)
        }

        return(.set_genome_experiment(x, genome_experiment))
    }

    SummarizedExperiment::assays(genome_experiment, withDimnames = FALSE) <- value
    .set_genome_experiment(x, genome_experiment)
})

#' Genome-Level Metadata
#'
#' `genomeData()` returns the genome-level metadata stored in
#' `rowData(genomeExperiment(x))`.
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
    genome_experiment <- genomeExperiment(x)

    if (is.null(genome_experiment)) {
        return(S4Vectors::DataFrame())
    }

    SummarizedExperiment::rowData(genome_experiment)
})

#' @rdname genomeData
#' @export
methods::setReplaceMethod("genomeData", "MTTKExperiment", function(x, value) {
    genome_experiment <- genomeExperiment(x)

    if (is.null(value) && !is.null(genome_experiment)) {
        value <- S4Vectors::DataFrame(row.names = rownames(genome_experiment))
    }

    value <- .normalize_genome_data(value)

    if (is.null(genome_experiment)) {
        genome_experiment <- .empty_genome_experiment(
            col_data = SummarizedExperiment::colData(x),
            genome_data = value
        )
    } else {
        if (nrow(genome_experiment) != nrow(value)) {
            stop(
                "'genomeData' must describe the same number of genomes as 'genomeExperiment(x)'.",
                call. = FALSE
            )
        }

        SummarizedExperiment::rowData(genome_experiment) <- value
    }

    .set_genome_experiment(x, genome_experiment)
})

#' Explicit Mapping Tables
#'
#' `links()` returns the named collection of mapping tables stored in an
#' `MTTKExperiment`. These tables can represent relationships such as
#' gene-to-annotation or annotation-to-parent mappings. The core nesting of
#' genes within genomes is defined by `rowData(x)$genome_id`; a stored
#' `gene_to_genome` table is optional and, when present, must mirror that
#' column exactly.
#'
#' The first two columns of each link table are treated as the source and target
#' identifier columns. They should therefore use explicit identifier names such
#' as `gene_id`, `ko_id`, `module_id`, or `pathway_id`. A non-core link table
#' is valid only when its source identifiers match canonical IDs already stored
#' in the object or target identifiers introduced by another link table.
#'
#' When an `MTTKExperiment` is subset by rows, link tables whose source IDs are
#' downstream of the retained genes are pruned so the remaining mappings stay
#' consistent with the subsetted object.
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

#' Gene-Level RNA Counts
#'
#' `rnaGeneCounts()` provides a convenient accessor for the gene-level assay
#' named `"rna_gene_counts"`.
#'
#' @param x An `MTTKExperiment`.
#' @param withDimnames Logical; should row and column names be included in the
#'   returned assay.
#' @param value For replacement methods, a matrix-like object or `NULL`.
#'
#' @return `rnaGeneCounts()` returns the `"rna_gene_counts"` assay when
#'   present, or `NULL` otherwise. The replacement method returns an updated
#'   `MTTKExperiment`.
#'
#' @rdname rnaGeneCounts
#' @export
methods::setMethod("rnaGeneCounts", "MTTKExperiment", function(x, withDimnames = TRUE) {
    .named_assay(x, assay_name = "rna_gene_counts", withDimnames = withDimnames)
})

#' @rdname rnaGeneCounts
#' @export
methods::setReplaceMethod("rnaGeneCounts", "MTTKExperiment", function(x, value) {
    .set_named_assay(x, assay_name = "rna_gene_counts", value = value)
})

#' Genome-Level RNA Counts
#'
#' `rnaGenomeCounts()` provides a convenient accessor for the genome-level assay
#' named `"rna_genome_counts"` stored in `genomeExperiment(x)`.
#'
#' @param x An `MTTKExperiment`.
#' @param withDimnames Logical; should row and column names be included in the
#'   returned assay.
#' @param value For replacement methods, a matrix-like object or `NULL`.
#'
#' @return `rnaGenomeCounts()` returns the `"rna_genome_counts"` genome assay
#'   when present, or `NULL` otherwise. The replacement method returns an
#'   updated `MTTKExperiment`.
#'
#' @rdname rnaGenomeCounts
#' @export
methods::setMethod("rnaGenomeCounts", "MTTKExperiment", function(x, withDimnames = TRUE) {
    .named_genome_assay(x, assay_name = "rna_genome_counts", withDimnames = withDimnames)
})

#' @rdname rnaGenomeCounts
#' @export
methods::setReplaceMethod("rnaGenomeCounts", "MTTKExperiment", function(x, value) {
    .set_named_genome_assay(x, assay_name = "rna_genome_counts", value = value)
})

#' Genome-Level DNA Counts
#'
#' `dnaGenomeCounts()` provides a convenient accessor for the genome-level assay
#' named `"dna_genome_counts"` stored in `genomeExperiment(x)`.
#'
#' @param x An `MTTKExperiment`.
#' @param withDimnames Logical; should row and column names be included in the
#'   returned assay.
#' @param value For replacement methods, a matrix-like object or `NULL`.
#'
#' @return `dnaGenomeCounts()` returns the `"dna_genome_counts"` genome assay
#'   when present, or `NULL` otherwise. The replacement method returns an
#'   updated `MTTKExperiment`.
#'
#' @rdname dnaGenomeCounts
#' @export
methods::setMethod("dnaGenomeCounts", "MTTKExperiment", function(x, withDimnames = TRUE) {
    .named_genome_assay(x, assay_name = "dna_genome_counts", withDimnames = withDimnames)
})

#' @rdname dnaGenomeCounts
#' @export
methods::setReplaceMethod("dnaGenomeCounts", "MTTKExperiment", function(x, value) {
    .set_named_genome_assay(x, assay_name = "dna_genome_counts", value = value)
})

.normalize_analysis_assays <- function(x, assays) {
    assay_names <- names(SummarizedExperiment::assays(x, withDimnames = FALSE))

    if (length(assay_names) == 0L) {
        stop("The object does not contain any assays to analyze.", call. = FALSE)
    }

    if (is.null(assays)) {
        return(assay_names)
    }

    assays <- as.character(assays)
    missing_assays <- setdiff(assays, assay_names)

    if (length(missing_assays) > 0L) {
        stop(
            "Unknown assay name(s): ",
            paste(missing_assays, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    assays
}

.normalize_genome_analysis_assay <- function(x, assay) {
    genome_experiment <- genomeExperiment(x)

    if (is.null(genome_experiment)) {
        stop(
            "A genome-level assay is required in 'genomeExperiment(x)' before this analysis can be run.",
            call. = FALSE
        )
    }

    assay_names <- names(
        SummarizedExperiment::assays(genome_experiment, withDimnames = FALSE)
    )

    if (length(assay_names) == 0L) {
        stop(
            "The genome-level experiment does not contain any assays to analyze.",
            call. = FALSE
        )
    }

    assay <- as.character(assay)

    if (length(assay) != 1L || is.na(assay) || assay == "") {
        stop("'genomeAssay' must be a single non-empty assay name.", call. = FALSE)
    }

    if (!(assay %in% assay_names)) {
        stop("Unknown genome assay name: ", assay, ".", call. = FALSE)
    }

    assay
}

.aligned_genome_assay <- function(x, assay_name, genome_ids, output_rownames) {
    genome_experiment <- genomeExperiment(x)
    assay_mat <- SummarizedExperiment::assay(
        genome_experiment,
        assay_name,
        withDimnames = TRUE
    )
    matched <- match(genome_ids, rownames(assay_mat))

    if (anyNA(matched)) {
        missing_ids <- unique(genome_ids[is.na(matched)])
        stop(
            "All requested genomes must be present in genome assay '",
            assay_name,
            "'. Missing: ",
            paste(missing_ids, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    aligned <- assay_mat[matched, , drop = FALSE]
    rownames(aligned) <- output_rownames
    aligned
}

.resolve_link_path <- function(x, path) {
    path <- as.character(path)

    if (length(path) == 0L || anyNA(path) || any(path == "")) {
        stop("'path' must contain one or more named link tables.", call. = FALSE)
    }

    feature_ids <- rownames(x)
    if (is.null(feature_ids) || anyNA(feature_ids) || any(feature_ids == "")) {
        stop(
            "Row names must be present on 'x' before link-based aggregation can be used.",
            call. = FALSE
        )
    }

    link_list <- links(x)
    available_links <- union(names(link_list), "gene_to_genome")
    missing_links <- setdiff(path, available_links)

    if (length(missing_links) > 0L) {
        stop(
            "Unknown link table(s): ",
            paste(missing_links, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    mapping <- data.frame(
        feature_id = feature_ids,
        current_id = feature_ids,
        feature_index = seq_along(feature_ids),
        stringsAsFactors = FALSE
    )

    target_column <- "group_id"

    for (link_name in path) {
        link_table <- if (identical(link_name, "gene_to_genome")) {
            .gene_to_genome_link_df(x)
        } else {
            as.data.frame(link_list[[link_name]])
        }

        if (ncol(link_table) < 2L) {
            stop(
                "Link table '",
                link_name,
                "' must contain at least two columns.",
                call. = FALSE
            )
        }

        target_column <- names(link_table)[2L]
        link_table <- unique(link_table[, 1:2, drop = FALSE])
        names(link_table) <- c("source_id", "target_id")

        mapping <- merge(
            mapping,
            link_table,
            by.x = "current_id",
            by.y = "source_id",
            all = FALSE,
            sort = FALSE
        )

        if (nrow(mapping) == 0L) {
            stop(
                "The link path did not map any features: ",
                paste(path, collapse = " -> "),
                ".",
                call. = FALSE
            )
        }

        mapping <- mapping[
            order(mapping$feature_index),
            c("feature_id", "feature_index", "target_id"),
            drop = FALSE
        ]
        names(mapping)[3L] <- "current_id"
    }

    mapping <- unique(mapping[, c("feature_id", "feature_index", "current_id"), drop = FALSE])
    mapping <- mapping[order(mapping$feature_index), , drop = FALSE]
    names(mapping)[3L] <- "group_id"

    list(
        mapping = S4Vectors::DataFrame(
            feature_id = mapping$feature_id,
            group_id = mapping$group_id
        ),
        target_column = target_column,
        path = path
    )
}

.aggregate_one_assay <- function(mat, mapping, fun) {
    mat <- as.matrix(mat)
    feature_index <- match(mapping$feature_id, rownames(mat))

    if (anyNA(feature_index)) {
        stop(
            "All mapped features must be present in the assay row names.",
            call. = FALSE
        )
    }

    grouped_mat <- mat[feature_index, , drop = FALSE]
    group_ids <- as.character(mapping$group_id)

    if (fun == "sum") {
        return(rowsum(grouped_mat, group = group_ids, reorder = FALSE))
    }

    sum_mat <- rowsum(grouped_mat, group = group_ids, reorder = FALSE)
    count_mat <- rowsum(
        matrix(1, nrow = nrow(grouped_mat), ncol = ncol(grouped_mat)),
        group = group_ids,
        reorder = FALSE
    )

    sum_mat / count_mat
}

.group_summary <- function(mapping, group_ids) {
    n_features <- vapply(
        group_ids,
        function(id) length(unique(mapping$feature_id[mapping$group_id == id])),
        integer(1)
    )
    n_links <- vapply(
        group_ids,
        function(id) sum(mapping$group_id == id),
        integer(1)
    )

    S4Vectors::DataFrame(
        n_features = as.integer(n_features),
        n_links = as.integer(n_links),
        row.names = group_ids
    )
}

.aggregation_row_data <- function(x, resolved, group_ids) {
    row_data <- .group_summary(resolved$mapping, group_ids)
    target_col <- resolved$target_column
    row_data[[target_col]] <- group_ids
    row_data <- row_data[, c(target_col, "n_features", "n_links")]
    rownames(row_data) <- group_ids

    if (identical(target_col, "genome_id")) {
        genome_data <- genomeData(x)
        matched <- match(group_ids, rownames(genome_data))

        if (all(!is.na(matched))) {
            row_data <- cbind(row_data, genome_data[matched, , drop = FALSE])
            rownames(row_data) <- group_ids
        }
    }

    row_data
}

#' Aggregate Assays Along a Link Path
#'
#' `aggregateByLink()` aggregates one or more gene-level assays in an
#' `MTTKExperiment` along one or more explicit link tables. This makes it
#' possible to summarize RNA measurements to genomes, KOs, modules, pathways,
#' or other mapped feature groups.
#'
#' The current implementation expects each link table in the path to be ordered
#' as a two-column mapping from source IDs to target IDs. Those first two
#' columns should use explicit identifier names such as `gene_id`, `ko_id`, or
#' `module_id`. For example, the packaged example data can be aggregated from
#' genes to modules with
#' `path = c("gene_to_ko", "ko_to_module")`. The special path
#' `"gene_to_genome"` is resolved from the authoritative
#' `rowData(x)$genome_id` mapping, so it does not require a stored
#' `links(x)[["gene_to_genome"]]` table.
#'
#' @param x An `MTTKExperiment`.
#' @param path Character vector naming one or more link tables stored in
#'   `links(x)`.
#' @param assays Character vector of assay names to aggregate. The default uses
#'   all assays in `x`.
#' @param fun Aggregation function. Supported values are `"sum"` and `"mean"`.
#'
#' @return A `SummarizedExperiment` with aggregated assays, copied `colData`,
#'   aggregated `rowData`, and provenance stored in `metadata()`.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#'
#' aggregateByLink(x, path = "gene_to_genome")
#' aggregateByLink(x, path = c("gene_to_ko", "ko_to_module"))
#'
#' @export
aggregateByLink <- function(x, path, assays = NULL, fun = c("sum", "mean")) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    fun <- match.arg(fun)
    assay_names <- .normalize_analysis_assays(x, assays)
    resolved <- .resolve_link_path(x, path)

    aggregated_assays <- lapply(assay_names, function(assay_name) {
        .aggregate_one_assay(
            SummarizedExperiment::assay(x, assay_name, withDimnames = TRUE),
            mapping = resolved$mapping,
            fun = fun
        )
    })
    names(aggregated_assays) <- assay_names

    group_ids <- rownames(aggregated_assays[[1L]])
    row_data <- .aggregation_row_data(x, resolved, group_ids)

    SummarizedExperiment::SummarizedExperiment(
        assays = aggregated_assays,
        rowData = row_data,
        colData = SummarizedExperiment::colData(x),
        metadata = c(
            S4Vectors::metadata(x),
            list(
                mttk_aggregation = list(
                    path = resolved$path,
                    assays = assay_names,
                    fun = fun,
                    target_column = resolved$target_column
                )
            )
        )
    )
}

#' Aggregate Assays to the Genome Level
#'
#' `aggregateToGenome()` is a convenience wrapper around `aggregateByLink()` for
#' the common case of aggregating gene-level measurements to genomes through the
#' core `rowData(x)$genome_id` mapping. When an input assay name contains `"_gene_"`,
#' the returned assay name is rewritten with `"_genome_"` so the level of the
#' aggregated data remains explicit.
#'
#' @param x An `MTTKExperiment`.
#' @param assays Character vector of assay names to aggregate. The default uses
#'   all assays in `x`.
#' @param fun Aggregation function. Supported values are `"sum"` and `"mean"`.
#'
#' @return A `SummarizedExperiment` with one row per genome.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' aggregateToGenome(x)
#'
#' @export
aggregateToGenome <- function(x, assays = NULL, fun = c("sum", "mean")) {
    aggregated <- aggregateByLink(
        x = x,
        path = "gene_to_genome",
        assays = assays,
        fun = match.arg(fun)
    )

    assay_list <- SummarizedExperiment::assays(aggregated, withDimnames = FALSE)
    renamed_assays <- stats::setNames(
        as.list(assay_list),
        vapply(names(assay_list), .to_genome_assay_name, character(1))
    )
    SummarizedExperiment::assays(aggregated, withDimnames = FALSE) <- renamed_assays

    S4Vectors::metadata(aggregated)$mttk_aggregation$output_assays <- names(
        SummarizedExperiment::assays(aggregated, withDimnames = FALSE)
    )

    aggregated
}

#' Summarize Gene or Genome Activity
#'
#' `summarizeActivity()` computes a simple RNA-over-DNA activity summary from a
#' gene-level RNA assay and a genome-level DNA assay. The summary can be
#' computed for each gene by pairing every gene with the DNA abundance of its
#' parent genome, or at the genome level after aggregating gene-level RNA to
#' genomes.
#'
#' This function is intended as a first, transparent analysis helper for the
#' packaged example data. It does not replace a formal statistical model, but it
#' provides an immediate way to compare relative transcriptional activity across
#' samples.
#'
#' @param x An `MTTKExperiment`.
#' @param by Level at which activity should be summarized. `"gene"` uses the
#'   gene-level RNA assay together with the genome-level DNA assay, while
#'   `"genome"` first aggregates the gene-level RNA assay with
#'   `aggregateToGenome()`.
#' @param numeratorAssay Gene-level assay used as the numerator.
#' @param genomeAssay Genome-level assay used as the denominator.
#' @param pseudocount Numeric pseudocount added to both assays before division.
#' @param transform Activity scale. `"log2_ratio"` returns
#'   `log2((numerator + pseudocount) / (denominator + pseudocount))`, while
#'   `"ratio"` returns the raw ratio.
#'
#' @return A `SummarizedExperiment` containing a single assay named
#'   `"log2_activity"` or `"activity"`.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#'
#' summarizeActivity(x, by = "gene")
#' summarizeActivity(x, by = "genome")
#'
#' @export
summarizeActivity <- function(
    x,
    by = c("gene", "genome"),
    numeratorAssay = "rna_gene_counts",
    genomeAssay = "dna_genome_counts",
    pseudocount = 1,
    transform = c("log2_ratio", "ratio")
) {
    if (!methods::is(x, "MTTKExperiment")) {
        stop("'x' must be an MTTKExperiment.", call. = FALSE)
    }

    by <- match.arg(by)
    transform <- match.arg(transform)
    numerator_assay <- .normalize_analysis_assays(x, numeratorAssay)
    genome_assay <- .normalize_genome_analysis_assay(x, genomeAssay)

    if (!is.numeric(pseudocount) || length(pseudocount) != 1L || is.na(pseudocount) ||
        pseudocount < 0) {
        stop("'pseudocount' must be a single non-negative numeric value.", call. = FALSE)
    }

    if (by == "gene") {
        numerator <- SummarizedExperiment::assay(
            x,
            numerator_assay,
            withDimnames = TRUE
        )
        gene_to_genome <- .gene_to_genome_map(x)
        denominator <- .aligned_genome_assay(
            x = x,
            assay_name = genome_assay,
            genome_ids = unname(gene_to_genome),
            output_rownames = names(gene_to_genome)
        )
        row_data <- SummarizedExperiment::rowData(x)
        col_data <- SummarizedExperiment::colData(x)
    } else {
        aggregated <- aggregateToGenome(
            x,
            assays = numerator_assay,
            fun = "sum"
        )
        numerator <- SummarizedExperiment::assay(
            aggregated,
            .to_genome_assay_name(numerator_assay),
            withDimnames = TRUE
        )
        denominator <- .aligned_genome_assay(
            x = x,
            assay_name = genome_assay,
            genome_ids = rownames(aggregated),
            output_rownames = rownames(aggregated)
        )
        row_data <- SummarizedExperiment::rowData(aggregated)
        col_data <- SummarizedExperiment::colData(aggregated)
    }

    ratio <- (as.matrix(numerator) + pseudocount) / (as.matrix(denominator) + pseudocount)
    assay_name <- if (transform == "log2_ratio") "log2_activity" else "activity"
    activity <- if (transform == "log2_ratio") log2(ratio) else ratio

    SummarizedExperiment::SummarizedExperiment(
        assays = stats::setNames(list(activity), assay_name),
        rowData = row_data,
        colData = col_data,
        metadata = c(
            S4Vectors::metadata(x),
            list(
                mttk_activity = list(
                    by = by,
                    numeratorAssay = numerator_assay,
                    genomeAssay = genome_assay,
                    pseudocount = pseudocount,
                    transform = transform
                )
            )
        )
    )
}

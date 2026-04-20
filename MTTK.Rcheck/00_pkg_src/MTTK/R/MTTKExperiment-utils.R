.mttk_genome_alt_name <- "genomes"

.empty_mttk_state <- function() {
    list(
        links = S4Vectors::SimpleList(),
        activeHierarchies = character()
    )
}

.normalize_genome_data <- function(x) {
    if (is.null(x)) {
        return(S4Vectors::DataFrame())
    }

    if (is.data.frame(x) && !methods::is(x, "DataFrame")) {
        x <- S4Vectors::DataFrame(x, check.names = FALSE)
    }

    if (!methods::is(x, "DataFrame")) {
        stop(
            "'genomeData' must be an S4Vectors::DataFrame, a data.frame, or NULL.",
            call. = FALSE
        )
    }

    genome_ids <- rownames(x)
    if (nrow(x) > 0L &&
        (is.null(genome_ids) || anyNA(genome_ids) || any(genome_ids == ""))) {
        stop(
            "'genomeData' must have non-empty row names when genomes are present.",
            call. = FALSE
        )
    }

    if (!is.null(genome_ids) && anyDuplicated(genome_ids)) {
        stop(
            "Row names of 'genomeData' must be unique when present.",
            call. = FALSE
        )
    }

    x
}

.normalize_single_link <- function(x) {
    if (is.data.frame(x) && !methods::is(x, "DataFrame")) {
        x <- S4Vectors::DataFrame(x, check.names = FALSE)
    }

    if (!methods::is(x, "DataFrame")) {
        stop(
            "Each entry in 'links' must be an S4Vectors::DataFrame or data.frame.",
            call. = FALSE
        )
    }

    x
}

.normalize_links <- function(x) {
    if (is.null(x)) {
        return(S4Vectors::SimpleList())
    }

    if (methods::is(x, "SimpleList")) {
        x <- as.list(x)
    } else if (!is.list(x)) {
        stop(
            "'links' must be an S4Vectors::SimpleList, a named list, or NULL.",
            call. = FALSE
        )
    }

    if (length(x) == 0L) {
        return(S4Vectors::SimpleList())
    }

    link_names <- names(x)
    if (is.null(link_names) || anyNA(link_names) || any(link_names == "")) {
        stop("'links' must be a named collection of mapping tables.", call. = FALSE)
    }

    x <- lapply(x, .normalize_single_link)
    methods::as(x, "SimpleList")
}

.normalize_active_hierarchies <- function(x) {
    if (is.null(x)) {
        return(character())
    }

    x <- as.character(x)
    x <- unique(x)

    if (anyNA(x) || any(x == "")) {
        stop(
            "'activeHierarchies' must contain non-missing, non-empty names.",
            call. = FALSE
        )
    }

    x
}

.normalize_assay_collection <- function(x, arg) {
    if (is.null(x)) {
        return(list())
    }

    if (methods::is(x, "SimpleList")) {
        x <- as.list(x)
    }

    if (!is.list(x)) {
        stop("'", arg, "' must be a named list of assays or NULL.", call. = FALSE)
    }

    if (length(x) == 0L) {
        return(list())
    }

    assay_names <- names(x)
    if (is.null(assay_names) || anyNA(assay_names) || any(assay_names == "")) {
        stop("'", arg, "' must be a named list of assays.", call. = FALSE)
    }

    x
}

.normalize_gene_experiment <- function(x) {
    if (methods::is(x, "MTTKExperiment")) {
        return(.as_tree_summarized_experiment(x))
    }

    if (methods::is(x, "TreeSummarizedExperiment")) {
        return(x)
    }

    if (!methods::is(x, "SummarizedExperiment")) {
        stop(
            "'geneExperiment' must be a TreeSummarizedExperiment, SummarizedExperiment, or MTTKExperiment.",
            call. = FALSE
        )
    }

    TreeSummarizedExperiment::TreeSummarizedExperiment(
        assays = SummarizedExperiment::assays(x, withDimnames = FALSE),
        rowData = SummarizedExperiment::rowData(x),
        colData = SummarizedExperiment::colData(x),
        metadata = S4Vectors::metadata(x)
    )
}

.default_genome_data <- function(genome_assays) {
    if (length(genome_assays) == 0L) {
        return(S4Vectors::DataFrame())
    }

    genome_ids <- rownames(genome_assays[[1L]])
    if (is.null(genome_ids) || anyNA(genome_ids) || any(genome_ids == "")) {
        stop(
            "Genome-level assays must have non-empty row names when 'genomeData' is not supplied.",
            call. = FALSE
        )
    }

    if (anyDuplicated(genome_ids)) {
        stop(
            "Genome-level assay row names must be unique when 'genomeData' is not supplied.",
            call. = FALSE
        )
    }

    S4Vectors::DataFrame(row.names = genome_ids)
}

.empty_genome_experiment <- function(col_data, genome_data = S4Vectors::DataFrame()) {
    SummarizedExperiment::SummarizedExperiment(
        rowData = genome_data,
        colData = col_data
    )
}

.as_tree_summarized_experiment <- function(x) {
    if (methods::is(x, "TreeSummarizedExperiment") &&
        !methods::is(x, "MTTKExperiment")) {
        return(x)
    }

    out <- methods::new("TreeSummarizedExperiment")

    for (slot_name in methods::slotNames("TreeSummarizedExperiment")) {
        methods::slot(out, slot_name, check = FALSE) <-
            methods::slot(x, slot_name)
    }

    out
}

.as_mttk_experiment <- function(x) {
    out <- methods::new("MTTKExperiment")

    for (slot_name in methods::slotNames("TreeSummarizedExperiment")) {
        methods::slot(out, slot_name, check = FALSE) <-
            methods::slot(x, slot_name)
    }

    out
}

.normalize_genome_experiment <- function(x, col_data) {
    if (is.null(x)) {
        return(NULL)
    }

    if (!methods::is(x, "SummarizedExperiment")) {
        stop(
            "'genomeExperiment' must be a SummarizedExperiment or NULL.",
            call. = FALSE
        )
    }

    if (ncol(x) != nrow(col_data)) {
        stop(
            "'genomeExperiment' must have one column per sample in 'x'.",
            call. = FALSE
        )
    }

    sample_ids <- rownames(col_data)
    if (is.null(colnames(x))) {
        colnames(x) <- sample_ids
    } else if (!identical(colnames(x), sample_ids)) {
        stop(
            "Column names of 'genomeExperiment' must match the sample names of 'x'.",
            call. = FALSE
        )
    }

    genome_ids <- rownames(x)
    if (nrow(x) > 0L &&
        (is.null(genome_ids) || anyNA(genome_ids) || any(genome_ids == ""))) {
        stop(
            "'genomeExperiment' must have non-empty row names when genomes are present.",
            call. = FALSE
        )
    }

    if (!is.null(genome_ids) && anyDuplicated(genome_ids)) {
        stop(
            "Row names of 'genomeExperiment' must be unique when present.",
            call. = FALSE
        )
    }

    x
}

.build_genome_experiment <- function(
    col_data,
    genomeExperiment = NULL,
    genomeAssays = NULL,
    genomeData = NULL
) {
    if (!is.null(genomeExperiment) && !is.null(genomeAssays)) {
        stop(
            "Use either 'genomeExperiment' or 'genomeAssays', not both.",
            call. = FALSE
        )
    }

    normalized_genome_data <- if (!is.null(genomeData)) {
        .normalize_genome_data(genomeData)
    } else {
        NULL
    }

    if (!is.null(genomeExperiment)) {
        genomeExperiment <- .normalize_genome_experiment(genomeExperiment, col_data)

        if (!is.null(normalized_genome_data)) {
            if (nrow(genomeExperiment) != nrow(normalized_genome_data)) {
                stop(
                    "'genomeData' and 'genomeExperiment' must describe the same number of genomes.",
                    call. = FALSE
                )
            }

            SummarizedExperiment::rowData(genomeExperiment) <- normalized_genome_data
        }

        return(.normalize_genome_experiment(genomeExperiment, col_data))
    }

    normalized_genome_assays <- .normalize_assay_collection(genomeAssays, "genomeAssays")

    if (length(normalized_genome_assays) == 0L && is.null(normalized_genome_data)) {
        return(NULL)
    }

    genomeExperiment <- SummarizedExperiment::SummarizedExperiment(
        assays = normalized_genome_assays,
        rowData = if (is.null(normalized_genome_data)) {
            .default_genome_data(normalized_genome_assays)
        } else {
            normalized_genome_data
        },
        colData = col_data
    )

    .normalize_genome_experiment(genomeExperiment, col_data)
}

.get_genome_experiment <- function(x) {
    parent <- .as_tree_summarized_experiment(x)
    alt_names <- SingleCellExperiment::altExpNames(parent)

    if (!(.mttk_genome_alt_name %in% alt_names)) {
        return(NULL)
    }

    SingleCellExperiment::altExp(parent, .mttk_genome_alt_name)
}

.set_genome_experiment <- function(x, value) {
    parent <- .as_tree_summarized_experiment(x)

    if (is.null(value)) {
        alt_names <- SingleCellExperiment::altExpNames(parent)

        if (.mttk_genome_alt_name %in% alt_names) {
            keep <- setdiff(alt_names, .mttk_genome_alt_name)
            SingleCellExperiment::altExps(parent) <-
                SingleCellExperiment::altExps(parent)[keep]
        }

        out <- .as_mttk_experiment(parent)
        methods::validObject(out)
        return(out)
    }

    value <- .normalize_genome_experiment(
        value,
        SummarizedExperiment::colData(x)
    )
    SingleCellExperiment::altExp(parent, .mttk_genome_alt_name) <- value
    out <- .as_mttk_experiment(parent)
    methods::validObject(out)
    out
}

.set_gene_experiment <- function(x, value) {
    value <- .normalize_gene_experiment(value)
    out <- .as_mttk_experiment(value)
    out <- .set_mttk_state(out, .mttk_state(x))

    current_genome_experiment <- .get_genome_experiment(x)
    if (!is.null(current_genome_experiment)) {
        out <- .set_genome_experiment(out, current_genome_experiment)
    } else {
        methods::validObject(out)
    }

    out
}

.mttk_state <- function(x) {
    defaults <- .empty_mttk_state()
    mttk <- S4Vectors::metadata(x)[["mttk"]]

    if (is.null(mttk)) {
        return(defaults)
    }

    if (!is.list(mttk)) {
        stop("'metadata(x)$mttk' must be a list.", call. = FALSE)
    }

    defaults$links <- .normalize_links(mttk[["links"]])
    defaults$activeHierarchies <- .normalize_active_hierarchies(
        mttk[["activeHierarchies"]]
    )

    defaults
}

.set_mttk_state <- function(x, state) {
    metadata_list <- S4Vectors::metadata(x)
    metadata_list[["mttk"]] <- state
    S4Vectors::metadata(x) <- metadata_list
    x
}

.replace_mttk_component <- function(x, name, value) {
    state <- .mttk_state(x)
    state[[name]] <- switch(
        name,
        links = .normalize_links(value),
        activeHierarchies = .normalize_active_hierarchies(value),
        stop("Unknown MTTK component: ", name, call. = FALSE)
    )

    x <- .set_mttk_state(x, state)
    methods::validObject(x)
    x
}

.update_mttk_state <- function(
    x,
    links = NULL,
    activeHierarchies = NULL
) {
    state <- .mttk_state(x)

    if (!is.null(links)) {
        state$links <- .normalize_links(links)
    }

    if (!is.null(activeHierarchies)) {
        state$activeHierarchies <- .normalize_active_hierarchies(activeHierarchies)
    }

    x <- .set_mttk_state(x, state)
    methods::validObject(x)
    x
}

.format_preview <- function(x, max_n = 4L) {
    if (length(x) == 0L) {
        return("none")
    }

    preview <- utils::head(x, n = max_n)
    label <- paste(preview, collapse = ", ")

    if (length(x) > max_n) {
        label <- paste0(label, ", ...")
    }

    label
}

.named_assay <- function(x, assay_name, withDimnames = TRUE) {
    assay_names <- names(SummarizedExperiment::assays(x, withDimnames = FALSE))

    if (!(assay_name %in% assay_names)) {
        return(NULL)
    }

    SummarizedExperiment::assay(
        x,
        i = assay_name,
        withDimnames = withDimnames
    )
}

.set_named_assay <- function(x, assay_name, value) {
    assay_list <- SummarizedExperiment::assays(x, withDimnames = FALSE)
    assay_list[[assay_name]] <- value
    SummarizedExperiment::assays(x, withDimnames = FALSE) <- assay_list
    methods::validObject(x)
    x
}

.empty_assay_list <- function() {
    methods::as(list(), "SimpleList")
}

.named_genome_assay <- function(x, assay_name, withDimnames = TRUE) {
    genome_experiment <- .get_genome_experiment(x)

    if (is.null(genome_experiment)) {
        return(NULL)
    }

    assay_names <- names(
        SummarizedExperiment::assays(genome_experiment, withDimnames = FALSE)
    )

    if (!(assay_name %in% assay_names)) {
        return(NULL)
    }

    SummarizedExperiment::assay(
        genome_experiment,
        i = assay_name,
        withDimnames = withDimnames
    )
}

.set_named_genome_assay <- function(x, assay_name, value) {
    genome_experiment <- .get_genome_experiment(x)

    if (is.null(genome_experiment)) {
        genome_data <- genomeData(x)

        if (nrow(genome_data) == 0L) {
            genome_ids <- rownames(value)

            if (is.null(genome_ids) || anyNA(genome_ids) || any(genome_ids == "")) {
                stop(
                    "Set 'genomeData' or provide row names on the genome-level assay before adding a genome-level assay.",
                    call. = FALSE
                )
            }

            genome_data <- S4Vectors::DataFrame(row.names = genome_ids)
        }

        genome_experiment <- .build_genome_experiment(
            col_data = SummarizedExperiment::colData(x),
            genomeAssays = stats::setNames(list(value), assay_name),
            genomeData = genome_data
        )

        return(.set_genome_experiment(x, genome_experiment))
    }

    assay_list <- SummarizedExperiment::assays(genome_experiment, withDimnames = FALSE)
    assay_list[[assay_name]] <- value
    SummarizedExperiment::assays(genome_experiment, withDimnames = FALSE) <- assay_list
    x <- .set_genome_experiment(x, genome_experiment)
    methods::validObject(x)
    x
}

.gene_to_genome_map <- function(x) {
    feature_ids <- rownames(x)

    if (is.null(feature_ids) || anyNA(feature_ids) || any(feature_ids == "")) {
        stop(
            "Gene-level row names must be present before gene-to-genome mappings can be used.",
            call. = FALSE
        )
    }

    row_data <- SummarizedExperiment::rowData(x)
    if ("genome_id" %in% names(row_data)) {
        genome_ids <- as.character(row_data$genome_id)

        if (anyNA(genome_ids) || any(genome_ids == "")) {
            stop(
                "'rowData(x)$genome_id' must contain non-empty genome identifiers.",
                call. = FALSE
            )
        }

        stats::setNames(genome_ids, feature_ids)
    } else {
        link_list <- links(x)

        if (!("gene_to_genome" %in% names(link_list))) {
            stop(
                "A gene-to-genome mapping is required in 'rowData(x)$genome_id' or links(x)[['gene_to_genome']].",
                call. = FALSE
            )
        }

        link_table <- as.data.frame(link_list[["gene_to_genome"]])

        if (ncol(link_table) < 2L) {
            stop(
                "The 'gene_to_genome' link table must contain at least two columns.",
                call. = FALSE
            )
        }

        link_table <- unique(link_table[, 1:2, drop = FALSE])
        names(link_table) <- c("gene_id", "genome_id")

        if (anyDuplicated(link_table$gene_id)) {
            stop(
                "The 'gene_to_genome' link table must map each gene to exactly one genome.",
                call. = FALSE
            )
        }

        matched <- match(feature_ids, link_table$gene_id)

        if (anyNA(matched)) {
            stop(
                "Every gene must map to a genome before genome-aware analyses can be run.",
                call. = FALSE
            )
        }

        stats::setNames(as.character(link_table$genome_id[matched]), feature_ids)
    }
}

.to_genome_assay_name <- function(assay_name) {
    if (grepl("_gene_", assay_name, fixed = TRUE)) {
        return(sub("_gene_", "_genome_", assay_name, fixed = TRUE))
    }

    assay_name
}

.validate_mttk_experiment <- function(object) {
    problems <- character()
    mttk <- S4Vectors::metadata(object)[["mttk"]]

    if (!is.null(mttk) && !is.list(mttk)) {
        return("'metadata(x)$mttk' must be a list.")
    }

    assay_names <- names(SummarizedExperiment::assays(object, withDimnames = FALSE))
    invalid_gene_assays <- intersect(
        assay_names,
        c("rna_counts", "dna_counts", "rna_genome_counts", "dna_genome_counts")
    )

    if (length(invalid_gene_assays) > 0L) {
        problems <- c(
            problems,
            paste0(
                "Gene-level assays must use explicit gene-level names. Invalid assay name(s): ",
                paste(invalid_gene_assays, collapse = ", "),
                ". Use names such as 'rna_gene_counts'."
            )
        )
    }

    genome_experiment <- .get_genome_experiment(object)
    if (!is.null(genome_experiment)) {
        tryCatch(
            .normalize_genome_experiment(
                genome_experiment,
                SummarizedExperiment::colData(object)
            ),
            error = function(e) {
                problems <<- c(problems, conditionMessage(e))
            }
        )

        genome_assay_names <- names(
            SummarizedExperiment::assays(genome_experiment, withDimnames = FALSE)
        )
        invalid_genome_assays <- intersect(
            genome_assay_names,
            c("rna_counts", "dna_counts", "rna_gene_counts")
        )

        if (length(invalid_genome_assays) > 0L) {
            problems <- c(
                problems,
                paste0(
                    "Genome-level assays must use explicit genome-level names. Invalid assay name(s): ",
                    paste(invalid_genome_assays, collapse = ", "),
                    ". Use names such as 'rna_genome_counts' or 'dna_genome_counts'."
                )
            )
        }
    }

    link_tables <- if (is.null(mttk)) NULL else mttk[["links"]]
    if (!is.null(link_tables)) {
        if (!methods::is(link_tables, "SimpleList")) {
            problems <- c(
                problems,
                "'metadata(x)$mttk$links' must be an S4Vectors::SimpleList."
            )
        } else {
            link_names <- names(link_tables)

            if (length(link_tables) > 0L &&
                (is.null(link_names) ||
                 anyNA(link_names) ||
                 any(link_names == "") ||
                 anyDuplicated(link_names))) {
                problems <- c(
                    problems,
                    "'links' must have unique, non-empty names."
                )
            }

            if (length(link_tables) > 0L) {
                are_data_frames <- vapply(
                    as.list(link_tables),
                    function(x) methods::is(x, "DataFrame"),
                    logical(1)
                )

                if (!all(are_data_frames)) {
                    problems <- c(
                        problems,
                        "Each entry in 'links' must be an S4Vectors::DataFrame."
                    )
                }
            }
        }
    }

    hierarchy_names <- if (is.null(mttk)) NULL else mttk[["activeHierarchies"]]
    if (!is.null(hierarchy_names)) {
        if (!is.character(hierarchy_names)) {
            problems <- c(
                problems,
                "'metadata(x)$mttk$activeHierarchies' must be a character vector."
            )
        } else if (anyNA(hierarchy_names) ||
                   any(hierarchy_names == "") ||
                   anyDuplicated(hierarchy_names)) {
            problems <- c(
                problems,
                "'activeHierarchies' must contain unique, non-empty names."
            )
        }
    }

    if (!is.null(genome_experiment)) {
        genome_ids <- rownames(genome_experiment)
        gene_to_genome <- tryCatch(
            .gene_to_genome_map(object),
            error = function(e) NULL
        )

        if (!is.null(gene_to_genome)) {
            missing_genomes <- setdiff(unique(unname(gene_to_genome)), genome_ids)

            if (length(missing_genomes) > 0L) {
                problems <- c(
                    problems,
                    paste0(
                        "All mapped genomes must be present in 'genomeExperiment(x)'. Missing: ",
                        paste(missing_genomes, collapse = ", "),
                        "."
                    )
                )
            }
        }
    }

    if (length(problems) == 0L) {
        TRUE
    } else {
        problems
    }
}

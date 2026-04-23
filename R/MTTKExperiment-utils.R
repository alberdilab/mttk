.mttk_genome_alt_name <- "genomes"
.mttk_genome_tree_metadata_name <- "mttk_genome_tree"

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

.normalize_feature_ids <- function(x) {
    feature_ids <- rownames(x)

    if (nrow(x) > 0L &&
        (is.null(feature_ids) || anyNA(feature_ids) || any(feature_ids == ""))) {
        stop(
            "Gene-level row names must be present, non-missing, and non-empty.",
            call. = FALSE
        )
    }

    if (!is.null(feature_ids) && anyDuplicated(feature_ids)) {
        stop("Gene-level row names must be unique.", call. = FALSE)
    }

    feature_ids
}

.normalize_sample_ids <- function(x) {
    sample_ids <- colnames(x)

    if (ncol(x) > 0L &&
        (is.null(sample_ids) || anyNA(sample_ids) || any(sample_ids == ""))) {
        stop(
            "Sample column names must be present, non-missing, and non-empty.",
            call. = FALSE
        )
    }

    if (!is.null(sample_ids) && anyDuplicated(sample_ids)) {
        stop("Sample column names must be unique.", call. = FALSE)
    }

    sample_ids
}

.validate_id_mirror <- function(ids, data, column, mirror_label, canonical_label) {
    if (!(column %in% names(data))) {
        return(NULL)
    }

    mirrored_ids <- as.character(data[[column]])

    if (length(mirrored_ids) != length(ids) ||
        anyNA(mirrored_ids) ||
        any(mirrored_ids == "")) {
        return(
            paste0(
                "'",
                mirror_label,
                "' must contain one non-empty identifier for every entry when present."
            )
        )
    }

    if (!identical(mirrored_ids, as.character(ids))) {
        return(
            paste0(
                "'",
                mirror_label,
                "' must match '",
                canonical_label,
                "' exactly when present."
            )
        )
    }

    NULL
}

.core_gene_to_genome_link <- function(x) {
    feature_ids <- .normalize_feature_ids(x)
    row_data <- SummarizedExperiment::rowData(x)

    if (!("genome_id" %in% names(row_data))) {
        stop(
            "rowData(x) must contain a 'genome_id' column for the core gene-to-genome nesting.",
            call. = FALSE
        )
    }

    genome_ids <- as.character(row_data$genome_id)

    if (length(genome_ids) != nrow(x) ||
        anyNA(genome_ids) ||
        any(genome_ids == "")) {
        stop(
            "'rowData(x)$genome_id' must contain one non-empty genome identifier for every gene.",
            call. = FALSE
        )
    }

    S4Vectors::DataFrame(
        gene_id = feature_ids,
        genome_id = genome_ids,
        row.names = feature_ids
    )
}

.gene_to_genome_link_df <- function(x) {
    as.data.frame(.core_gene_to_genome_link(x))
}

.link_columns <- function(link_table) {
    column_names <- names(link_table)

    if (length(column_names) < 2L) {
        return(NULL)
    }

    stats::setNames(column_names[1:2], c("source", "target"))
}

.link_spec <- function(link_table, link_name) {
    link_cols <- .link_columns(link_table)

    if (is.null(link_cols)) {
        return(list(
            problems = paste0(
                "Link table '",
                link_name,
                "' must contain at least two columns."
            )
        ))
    }

    source_col <- link_cols[["source"]]
    target_col <- link_cols[["target"]]
    problems <- character()

    if (anyNA(c(source_col, target_col)) || any(c(source_col, target_col) == "")) {
        problems <- c(
            problems,
            paste0(
                "The first two columns of link table '",
                link_name,
                "' must be named."
            )
        )
    }

    if (!all(grepl("_id$", c(source_col, target_col)))) {
        problems <- c(
            problems,
            paste0(
                "The first two columns of link table '",
                link_name,
                "' must be identifier columns ending in '_id'."
            )
        )
    }

    source_ids <- as.character(link_table[[source_col]])
    target_ids <- as.character(link_table[[target_col]])

    if (length(source_ids) != nrow(link_table) ||
        anyNA(source_ids) ||
        any(source_ids == "")) {
        problems <- c(
            problems,
            paste0(
                "Column '",
                source_col,
                "' in link table '",
                link_name,
                "' must contain one non-empty identifier for every row."
            )
        )
    }

    if (length(target_ids) != nrow(link_table) ||
        anyNA(target_ids) ||
        any(target_ids == "")) {
        problems <- c(
            problems,
            paste0(
                "Column '",
                target_col,
                "' in link table '",
                link_name,
                "' must contain one non-empty identifier for every row."
            )
        )
    }

    if (length(problems) > 0L) {
        return(list(problems = problems))
    }

    list(
        problems = character(),
        source_col = source_col,
        target_col = target_col,
        source_ids = source_ids,
        target_ids = target_ids
    )
}

.validate_link_referential_integrity <- function(link_tables, known_ids) {
    if (is.null(link_tables) ||
        !methods::is(link_tables, "SimpleList") ||
        length(link_tables) == 0L) {
        return(character())
    }

    link_specs <- list()
    problems <- character()

    for (link_name in names(link_tables)) {
        if (identical(link_name, "gene_to_genome")) {
            next
        }

        spec <- .link_spec(link_tables[[link_name]], link_name)
        problems <- c(problems, spec$problems)

        if (length(spec$problems) == 0L) {
            link_specs[[link_name]] <- spec
        }
    }

    if (length(link_specs) == 0L) {
        return(problems)
    }

    remaining <- link_specs
    progressed <- TRUE

    while (progressed && length(remaining) > 0L) {
        progressed <- FALSE
        next_remaining <- list()

        for (link_name in names(remaining)) {
            spec <- remaining[[link_name]]
            source_pool <- known_ids[[spec$source_col]]

            if (is.null(source_pool)) {
                next_remaining[[link_name]] <- spec
                next
            }

            source_pool <- unique(as.character(source_pool))
            source_pool <- source_pool[!is.na(source_pool) & source_pool != ""]
            source_ids <- unique(spec$source_ids)
            missing_ids <- setdiff(source_ids, source_pool)

            if (length(missing_ids) > 0L) {
                problems <- c(
                    problems,
                    paste0(
                        "Link table '",
                        link_name,
                        "' contains ",
                        spec$source_col,
                        " values that are not defined in the object or upstream links: ",
                        .format_preview(missing_ids),
                        "."
                    )
                )
            }

            matched_rows <- spec$source_ids %in% source_pool
            reached_targets <- unique(spec$target_ids[matched_rows])
            reached_targets <- reached_targets[!is.na(reached_targets) & reached_targets != ""]

            old_targets <- known_ids[[spec$target_col]]
            new_targets <- if (is.null(old_targets)) {
                reached_targets
            } else {
                unique(c(as.character(old_targets), reached_targets))
            }

            if (is.null(old_targets) || !identical(old_targets, new_targets)) {
                known_ids[[spec$target_col]] <- new_targets
                progressed <- TRUE
            }
        }

        remaining <- next_remaining
    }

    if (length(remaining) > 0L) {
        for (link_name in names(remaining)) {
            spec <- remaining[[link_name]]
            problems <- c(
                problems,
                paste0(
                    "Link table '",
                    link_name,
                    "' has an unresolved source namespace '",
                    spec$source_col,
                    "'. Its source IDs must match canonical IDs in the object or the target IDs produced by another link table."
                )
            )
        }
    }

    unique(problems)
}

.subset_links_after_rows <- function(x, link_list) {
    link_list <- .normalize_links(link_list)

    if (length(link_list) == 0L) {
        return(link_list)
    }

    keep_gene_to_genome <- "gene_to_genome" %in% names(link_list)
    working_links <- as.list(link_list)

    if (keep_gene_to_genome) {
        working_links[["gene_to_genome"]] <- .core_gene_to_genome_link(x)
    }

    retained_ids <- list(
        gene_id = .normalize_feature_ids(x),
        genome_id = unique(as.character(.core_gene_to_genome_link(x)$genome_id))
    )

    changed <- TRUE
    while (changed) {
        changed <- FALSE

        for (link_name in names(working_links)) {
            link_table <- working_links[[link_name]]
            link_cols <- .link_columns(link_table)

            if (is.null(link_cols)) {
                next
            }

            source_ids <- retained_ids[[link_cols[["source"]]]]
            if (is.null(source_ids)) {
                next
            }

            keep <- as.character(link_table[[link_cols[["source"]]]]) %in% source_ids
            target_ids <- unique(as.character(link_table[[link_cols[["target"]]]][keep]))
            target_ids <- target_ids[!is.na(target_ids) & target_ids != ""]

            old_target_ids <- retained_ids[[link_cols[["target"]]]]
            new_target_ids <- if (is.null(old_target_ids)) {
                target_ids
            } else {
                unique(c(old_target_ids, target_ids))
            }

            if (is.null(old_target_ids) || !identical(old_target_ids, new_target_ids)) {
                retained_ids[[link_cols[["target"]]]] <- new_target_ids
                changed <- TRUE
            }
        }
    }

    pruned_links <- lapply(names(working_links), function(link_name) {
        link_table <- working_links[[link_name]]
        link_cols <- .link_columns(link_table)

        if (is.null(link_cols)) {
            return(link_table)
        }

        source_ids <- retained_ids[[link_cols[["source"]]]]
        if (is.null(source_ids)) {
            return(link_table)
        }

        keep <- as.character(link_table[[link_cols[["source"]]]]) %in% source_ids
        link_table[keep, , drop = FALSE]
    })
    names(pruned_links) <- names(working_links)

    methods::as(pruned_links, "SimpleList")
}

.prune_genome_experiment_after_rows <- function(x) {
    genome_experiment <- .get_genome_experiment(x)

    if (is.null(genome_experiment)) {
        return(x)
    }

    kept_genomes <- unique(as.character(.core_gene_to_genome_link(x)$genome_id))
    keep <- rownames(genome_experiment) %in% kept_genomes
    genome_experiment <- genome_experiment[keep, , drop = FALSE]
    genome_experiment <- .set_stored_genome_tree(
        genome_experiment,
        .prune_genome_tree(.stored_genome_tree(genome_experiment), rownames(genome_experiment))
    )

    .set_genome_experiment(x, genome_experiment)
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

.normalize_genome_tree <- function(tree, genome_ids = NULL) {
    if (is.null(tree)) {
        return(NULL)
    }

    if (!inherits(tree, "phylo")) {
        stop("'genomeTree' must be an ape::phylo object or NULL.", call. = FALSE)
    }

    tip_labels <- tree$tip.label
    if (is.null(tip_labels) || anyNA(tip_labels) || any(tip_labels == "")) {
        stop(
            "'genomeTree' must have non-empty, non-missing tip labels.",
            call. = FALSE
        )
    }

    tip_labels <- as.character(tip_labels)
    if (anyDuplicated(tip_labels)) {
        stop("'genomeTree' tip labels must be unique.", call. = FALSE)
    }

    tree$tip.label <- tip_labels

    if (!is.null(genome_ids)) {
        genome_ids <- as.character(genome_ids)
        genome_ids <- genome_ids[!is.na(genome_ids) & genome_ids != ""]

        if (!setequal(tip_labels, genome_ids)) {
            stop(
                paste(
                    "'genomeTree' tip labels must match",
                    "'rownames(genomeExperiment(x))' exactly."
                ),
                call. = FALSE
            )
        }
    }

    tree
}

.stored_genome_tree <- function(genome_experiment) {
    if (is.null(genome_experiment)) {
        return(NULL)
    }

    S4Vectors::metadata(genome_experiment)[[.mttk_genome_tree_metadata_name]]
}

.set_stored_genome_tree <- function(genome_experiment, tree) {
    metadata_list <- S4Vectors::metadata(genome_experiment)
    metadata_list[[.mttk_genome_tree_metadata_name]] <- tree
    S4Vectors::metadata(genome_experiment) <- metadata_list
    genome_experiment
}

.get_genome_tree <- function(x) {
    genome_experiment <- .get_genome_experiment(x)

    if (is.null(genome_experiment)) {
        return(NULL)
    }

    .stored_genome_tree(genome_experiment)
}

.empty_genome_experiment_from_mapping <- function(x) {
    genome_ids <- unique(as.character(.core_gene_to_genome_link(x)$genome_id))
    genome_ids <- genome_ids[!is.na(genome_ids) & genome_ids != ""]

    if (length(genome_ids) == 0L) {
        stop(
            "No genome identifiers are available to store a genome phylogeny.",
            call. = FALSE
        )
    }

    .empty_genome_experiment(
        col_data = SummarizedExperiment::colData(x),
        genome_data = S4Vectors::DataFrame(row.names = genome_ids)
    )
}

.set_genome_tree <- function(x, value) {
    genome_experiment <- .get_genome_experiment(x)

    if (is.null(genome_experiment)) {
        if (is.null(value)) {
            methods::validObject(x)
            return(x)
        }

        genome_experiment <- .empty_genome_experiment_from_mapping(x)
    }

    value <- .normalize_genome_tree(
        value,
        genome_ids = rownames(genome_experiment)
    )
    genome_experiment <- .set_stored_genome_tree(genome_experiment, value)
    .set_genome_experiment(x, genome_experiment)
}

.prune_genome_tree <- function(tree, genome_ids) {
    if (is.null(tree)) {
        return(NULL)
    }

    genome_ids <- unique(as.character(genome_ids))
    genome_ids <- genome_ids[!is.na(genome_ids) & genome_ids != ""]

    if (length(genome_ids) == 0L) {
        return(NULL)
    }

    drop_ids <- setdiff(tree$tip.label, genome_ids)
    if (length(drop_ids) > 0L) {
        tree <- ape::drop.tip(tree, tip = drop_ids)
    }

    .normalize_genome_tree(tree, genome_ids = genome_ids)
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
        stop("Unknown mttk component: ", name, call. = FALSE)
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
    link_table <- .core_gene_to_genome_link(x)
    stats::setNames(as.character(link_table$genome_id), link_table$gene_id)
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

    gene_ids <- tryCatch(
        .normalize_feature_ids(object),
        error = function(e) {
            problems <<- c(problems, conditionMessage(e))
            NULL
        }
    )
    sample_ids <- tryCatch(
        .normalize_sample_ids(object),
        error = function(e) {
            problems <<- c(problems, conditionMessage(e))
            NULL
        }
    )

    if (!is.null(gene_ids)) {
        gene_id_problem <- .validate_id_mirror(
            ids = gene_ids,
            data = SummarizedExperiment::rowData(object),
            column = "gene_id",
            mirror_label = "rowData(x)$gene_id",
            canonical_label = "rownames(x)"
        )

        if (!is.null(gene_id_problem)) {
            problems <- c(problems, gene_id_problem)
        }
    }

    if (!is.null(sample_ids)) {
        sample_id_problem <- .validate_id_mirror(
            ids = sample_ids,
            data = SummarizedExperiment::colData(object),
            column = "sample_id",
            mirror_label = "colData(x)$sample_id",
            canonical_label = "colnames(x)"
        )

        if (!is.null(sample_id_problem)) {
            problems <- c(problems, sample_id_problem)
        }
    }

    core_gene_to_genome <- tryCatch(
        .core_gene_to_genome_link(object),
        error = function(e) {
            problems <<- c(problems, conditionMessage(e))
            NULL
        }
    )

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

        genome_id_problem <- .validate_id_mirror(
            ids = rownames(genome_experiment),
            data = SummarizedExperiment::rowData(genome_experiment),
            column = "genome_id",
            mirror_label = "genomeData(x)$genome_id",
            canonical_label = "rownames(genomeExperiment(x))"
        )

        if (!is.null(genome_id_problem)) {
            problems <- c(problems, genome_id_problem)
        }

        genome_tree_problem <- tryCatch(
            {
                .normalize_genome_tree(
                    .stored_genome_tree(genome_experiment),
                    genome_ids = rownames(genome_experiment)
                )
                NULL
            },
            error = function(e) conditionMessage(e)
        )

        if (!is.null(genome_tree_problem)) {
            problems <- c(problems, genome_tree_problem)
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

    if (!is.null(core_gene_to_genome) &&
        !is.null(link_tables) &&
        methods::is(link_tables, "SimpleList") &&
        "gene_to_genome" %in% names(link_tables)) {
        gene_to_genome_link <- link_tables[["gene_to_genome"]]

        if (!all(c("gene_id", "genome_id") %in% names(gene_to_genome_link))) {
            problems <- c(
                problems,
                "The 'gene_to_genome' link table must contain 'gene_id' and 'genome_id' columns."
            )
        } else {
            link_df <- as.data.frame(
                gene_to_genome_link[, c("gene_id", "genome_id"), drop = FALSE]
            )

            if (anyDuplicated(link_df$gene_id)) {
                problems <- c(
                    problems,
                    "The 'gene_to_genome' link table must map each gene to exactly one genome."
                )
            } else if (!setequal(link_df$gene_id, core_gene_to_genome$gene_id)) {
                problems <- c(
                    problems,
                    "The 'gene_to_genome' link table must describe the same genes as rowData(x)$genome_id."
                )
            } else {
                matched <- match(core_gene_to_genome$gene_id, link_df$gene_id)
                link_genome_ids <- as.character(link_df$genome_id[matched])

                if (!identical(link_genome_ids, as.character(core_gene_to_genome$genome_id))) {
                    problems <- c(
                        problems,
                        "The 'gene_to_genome' link table must match rowData(x)$genome_id exactly."
                    )
                }
            }
        }
    }

    if (!is.null(link_tables) &&
        methods::is(link_tables, "SimpleList")) {
        known_ids <- list(
            gene_id = gene_ids,
            sample_id = sample_ids
        )

        if (!is.null(core_gene_to_genome)) {
            known_ids[["genome_id"]] <- unique(as.character(core_gene_to_genome$genome_id))
        }

        problems <- c(
            problems,
            .validate_link_referential_integrity(link_tables, known_ids)
        )
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
        if (!is.null(core_gene_to_genome)) {
            gene_to_genome <- stats::setNames(
                as.character(core_gene_to_genome$genome_id),
                core_gene_to_genome$gene_id
            )
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

.empty_mttk_state <- function() {
    list(
        genomeData = S4Vectors::DataFrame(),
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

.mttk_state <- function(x) {
    defaults <- .empty_mttk_state()
    mttk <- S4Vectors::metadata(x)[["mttk"]]

    if (is.null(mttk)) {
        return(defaults)
    }

    if (!is.list(mttk)) {
        stop("'metadata(x)$mttk' must be a list.", call. = FALSE)
    }

    defaults$genomeData <- .normalize_genome_data(mttk[["genomeData"]])
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
        genomeData = .normalize_genome_data(value),
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
    genomeData = NULL,
    links = NULL,
    activeHierarchies = NULL
) {
    state <- .mttk_state(x)

    if (!is.null(genomeData)) {
        state$genomeData <- .normalize_genome_data(genomeData)
    }

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

.validate_mttk_experiment <- function(object) {
    problems <- character()
    mttk <- S4Vectors::metadata(object)[["mttk"]]

    if (is.null(mttk)) {
        return(TRUE)
    }

    if (!is.list(mttk)) {
        return("'metadata(x)$mttk' must be a list.")
    }

    genome_data <- mttk[["genomeData"]]
    if (!is.null(genome_data)) {
        if (!methods::is(genome_data, "DataFrame")) {
            problems <- c(
                problems,
                "'metadata(x)$mttk$genomeData' must be an S4Vectors::DataFrame."
            )
        } else {
            genome_ids <- rownames(genome_data)
            if (!is.null(genome_ids) && anyDuplicated(genome_ids)) {
                problems <- c(
                    problems,
                    "Row names of 'genomeData' must be unique when present."
                )
            }
        }
    }

    link_tables <- mttk[["links"]]
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

    hierarchy_names <- mttk[["activeHierarchies"]]
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

    if (length(problems) == 0L) {
        TRUE
    } else {
        problems
    }
}

#' Construct an `MTTKFit`
#'
#' `MTTKFit()` builds a compact result container for MTTK modeling workflows.
#' The object behaves like an `S4Vectors::DataFrame`, with one row per modeled
#' feature and columns for coefficient summaries, fit diagnostics, and feature
#' metadata.
#'
#' @param results A `S4Vectors::DataFrame`, `data.frame`, or `NULL`.
#' @param info A named list describing the model specification and provenance.
#' @param models An optional named list of backend model objects.
#'
#' @return A valid `MTTKFit`.
#'
#' @examples
#' fit <- MTTKFit(
#'     results = data.frame(
#'         ko_id = "K00001",
#'         estimate = 0.5,
#'         p_value = 0.01,
#'         row.names = "K00001"
#'     ),
#'     info = list(backend = "example")
#' )
#'
#' fit
#' fitInfo(fit)$backend
#'
#' @export
MTTKFit <- function(results = NULL, info = list(), models = list()) {
    if (is.null(results)) {
        results <- S4Vectors::DataFrame()
    } else if (is.data.frame(results) && !methods::is(results, "DataFrame")) {
        results <- S4Vectors::DataFrame(results, check.names = FALSE)
    }

    if (!methods::is(results, "DataFrame")) {
        stop("'results' must be an S4Vectors::DataFrame, a data.frame, or NULL.", call. = FALSE)
    }

    if (is.null(info)) {
        info <- list()
    }
    if (!is.list(info)) {
        stop("'info' must be a list.", call. = FALSE)
    }

    if (is.null(models)) {
        models <- list()
    }
    if (!is.list(models)) {
        stop("'models' must be a list.", call. = FALSE)
    }

    out <- methods::as(results, "MTTKFit")
    metadata_list <- S4Vectors::metadata(out)
    metadata_list$mttk_fit <- list(
        info = info,
        models = models
    )
    S4Vectors::metadata(out) <- metadata_list
    methods::validObject(out)
    out
}

.fit_feature_ids <- function(x) {
    info <- fitInfo(x)
    feature_id_column <- if (!is.null(info$featureIdColumn)) {
        as.character(info$featureIdColumn)[1L]
    } else {
        NA_character_
    }

    if (!is.na(feature_id_column) &&
        feature_id_column != "" &&
        feature_id_column %in% names(x)) {
        ids <- as.character(x[[feature_id_column]])
    } else {
        id_columns <- c("ko_id", "gene_id", "genome_id", "module_id", "pathway_id")
        present_id_columns <- intersect(id_columns, names(x))

        if (length(present_id_columns) > 0L) {
            ids <- as.character(x[[present_id_columns[[1L]]]])
        } else {
            ids <- rownames(x)
        }
    }

    if (is.null(ids) || anyNA(ids) || any(ids == "")) {
        stop(
            paste(
                "The fit must contain canonical feature identifiers in",
                "'ko_id', 'gene_id', 'genome_id', 'module_id', 'pathway_id',",
                "the configured featureIdColumn, or row names."
            ),
            call. = FALSE
        )
    }

    ids
}

.subset_named_models <- function(models, feature_ids) {
    if (!is.list(models) || length(models) == 0L) {
        return(list())
    }

    matched <- match(feature_ids, names(models))
    kept <- matched[!is.na(matched)]
    out <- models[kept]
    names(out) <- feature_ids[!is.na(matched)]
    out
}

.subset_group_effects <- function(group_effects, feature_ids, feature_id_column) {
    if (is.null(group_effects)) {
        return(NULL)
    }

    if (is.data.frame(group_effects) && !methods::is(group_effects, "DataFrame")) {
        group_effects <- S4Vectors::DataFrame(group_effects, check.names = FALSE)
    }

    if (!methods::is(group_effects, "DataFrame")) {
        return(group_effects)
    }

    if (is.null(feature_id_column) ||
        length(feature_id_column) != 1L ||
        is.na(feature_id_column) ||
        feature_id_column == "" ||
        !(feature_id_column %in% names(group_effects))) {
        return(group_effects)
    }

    keep <- as.character(group_effects[[feature_id_column]]) %in% as.character(feature_ids)
    group_effects[keep, , drop = FALSE]
}

.sort_fit_indices <- function(x, sortBy, decreasing) {
    if (is.null(sortBy)) {
        return(seq_len(nrow(x)))
    }

    if (!is.character(sortBy) || length(sortBy) != 1L || is.na(sortBy) || sortBy == "") {
        stop("'sortBy' must be NULL or a single non-empty sort key.", call. = FALSE)
    }

    if (identical(sortBy, "abs_estimate")) {
        if (!("estimate" %in% names(x))) {
            stop("The fit does not contain an 'estimate' column.", call. = FALSE)
        }
        key <- abs(as.numeric(x$estimate))
    } else {
        if (!(sortBy %in% names(x))) {
            stop("Unknown fit result column: ", sortBy, ".", call. = FALSE)
        }
        key <- x[[sortBy]]
    }

    if (missing(decreasing) || is.null(decreasing)) {
        decreasing <- !(sortBy %in% c("p_value", "q_value"))
    }

    if (!is.logical(decreasing) || length(decreasing) != 1L || is.na(decreasing)) {
        stop("'decreasing' must be TRUE, FALSE, or NULL.", call. = FALSE)
    }

    key <- S4Vectors::decode(key)

    if (is.numeric(key) || is.integer(key)) {
        if (decreasing) {
            return(order(is.na(key), -key, seq_along(key)))
        }
        return(order(is.na(key), key, seq_along(key)))
    }

    if (decreasing) {
        key <- base::xtfrm(key)
        return(order(is.na(key), -key, seq_along(key)))
    }

    order(is.na(key), key, seq_along(key))
}

.join_fit_annotations <- function(feature_ids, link_list, path, collapse) {
    path <- as.character(path)

    if (length(path) == 0L || anyNA(path) || any(path == "")) {
        stop("'path' must contain one or more named link tables.", call. = FALSE)
    }

    if (!is.character(collapse) || length(collapse) != 1L || is.na(collapse)) {
        stop("'collapse' must be a single character string.", call. = FALSE)
    }

    missing_links <- setdiff(path, names(link_list))
    if (length(missing_links) > 0L) {
        stop(
            "Unknown link table(s): ",
            paste(missing_links, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    mapping <- data.frame(
        feature_id = as.character(feature_ids),
        current_id = as.character(feature_ids),
        stringsAsFactors = FALSE
    )
    annotation_cols <- list()

    for (link_name in path) {
        link_table <- as.data.frame(link_list[[link_name]])

        if (ncol(link_table) < 2L) {
            stop(
                "Link table '",
                link_name,
                "' must contain at least two columns.",
                call. = FALSE
            )
        }

        source_col <- names(link_table)[1L]
        target_col <- names(link_table)[2L]
        names(link_table)[1:2] <- c("source_id", "target_id")
        link_table <- unique(link_table[, c("source_id", "target_id"), drop = FALSE])

        mapping <- merge(
            mapping,
            link_table,
            by.x = "current_id",
            by.y = "source_id",
            all.x = TRUE,
            sort = FALSE
        )

        split_values <- split(mapping$target_id, mapping$feature_id)
        values <- rep(NA_character_, length(feature_ids))
        names(values) <- feature_ids
        values[names(split_values)] <- vapply(
            split_values,
            function(x) {
                x <- unique(as.character(x))
                x <- x[!is.na(x) & x != ""]
                if (length(x) == 0L) {
                    NA_character_
                } else {
                    paste(x, collapse = collapse)
                }
            },
            character(1)
        )
        annotation_cols[[target_col]] <- values[feature_ids]

        mapping <- mapping[!is.na(mapping$target_id) & mapping$target_id != "", , drop = FALSE]
        if (nrow(mapping) == 0L) {
            break
        }

        mapping <- mapping[, c("feature_id", "target_id"), drop = FALSE]
        names(mapping)[2L] <- "current_id"
    }

    if (length(annotation_cols) == 0L) {
        return(S4Vectors::DataFrame(row.names = feature_ids))
    }

    out <- do.call(S4Vectors::DataFrame, c(annotation_cols, list(row.names = feature_ids)))
    rownames(out) <- feature_ids
    out
}

.normalize_annotation_paths <- function(path) {
    if (is.character(path)) {
        if (length(path) == 0L || anyNA(path) || any(path == "")) {
            stop("'path' must contain one or more named link tables.", call. = FALSE)
        }
        return(list(path))
    }

    if (!is.list(path) || length(path) == 0L) {
        stop(
            "'path' must be a character vector or a list of character vectors naming link tables.",
            call. = FALSE
        )
    }

    out <- lapply(path, function(one_path) {
        if (!is.character(one_path) || length(one_path) == 0L || anyNA(one_path) || any(one_path == "")) {
            stop(
                "Each element of 'path' must be a non-empty character vector naming link tables.",
                call. = FALSE
            )
        }

        one_path
    })

    if (is.null(names(out))) {
        names(out) <- paste0("path_", seq_along(out))
    }

    out
}

.normalize_module_table <- function(module_table) {
    if (is.null(module_table)) {
        return(NULL)
    }

    if (is.character(module_table) && length(module_table) == 1L && !is.na(module_table)) {
        return(readKEGGModuleTable(module_table))
    }

    if (is.data.frame(module_table) && !methods::is(module_table, "DataFrame")) {
        module_table <- S4Vectors::DataFrame(module_table, check.names = FALSE)
    }

    if (!methods::is(module_table, "DataFrame")) {
        stop(
            "'moduleTable' must be NULL, a file path, an S4Vectors::DataFrame, or a data.frame.",
            call. = FALSE
        )
    }

    normalized <- module_table
    current_names <- names(normalized)

    rename_map <- c(
        module = "module_id",
        definition = "module_definition",
        name = "module_name",
        class = "module_class"
    )

    for (from in names(rename_map)) {
        to <- rename_map[[from]]
        if (from %in% current_names && !(to %in% current_names)) {
            names(normalized)[names(normalized) == from] <- to
        }
    }

    if (!("module_id" %in% names(normalized))) {
        stop("'moduleTable' must contain a 'module_id' column.", call. = FALSE)
    }

    keep_cols <- c(
        "module_id",
        intersect(
            c("module_name", "module_class", "module_definition"),
            names(normalized)
        )
    )
    normalized <- normalized[, keep_cols, drop = FALSE]
    normalized$module_id <- as.character(normalized$module_id)
    normalized <- normalized[!is.na(normalized$module_id) & normalized$module_id != "", , drop = FALSE]

    if (nrow(normalized) == 0L) {
        stop("'moduleTable' does not contain any non-empty module identifiers.", call. = FALSE)
    }

    keep_rows <- !duplicated(normalized$module_id)
    normalized <- normalized[keep_rows, , drop = FALSE]
    rownames(normalized) <- normalized$module_id
    normalized
}

.normalize_kegg_module_ids <- function(module_ids) {
    if (is.null(module_ids)) {
        return(NULL)
    }

    if (!is.character(module_ids)) {
        stop(
            "'moduleIds' must be NULL or a character vector of KEGG module identifiers.",
            call. = FALSE
        )
    }

    module_ids <- trimws(as.character(module_ids))
    module_ids <- module_ids[!is.na(module_ids) & module_ids != ""]

    if (length(module_ids) == 0L) {
        stop("'moduleIds' does not contain any non-empty module identifiers.", call. = FALSE)
    }

    module_ids <- unlist(strsplit(module_ids, ";", fixed = TRUE), use.names = FALSE)
    module_ids <- trimws(module_ids)
    module_ids <- module_ids[module_ids != ""]
    module_ids <- sub("^md:", "", module_ids)
    unique(module_ids)
}

.collapse_kegg_text <- function(values, collapse = "; ") {
    values <- trimws(as.character(values))
    values <- unique(values[!is.na(values) & values != ""])

    if (length(values) == 0L) {
        return(NA_character_)
    }

    paste(values, collapse = collapse)
}

.kegg_entry_id <- function(entry) {
    entry_value <- entry$ENTRY
    if (length(entry_value) == 0L) {
        return(NA_character_)
    }

    entry_value <- trimws(as.character(entry_value)[1L])
    entry_value <- sub("\\s.*$", "", entry_value)
    sub("^md:", "", entry_value)
}

.fetch_kegg_module_listing <- function() {
    tryCatch(
        KEGGREST::keggList("module"),
        error = function(e) {
            stop(
                "KEGGREST failed to retrieve the KEGG module list: ",
                conditionMessage(e),
                call. = FALSE
            )
        }
    )
}

.fetch_kegg_module_entries <- function(module_ids, batch_size, sleep) {
    batches <- split(
        paste0("md:", module_ids),
        ceiling(seq_along(module_ids) / batch_size)
    )

    entries <- list()

    for (i in seq_along(batches)) {
        batch_entries <- tryCatch(
            KEGGREST::keggGet(batches[[i]]),
            error = function(e) {
                stop(
                    "KEGGREST failed to retrieve KEGG module details: ",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        )

        entries <- c(entries, batch_entries)

        if (sleep > 0 && i < length(batches)) {
            Sys.sleep(sleep)
        }
    }

    entries
}

.collapse_unique_strings <- function(values, collapse) {
    values <- unique(as.character(values))
    values <- values[!is.na(values) & values != ""]

    if (length(values) == 0L) {
        return(NA_character_)
    }

    paste(values, collapse = collapse)
}

.split_collapsed_ids <- function(value, collapse) {
    if (is.na(value) || value == "") {
        return(character())
    }

    values <- strsplit(as.character(value), collapse, fixed = TRUE)[[1L]]
    values <- trimws(values)
    values[values != ""]
}

.expand_module_annotations <- function(module_ids, module_table, collapse) {
    module_table <- .normalize_module_table(module_table)

    if (is.null(module_table)) {
        return(S4Vectors::DataFrame(row.names = names(module_ids)))
    }

    extra_cols <- setdiff(names(module_table), "module_id")
    if (length(extra_cols) == 0L) {
        return(S4Vectors::DataFrame(row.names = names(module_ids)))
    }

    lookup <- as.data.frame(module_table)
    lookup_idx <- stats::setNames(seq_len(nrow(lookup)), lookup$module_id)

    annotation_cols <- lapply(extra_cols, function(col_name) {
        values <- vapply(
            module_ids,
            function(one_entry) {
                ids <- .split_collapsed_ids(one_entry, collapse)
                idx <- unname(lookup_idx[ids])
                idx <- idx[!is.na(idx)]

                if (length(idx) == 0L) {
                    return(NA_character_)
                }

                .collapse_unique_strings(lookup[[col_name]][idx], collapse = collapse)
            },
            character(1)
        )

        values
    })
    names(annotation_cols) <- extra_cols

    out <- do.call(
        S4Vectors::DataFrame,
        c(annotation_cols, list(row.names = names(module_ids)))
    )
    rownames(out) <- names(module_ids)
    out
}

.join_fit_annotations_multi <- function(feature_ids, link_list, path, collapse) {
    paths <- .normalize_annotation_paths(path)
    annotations <- S4Vectors::DataFrame(row.names = feature_ids)

    for (one_path in paths) {
        joined <- .join_fit_annotations(
            feature_ids = feature_ids,
            link_list = link_list,
            path = one_path,
            collapse = collapse
        )

        duplicated_cols <- intersect(names(annotations), names(joined))
        if (length(duplicated_cols) > 0L) {
            stop(
                "Annotation paths produce duplicate output column(s): ",
                paste(duplicated_cols, collapse = ", "),
                ".",
                call. = FALSE
            )
        }

        annotations <- cbind(annotations, joined)
    }

    annotations
}

#' Fit Metadata
#'
#' `fitInfo()` returns the metadata list stored alongside an `MTTKFit`.
#'
#' @param x An `MTTKFit`.
#'
#' @return A named list.
#'
#' @rdname fitInfo
#' @export
methods::setMethod("fitInfo", "MTTKFit", function(x) {
    fit_state <- S4Vectors::metadata(x)$mttk_fit

    if (is.null(fit_state) || is.null(fit_state$info)) {
        return(list())
    }

    fit_state$info
})

#' Stored Model Objects
#'
#' `modelObjects()` returns the optional backend model objects stored in an
#' `MTTKFit`.
#'
#' @param x An `MTTKFit`.
#'
#' @return A named list.
#'
#' @rdname modelObjects
#' @export
methods::setMethod("modelObjects", "MTTKFit", function(x) {
    fit_state <- S4Vectors::metadata(x)$mttk_fit

    if (is.null(fit_state) || is.null(fit_state$models)) {
        return(list())
    }

    fit_state$models
})

#' Extract KO-by-Genome Effects From a Random-Slope KO Fit
#'
#' `koGenomeEffects()` returns the KO-by-genome conditional coefficients stored
#' alongside a KO random-slope model fit. Each row corresponds to one KO/genome
#' combination and contains the conditional KO effect estimated for that genome.
#'
#' This workflow is intended as the first bridge between KO-level and
#' genome-level analyses in MTTK. The returned table can be used to ask whether
#' the same KO responds similarly or differently across genomes, and can later
#' be joined to `genomeData(x)` for taxonomy-aware summaries.
#'
#' @param x An `MTTKFit` returned by [fitKORandomSlopeModel()].
#'
#' @return An `S4Vectors::DataFrame` with one row per KO/genome combination.
#'
#' @examples
#' fit <- MTTKFit(
#'     results = data.frame(
#'         ko_id = "K00001",
#'         estimate = 0.5,
#'         row.names = "K00001"
#'     ),
#'     info = list(
#'         model = "ko_random_slope_model",
#'         featureIdColumn = "ko_id",
#'         groupEffectColumn = "genome_id"
#'     )
#' )
#' S4Vectors::metadata(fit)$mttk_fit$groupEffects <- S4Vectors::DataFrame(
#'     ko_id = "K00001",
#'     genome_id = "genome_1",
#'     conditional_effect_estimate = 0.6,
#'     row.names = "K00001::genome_1"
#' )
#' koGenomeEffects(fit)
#'
#' @export
koGenomeEffects <- function(x) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    info <- fitInfo(x)
    if (!identical(info$model, "ko_random_slope_model")) {
        stop(
            "'x' must be an MTTKFit returned by fitKORandomSlopeModel().",
            call. = FALSE
        )
    }

    .stored_group_effects(x)
}

.stored_group_effects <- function(x) {
    info <- fitInfo(x)
    fit_state <- S4Vectors::metadata(x)$mttk_fit
    if (!is.null(fit_state$groupEffects)) {
        out <- fit_state$groupEffects
        if (is.data.frame(out) && !methods::is(out, "DataFrame")) {
            out <- S4Vectors::DataFrame(out, check.names = FALSE)
        }
        return(out)
    }

    models <- modelObjects(x)
    if (length(models) == 0L) {
        stop(
            paste(
                "This fit does not contain stored genome-specific conditional effects or backend model objects.",
                "Refit with keepFits = TRUE or use the current random-slope fit output."
            ),
            call. = FALSE
        )
    }

    stop(
        paste(
            "The genome-specific conditional effects are not stored in",
            "metadata(x)$mttk_fit$groupEffects.",
            "Please refit with the current version of",
            if (!is.null(info$model)) info$model else "the random-slope workflow",
            "."
        ),
        call. = FALSE
    )
}

#' Extract Module-by-Genome Effects From a Random-Slope Module Fit
#'
#' `moduleGenomeEffects()` returns the module-by-genome conditional coefficients
#' stored alongside a module random-slope model fit. Each row corresponds to
#' one module/genome combination and contains the conditional module effect
#' estimated for that genome.
#'
#' @param x An `MTTKFit` returned by [fitModuleRandomSlopeModel()].
#'
#' @return An `S4Vectors::DataFrame` with one row per module/genome
#'   combination.
#'
#' @examples
#' fit <- MTTKFit(
#'     results = data.frame(
#'         module_id = "M00001",
#'         estimate = 0.5,
#'         row.names = "M00001"
#'     ),
#'     info = list(
#'         model = "module_random_slope_model",
#'         featureIdColumn = "module_id",
#'         groupEffectColumn = "genome_id"
#'     )
#' )
#' S4Vectors::metadata(fit)$mttk_fit$groupEffects <- S4Vectors::DataFrame(
#'     module_id = "M00001",
#'     genome_id = "genome_1",
#'     conditional_effect_estimate = 0.6,
#'     row.names = "M00001::genome_1"
#' )
#' moduleGenomeEffects(fit)
#'
#' @export
moduleGenomeEffects <- function(x) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    info <- fitInfo(x)
    if (!identical(info$model, "module_random_slope_model")) {
        stop(
            "'x' must be an MTTKFit returned by fitModuleRandomSlopeModel().",
            call. = FALSE
        )
    }

    .stored_group_effects(x)
}

#' Extract Pathway-by-Genome Effects From a Random-Slope Pathway Fit
#'
#' `pathwayGenomeEffects()` returns the pathway-by-genome conditional
#' coefficients stored alongside a pathway random-slope model fit. Each row
#' corresponds to one pathway/genome combination and contains the conditional
#' pathway effect estimated for that genome.
#'
#' @param x An `MTTKFit` returned by [fitPathwayRandomSlopeModel()].
#'
#' @return An `S4Vectors::DataFrame` with one row per pathway/genome
#'   combination.
#'
#' @examples
#' fit <- MTTKFit(
#'     results = data.frame(
#'         pathway_id = "map00010",
#'         estimate = 0.5,
#'         row.names = "map00010"
#'     ),
#'     info = list(
#'         model = "pathway_random_slope_model",
#'         featureIdColumn = "pathway_id",
#'         groupEffectColumn = "genome_id"
#'     )
#' )
#' S4Vectors::metadata(fit)$mttk_fit$groupEffects <- S4Vectors::DataFrame(
#'     pathway_id = "map00010",
#'     genome_id = "genome_1",
#'     conditional_effect_estimate = 0.6,
#'     row.names = "map00010::genome_1"
#' )
#' pathwayGenomeEffects(fit)
#'
#' @export
pathwayGenomeEffects <- function(x) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    info <- fitInfo(x)
    if (!identical(info$model, "pathway_random_slope_model")) {
        stop(
            "'x' must be an MTTKFit returned by fitPathwayRandomSlopeModel().",
            call. = FALSE
        )
    }

    .stored_group_effects(x)
}

#' Subset an MTTKFit
#'
#' `x[i, j]` subsets an `MTTKFit` by rows and columns while preserving fit
#' metadata. When rows are subset, stored backend model objects are pruned to
#' the retained features.
#'
#' @param x An `MTTKFit`.
#' @param i,j Row and column indices.
#' @param ... Additional arguments passed to the inherited subsetting method.
#' @param drop Included for compatibility and ignored for `MTTKFit`.
#'
#' @return A valid subsetted `MTTKFit`.
#'
#' @rdname subset-MTTKFit
#' @name subset-MTTKFit
#' @aliases [,MTTKFit,ANY,ANY,ANY-method
#' @exportMethod [
methods::setMethod(
    "[",
    c("MTTKFit", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE) {
        out <- methods::as(methods::callNextMethod(), "MTTKFit")
        fit_state <- S4Vectors::metadata(out)$mttk_fit

        if (!is.null(fit_state)) {
            if (!is.null(fit_state$models)) {
                fit_state$models <- .subset_named_models(
                    fit_state$models,
                    feature_ids = rownames(out)
                )
            }
            fit_state$groupEffects <- .subset_group_effects(
                fit_state$groupEffects,
                feature_ids = rownames(out),
                feature_id_column = if (!is.null(fit_state$info$featureIdColumn)) {
                    as.character(fit_state$info$featureIdColumn)[1L]
                } else {
                    NA_character_
                }
            )
            metadata_list <- S4Vectors::metadata(out)
            metadata_list$mttk_fit <- fit_state
            S4Vectors::metadata(out) <- metadata_list
        }

        methods::validObject(out)
        out
    }
)

#' Extract a Sorted Fit Result Table
#'
#' `fitTable()` returns an ordered `S4Vectors::DataFrame` view of an
#' `MTTKFit`. This is useful when the fit should be inspected as a regular
#' result table, optionally filtered to a subset of fit statuses and sorted by
#' p-values, q-values, effect size, or another result column.
#'
#' @param x An `MTTKFit`.
#' @param status Optional character vector of fit statuses to keep, such as
#'   `"ok"` or `c("ok", "skipped")`. Use `NULL` to keep all rows.
#' @param sortBy Optional sort key. This can be any result column name or the
#'   special value `"abs_estimate"`.
#' @param decreasing Logical or `NULL`. When `NULL`, p-values and q-values sort
#'   in ascending order and other keys sort in descending order.
#' @param n Optional maximum number of rows to return after filtering and
#'   sorting.
#'
#' @return An `S4Vectors::DataFrame`.
#'
#' @examples
#' fit <- MTTKFit(
#'     results = data.frame(
#'         ko_id = c("K1", "K2"),
#'         estimate = c(0.3, -1.2),
#'         q_value = c(0.2, 0.01),
#'         status = c("ok", "ok"),
#'         row.names = c("K1", "K2")
#'     )
#' )
#'
#' fitTable(fit, sortBy = "q_value")
#'
#' @export
fitTable <- function(x, status = NULL, sortBy = NULL, decreasing = NULL, n = NULL) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    out <- x

    if (!is.null(status)) {
        if (!("status" %in% names(out))) {
            stop("The fit does not contain a 'status' column.", call. = FALSE)
        }

        keep <- as.character(out$status) %in% as.character(status)
        out <- out[keep, , drop = FALSE]
    }

    if (!is.null(n)) {
        if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0) {
            stop("'n' must be NULL or a single non-negative number.", call. = FALSE)
        }
        n <- as.integer(n)
    }

    ord <- .sort_fit_indices(out, sortBy = sortBy, decreasing = decreasing)
    out <- out[ord, , drop = FALSE]

    if (!is.null(n) && n < nrow(out)) {
        out <- out[seq_len(n), , drop = FALSE]
    }

    S4Vectors::DataFrame(out, check.names = FALSE)
}

#' Filter Significant Fit Results
#'
#' `significantResults()` returns the subset of an `MTTKFit` that passes a
#' significance threshold. By default it uses `q_value <= 0.05` and keeps only
#' rows with `status == "ok"`.
#'
#' The returned object is still an `MTTKFit`, so stored model objects and fit
#' metadata remain aligned to the retained rows.
#'
#' @param x An `MTTKFit`.
#' @param alpha Significance threshold.
#' @param value Which significance column to use, typically `"q_value"` or
#'   `"p_value"`.
#' @param status Optional fit status to keep before significance filtering.
#' @param direction Optional effect-direction filter. `"up"` keeps positive
#'   estimates, `"down"` keeps negative estimates, and `"both"` keeps both.
#' @param n Optional maximum number of rows to return after filtering and
#'   sorting.
#' @param sortBy Optional sort key. When `NULL`, the significance column in
#'   `value` is used.
#' @param decreasing Logical or `NULL`, passed to the sorting step.
#'
#' @return An `MTTKFit`.
#'
#' @examples
#' fit <- MTTKFit(
#'     results = data.frame(
#'         ko_id = c("K1", "K2"),
#'         estimate = c(0.3, -1.2),
#'         q_value = c(0.2, 0.01),
#'         status = c("ok", "ok"),
#'         row.names = c("K1", "K2")
#'     )
#' )
#'
#' significantResults(fit)
#'
#' @export
significantResults <- function(
    x,
    alpha = 0.05,
    value = c("q_value", "p_value"),
    status = "ok",
    direction = c("both", "up", "down"),
    n = NULL,
    sortBy = NULL,
    decreasing = NULL
) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    value <- match.arg(value)
    direction <- match.arg(direction)

    if (!(value %in% names(x))) {
        stop("The fit does not contain a '", value, "' column.", call. = FALSE)
    }

    if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha < 0 || alpha > 1) {
        stop("'alpha' must be a single number between 0 and 1.", call. = FALSE)
    }

    out <- x
    if (!is.null(status)) {
        if (!("status" %in% names(out))) {
            stop("The fit does not contain a 'status' column.", call. = FALSE)
        }
        out <- out[as.character(out$status) %in% as.character(status), , drop = FALSE]
    }

    sig_values <- as.numeric(out[[value]])
    keep <- !is.na(sig_values) & sig_values <= alpha

    if (!identical(direction, "both")) {
        if (!("estimate" %in% names(out))) {
            stop("The fit does not contain an 'estimate' column.", call. = FALSE)
        }

        estimate <- as.numeric(out$estimate)
        if (identical(direction, "up")) {
            keep <- keep & !is.na(estimate) & estimate > 0
        } else {
            keep <- keep & !is.na(estimate) & estimate < 0
        }
    }

    out <- out[keep, , drop = FALSE]
    ord <- .sort_fit_indices(
        out,
        sortBy = if (is.null(sortBy)) value else sortBy,
        decreasing = decreasing
    )
    out <- out[ord, , drop = FALSE]

    if (!is.null(n)) {
        if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0) {
            stop("'n' must be NULL or a single non-negative number.", call. = FALSE)
        }
        n <- as.integer(n)
        if (n < nrow(out)) {
            out <- out[seq_len(n), , drop = FALSE]
        }
    }

    out
}

#' Join KO Annotations Onto a Fit
#'
#' `annotateKOFit()` appends KO-centered annotation columns from an
#' `MTTKExperiment` onto an `MTTKFit`. This is intended for fits whose rows
#' represent KO identifiers, such as the first KO-level mixed-model workflow in
#' MTTK.
#'
#' By default the function follows two direct KO-centered mappings,
#' `ko_to_module` and `ko_to_pathway`, and adds `module_id` and `pathway_id`
#' columns. When a KO maps to multiple downstream identifiers, the unique
#' values are collapsed into a single string.
#'
#' @param x An `MTTKFit` whose rows represent KO identifiers.
#' @param object An `MTTKExperiment` supplying the annotation link tables.
#' @param path A character vector naming one sequential annotation path, or a
#'   list of character vectors naming multiple paths to traverse in parallel.
#'   The default adds direct `ko_to_module` and `ko_to_pathway` annotations.
#' @param collapse String used to combine multiple target identifiers per KO.
#' @param moduleTable Optional KEGG module annotation table used to append
#'   module description columns when `module_id` annotations are present. This
#'   can be supplied as a file path, `data.frame`, or `S4Vectors::DataFrame`.
#'   [fetchKEGGModuleTable()] can be used to retrieve this table directly from
#'   KEGG through `KEGGREST`.
#'
#' @return An `MTTKFit` with appended annotation columns.
#'
#' @examples
#' fit <- MTTKFit(
#'     results = data.frame(
#'         ko_id = c("K03043", "K02111"),
#'         estimate = c(0.3, -1.2),
#'         row.names = c("K03043", "K02111")
#'     )
#' )
#' x <- makeExampleMTTKExperiment()
#' annotateKOFit(fit, x)
#'
#' @export
annotateKOFit <- function(
    x,
    object,
    path = list(
        module = "ko_to_module",
        pathway = "ko_to_pathway"
    ),
    collapse = ";",
    moduleTable = NULL
) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    if (!methods::is(object, "MTTKExperiment")) {
        stop("'object' must be an MTTKExperiment.", call. = FALSE)
    }

    feature_ids <- .fit_feature_ids(x)
    annotations <- .join_fit_annotations_multi(
        feature_ids = feature_ids,
        link_list = links(object),
        path = path,
        collapse = collapse
    )

    if ("module_id" %in% names(annotations) && !is.null(moduleTable)) {
        module_annotations <- .expand_module_annotations(
            module_ids = stats::setNames(as.character(annotations$module_id), rownames(annotations)),
            module_table = moduleTable,
            collapse = collapse
        )

        duplicated_cols <- intersect(names(annotations), names(module_annotations))
        if (length(duplicated_cols) > 0L) {
            stop(
                "Module annotations produce duplicate output column(s): ",
                paste(duplicated_cols, collapse = ", "),
                ".",
                call. = FALSE
            )
        }

        annotations <- cbind(annotations, module_annotations[rownames(annotations), , drop = FALSE])
    }

    conflicting <- intersect(colnames(annotations), colnames(x))
    if (length(conflicting) > 0L) {
        stop(
            "The fit already contains column(s): ",
            paste(conflicting, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    results <- cbind(
        S4Vectors::DataFrame(x, check.names = FALSE),
        annotations[rownames(x), , drop = FALSE]
    )

    MTTKFit(
        results = results,
        info = fitInfo(x),
        models = modelObjects(x)
    )
}

#' Fetch KEGG Module Annotations with KEGGREST
#'
#' `fetchKEGGModuleTable()` retrieves KEGG module metadata directly from the
#' KEGG MODULE database through `KEGGREST` and returns it in the standardized
#' format used by MTTK.
#'
#' The resulting table contains `module_id`, `module_name`, `module_class`, and
#' `module_definition` columns and can be supplied directly to
#' [annotateKOFit()] through its `moduleTable` argument.
#'
#' By default the function requests the full KEGG module catalog. To reduce
#' runtime and network traffic, pass a subset of `moduleIds` when only a small
#' number of modules are needed.
#'
#' Access to the KEGG REST API is provided by `KEGGREST` and subject to KEGG's
#' academic-use terms.
#'
#' @param moduleIds Optional character vector of KEGG module identifiers. Bare
#'   IDs such as `M00001` and prefixed IDs such as `md:M00001` are both
#'   accepted. Values collapsed with `";"` are also split automatically.
#' @param batchSize Number of module entries to request per `keggGet()` call.
#'   KEGG limits `keggGet()` requests to at most 10 entries.
#' @param sleep Optional number of seconds to wait between successive KEGG REST
#'   requests.
#'
#' @return An `S4Vectors::DataFrame`.
#'
#' @examples
#' \dontrun{
#' module_table <- fetchKEGGModuleTable(c("M00001", "M00002"))
#' }
#'
#' @export
fetchKEGGModuleTable <- function(moduleIds = NULL, batchSize = 10L, sleep = 0) {
    module_ids <- .normalize_kegg_module_ids(moduleIds)

    if (!is.numeric(batchSize) || length(batchSize) != 1L || is.na(batchSize)) {
        stop("'batchSize' must be a single number between 1 and 10.", call. = FALSE)
    }

    batchSize <- as.integer(batchSize)
    if (batchSize < 1L || batchSize > 10L) {
        stop("'batchSize' must be between 1 and 10.", call. = FALSE)
    }

    if (!is.numeric(sleep) || length(sleep) != 1L || is.na(sleep) || sleep < 0) {
        stop("'sleep' must be a single non-negative number.", call. = FALSE)
    }

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop(
            "Package 'KEGGREST' is required for fetchKEGGModuleTable(). ",
            "Install it with BiocManager::install(\"KEGGREST\").",
            call. = FALSE
        )
    }

    module_listing <- .fetch_kegg_module_listing()
    listing_ids <- sub("^md:", "", names(module_listing))
    listing_names <- unname(as.character(module_listing))

    listing <- S4Vectors::DataFrame(
        module_id = listing_ids,
        module_name = listing_names,
        row.names = listing_ids
    )
    listing <- listing[!duplicated(listing$module_id), , drop = FALSE]

    if (is.null(module_ids)) {
        module_ids <- as.character(listing$module_id)
    } else {
        missing_ids <- setdiff(module_ids, as.character(listing$module_id))
        if (length(missing_ids) > 0L) {
            warning(
                "The following KEGG module IDs were not found in the KEGG module list and will be skipped: ",
                paste(missing_ids, collapse = ", "),
                call. = FALSE
            )
        }

        module_ids <- intersect(module_ids, as.character(listing$module_id))
        if (length(module_ids) == 0L) {
            stop("None of the requested 'moduleIds' were found in the KEGG module list.", call. = FALSE)
        }
    }

    listing <- listing[module_ids, , drop = FALSE]
    module_entries <- .fetch_kegg_module_entries(module_ids, batchSize, sleep)

    result <- S4Vectors::DataFrame(
        module_id = as.character(listing$module_id),
        module_name = as.character(listing$module_name),
        module_class = rep(NA_character_, nrow(listing)),
        module_definition = rep(NA_character_, nrow(listing)),
        row.names = as.character(listing$module_id)
    )

    for (entry in module_entries) {
        entry_id <- .kegg_entry_id(entry)
        if (is.na(entry_id) || !(entry_id %in% rownames(result))) {
            next
        }

        if (length(entry$NAME) > 0L) {
            result[entry_id, "module_name"] <- .collapse_kegg_text(entry$NAME)
        }

        if (length(entry$CLASS) > 0L) {
            result[entry_id, "module_class"] <- .collapse_kegg_text(entry$CLASS)
        }

        if (length(entry$DEFINITION) > 0L) {
            result[entry_id, "module_definition"] <- .collapse_kegg_text(entry$DEFINITION)
        }
    }

    .normalize_module_table(result)
}

#' Read a KEGG Module Annotation Table
#'
#' `readKEGGModuleTable()` imports a tab-separated KEGG module table and
#' standardizes the main columns used by MTTK. The returned table uses
#' `module_id`, `module_name`, `module_class`, and `module_definition` column
#' names when those fields are available.
#'
#' This is intended for tables exported from KEGG-derived resources, where the
#' original columns are named `module`,
#' `name`, `class`, and `definition`.
#'
#' @param file Path to a tab-separated module table.
#'
#' @return An `S4Vectors::DataFrame`.
#'
#' @examples
#' tf <- tempfile(fileext = ".tsv")
#' writeLines(
#'     c(
#'         "module\tdefinition\tname\tclass",
#'         "M00001\tK00001 K00002\tExample module\tPathway modules; Example"
#'     ),
#'     tf
#' )
#' readKEGGModuleTable(tf)
#'
#' @export
readKEGGModuleTable <- function(file) {
    if (!is.character(file) || length(file) != 1L || is.na(file) || file == "") {
        stop("'file' must be a single non-empty file path.", call. = FALSE)
    }

    if (!file.exists(file)) {
        stop("The module table file does not exist: ", file, ".", call. = FALSE)
    }

    module_table <- utils::read.delim(
        file,
        header = TRUE,
        sep = "\t",
        quote = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    .normalize_module_table(module_table)
}

methods::setMethod("show", "MTTKFit", function(object) {
    info <- fitInfo(object)
    backend <- info$backend
    if (is.null(backend) || length(backend) != 1L || is.na(backend) || backend == "") {
        backend <- "unspecified"
    }

    model_name <- info$model
    if (is.null(model_name) || length(model_name) != 1L || is.na(model_name) || model_name == "") {
        model_name <- "unspecified"
    }

    effect_label <- info$effectLabel
    if (is.null(effect_label) || length(effect_label) != 1L || is.na(effect_label) || effect_label == "") {
        effect_label <- NULL
    }

    status_counts <- if ("status" %in% names(object)) {
        stats::setNames(as.integer(table(as.character(object$status))), names(table(as.character(object$status))))
    } else {
        integer()
    }

    cat("class: MTTKFit\n", sep = "")
    cat("backend: ", backend, "\n", sep = "")
    cat("model: ", model_name, "\n", sep = "")
    if (!is.null(effect_label)) {
        cat("effect: ", effect_label, "\n", sep = "")
    }
    cat("features(", nrow(object), "): ", .format_preview(rownames(object)), "\n", sep = "")

    if (length(status_counts) > 0L) {
        status_label <- paste(
            paste(names(status_counts), status_counts, sep = "="),
            collapse = ", "
        )
        cat("status: ", status_label, "\n", sep = "")
    }
})

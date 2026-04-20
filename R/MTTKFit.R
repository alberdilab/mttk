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
    if ("ko_id" %in% names(x)) {
        ids <- as.character(x$ko_id)
    } else {
        ids <- rownames(x)
    }

    if (is.null(ids) || anyNA(ids) || any(ids == "")) {
        stop(
            "The fit must contain canonical feature identifiers in 'ko_id' or row names.",
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

        if (!is.null(fit_state) && !is.null(fit_state$models)) {
            fit_state$models <- .subset_named_models(
                fit_state$models,
                feature_ids = rownames(out)
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
#' By default the function follows `ko_to_module -> module_to_pathway` and adds
#' `module_id` and `pathway_id` columns. When a KO maps to multiple downstream
#' identifiers, the unique values are collapsed into a single string.
#'
#' @param x An `MTTKFit` whose rows represent KO identifiers.
#' @param object An `MTTKExperiment` supplying the annotation link tables.
#' @param path Character vector naming one or more annotation link tables to
#'   traverse.
#' @param collapse String used to combine multiple target identifiers per KO.
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
    path = c("ko_to_module", "module_to_pathway"),
    collapse = ";"
) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    if (!methods::is(object, "MTTKExperiment")) {
        stop("'object' must be an MTTKExperiment.", call. = FALSE)
    }

    feature_ids <- .fit_feature_ids(x)
    annotations <- .join_fit_annotations(
        feature_ids = feature_ids,
        link_list = links(object),
        path = path,
        collapse = collapse
    )

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

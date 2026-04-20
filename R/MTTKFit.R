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

    status_counts <- if ("status" %in% names(object)) {
        stats::setNames(as.integer(table(as.character(object$status))), names(table(as.character(object$status))))
    } else {
        integer()
    }

    cat("class: MTTKFit\n", sep = "")
    cat("backend: ", backend, "\n", sep = "")
    cat("model: ", model_name, "\n", sep = "")
    cat("features(", nrow(object), "): ", .format_preview(rownames(object)), "\n", sep = "")

    if (length(status_counts) > 0L) {
        status_label <- paste(
            paste(names(status_counts), status_counts, sep = "="),
            collapse = ", "
        )
        cat("status: ", status_label, "\n", sep = "")
    }
})

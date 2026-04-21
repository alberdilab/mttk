#' MTTKFit Class
#'
#' `MTTKFit` stores feature-level model results as a Bioconductor-compatible
#' table. Each row represents one modeled feature, such as a KO, and the
#' columns store coefficient summaries, fit diagnostics, and feature-level
#' metadata.
#'
#' The class extends [S4Vectors::DFrame-class], so it behaves like a regular
#' `S4Vectors::DataFrame`. Additional fit metadata is stored in
#' `metadata(x)$mttk_fit`.
#'
#' @section Canonical identifiers:
#' `rownames(x)` are the canonical identifiers of the modeled features. When
#' the modeled feature ID is also stored in a result column such as `ko_id`,
#' `gene_id`, `genome_id`, `module_id`, or `pathway_id`, that canonical ID
#' column must match `rownames(x)` exactly.
#'
#' @section Stored metadata:
#' `metadata(x)$mttk_fit` is reserved for MTTK-specific fit metadata and may
#' contain:
#'
#' - `info`: a named list describing the model specification, backend, and
#'   provenance.
#' - `models`: an optional named list of backend model objects.
#'
#' @name MTTKFit-class
#' @rdname MTTKFit-class
#' @aliases MTTKFit-class
#' @exportClass MTTKFit
NULL

.validate_mttk_fit <- function(object) {
    problems <- character()
    feature_ids <- rownames(object)

    if (nrow(object) > 0L &&
        (is.null(feature_ids) || anyNA(feature_ids) || any(feature_ids == ""))) {
        problems <- c(
            problems,
            "Modeled feature row names must be present, non-missing, and non-empty."
        )
    }

    if (!is.null(feature_ids) && anyDuplicated(feature_ids)) {
        problems <- c(problems, "Modeled feature row names must be unique.")
    }

    if ("p_value" %in% names(object)) {
        p_value <- as.numeric(object$p_value)
        if (any(!is.na(p_value) & (p_value < 0 | p_value > 1))) {
            problems <- c(problems, "'p_value' must lie between 0 and 1 when present.")
        }
    }

    if ("q_value" %in% names(object)) {
        q_value <- as.numeric(object$q_value)
        if (any(!is.na(q_value) & (q_value < 0 | q_value > 1))) {
            problems <- c(problems, "'q_value' must lie between 0 and 1 when present.")
        }
    }

    fit_state <- S4Vectors::metadata(object)$mttk_fit
    if (!is.null(fit_state)) {
        if (!is.list(fit_state)) {
            problems <- c(problems, "'metadata(x)$mttk_fit' must be a list when present.")
        } else {
            if (!is.null(fit_state$info) && !is.list(fit_state$info)) {
                problems <- c(problems, "'metadata(x)$mttk_fit$info' must be a list.")
            }
            if (!is.null(fit_state$models) && !is.list(fit_state$models)) {
                problems <- c(problems, "'metadata(x)$mttk_fit$models' must be a list.")
            }
            if (!is.null(fit_state$groupEffects) &&
                !methods::is(fit_state$groupEffects, "DataFrame") &&
                !is.data.frame(fit_state$groupEffects)) {
                problems <- c(
                    problems,
                    "'metadata(x)$mttk_fit$groupEffects' must be a DataFrame or data.frame."
                )
            }
        }
    }

    id_columns <- c("ko_id", "gene_id", "genome_id", "module_id", "pathway_id")
    canonical_id_column <- NULL

    if (!is.null(fit_state) &&
        is.list(fit_state) &&
        !is.null(fit_state$info) &&
        is.list(fit_state$info) &&
        !is.null(fit_state$info$featureIdColumn)) {
        canonical_id_column <- as.character(fit_state$info$featureIdColumn)[1L]
    }

    if (is.null(canonical_id_column) || is.na(canonical_id_column) || canonical_id_column == "") {
        present_id_columns <- intersect(id_columns, names(object))
        if (length(present_id_columns) > 0L) {
            canonical_id_column <- present_id_columns[[1L]]
        }
    }

    if (!is.null(canonical_id_column) &&
        canonical_id_column %in% names(object) &&
        !identical(as.character(object[[canonical_id_column]]), as.character(feature_ids))) {
        problems <- c(
            problems,
            paste0(
                "'",
                canonical_id_column,
                "' must match 'rownames(x)' exactly when present."
            )
        )
    }

    if (length(problems) == 0L) {
        TRUE
    } else {
        problems
    }
}

methods::setClass(
    "MTTKFit",
    contains = "DFrame"
)

methods::setValidity("MTTKFit", .validate_mttk_fit)

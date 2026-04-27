.extract_variance_components <- function(fit) {
    vc <- tryCatch(
        glmmTMB::VarCorr(fit)$cond,
        error = function(e) NULL
    )

    genome_var <- tryCatch({
        v <- vc[["genome_id"]]
        if (is.null(v)) NA_real_ else as.numeric(v[1L, 1L])
    }, error = function(e) NA_real_)

    block_var <- tryCatch({
        block_names <- setdiff(names(vc), "genome_id")
        if (length(block_names) == 0L) {
            NA_real_
        } else {
            v <- vc[[block_names[[1L]]]]
            if (is.null(v)) NA_real_ else as.numeric(v[1L, 1L])
        }
    }, error = function(e) NA_real_)

    # phi (NB2 size parameter): larger phi = less overdispersion
    # glmmTMB stores log(phi) as the intercept of dispformula
    phi <- tryCatch({
        disp_tab <- summary(fit)$coefficients$disp
        if (is.null(disp_tab) || !("(Intercept)" %in% rownames(disp_tab))) {
            NA_real_
        } else {
            exp(as.numeric(disp_tab["(Intercept)", "Estimate"]))
        }
    }, error = function(e) NA_real_)

    # Latent-scale residual variance for NB2 on the log link:
    # sigma^2_epsilon approx log(1 + 1/phi)  (Nakagawa et al. 2017)
    residual_var <- if (!is.na(phi) && phi > 0) log(1 + 1 / phi) else NA_real_

    list(
        genome_var   = genome_var,
        block_var    = block_var,
        phi          = phi,
        residual_var = residual_var
    )
}

#' Decompose Variance from a Fitted Mixed Model
#'
#' `varianceDecomposition()` extracts per-feature variance components from the
#' backend `glmmTMB` model objects stored in an `MTTKFit`. It reports the genome
#' random-effect variance, optional sample-block variance, the NB2 dispersion
#' parameter `phi`, an approximated latent residual variance, and intraclass
#' correlation coefficients (ICC) attributable to genome identity.
#'
#' The function requires that the fit was created with `keepFits = TRUE` so that
#' backend model objects are available.
#'
#' **Variance on the log (latent) scale** — all variance components are on the
#' log link scale, making them comparable across features with different baseline
#' expression levels:
#'
#' - `genome_var`: variance of the genome random intercept `(1 | genome_id)`.
#'   Captures how much baseline expression differs between genomes.
#' - `block_var`: variance of the sample-block random intercept, if a
#'   `sampleBlock` was used during model fitting. `NA` when absent.
#' - `phi`: NB2 size parameter (larger = less overdispersion, more
#'   Poisson-like). Extracted from the fitted dispersion sub-model.
#' - `residual_var`: approximated distribution-specific residual variance on
#'   the log scale, `log(1 + 1/phi)` (Nakagawa *et al.* 2017). Smaller when
#'   expression is more consistent across replicates.
#' - `total_var`: sum of all variance components.
#' - `icc_genome`: proportion of total variance attributable to genome identity.
#'   High ICC means expression is genome-driven; low ICC means it is mostly
#'   condition- or replicate-driven.
#' - `icc_block`: proportion of total variance attributable to the sample block,
#'   when applicable.
#'
#' @param x An `MTTKFit` created with `keepFits = TRUE`.
#'
#' @return An `S4Vectors::DataFrame` with one row per feature, containing
#'   columns `genome_var`, `block_var`, `phi`, `residual_var`, `total_var`,
#'   `icc_genome`, and `icc_block`.
#'
#' @references
#' Nakagawa S, Johnson PCD, Schielzeth H (2017). The coefficient of
#' determination R² and intra-class correlation coefficient from generalized
#' linear mixed-effects models revisited and expanded. *Journal of the Royal
#' Society Interface*, 14(134), 20170213.
#'
#' @seealso [fitKOMixedModel()], [fitKORandomSlopeModel()],
#'   [fitKODispersion()] for the dispersion-focused alternative.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     fit <- fitKOMixedModel(x, variable = "condition", keepFits = TRUE)
#'     varianceDecomposition(fit)
#' }
#'
#' @export
varianceDecomposition <- function(x) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    models <- modelObjects(x)
    if (length(models) == 0L) {
        stop(
            "No backend model objects are stored in this MTTKFit. ",
            "Refit the model with 'keepFits = TRUE'.",
            call. = FALSE
        )
    }
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop(
            "The 'glmmTMB' package must be installed to use varianceDecomposition().",
            call. = FALSE
        )
    }

    feature_ids <- rownames(x)
    info        <- fitInfo(x)
    feature_id_col <- if (!is.null(info$featureIdColumn) &&
                          !is.na(info$featureIdColumn) &&
                          info$featureIdColumn != "") {
        info$featureIdColumn
    } else {
        "feature_id"
    }

    rows <- lapply(feature_ids, function(fid) {
        fit <- models[[fid]]
        if (is.null(fit)) {
            return(S4Vectors::DataFrame(
                genome_var   = NA_real_,
                block_var    = NA_real_,
                phi          = NA_real_,
                residual_var = NA_real_,
                total_var    = NA_real_,
                icc_genome   = NA_real_,
                icc_block    = NA_real_,
                row.names    = fid
            ))
        }

        vc <- .extract_variance_components(fit)

        genome_var   <- vc$genome_var
        block_var    <- vc$block_var
        phi          <- vc$phi
        residual_var <- vc$residual_var

        block_contrib <- if (!is.na(block_var)) block_var else 0
        genome_contrib <- if (!is.na(genome_var)) genome_var else 0
        resid_contrib  <- if (!is.na(residual_var)) residual_var else 0

        total_var <- genome_contrib + block_contrib + resid_contrib
        icc_genome <- if (total_var > 0) genome_contrib / total_var else NA_real_
        icc_block  <- if (!is.na(block_var) && total_var > 0) {
            block_contrib / total_var
        } else NA_real_

        S4Vectors::DataFrame(
            genome_var   = genome_var,
            block_var    = block_var,
            phi          = phi,
            residual_var = residual_var,
            total_var    = if (total_var > 0) total_var else NA_real_,
            icc_genome   = icc_genome,
            icc_block    = icc_block,
            row.names    = fid
        )
    })

    out <- do.call(rbind, rows)
    feature_col <- as.data.frame(x)[[feature_id_col]]
    if (!is.null(feature_col)) {
        out <- cbind(
            S4Vectors::DataFrame(stats::setNames(list(feature_col), feature_id_col)),
            out
        )
    }
    S4Vectors::DataFrame(out, check.names = FALSE)
}

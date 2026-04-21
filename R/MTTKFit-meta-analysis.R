.normalize_meta_method <- function(method) {
    match.arg(method, c("random", "fixed"))
}

.normalize_meta_min_kos <- function(min_kos) {
    if (!is.numeric(min_kos) ||
        length(min_kos) != 1L ||
        is.na(min_kos) ||
        min_kos < 1) {
        stop("'minKOs' must be a single integer greater than or equal to 1.", call. = FALSE)
    }

    as.integer(min_kos)
}

.normalize_meta_status_filter <- function(status) {
    if (is.null(status)) {
        return(NULL)
    }

    status <- as.character(status)
    status <- status[!is.na(status) & status != ""]

    if (length(status) == 0L) {
        stop("'status' must be NULL or a non-empty character vector.", call. = FALSE)
    }

    unique(status)
}

.normalize_meta_collapse <- function(collapse) {
    if (!is.character(collapse) || length(collapse) != 1L || is.na(collapse)) {
        stop("'collapse' must be a single character string.", call. = FALSE)
    }

    collapse
}

.normalize_meta_min_genomes <- function(min_genomes) {
    if (!is.numeric(min_genomes) ||
        length(min_genomes) != 1L ||
        is.na(min_genomes) ||
        min_genomes < 1) {
        stop("'minGenomes' must be a single integer greater than or equal to 1.", call. = FALSE)
    }

    as.integer(min_genomes)
}

.normalize_meta_membership_mode <- function(membership_mode) {
    match.arg(membership_mode, c("duplicate", "split", "exclusive"))
}

.require_ko_fit <- function(x) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    info <- fitInfo(x)
    feature_id_column <- info$featureIdColumn

    if (is.null(feature_id_column) || length(feature_id_column) != 1L ||
        is.na(feature_id_column) || feature_id_column == "") {
        if ("ko_id" %in% names(x)) {
            feature_id_column <- "ko_id"
        } else if (
            identical(info$model, "ko_mixed_model") ||
                identical(info$featureType, "KO") ||
                identical(info$sourceFeatureType, "KO")
        ) {
            feature_id_column <- "ko_id"
            x$ko_id <- rownames(x)
        } else {
            feature_id_column <- NA_character_
        }
    }

    if (!identical(feature_id_column, "ko_id")) {
        stop(
            "'x' must represent KO-level results, with canonical identifiers stored as 'ko_id'.",
            call. = FALSE
        )
    }

    required_cols <- c("estimate", "std_error")
    missing_cols <- setdiff(required_cols, names(x))
    if (length(missing_cols) > 0L) {
        stop(
            "The KO fit is missing required column(s): ",
            paste(missing_cols, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    S4Vectors::DataFrame(x, check.names = FALSE)
}

.require_genome_fit <- function(x) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    info <- fitInfo(x)
    feature_id_column <- info$featureIdColumn

    if (is.null(feature_id_column) || length(feature_id_column) != 1L ||
        is.na(feature_id_column) || feature_id_column == "") {
        if ("genome_id" %in% names(x)) {
            feature_id_column <- "genome_id"
        } else if (
            identical(info$model, "genome_model") ||
                identical(info$featureType, "genome")
        ) {
            feature_id_column <- "genome_id"
            x$genome_id <- rownames(x)
        } else {
            feature_id_column <- NA_character_
        }
    }

    if (!identical(feature_id_column, "genome_id")) {
        stop(
            "'x' must represent genome-level results, with canonical identifiers stored as 'genome_id'.",
            call. = FALSE
        )
    }

    required_cols <- c("estimate", "std_error")
    missing_cols <- setdiff(required_cols, names(x))
    if (length(missing_cols) > 0L) {
        stop(
            "The genome fit is missing required column(s): ",
            paste(missing_cols, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    S4Vectors::DataFrame(x, check.names = FALSE)
}

.ko_membership_link_df <- function(object, link_name, target_column, membership_mode) {
    membership_mode <- .normalize_meta_membership_mode(membership_mode)

    if (!methods::is(object, "MTTKExperiment")) {
        stop("'object' must be an MTTKExperiment.", call. = FALSE)
    }

    link_list <- links(object)
    if (!(link_name %in% names(link_list))) {
        stop("Unknown link table: ", link_name, ".", call. = FALSE)
    }

    mapping <- as.data.frame(link_list[[link_name]])
    required_cols <- c("ko_id", target_column)

    if (!all(required_cols %in% names(mapping))) {
        stop(
            "Link table '",
            link_name,
            "' must contain '",
            paste(required_cols, collapse = "' and '"),
            "' columns.",
            call. = FALSE
        )
    }

    mapping <- unique(mapping[, required_cols, drop = FALSE])
    mapping$ko_id <- as.character(mapping$ko_id)
    mapping[[target_column]] <- as.character(mapping[[target_column]])
    mapping <- mapping[
        !is.na(mapping$ko_id) & mapping$ko_id != "" &
            !is.na(mapping[[target_column]]) & mapping[[target_column]] != "",
        ,
        drop = FALSE
    ]

    if (nrow(mapping) == 0L) {
        stop("The link table '", link_name, "' does not contain any usable KO memberships.", call. = FALSE)
    }

    membership_count <- vapply(
        split(mapping[[target_column]], mapping$ko_id),
        function(ids) length(unique(ids)),
        integer(1)
    )
    mapping$membership_count <- as.integer(membership_count[mapping$ko_id])

    if (identical(membership_mode, "duplicate")) {
        mapping$membership_weight <- 1
    } else if (identical(membership_mode, "split")) {
        mapping$membership_weight <- 1 / mapping$membership_count
    } else {
        mapping <- mapping[mapping$membership_count == 1L, , drop = FALSE]
    }

    if (nrow(mapping) == 0L) {
        stop(
            "No usable KO memberships remained after applying membershipMode = '",
            membership_mode,
            "'.",
            call. = FALSE
        )
    }

    if (!("membership_weight" %in% names(mapping))) {
        mapping$membership_weight <- 1
    }

    mapping
}

.collapse_unique_ids <- function(values, collapse) {
    values <- unique(as.character(values))
    values <- values[!is.na(values) & values != ""]

    if (length(values) == 0L) {
        return(NA_character_)
    }

    paste(values, collapse = collapse)
}

.meta_tested_term <- function(ko_results) {
    if (!("tested_term" %in% names(ko_results))) {
        return(NA_character_)
    }

    tested_terms <- unique(as.character(ko_results$tested_term))
    tested_terms <- tested_terms[!is.na(tested_terms) & tested_terms != ""]

    if (length(tested_terms) == 0L) {
        return(NA_character_)
    }

    if (length(tested_terms) == 1L) {
        return(tested_terms)
    }

    paste(tested_terms, collapse = ";")
}

.meta_variable_info <- function(ko_fit) {
    info <- fitInfo(ko_fit)

    list(
        testedTerm = .meta_tested_term(S4Vectors::DataFrame(ko_fit, check.names = FALSE)),
        type = if (!is.null(info$variableType)) as.character(info$variableType) else NA_character_,
        referenceLevel = if (!is.null(info$referenceLevel)) as.character(info$referenceLevel) else NA_character_,
        contrastLevel = if (!is.null(info$contrastLevel)) as.character(info$contrastLevel) else NA_character_,
        effectLabel = if (!is.null(info$effectLabel)) as.character(info$effectLabel) else NA_character_,
        variable = if (!is.null(info$variable)) as.character(info$variable) else NA_character_
    )
}

.empty_meta_row <- function(
    feature_id,
    feature_id_column,
    variable_info,
    meta_method,
    membership_mode,
    mapped_membership,
    collapse
) {
    mapped_kos <- mapped_membership$ko_id
    row <- S4Vectors::DataFrame(
        tested_term = variable_info$testedTerm,
        variable_type = variable_info$type,
        reference_level = variable_info$referenceLevel,
        contrast_level = variable_info$contrastLevel,
        effect_label = variable_info$effectLabel,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        q_value = NA_real_,
        n_kos_mapped = as.integer(length(mapped_kos)),
        n_kos_available = NA_integer_,
        n_kos_tested = NA_integer_,
        member_ko_ids = .collapse_unique_ids(mapped_kos, collapse = collapse),
        tested_ko_ids = NA_character_,
        mapped_weight_sum = sum(as.numeric(mapped_membership$membership_weight)),
        tested_weight_sum = NA_real_,
        n_positive_kos = NA_integer_,
        n_negative_kos = NA_integer_,
        Q = NA_real_,
        Q_df = NA_integer_,
        Q_p_value = NA_real_,
        I2 = NA_real_,
        tau2 = NA_real_,
        meta_method = meta_method,
        membership_mode = membership_mode,
        warning_message = NA_character_,
        error_message = NA_character_,
        status = NA_character_,
        row.names = feature_id
    )
    row[[feature_id_column]] <- feature_id
    row <- row[, c(
        feature_id_column,
        setdiff(names(row), feature_id_column)
    )]
    row
}

.compute_inverse_variance_meta <- function(estimates, standard_errors, method, membership_weights) {
    weights_fixed <- membership_weights / (standard_errors^2)
    fixed_estimate <- sum(weights_fixed * estimates) / sum(weights_fixed)
    Q <- sum(weights_fixed * (estimates - fixed_estimate)^2)
    q_df <- length(estimates) - 1L

    if (identical(method, "random") && length(estimates) > 1L) {
        c_val <- sum(weights_fixed) - (sum(weights_fixed^2) / sum(weights_fixed))
        tau2 <- if (c_val > 0) max((Q - q_df) / c_val, 0) else 0
        weights <- membership_weights / (standard_errors^2 + tau2)
    } else {
        tau2 <- 0
        weights <- weights_fixed
    }

    estimate <- sum(weights * estimates) / sum(weights)
    std_error <- sqrt(1 / sum(weights))
    statistic <- estimate / std_error
    p_value <- 2 * stats::pnorm(-abs(statistic))

    list(
        estimate = estimate,
        std_error = std_error,
        statistic = statistic,
        p_value = p_value,
        Q = Q,
        Q_df = as.integer(q_df),
        Q_p_value = if (q_df > 0L) stats::pchisq(Q, df = q_df, lower.tail = FALSE) else NA_real_,
        I2 = if (q_df > 0L && Q > 0) max((Q - q_df) / Q, 0) * 100 else NA_real_,
        tau2 = tau2
    )
}

.fit_one_ko_meta_group <- function(
    feature_id,
    mapped_membership,
    ko_results,
    feature_id_column,
    variable_info,
    meta_method,
    membership_mode,
    min_kos,
    status_filter,
    collapse
) {
    row <- .empty_meta_row(
        feature_id = feature_id,
        feature_id_column = feature_id_column,
        variable_info = variable_info,
        meta_method = meta_method,
        membership_mode = membership_mode,
        mapped_membership = mapped_membership,
        collapse = collapse
    )

    matched <- match(mapped_membership$ko_id, rownames(ko_results))
    available <- ko_results[matched[!is.na(matched)], , drop = FALSE]
    available$membership_weight <- as.numeric(mapped_membership$membership_weight[!is.na(matched)])
    row$n_kos_available <- as.integer(nrow(available))

    if (!is.null(status_filter)) {
        if (!("status" %in% names(available))) {
            stop("The KO fit does not contain a 'status' column.", call. = FALSE)
        }

        available <- available[as.character(available$status) %in% status_filter, , drop = FALSE]
    }

    valid <- available[
        !is.na(available$estimate) &
            !is.na(available$std_error) &
            is.finite(available$estimate) &
            is.finite(available$std_error) &
            as.numeric(available$std_error) > 0,
        ,
        drop = FALSE
    ]

    row$n_kos_tested <- as.integer(nrow(valid))
    row$tested_ko_ids <- .collapse_unique_ids(rownames(valid), collapse = collapse)
    row$tested_weight_sum <- sum(as.numeric(valid$membership_weight))
    row$n_positive_kos <- as.integer(sum(valid$estimate > 0, na.rm = TRUE))
    row$n_negative_kos <- as.integer(sum(valid$estimate < 0, na.rm = TRUE))

    dropped <- row$n_kos_available - row$n_kos_tested
    if (dropped > 0L) {
        row$warning_message <- paste0(
            dropped,
            " KO(s) were excluded because they did not pass the status or standard-error filters."
        )
    }

    if (row$n_kos_available == 0L) {
        row$status <- "skipped"
        row$error_message <- "No KO-level fit rows were available for this feature."
        return(row)
    }

    if (row$n_kos_tested < min_kos) {
        row$status <- "skipped"
        row$error_message <- paste0(
            "At least ",
            min_kos,
            " KO(s) with finite estimates and positive standard errors are required."
        )
        return(row)
    }

    meta <- .compute_inverse_variance_meta(
        estimates = as.numeric(valid$estimate),
        standard_errors = as.numeric(valid$std_error),
        method = meta_method,
        membership_weights = as.numeric(valid$membership_weight)
    )

    row$estimate <- meta$estimate
    row$std_error <- meta$std_error
    row$statistic <- meta$statistic
    row$p_value <- meta$p_value
    row$Q <- meta$Q
    row$Q_df <- meta$Q_df
    row$Q_p_value <- meta$Q_p_value
    row$I2 <- meta$I2
    row$tau2 <- meta$tau2
    row$status <- "ok"
    row$error_message <- NA_character_
    row
}

.fit_ko_meta_analysis <- function(
    x,
    object,
    link_name,
    feature_id_column,
    feature_label,
    method,
    membershipMode,
    minKOs,
    status,
    collapse
) {
    ko_results <- .require_ko_fit(x)
    meta_method <- .normalize_meta_method(method)
    membership_mode <- .normalize_meta_membership_mode(membershipMode)
    min_kos <- .normalize_meta_min_kos(minKOs)
    status_filter <- .normalize_meta_status_filter(status)
    collapse <- .normalize_meta_collapse(collapse)
    mapping <- .ko_membership_link_df(
        object,
        link_name = link_name,
        target_column = feature_id_column,
        membership_mode = membership_mode
    )

    memberships <- split(mapping, mapping[[feature_id_column]])
    feature_ids <- names(memberships)
    variable_info <- .meta_variable_info(x)

    results <- do.call(
        rbind,
        lapply(feature_ids, function(feature_id) {
            .fit_one_ko_meta_group(
                feature_id = feature_id,
                mapped_membership = memberships[[feature_id]],
                ko_results = ko_results,
                feature_id_column = feature_id_column,
                variable_info = variable_info,
                meta_method = meta_method,
                membership_mode = membership_mode,
                min_kos = min_kos,
                status_filter = status_filter,
                collapse = collapse
            )
        })
    )
    rownames(results) <- feature_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    source_info <- fitInfo(x)
    info <- list(
        backend = "meta_analysis",
        model = paste0(feature_label, "_meta_analysis"),
        featureType = feature_label,
        featureIdColumn = feature_id_column,
        sourceFeatureType = "KO",
        sourceFeatureIdColumn = "ko_id",
        sourceFitModel = if (!is.null(source_info$model)) source_info$model else NA_character_,
        sourceBackend = if (!is.null(source_info$backend)) source_info$backend else NA_character_,
        sourceRandomEffects = if (!is.null(source_info$randomEffects)) {
            source_info$randomEffects
        } else {
            NA_character_
        },
        sourceGroupEffectColumn = if (!is.null(source_info$groupEffectColumn)) {
            source_info$groupEffectColumn
        } else {
            NA_character_
        },
        variable = variable_info$variable,
        variableType = variable_info$type,
        referenceLevel = variable_info$referenceLevel,
        contrastLevel = variable_info$contrastLevel,
        effectLabel = variable_info$effectLabel,
        sourceFormula = if (!is.null(source_info$formula)) source_info$formula else NA_character_,
        membershipLink = link_name,
        metaMethod = meta_method,
        membershipMode = membership_mode,
        tau2Estimator = if (identical(meta_method, "random")) "DerSimonian-Laird" else NA_character_,
        minKOs = min_kos,
        statusFilter = if (is.null(status_filter)) NA_character_ else status_filter,
        n_features = nrow(results),
        n_ok = sum(results$status == "ok"),
        n_skipped = sum(results$status == "skipped"),
        n_error = sum(results$status == "error")
    )

    MTTKFit(
        results = results,
        info = info,
        models = list()
    )
}

.genome_group_mapping_df <- function(object, group_column) {
    if (!methods::is(object, "MTTKExperiment")) {
        stop("'object' must be an MTTKExperiment.", call. = FALSE)
    }

    if (!is.character(group_column) ||
        length(group_column) != 1L ||
        is.na(group_column) ||
        group_column == "") {
        stop("'group' must be a single non-empty genomeData column name.", call. = FALSE)
    }

    genome_data <- genomeData(object)
    if (!(group_column %in% names(genome_data))) {
        stop("Unknown genome grouping column: ", group_column, ".", call. = FALSE)
    }

    mapping <- data.frame(
        genome_id = rownames(genome_data),
        group_id = as.character(S4Vectors::decode(genome_data[[group_column]])),
        stringsAsFactors = FALSE
    )
    mapping <- mapping[
        !is.na(mapping$genome_id) & mapping$genome_id != "" &
            !is.na(mapping$group_id) & mapping$group_id != "",
        ,
        drop = FALSE
    ]

    if (nrow(mapping) == 0L) {
        stop(
            "No usable genome group memberships were found in genomeData(x)$",
            group_column,
            ".",
            call. = FALSE
        )
    }

    mapping
}

.empty_genome_group_meta_row <- function(
    group_id,
    group_column,
    variable_info,
    meta_method,
    mapped_membership,
    collapse
) {
    mapped_genomes <- mapped_membership$genome_id
    row <- S4Vectors::DataFrame(
        tested_term = variable_info$testedTerm,
        variable_type = variable_info$type,
        reference_level = variable_info$referenceLevel,
        contrast_level = variable_info$contrastLevel,
        effect_label = variable_info$effectLabel,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        q_value = NA_real_,
        n_genomes_mapped = as.integer(length(mapped_genomes)),
        n_genomes_available = NA_integer_,
        n_genomes_tested = NA_integer_,
        member_genome_ids = .collapse_unique_ids(mapped_genomes, collapse = collapse),
        tested_genome_ids = NA_character_,
        n_positive_genomes = NA_integer_,
        n_negative_genomes = NA_integer_,
        Q = NA_real_,
        Q_df = NA_integer_,
        Q_p_value = NA_real_,
        I2 = NA_real_,
        tau2 = NA_real_,
        meta_method = meta_method,
        warning_message = NA_character_,
        error_message = NA_character_,
        status = NA_character_,
        row.names = group_id
    )
    row[[group_column]] <- group_id
    row <- row[, c(group_column, setdiff(names(row), group_column))]
    row
}

.fit_one_genome_meta_group <- function(
    group_id,
    mapped_membership,
    genome_results,
    group_column,
    variable_info,
    meta_method,
    min_genomes,
    status_filter,
    collapse
) {
    row <- .empty_genome_group_meta_row(
        group_id = group_id,
        group_column = group_column,
        variable_info = variable_info,
        meta_method = meta_method,
        mapped_membership = mapped_membership,
        collapse = collapse
    )

    matched <- match(mapped_membership$genome_id, rownames(genome_results))
    available <- genome_results[matched[!is.na(matched)], , drop = FALSE]
    row$n_genomes_available <- as.integer(nrow(available))

    if (!is.null(status_filter)) {
        if (!("status" %in% names(available))) {
            stop("The genome fit does not contain a 'status' column.", call. = FALSE)
        }

        available <- available[as.character(available$status) %in% status_filter, , drop = FALSE]
    }

    valid <- available[
        !is.na(available$estimate) &
            !is.na(available$std_error) &
            is.finite(available$estimate) &
            is.finite(available$std_error) &
            as.numeric(available$std_error) > 0,
        ,
        drop = FALSE
    ]

    row$n_genomes_tested <- as.integer(nrow(valid))
    row$tested_genome_ids <- .collapse_unique_ids(rownames(valid), collapse = collapse)
    row$n_positive_genomes <- as.integer(sum(valid$estimate > 0, na.rm = TRUE))
    row$n_negative_genomes <- as.integer(sum(valid$estimate < 0, na.rm = TRUE))

    dropped <- row$n_genomes_available - row$n_genomes_tested
    if (dropped > 0L) {
        row$warning_message <- paste0(
            dropped,
            " genome(s) were excluded because they did not pass the status or standard-error filters."
        )
    }

    if (row$n_genomes_available == 0L) {
        row$status <- "skipped"
        row$error_message <- "No genome-level fit rows were available for this group."
        return(row)
    }

    if (row$n_genomes_tested < min_genomes) {
        row$status <- "skipped"
        row$error_message <- paste0(
            "At least ",
            min_genomes,
            " genome(s) with finite estimates and positive standard errors are required."
        )
        return(row)
    }

    meta <- .compute_inverse_variance_meta(
        estimates = as.numeric(valid$estimate),
        standard_errors = as.numeric(valid$std_error),
        method = meta_method,
        membership_weights = rep(1, nrow(valid))
    )

    row$estimate <- meta$estimate
    row$std_error <- meta$std_error
    row$statistic <- meta$statistic
    row$p_value <- meta$p_value
    row$Q <- meta$Q
    row$Q_df <- meta$Q_df
    row$Q_p_value <- meta$Q_p_value
    row$I2 <- meta$I2
    row$tau2 <- meta$tau2
    row$status <- "ok"
    row$error_message <- NA_character_
    row
}

#' Meta-Analyze KO Effects to the Module Level
#'
#' `fitModuleMetaAnalysis()` summarizes KO-level effects to module-level
#' associations without refitting the count model. This workflow is intended for
#' the question: do the KOs assigned to a KEGG module show a coherent
#' condition-associated shift?
#'
#' The input `x` must be a KO-level `MTTKFit`, typically returned by
#' [fitKOMixedModel()] or [fitKORandomSlopeModel()]. KO estimates and standard
#' errors are combined within each module using inverse-variance meta-analysis.
#' By default a random-effects model with a DerSimonian-Laird `tau^2` estimate
#' is used.
#'
#' Because KEGG module membership is many-to-many, the same KO can contribute
#' its KO-level effect to multiple modules. This workflow therefore answers a
#' different question from [fitModuleMixedModel()], which re-aggregates counts
#' and tests shifts in total module activity.
#'
#' `membershipMode = "duplicate"` counts each KO fully in every mapped module
#' and is appropriate when the question is annotation-centric: if a KO belongs
#' to a module, should it contribute fully to that module's synthesis?
#' `membershipMode = "split"` divides each KO's contribution across its
#' memberships, which is useful when broad or promiscuous KOs would otherwise
#' dominate overlapping sets.
#' `membershipMode = "exclusive"` keeps only uniquely assigned KOs, which is
#' the most specific but also the most conservative option.
#'
#' The current implementation treats KO-level effects as approximately
#' independent, so higher-level p-values should be interpreted as an
#' approximate synthesis layer rather than as a fully correlation-aware test.
#'
#' @param x A KO-level `MTTKFit`, typically from [fitKOMixedModel()] or
#'   [fitKORandomSlopeModel()].
#' @param object An `MTTKExperiment` supplying `ko_to_module` memberships.
#' @param method Meta-analysis method. `"random"` uses inverse-variance
#'   random-effects meta-analysis with a DerSimonian-Laird `tau^2` estimate,
#'   while `"fixed"` uses fixed-effect inverse-variance pooling.
#' @param membershipMode How to handle many-to-many KO memberships.
#'   `"duplicate"` carries each KO fully into every mapped set,
#'   `"split"` divides each KO's contribution across its memberships, and
#'   `"exclusive"` keeps only uniquely assigned KOs.
#' @param minKOs Minimum number of tested KO-level effects required for a module
#'   to receive an `"ok"` result. Modules with fewer tested KOs are returned
#'   with `status = "skipped"`.
#' @param status Optional KO-fit status filter applied before the meta-analysis.
#'   The default keeps only KO rows with `status == "ok"`. Use `NULL` to keep
#'   all KO rows.
#' @param collapse String used to combine KO identifiers in the
#'   `member_ko_ids` and `tested_ko_ids` columns.
#'
#' @return An `MTTKFit` with one row per module.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     ko_fit <- fitKOMixedModel(x, variable = "condition")
#'     module_fit <- fitModuleMetaAnalysis(ko_fit, x)
#'     as.data.frame(module_fit)[, c("module_id", "estimate", "q_value", "n_kos_tested")]
#' }
#'
#' @export
fitModuleMetaAnalysis <- function(
    x,
    object,
    method = c("random", "fixed"),
    membershipMode = c("duplicate", "split", "exclusive"),
    minKOs = 1L,
    status = "ok",
    collapse = ";"
) {
    .fit_ko_meta_analysis(
        x = x,
        object = object,
        link_name = "ko_to_module",
        feature_id_column = "module_id",
        feature_label = "module",
        method = match.arg(method),
        membershipMode = match.arg(membershipMode),
        minKOs = minKOs,
        status = status,
        collapse = collapse
    )
}

#' Meta-Analyze KO Effects to the Pathway Level
#'
#' `fitPathwayMetaAnalysis()` summarizes KO-level effects to pathway-level
#' associations without refitting the count model. This workflow is intended for
#' the question: do the KOs assigned to a KEGG pathway show a coherent
#' condition-associated shift?
#'
#' The input `x` must be a KO-level `MTTKFit`, typically returned by
#' [fitKOMixedModel()] or [fitKORandomSlopeModel()]. KO estimates and standard
#' errors are combined within each pathway using inverse-variance meta-analysis.
#' By default a random-effects model with a DerSimonian-Laird `tau^2` estimate
#' is used.
#'
#' Because KEGG pathway membership is many-to-many, the same KO can contribute
#' its KO-level effect to multiple pathways. This workflow therefore answers a
#' different question from [fitPathwayMixedModel()], which re-aggregates counts
#' and tests shifts in total pathway activity.
#'
#' The interpretation of `membershipMode` is the same as in
#' [fitModuleMetaAnalysis()]: use `"duplicate"` for annotation-centric
#' synthesis, `"split"` to divide each KO contribution across memberships, and
#' `"exclusive"` when only uniquely assigned KOs should contribute.
#'
#' The current implementation treats KO-level effects as approximately
#' independent, so higher-level p-values should be interpreted as an
#' approximate synthesis layer rather than as a fully correlation-aware test.
#'
#' @inheritParams fitModuleMetaAnalysis
#'
#' @return An `MTTKFit` with one row per pathway.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     ko_fit <- fitKOMixedModel(x, variable = "condition")
#'     pathway_fit <- fitPathwayMetaAnalysis(ko_fit, x)
#'     as.data.frame(pathway_fit)[, c("pathway_id", "estimate", "q_value", "n_kos_tested")]
#' }
#'
#' @export
fitPathwayMetaAnalysis <- function(
    x,
    object,
    method = c("random", "fixed"),
    membershipMode = c("duplicate", "split", "exclusive"),
    minKOs = 1L,
    status = "ok",
    collapse = ";"
) {
    .fit_ko_meta_analysis(
        x = x,
        object = object,
        link_name = "ko_to_pathway",
        feature_id_column = "pathway_id",
        feature_label = "pathway",
        method = match.arg(method),
        membershipMode = match.arg(membershipMode),
        minKOs = minKOs,
        status = status,
        collapse = collapse
    )
}

#' Meta-Analyze Genome Effects by a Genome Metadata Group
#'
#' `fitGenomeGroupMetaAnalysis()` summarizes genome-level effects to a grouping
#' variable stored in `genomeData(x)`, such as `clade`, `domain`, or another
#' taxonomic/phylogenetic grouping column.
#'
#' The input `x` must be a genome-level `MTTKFit`, typically returned by
#' [fitGenomeModel()]. Genome estimates and standard errors are combined within
#' each selected group using inverse-variance meta-analysis. By default a
#' random-effects model with a DerSimonian-Laird `tau^2` estimate is used.
#'
#' This workflow is intended for the question: do the genomes assigned to a
#' given clade or taxonomic group show a coherent response, and how
#' heterogeneous are those genome-level responses?
#'
#' @param x A genome-level `MTTKFit`, typically from [fitGenomeModel()].
#' @param object An `MTTKExperiment` supplying genome metadata.
#' @param group A single column name from `genomeData(object)` used to group
#'   genomes, such as `"clade"` or `"domain"`.
#' @param method Meta-analysis method. `"random"` uses inverse-variance
#'   random-effects meta-analysis with a DerSimonian-Laird `tau^2` estimate,
#'   while `"fixed"` uses fixed-effect inverse-variance pooling.
#' @param minGenomes Minimum number of tested genome-level effects required for
#'   a group to receive an `"ok"` result. Groups with fewer tested genomes are
#'   returned with `status = "skipped"`.
#' @param status Optional genome-fit status filter applied before the
#'   meta-analysis. The default keeps only genome rows with `status == "ok"`.
#'   Use `NULL` to keep all genome rows.
#' @param collapse String used to combine genome identifiers in the
#'   `member_genome_ids` and `tested_genome_ids` columns.
#'
#' @return An `MTTKFit` with one row per selected genome group.
#'
#' @examples
#' fit <- MTTKFit(
#'     results = data.frame(
#'         genome_id = c("genome_1", "genome_2", "genome_3"),
#'         estimate = c(0.5, 0.4, -0.2),
#'         std_error = c(0.2, 0.25, 0.3),
#'         status = c("ok", "ok", "ok"),
#'         row.names = c("genome_1", "genome_2", "genome_3")
#'     ),
#'     info = list(
#'         model = "genome_model",
#'         featureIdColumn = "genome_id",
#'         variable = "condition",
#'         variableType = "two_level_factor",
#'         referenceLevel = "control",
#'         contrastLevel = "treated",
#'         effectLabel = "condition: treated vs control"
#'     )
#' )
#'
#' x <- makeExampleMTTKExperiment()
#' fitGenomeGroupMetaAnalysis(fit, x, group = "domain")
#'
#' @export
fitGenomeGroupMetaAnalysis <- function(
    x,
    object,
    group,
    method = c("random", "fixed"),
    minGenomes = 1L,
    status = "ok",
    collapse = ";"
) {
    genome_results <- .require_genome_fit(x)
    meta_method <- .normalize_meta_method(method)
    min_genomes <- .normalize_meta_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    collapse <- .normalize_meta_collapse(collapse)
    mapping <- .genome_group_mapping_df(object, group_column = group)

    memberships <- split(mapping, mapping$group_id)
    group_ids <- names(memberships)
    variable_info <- .meta_variable_info(x)

    results <- do.call(
        rbind,
        lapply(group_ids, function(group_id) {
            .fit_one_genome_meta_group(
                group_id = group_id,
                mapped_membership = memberships[[group_id]],
                genome_results = genome_results,
                group_column = group,
                variable_info = variable_info,
                meta_method = meta_method,
                min_genomes = min_genomes,
                status_filter = status_filter,
                collapse = collapse
            )
        })
    )
    rownames(results) <- group_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    source_info <- fitInfo(x)
    info <- list(
        backend = "meta_analysis",
        model = "genome_group_meta_analysis",
        featureType = paste0("genome_group:", group),
        featureIdColumn = group,
        sourceFeatureType = "genome",
        sourceFeatureIdColumn = "genome_id",
        sourceFitModel = if (!is.null(source_info$model)) source_info$model else NA_character_,
        sourceBackend = if (!is.null(source_info$backend)) source_info$backend else NA_character_,
        variable = variable_info$variable,
        variableType = variable_info$type,
        referenceLevel = variable_info$referenceLevel,
        contrastLevel = variable_info$contrastLevel,
        effectLabel = variable_info$effectLabel,
        sourceFormula = if (!is.null(source_info$formula)) source_info$formula else NA_character_,
        groupColumn = group,
        metaMethod = meta_method,
        tau2Estimator = if (identical(meta_method, "random")) "DerSimonian-Laird" else NA_character_,
        minGenomes = min_genomes,
        statusFilter = if (is.null(status_filter)) NA_character_ else status_filter,
        n_features = nrow(results),
        n_ok = sum(results$status == "ok"),
        n_skipped = sum(results$status == "skipped"),
        n_error = sum(results$status == "error")
    )

    MTTKFit(
        results = results,
        info = info,
        models = list()
    )
}

.require_ko_random_slope_fit <- function(x) {
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

    ko_results <- .require_ko_fit(x)
    group_effects <- koGenomeEffects(x)
    if (is.data.frame(group_effects) && !methods::is(group_effects, "DataFrame")) {
        group_effects <- S4Vectors::DataFrame(group_effects, check.names = FALSE)
    }

    required_cols <- c(
        "ko_id",
        "genome_id",
        "conditional_effect_estimate",
        "conditional_effect_std_error"
    )
    missing_cols <- setdiff(required_cols, names(group_effects))
    if (length(missing_cols) > 0L) {
        stop(
            "The KO random-slope fit is missing required group-effect column(s): ",
            paste(missing_cols, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    list(
        koResults = ko_results,
        groupEffects = group_effects
    )
}

.empty_ko_genome_group_meta_row <- function(
    ko_id,
    group_id,
    group_column,
    variable_info,
    meta_method,
    collapse
) {
    ko_group_id <- paste(ko_id, group_id, sep = "::")
    row <- S4Vectors::DataFrame(
        tested_term = variable_info$testedTerm,
        variable_type = variable_info$type,
        reference_level = variable_info$referenceLevel,
        contrast_level = variable_info$contrastLevel,
        effect_label = variable_info$effectLabel,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        q_value = NA_real_,
        n_genomes_mapped = NA_integer_,
        n_genomes_available = NA_integer_,
        n_genomes_tested = NA_integer_,
        member_genome_ids = NA_character_,
        tested_genome_ids = NA_character_,
        n_positive_genomes = NA_integer_,
        n_negative_genomes = NA_integer_,
        Q = NA_real_,
        Q_df = NA_integer_,
        Q_p_value = NA_real_,
        I2 = NA_real_,
        tau2 = NA_real_,
        meta_method = meta_method,
        effect_source = "conditional_effect_estimate",
        std_error_source = "conditional_effect_std_error",
        warning_message = NA_character_,
        error_message = NA_character_,
        status = NA_character_,
        row.names = ko_group_id
    )
    row$ko_group_id <- ko_group_id
    row$ko_id <- ko_id
    row[[group_column]] <- group_id
    row <- row[, c(
        "ko_group_id",
        "ko_id",
        group_column,
        setdiff(names(row), c("ko_group_id", "ko_id", group_column))
    )]
    row
}

.fit_one_ko_genome_group <- function(
    ko_id,
    group_id,
    data,
    group_column,
    variable_info,
    meta_method,
    min_genomes,
    status_filter,
    collapse
) {
    row <- .empty_ko_genome_group_meta_row(
        ko_id = ko_id,
        group_id = group_id,
        group_column = group_column,
        variable_info = variable_info,
        meta_method = meta_method,
        collapse = collapse
    )

    row$n_genomes_mapped <- as.integer(length(unique(as.character(data$genome_id))))
    row$member_genome_ids <- .collapse_unique_ids(as.character(data$genome_id), collapse = collapse)
    row$n_genomes_available <- as.integer(row$n_genomes_mapped)

    if (!is.null(status_filter)) {
        if (!("source_fit_status" %in% names(data))) {
            stop(
                "The KO genome-effect table does not contain a 'source_fit_status' column.",
                call. = FALSE
            )
        }

        data <- data[as.character(data$source_fit_status) %in% status_filter, , drop = FALSE]
    }

    valid <- data[
        !is.na(data$conditional_effect_estimate) &
            !is.na(data$conditional_effect_std_error) &
            is.finite(data$conditional_effect_estimate) &
            is.finite(data$conditional_effect_std_error) &
            as.numeric(data$conditional_effect_std_error) > 0,
        ,
        drop = FALSE
    ]

    row$n_genomes_tested <- as.integer(length(unique(as.character(valid$genome_id))))
    row$tested_genome_ids <- .collapse_unique_ids(as.character(valid$genome_id), collapse = collapse)
    row$n_positive_genomes <- as.integer(sum(valid$conditional_effect_estimate > 0, na.rm = TRUE))
    row$n_negative_genomes <- as.integer(sum(valid$conditional_effect_estimate < 0, na.rm = TRUE))

    dropped <- row$n_genomes_available - row$n_genomes_tested
    if (dropped > 0L) {
        row$warning_message <- paste0(
            dropped,
            " genome(s) were excluded because they did not pass the status or standard-error filters."
        )
    }

    if (row$n_genomes_available == 0L) {
        row$status <- "skipped"
        row$error_message <- "No genome-specific KO effects were available for this KO/group combination."
        return(row)
    }

    if (row$n_genomes_tested < min_genomes) {
        row$status <- "skipped"
        row$error_message <- paste0(
            "At least ",
            min_genomes,
            " genome-specific KO effect(s) with finite estimates and positive standard errors are required."
        )
        return(row)
    }

    meta <- .compute_inverse_variance_meta(
        estimates = as.numeric(valid$conditional_effect_estimate),
        standard_errors = as.numeric(valid$conditional_effect_std_error),
        method = meta_method,
        membership_weights = rep(1, nrow(valid))
    )

    row$estimate <- meta$estimate
    row$std_error <- meta$std_error
    row$statistic <- meta$statistic
    row$p_value <- meta$p_value
    row$Q <- meta$Q
    row$Q_df <- meta$Q_df
    row$Q_p_value <- meta$Q_p_value
    row$I2 <- meta$I2
    row$tau2 <- meta$tau2
    row$status <- "ok"
    row$error_message <- NA_character_
    row
}

#' Meta-Analyze Genome-Specific KO Effects by a Genome Metadata Group
#'
#' `fitKOGenomeGroupMetaAnalysis()` follows a KO random-slope model with a
#' group-level synthesis step. It starts from the genome-specific KO effects
#' extracted by [koGenomeEffects()] and summarizes those effects within a genome
#' grouping variable stored in `genomeData(x)`, such as `domain`, `clade`, or
#' another taxonomy/phylogeny-derived label.
#'
#' This workflow is intended for the question: for a given KO, do genomes in a
#' taxonomic or clade-defined group show a coherent response, and how
#' heterogeneous are those genome-specific KO responses within the group?
#'
#' The input `x` must be an `MTTKFit` returned by [fitKORandomSlopeModel()].
#' The meta-analysis operates on the genome-specific conditional KO effects
#' stored in that fit, not on the KO-wide fixed effect. By default a
#' random-effects inverse-variance meta-analysis with a DerSimonian-Laird
#' `tau^2` estimate is used.
#'
#' The per-genome standard errors are approximate. MTTK derives them from the
#' fixed-effect variance and the conditional random-slope variance, without
#' modeling the full cross-covariance between those components. As a result,
#' this workflow should be interpreted as an approximate synthesis layer for
#' taxonomy-aware KO follow-up rather than as fully exact hierarchical
#' inference.
#'
#' @param x A KO random-slope `MTTKFit`, from [fitKORandomSlopeModel()].
#' @param object An `MTTKExperiment` supplying genome metadata.
#' @param group A single column name from `genomeData(object)` used to group
#'   genomes, such as `"domain"` or `"clade"`.
#' @param method Meta-analysis method. `"random"` uses inverse-variance
#'   random-effects meta-analysis with a DerSimonian-Laird `tau^2` estimate,
#'   while `"fixed"` uses fixed-effect inverse-variance pooling.
#' @param minGenomes Minimum number of tested genome-specific KO effects
#'   required for a KO/group combination to receive an `"ok"` result.
#' @param status Optional KO-fit status filter applied before the meta-analysis.
#'   The default keeps only KO rows with `status == "ok"`. Use `NULL` to keep
#'   all KO rows.
#' @param collapse String used to combine genome identifiers in the
#'   `member_genome_ids` and `tested_genome_ids` columns.
#'
#' @return An `MTTKFit` with one row per KO/group combination.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     ko_fit <- fitKORandomSlopeModel(x, variable = "condition")
#'     ko_domain_fit <- fitKOGenomeGroupMetaAnalysis(ko_fit, x, group = "domain")
#'     fitTable(ko_domain_fit, sortBy = "q_value")
#' }
#'
#' @export
fitKOGenomeGroupMetaAnalysis <- function(
    x,
    object,
    group,
    method = c("random", "fixed"),
    minGenomes = 1L,
    status = "ok",
    collapse = ";"
) {
    parsed <- .require_ko_random_slope_fit(x)
    meta_method <- .normalize_meta_method(method)
    min_genomes <- .normalize_meta_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    collapse <- .normalize_meta_collapse(collapse)
    group_mapping <- .genome_group_mapping_df(object, group_column = group)
    variable_info <- .meta_variable_info(x)

    group_effects <- as.data.frame(parsed$groupEffects, stringsAsFactors = FALSE)
    group_effects$ko_id <- as.character(group_effects$ko_id)
    group_effects$genome_id <- as.character(group_effects$genome_id)
    ko_status <- data.frame(
        ko_id = rownames(parsed$koResults),
        source_fit_status = if ("status" %in% names(parsed$koResults)) {
            as.character(parsed$koResults$status)
        } else {
            NA_character_
        },
        stringsAsFactors = FALSE
    )

    group_effects <- merge(
        group_effects,
        ko_status,
        by = "ko_id",
        all.x = TRUE,
        sort = FALSE
    )
    group_effects <- merge(
        group_effects,
        group_mapping,
        by = "genome_id",
        all.x = FALSE,
        sort = FALSE
    )

    if (nrow(group_effects) == 0L) {
        stop(
            "No usable KO genome-effect rows could be aligned to genomeData(object)$",
            group,
            ".",
            call. = FALSE
        )
    }

    split_key <- paste(group_effects$ko_id, group_effects$group_id, sep = "::")
    memberships <- split(group_effects, split_key)
    feature_ids <- names(memberships)

    results <- do.call(
        rbind,
        lapply(feature_ids, function(feature_id) {
            ko_id <- memberships[[feature_id]]$ko_id[[1L]]
            group_id <- memberships[[feature_id]]$group_id[[1L]]
            .fit_one_ko_genome_group(
                ko_id = ko_id,
                group_id = group_id,
                data = memberships[[feature_id]],
                group_column = group,
                variable_info = variable_info,
                meta_method = meta_method,
                min_genomes = min_genomes,
                status_filter = status_filter,
                collapse = collapse
            )
        })
    )
    rownames(results) <- feature_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    source_info <- fitInfo(x)
    info <- list(
        backend = "meta_analysis",
        model = "ko_genome_group_meta_analysis",
        featureType = paste0("ko_group:", group),
        featureIdColumn = "ko_group_id",
        sourceFeatureType = "KO_genome_effect",
        sourceFeatureIdColumn = "ko_id",
        sourceFitModel = if (!is.null(source_info$model)) source_info$model else NA_character_,
        sourceBackend = if (!is.null(source_info$backend)) source_info$backend else NA_character_,
        sourceRandomEffects = if (!is.null(source_info$randomEffects)) {
            source_info$randomEffects
        } else {
            NA_character_
        },
        sourceGroupEffectColumn = if (!is.null(source_info$groupEffectColumn)) {
            source_info$groupEffectColumn
        } else {
            NA_character_
        },
        variable = variable_info$variable,
        variableType = variable_info$type,
        referenceLevel = variable_info$referenceLevel,
        contrastLevel = variable_info$contrastLevel,
        effectLabel = variable_info$effectLabel,
        groupColumn = group,
        sourceFormula = if (!is.null(source_info$formula)) source_info$formula else NA_character_,
        metaMethod = meta_method,
        tau2Estimator = if (identical(meta_method, "random")) "DerSimonian-Laird" else NA_character_,
        effectSource = "conditional_effect_estimate",
        stdErrorSource = "conditional_effect_std_error",
        stdErrorApproximation = paste(
            "sqrt(fixed_effect_variance + conditional_random_slope_variance)",
            "without fixed/random cross-covariance"
        ),
        minGenomes = min_genomes,
        statusFilter = if (is.null(status_filter)) NA_character_ else status_filter,
        n_features = nrow(results),
        n_ok = sum(results$status == "ok"),
        n_skipped = sum(results$status == "skipped"),
        n_error = sum(results$status == "error")
    )

    MTTKFit(
        results = results,
        info = info,
        models = list()
    )
}

.require_ko_genome_group_fit <- function(x) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    info <- fitInfo(x)
    if (!identical(info$model, "ko_genome_group_meta_analysis")) {
        stop(
            "'x' must be an MTTKFit returned by fitKOGenomeGroupMetaAnalysis().",
            call. = FALSE
        )
    }

    group_column <- if (!is.null(info$groupColumn)) as.character(info$groupColumn)[1L] else NA_character_
    if (is.na(group_column) || group_column == "") {
        stop(
            "The KO genome-group fit does not declare the source grouping column in fitInfo(x)$groupColumn.",
            call. = FALSE
        )
    }

    required_cols <- c("ko_id", group_column, "estimate", "std_error")
    missing_cols <- setdiff(required_cols, names(x))
    if (length(missing_cols) > 0L) {
        stop(
            "The KO genome-group fit is missing required column(s): ",
            paste(missing_cols, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    list(
        results = S4Vectors::DataFrame(x, check.names = FALSE),
        groupColumn = group_column,
        info = info
    )
}

.normalize_group_comparison_pairs <- function(group_values, reference, contrast) {
    groups <- unique(as.character(group_values))
    groups <- groups[!is.na(groups) & groups != ""]

    if (length(groups) < 2L) {
        stop("At least two groups are required to perform a comparison.", call. = FALSE)
    }

    if (is.null(reference) && is.null(contrast)) {
        pair_matrix <- utils::combn(sort(groups), 2L)
        return(data.frame(
            reference_group = pair_matrix[1L, ],
            contrast_group = pair_matrix[2L, ],
            stringsAsFactors = FALSE
        ))
    }

    if (!is.character(reference) || length(reference) != 1L || is.na(reference) || reference == "") {
        stop(
            "'reference' must be NULL or a single non-empty group identifier.",
            call. = FALSE
        )
    }

    if (!(reference %in% groups)) {
        stop("Unknown reference group: ", reference, ".", call. = FALSE)
    }

    if (is.null(contrast)) {
        contrast_groups <- setdiff(sort(groups), reference)
        return(data.frame(
            reference_group = rep(reference, length(contrast_groups)),
            contrast_group = contrast_groups,
            stringsAsFactors = FALSE
        ))
    }

    if (!is.character(contrast) || length(contrast) != 1L || is.na(contrast) || contrast == "") {
        stop(
            "'contrast' must be NULL or a single non-empty group identifier.",
            call. = FALSE
        )
    }

    if (!(contrast %in% groups)) {
        stop("Unknown contrast group: ", contrast, ".", call. = FALSE)
    }

    if (identical(reference, contrast)) {
        stop("'reference' and 'contrast' must differ.", call. = FALSE)
    }

    data.frame(
        reference_group = reference,
        contrast_group = contrast,
        stringsAsFactors = FALSE
    )
}

.empty_ko_group_comparison_row <- function(
    ko_id,
    group_column,
    reference_group,
    contrast_group,
    variable_info
) {
    comparison_label <- paste0(contrast_group, " vs ", reference_group)
    comparison_id <- paste(ko_id, contrast_group, reference_group, sep = "::")

    row <- S4Vectors::DataFrame(
        tested_term = variable_info$testedTerm,
        variable_type = variable_info$type,
        reference_level = variable_info$referenceLevel,
        contrast_level = variable_info$contrastLevel,
        effect_label = variable_info$effectLabel,
        group_comparison = comparison_label,
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        q_value = NA_real_,
        reference_group_estimate = NA_real_,
        reference_group_std_error = NA_real_,
        contrast_group_estimate = NA_real_,
        contrast_group_std_error = NA_real_,
        reference_group_status = NA_character_,
        contrast_group_status = NA_character_,
        reference_n_genomes_tested = NA_integer_,
        contrast_n_genomes_tested = NA_integer_,
        warning_message = NA_character_,
        error_message = NA_character_,
        status = NA_character_,
        row.names = comparison_id
    )
    row$ko_group_comparison_id <- comparison_id
    row$ko_id <- ko_id
    row[[paste0("reference_", group_column)]] <- reference_group
    row[[paste0("contrast_", group_column)]] <- contrast_group
    row <- row[, c(
        "ko_group_comparison_id",
        "ko_id",
        paste0("reference_", group_column),
        paste0("contrast_", group_column),
        setdiff(
            names(row),
            c(
                "ko_group_comparison_id",
                "ko_id",
                paste0("reference_", group_column),
                paste0("contrast_", group_column)
            )
        )
    )]
    row
}

.fit_one_ko_group_comparison <- function(
    ko_id,
    reference_group,
    contrast_group,
    data,
    group_column,
    variable_info,
    status_filter
) {
    row <- .empty_ko_group_comparison_row(
        ko_id = ko_id,
        group_column = group_column,
        reference_group = reference_group,
        contrast_group = contrast_group,
        variable_info = variable_info
    )

    ref_row <- data[
        as.character(data[[group_column]]) == reference_group,
        ,
        drop = FALSE
    ]
    contrast_row <- data[
        as.character(data[[group_column]]) == contrast_group,
        ,
        drop = FALSE
    ]

    if (nrow(ref_row) > 1L || nrow(contrast_row) > 1L) {
        row$status <- "error"
        row$error_message <- "Each KO/group combination must occur at most once in the input fit."
        return(row)
    }

    if (nrow(ref_row) == 0L || nrow(contrast_row) == 0L) {
        row$status <- "skipped"
        row$error_message <- "At least one requested group was unavailable for this KO."
        return(row)
    }

    row$reference_group_status <- as.character(ref_row$status)
    row$contrast_group_status <- as.character(contrast_row$status)
    row$reference_n_genomes_tested <- as.integer(ref_row$n_genomes_tested)
    row$contrast_n_genomes_tested <- as.integer(contrast_row$n_genomes_tested)
    row$reference_group_estimate <- as.numeric(ref_row$estimate)
    row$reference_group_std_error <- as.numeric(ref_row$std_error)
    row$contrast_group_estimate <- as.numeric(contrast_row$estimate)
    row$contrast_group_std_error <- as.numeric(contrast_row$std_error)

    if (!is.null(status_filter)) {
        valid_status <- as.character(c(ref_row$status, contrast_row$status)) %in% status_filter
        if (!all(valid_status)) {
            row$status <- "skipped"
            row$error_message <- "At least one requested group did not pass the status filter."
            return(row)
        }
    }

    if (anyNA(c(
        row$reference_group_estimate,
        row$reference_group_std_error,
        row$contrast_group_estimate,
        row$contrast_group_std_error
    )) ||
        any(!is.finite(c(
            row$reference_group_estimate,
            row$reference_group_std_error,
            row$contrast_group_estimate,
            row$contrast_group_std_error
        ))) ||
        any(c(
            row$reference_group_std_error,
            row$contrast_group_std_error
        ) <= 0)) {
        row$status <- "skipped"
        row$error_message <- "Both groups must have finite estimates and positive standard errors."
        return(row)
    }

    row$estimate <- row$contrast_group_estimate - row$reference_group_estimate
    row$std_error <- sqrt(row$contrast_group_std_error^2 + row$reference_group_std_error^2)
    row$statistic <- row$estimate / row$std_error
    row$p_value <- 2 * stats::pnorm(-abs(row$statistic))
    row$status <- "ok"
    row$error_message <- NA_character_
    row
}

#' Compare KO Genome-Group Responses Between Groups
#'
#' `compareKOGenomeGroups()` follows [fitKOGenomeGroupMetaAnalysis()] with a
#' direct between-group comparison step. It compares the KO-level group
#' meta-estimates for two groups and tests whether the estimated KO response
#' differs between them.
#'
#' This workflow is intended for the question: for a given KO, is the
#' group-level response different between two clades, domains, or other
#' genome-derived groups?
#'
#' The input `x` must be an `MTTKFit` returned by
#' [fitKOGenomeGroupMetaAnalysis()]. The comparison estimate is calculated as
#' `contrast - reference`, using the group-level KO meta-estimates stored in the
#' input fit. The comparison standard error is computed as
#' `sqrt(se_contrast^2 + se_reference^2)`, which is appropriate when the two
#' compared groups contain disjoint genomes.
#'
#' Because the input group-level KO fits are themselves approximate summaries of
#' genome-specific conditional effects, the resulting between-group p-values
#' should also be interpreted as approximate.
#'
#' @param x An `MTTKFit` returned by [fitKOGenomeGroupMetaAnalysis()].
#' @param reference Optional single reference group identifier. When `NULL` and
#'   `contrast = NULL`, all pairwise group comparisons are returned. When
#'   `reference` is supplied and `contrast = NULL`, MTTK compares that
#'   reference group against every other available group.
#' @param contrast Optional single contrast group identifier. When supplied, a
#'   single directed contrast `contrast - reference` is returned.
#' @param status Optional fit-status filter applied to the input KO/group rows
#'   before comparison. The default keeps only rows with `status == "ok"`. Use
#'   `NULL` to keep all KO/group rows.
#'
#' @return An `MTTKFit` with one row per KO/group comparison.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeExampleMTTKExperiment()
#'     ko_fit <- fitKORandomSlopeModel(x, variable = "condition")
#'     ko_domain_fit <- fitKOGenomeGroupMetaAnalysis(ko_fit, x, group = "domain")
#'     compareKOGenomeGroups(
#'         ko_domain_fit,
#'         reference = "Bacteria",
#'         contrast = "Archaea"
#'     )
#' }
#'
#' @export
compareKOGenomeGroups <- function(
    x,
    reference = NULL,
    contrast = NULL,
    status = "ok"
) {
    parsed <- .require_ko_genome_group_fit(x)
    group_column <- parsed$groupColumn
    source_info <- parsed$info
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    comparison_pairs <- .normalize_group_comparison_pairs(
        group_values = parsed$results[[group_column]],
        reference = reference,
        contrast = contrast
    )

    ko_ids <- unique(as.character(parsed$results$ko_id))
    results <- do.call(
        rbind,
        lapply(seq_len(nrow(comparison_pairs)), function(i) {
            ref_group <- comparison_pairs$reference_group[[i]]
            contrast_group <- comparison_pairs$contrast_group[[i]]

            do.call(
                rbind,
                lapply(ko_ids, function(ko_id) {
                    ko_rows <- parsed$results[
                        as.character(parsed$results$ko_id) == ko_id,
                        ,
                        drop = FALSE
                    ]
                    .fit_one_ko_group_comparison(
                        ko_id = ko_id,
                        reference_group = ref_group,
                        contrast_group = contrast_group,
                        data = ko_rows,
                        group_column = group_column,
                        variable_info = variable_info,
                        status_filter = status_filter
                    )
                })
            )
        })
    )
    rownames(results) <- as.character(results$ko_group_comparison_id)

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    comparison_mode <- if (is.null(reference) && is.null(contrast)) {
        "all_pairwise"
    } else if (!is.null(reference) && is.null(contrast)) {
        "reference_vs_all"
    } else {
        "single_contrast"
    }

    info <- list(
        backend = "group_comparison",
        model = "ko_genome_group_comparison",
        featureType = paste0("ko_group_comparison:", group_column),
        featureIdColumn = "ko_group_comparison_id",
        sourceFeatureType = if (!is.null(source_info$featureType)) {
            source_info$featureType
        } else {
            paste0("ko_group:", group_column)
        },
        sourceFeatureIdColumn = if (!is.null(source_info$featureIdColumn)) {
            source_info$featureIdColumn
        } else {
            "ko_group_id"
        },
        sourceFitModel = if (!is.null(source_info$model)) source_info$model else NA_character_,
        sourceBackend = if (!is.null(source_info$backend)) source_info$backend else NA_character_,
        sourceGroupColumn = group_column,
        sourceEffectSource = if (!is.null(source_info$effectSource)) {
            source_info$effectSource
        } else {
            NA_character_
        },
        sourceStdErrorSource = if (!is.null(source_info$stdErrorSource)) {
            source_info$stdErrorSource
        } else {
            NA_character_
        },
        variable = variable_info$variable,
        variableType = variable_info$type,
        referenceLevel = variable_info$referenceLevel,
        contrastLevel = variable_info$contrastLevel,
        effectLabel = variable_info$effectLabel,
        groupColumn = group_column,
        comparisonMode = comparison_mode,
        requestedReference = if (is.null(reference)) NA_character_ else reference,
        requestedContrast = if (is.null(contrast)) NA_character_ else contrast,
        statusFilter = if (is.null(status_filter)) NA_character_ else status_filter,
        comparisonFormula = "contrast_group_estimate - reference_group_estimate",
        comparisonStdError = "sqrt(se_contrast^2 + se_reference^2)",
        n_features = nrow(results),
        n_ok = sum(results$status == "ok"),
        n_skipped = sum(results$status == "skipped"),
        n_error = sum(results$status == "error")
    )

    MTTKFit(
        results = results,
        info = info,
        models = list()
    )
}

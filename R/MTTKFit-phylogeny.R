.normalize_phylo_correlation_method <- function(method) {
    match.arg(method, c("spearman", "pearson"))
}

.normalize_phylo_alternative <- function(alternative) {
    match.arg(alternative, c("two.sided", "greater", "less"))
}

.normalize_phylo_permutations <- function(n_perm) {
    if (!is.numeric(n_perm) ||
        length(n_perm) != 1L ||
        is.na(n_perm) ||
        n_perm < 0) {
        stop("'nPerm' must be a single non-negative integer.", call. = FALSE)
    }

    as.integer(n_perm)
}

.normalize_tree_min_genomes <- function(min_genomes) {
    if (!is.numeric(min_genomes) ||
        length(min_genomes) != 1L ||
        is.na(min_genomes) ||
        min_genomes < 3) {
        stop("'minGenomes' must be a single integer greater than or equal to 3.", call. = FALSE)
    }

    as.integer(min_genomes)
}

.normalize_tree_min_tips <- function(min_tips) {
    if (!is.numeric(min_tips) ||
        length(min_tips) != 1L ||
        is.na(min_tips) ||
        min_tips < 2) {
        stop("'minTips' must be a single integer greater than or equal to 2.", call. = FALSE)
    }

    as.integer(min_tips)
}

.normalize_tree_max_tips <- function(max_tips, n_tips) {
    if (is.null(max_tips)) {
        return(NULL)
    }

    if (!is.numeric(max_tips) ||
        length(max_tips) != 1L ||
        is.na(max_tips) ||
        max_tips < 2) {
        stop("'maxTips' must be NULL or a single integer greater than or equal to 2.", call. = FALSE)
    }

    max_tips <- as.integer(max_tips)
    if (max_tips >= n_tips) {
        return(NULL)
    }

    max_tips
}

.normalize_tree_collapse <- function(collapse) {
    if (!is.character(collapse) || length(collapse) != 1L || is.na(collapse)) {
        stop("'collapse' must be a single character string.", call. = FALSE)
    }

    collapse
}

.require_genome_tree_object <- function(object, tree = NULL) {
    if (!methods::is(object, "MTTKExperiment")) {
        stop("'object' must be an MTTKExperiment.", call. = FALSE)
    }

    if (is.null(tree)) {
        tree <- genomeTree(object)
    }

    if (is.null(tree)) {
        stop(
            "No genome phylogeny is stored in genomeTree(object). Set genomeTree(x) before running tree-based analyses.",
            call. = FALSE
        )
    }

    .normalize_genome_tree(tree)
}

.ensure_tree_branch_lengths <- function(tree) {
    if (is.null(tree$edge.length)) {
        return(ape::compute.brlen(tree, method = 1))
    }

    tree
}

.align_tree_values <- function(tree, genome_ids, values) {
    genome_ids <- as.character(genome_ids)
    values <- as.numeric(values)

    keep <- !is.na(genome_ids) &
        genome_ids != "" &
        !is.na(values) &
        is.finite(values)
    genome_ids <- genome_ids[keep]
    values <- values[keep]

    if (length(genome_ids) == 0L) {
        return(list(
            tree = NULL,
            values = numeric(),
            genomeIds = character(),
            nAvailable = 0L,
            nDroppedFromTree = 0L,
            nDroppedFromValues = 0L
        ))
    }

    deduplicated <- !duplicated(genome_ids)
    genome_ids <- genome_ids[deduplicated]
    values <- values[deduplicated]

    available_ids <- intersect(tree$tip.label, genome_ids)
    aligned_tree <- .prune_genome_tree(tree, available_ids)
    aligned_values <- stats::setNames(values, genome_ids)[aligned_tree$tip.label]

    list(
        tree = aligned_tree,
        values = as.numeric(aligned_values),
        genomeIds = aligned_tree$tip.label,
        nAvailable = length(genome_ids),
        nDroppedFromTree = length(setdiff(tree$tip.label, available_ids)),
        nDroppedFromValues = length(setdiff(genome_ids, tree$tip.label))
    )
}

.align_tree_data_frame <- function(tree, data, genome_column = "genome_id", required_columns = character()) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)

    if (!(genome_column %in% names(data))) {
        stop("Missing genome identifier column: ", genome_column, ".", call. = FALSE)
    }

    keep <- !is.na(data[[genome_column]]) & as.character(data[[genome_column]]) != ""
    for (column in required_columns) {
        keep <- keep &
            !is.na(data[[column]]) &
            is.finite(as.numeric(data[[column]]))
    }

    data <- data[keep, , drop = FALSE]
    data[[genome_column]] <- as.character(data[[genome_column]])

    if (nrow(data) == 0L) {
        return(list(
            tree = NULL,
            data = data,
            genomeIds = character(),
            nAvailable = 0L
        ))
    }

    data <- data[!duplicated(data[[genome_column]]), , drop = FALSE]
    available_ids <- intersect(tree$tip.label, data[[genome_column]])
    aligned_tree <- .prune_genome_tree(tree, available_ids)
    data <- data[match(aligned_tree$tip.label, data[[genome_column]]), , drop = FALSE]
    rownames(data) <- data[[genome_column]]

    list(
        tree = aligned_tree,
        data = data,
        genomeIds = aligned_tree$tip.label,
        nAvailable = nrow(data)
    )
}

.normalize_phylo_gls_method <- function(method) {
    match.arg(method, c("REML", "ML"))
}

.require_feature_random_slope_fit <- function(
    x,
    model_name,
    feature_id_column,
    function_label
) {
    if (!methods::is(x, "MTTKFit")) {
        stop("'x' must be an MTTKFit.", call. = FALSE)
    }

    info <- fitInfo(x)
    if (!identical(info$model, model_name)) {
        stop(
            "'x' must be an MTTKFit returned by ",
            function_label,
            "().",
            call. = FALSE
        )
    }

    feature_results <- S4Vectors::DataFrame(x, check.names = FALSE)
    if (!(feature_id_column %in% names(feature_results))) {
        feature_results[[feature_id_column]] <- rownames(feature_results)
    }

    group_effects <- .stored_group_effects(x)
    if (is.data.frame(group_effects) && !methods::is(group_effects, "DataFrame")) {
        group_effects <- S4Vectors::DataFrame(group_effects, check.names = FALSE)
    }

    required_cols <- c(
        feature_id_column,
        "genome_id",
        "conditional_effect_estimate",
        "conditional_effect_std_error"
    )
    missing_cols <- setdiff(required_cols, names(group_effects))
    if (length(missing_cols) > 0L) {
        stop(
            "The random-slope fit is missing required group-effect column(s): ",
            paste(missing_cols, collapse = ", "),
            ".",
            call. = FALSE
        )
    }

    list(
        featureResults = feature_results,
        groupEffects = group_effects,
        info = info,
        featureIdColumn = feature_id_column
    )
}

.feature_status_map <- function(feature_results, feature_id_column) {
    feature_ids <- if (feature_id_column %in% names(feature_results)) {
        as.character(feature_results[[feature_id_column]])
    } else {
        rownames(feature_results)
    }

    stats::setNames(
        if ("status" %in% names(feature_results)) {
            as.character(feature_results$status)
        } else {
            rep(NA_character_, nrow(feature_results))
        },
        feature_ids
    )
}

.fit_feature_phylogenetic_signal_results <- function(
    parsed,
    tree,
    feature_id_column,
    correlation_method,
    alternative,
    n_perm,
    min_genomes,
    status_filter,
    variable_info
) {
    feature_status <- .feature_status_map(parsed$featureResults, feature_id_column)
    group_effects <- as.data.frame(parsed$groupEffects, stringsAsFactors = FALSE)
    group_effects[[feature_id_column]] <- as.character(group_effects[[feature_id_column]])
    group_effects$genome_id <- as.character(group_effects$genome_id)
    feature_split <- split(group_effects, group_effects[[feature_id_column]])
    feature_ids <- rownames(parsed$featureResults)

    results <- do.call(
        rbind,
        lapply(feature_ids, function(feature_id) {
            data <- feature_split[[feature_id]]

            if (is.null(data)) {
                return(.phylo_signal_row(
                    analysis_id = feature_id,
                    feature_id_column = feature_id_column,
                    variable_info = variable_info,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = 0L,
                    n_tested = 0L,
                    genome_ids = character(),
                    correlation_method = correlation_method,
                    alternative = alternative,
                    n_perm = n_perm,
                    status = "skipped",
                    error_message = "No genome-specific feature effects were available."
                ))
            }

            if (!is.null(status_filter) &&
                !(feature_status[[feature_id]] %in% status_filter)) {
                return(.phylo_signal_row(
                    analysis_id = feature_id,
                    feature_id_column = feature_id_column,
                    variable_info = variable_info,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = nrow(data),
                    n_tested = 0L,
                    genome_ids = character(),
                    correlation_method = correlation_method,
                    alternative = alternative,
                    n_perm = n_perm,
                    status = "skipped",
                    error_message = "The feature did not pass the requested status filter."
                ))
            }

            aligned <- .align_tree_values(
                tree = tree,
                genome_ids = data$genome_id,
                values = data$conditional_effect_estimate
            )

            if (length(aligned$values) < min_genomes) {
                return(.phylo_signal_row(
                    analysis_id = feature_id,
                    feature_id_column = feature_id_column,
                    variable_info = variable_info,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = aligned$nAvailable,
                    n_tested = length(aligned$values),
                    genome_ids = aligned$genomeIds,
                    correlation_method = correlation_method,
                    alternative = alternative,
                    n_perm = n_perm,
                    status = "skipped",
                    error_message = paste0(
                        "At least ",
                        min_genomes,
                        " genomes with finite feature effects are required."
                    )
                ))
            }

            statistic <- .phylogenetic_signal_statistic(
                tree = aligned$tree,
                values = aligned$values,
                method = correlation_method
            )
            permuted <- .permuted_signal_statistics(
                tree = aligned$tree,
                values = aligned$values,
                method = correlation_method,
                n_perm = n_perm
            )

            .phylo_signal_row(
                analysis_id = feature_id,
                feature_id_column = feature_id_column,
                variable_info = variable_info,
                statistic = statistic,
                p_value = .permutation_p_value(statistic, permuted, alternative),
                n_available = aligned$nAvailable,
                n_tested = length(aligned$values),
                genome_ids = aligned$genomeIds,
                correlation_method = correlation_method,
                alternative = alternative,
                n_perm = n_perm,
                status = if (is.finite(statistic)) "ok" else "skipped",
                error_message = if (is.finite(statistic)) NA_character_ else {
                    "The genome-specific feature effects did not vary enough to compute a distance correlation."
                }
            )
        })
    )
    rownames(results) <- feature_ids
    results
}

.scan_feature_clades_results <- function(
    parsed,
    tree,
    feature_id_column,
    result_id_column,
    alternative,
    n_perm,
    min_tips,
    max_tips,
    min_genomes,
    status_filter,
    variable_info
) {
    feature_status <- .feature_status_map(parsed$featureResults, feature_id_column)
    group_effects <- as.data.frame(parsed$groupEffects, stringsAsFactors = FALSE)
    group_effects[[feature_id_column]] <- as.character(group_effects[[feature_id_column]])
    group_effects$genome_id <- as.character(group_effects$genome_id)
    feature_split <- split(group_effects, group_effects[[feature_id_column]])
    feature_ids <- rownames(parsed$featureResults)

    results <- do.call(
        rbind,
        lapply(feature_ids, function(feature_id) {
            data <- feature_split[[feature_id]]

            if (is.null(data)) {
                return(.clade_scan_row(
                    row_id = paste0(feature_id, "::node_na"),
                    variable_info = variable_info,
                    estimate = NA_real_,
                    p_value = NA_real_,
                    clade_id = NA_character_,
                    clade_label = NA_character_,
                    node_id = NA_integer_,
                    feature_id_column = feature_id_column,
                    feature_id = feature_id,
                    member_ids = character(),
                    all_ids = character(),
                    n_perm = n_perm,
                    alternative = alternative,
                    status = "skipped",
                    error_message = "No genome-specific feature effects were available."
                ))
            }

            if (!is.null(status_filter) &&
                !(feature_status[[feature_id]] %in% status_filter)) {
                return(.clade_scan_row(
                    row_id = paste0(feature_id, "::node_na"),
                    variable_info = variable_info,
                    estimate = NA_real_,
                    p_value = NA_real_,
                    clade_id = NA_character_,
                    clade_label = NA_character_,
                    node_id = NA_integer_,
                    feature_id_column = feature_id_column,
                    feature_id = feature_id,
                    member_ids = character(),
                    all_ids = character(),
                    n_perm = n_perm,
                    alternative = alternative,
                    status = "skipped",
                    error_message = "The feature did not pass the requested status filter."
                ))
            }

            aligned <- .align_tree_values(
                tree = tree,
                genome_ids = data$genome_id,
                values = data$conditional_effect_estimate
            )

            if (length(aligned$values) < min_genomes) {
                return(.clade_scan_row(
                    row_id = paste0(feature_id, "::node_na"),
                    variable_info = variable_info,
                    estimate = NA_real_,
                    p_value = NA_real_,
                    clade_id = NA_character_,
                    clade_label = NA_character_,
                    node_id = NA_integer_,
                    feature_id_column = feature_id_column,
                    feature_id = feature_id,
                    member_ids = character(),
                    all_ids = aligned$genomeIds,
                    n_perm = n_perm,
                    alternative = alternative,
                    status = "skipped",
                    error_message = paste0(
                        "At least ",
                        min_genomes,
                        " genomes with finite feature effects are required."
                    )
                ))
            }

            values <- stats::setNames(aligned$values, aligned$genomeIds)
            clades <- .clade_memberships(
                tree = aligned$tree,
                min_tips = min_tips,
                max_tips = max_tips,
                collapse = ";"
            )

            if (length(clades) == 0L) {
                return(.clade_scan_row(
                    row_id = paste0(feature_id, "::node_na"),
                    variable_info = variable_info,
                    estimate = NA_real_,
                    p_value = NA_real_,
                    clade_id = "node_na",
                    clade_label = NA_character_,
                    node_id = NA_integer_,
                    feature_id_column = feature_id_column,
                    feature_id = feature_id,
                    member_ids = character(),
                    all_ids = aligned$genomeIds,
                    n_perm = n_perm,
                    alternative = alternative,
                    status = "skipped",
                    error_message = "No clades satisfied the requested subtree-size constraints."
                ))
            }

            do.call(
                rbind,
                lapply(clades, function(clade) {
                    estimate <- .clade_difference_statistic(values, member_ids = clade$memberGenomeIds)
                    permuted <- .permuted_clade_statistics(
                        values = values,
                        member_ids = clade$memberGenomeIds,
                        n_perm = n_perm
                    )
                    .clade_scan_row(
                        row_id = paste0(feature_id, "::", clade$cladeId),
                        variable_info = variable_info,
                        estimate = estimate,
                        p_value = .permutation_p_value(estimate, permuted, alternative),
                        clade_id = clade$cladeId,
                        clade_label = clade$cladeLabel,
                        node_id = clade$nodeId,
                        feature_id_column = feature_id_column,
                        feature_id = feature_id,
                        member_ids = clade$memberGenomeIds,
                        all_ids = aligned$genomeIds,
                        n_perm = n_perm,
                        alternative = alternative,
                        status = if (is.finite(estimate)) "ok" else "skipped",
                        error_message = if (is.finite(estimate)) NA_character_ else {
                            "The clade statistic could not be computed."
                        }
                    )
                })
            )
        })
    )

    results[[result_id_column]] <- rownames(results)
    results[, c(result_id_column, setdiff(names(results), result_id_column))]
}

.phylo_gls_row <- function(
    analysis_id,
    feature_id_column,
    variable_info,
    estimate,
    std_error,
    statistic,
    p_value,
    n_available,
    n_tested,
    genome_ids,
    gls_method,
    use_weights,
    status,
    error_message = NA_character_,
    warning_message = NA_character_
) {
    row <- S4Vectors::DataFrame(
        tested_term = "(Intercept)",
        variable_type = variable_info$type,
        reference_level = variable_info$referenceLevel,
        contrast_level = variable_info$contrastLevel,
        effect_label = variable_info$effectLabel,
        estimate = estimate,
        std_error = std_error,
        statistic = statistic,
        p_value = p_value,
        q_value = NA_real_,
        n_genomes_available = as.integer(n_available),
        n_genomes_tested = as.integer(n_tested),
        tested_genome_ids = if (length(genome_ids) == 0L) NA_character_ else {
            paste(genome_ids, collapse = ";")
        },
        gls_method = gls_method,
        correlation_model = "Brownian",
        use_measurement_error_weights = use_weights,
        warning_message = warning_message,
        error_message = error_message,
        status = status,
        row.names = analysis_id
    )
    row[[feature_id_column]] <- analysis_id
    row[, c(feature_id_column, setdiff(names(row), feature_id_column))]
}

.fit_brownian_gls <- function(aligned_data, tree, method, use_weights) {
    model_data <- data.frame(
        genome_id = factor(
            as.character(aligned_data$genome_id),
            levels = tree$tip.label
        ),
        response = as.numeric(aligned_data$response),
        effect_variance = if ("std_error" %in% names(aligned_data)) {
            pmax(as.numeric(aligned_data$std_error)^2, .Machine$double.eps)
        } else {
            1
        },
        stringsAsFactors = FALSE
    )

    gls_fit <- nlme::gls(
        response ~ 1,
        data = model_data,
        correlation = ape::corBrownian(phy = .ensure_tree_branch_lengths(tree), form = ~ genome_id),
        weights = if (use_weights && "std_error" %in% names(aligned_data)) {
            nlme::varFixed(~ effect_variance)
        } else {
            NULL
        },
        method = method
    )

    coef_row <- summary(gls_fit)$tTable["(Intercept)", , drop = FALSE]

    list(
        fit = gls_fit,
        estimate = as.numeric(coef_row[1L, "Value"]),
        std_error = as.numeric(coef_row[1L, "Std.Error"]),
        statistic = as.numeric(coef_row[1L, "t-value"]),
        p_value = as.numeric(coef_row[1L, "p-value"])
    )
}

.fit_feature_phylogenetic_gls_results <- function(
    parsed,
    tree,
    feature_id_column,
    gls_method,
    use_weights,
    min_genomes,
    status_filter,
    variable_info,
    keep_fits = FALSE
) {
    feature_status <- .feature_status_map(parsed$featureResults, feature_id_column)
    group_effects <- as.data.frame(parsed$groupEffects, stringsAsFactors = FALSE)
    group_effects[[feature_id_column]] <- as.character(group_effects[[feature_id_column]])
    group_effects$genome_id <- as.character(group_effects$genome_id)
    feature_split <- split(group_effects, group_effects[[feature_id_column]])
    feature_ids <- rownames(parsed$featureResults)
    fit_list <- vector("list", length(feature_ids))
    names(fit_list) <- feature_ids

    results <- do.call(
        rbind,
        lapply(feature_ids, function(feature_id) {
            data <- feature_split[[feature_id]]

            if (is.null(data)) {
                return(.phylo_gls_row(
                    analysis_id = feature_id,
                    feature_id_column = feature_id_column,
                    variable_info = variable_info,
                    estimate = NA_real_,
                    std_error = NA_real_,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = 0L,
                    n_tested = 0L,
                    genome_ids = character(),
                    gls_method = gls_method,
                    use_weights = use_weights,
                    status = "skipped",
                    error_message = "No genome-specific feature effects were available."
                ))
            }

            if (!is.null(status_filter) &&
                !(feature_status[[feature_id]] %in% status_filter)) {
                return(.phylo_gls_row(
                    analysis_id = feature_id,
                    feature_id_column = feature_id_column,
                    variable_info = variable_info,
                    estimate = NA_real_,
                    std_error = NA_real_,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = nrow(data),
                    n_tested = 0L,
                    genome_ids = character(),
                    gls_method = gls_method,
                    use_weights = use_weights,
                    status = "skipped",
                    error_message = "The feature did not pass the requested status filter."
                ))
            }

            aligned <- .align_tree_data_frame(
                tree = tree,
                data = data.frame(
                    genome_id = data$genome_id,
                    response = data$conditional_effect_estimate,
                    std_error = data$conditional_effect_std_error,
                    stringsAsFactors = FALSE
                ),
                genome_column = "genome_id",
                required_columns = c("response", "std_error")
            )

            if (nrow(aligned$data) < min_genomes) {
                return(.phylo_gls_row(
                    analysis_id = feature_id,
                    feature_id_column = feature_id_column,
                    variable_info = variable_info,
                    estimate = NA_real_,
                    std_error = NA_real_,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = aligned$nAvailable,
                    n_tested = nrow(aligned$data),
                    genome_ids = aligned$genomeIds,
                    gls_method = gls_method,
                    use_weights = use_weights,
                    status = "skipped",
                    error_message = paste0(
                        "At least ",
                        min_genomes,
                        " genomes with finite feature effects are required."
                    )
                ))
            }

            fitted <- tryCatch(
                .fit_brownian_gls(
                    aligned_data = aligned$data,
                    tree = aligned$tree,
                    method = gls_method,
                    use_weights = use_weights
                ),
                error = function(e) e
            )

            if (inherits(fitted, "error")) {
                return(.phylo_gls_row(
                    analysis_id = feature_id,
                    feature_id_column = feature_id_column,
                    variable_info = variable_info,
                    estimate = NA_real_,
                    std_error = NA_real_,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = aligned$nAvailable,
                    n_tested = nrow(aligned$data),
                    genome_ids = aligned$genomeIds,
                    gls_method = gls_method,
                    use_weights = use_weights,
                    status = "error",
                    error_message = conditionMessage(fitted)
                ))
            }

            if (keep_fits) {
                fit_list[[feature_id]] <<- fitted$fit
            }

            .phylo_gls_row(
                analysis_id = feature_id,
                feature_id_column = feature_id_column,
                variable_info = variable_info,
                estimate = fitted$estimate,
                std_error = fitted$std_error,
                statistic = fitted$statistic,
                p_value = fitted$p_value,
                n_available = aligned$nAvailable,
                n_tested = nrow(aligned$data),
                genome_ids = aligned$genomeIds,
                gls_method = gls_method,
                use_weights = use_weights,
                status = "ok"
            )
        })
    )
    rownames(results) <- feature_ids

    list(
        results = results,
        models = if (keep_fits) {
            fit_list[!vapply(fit_list, is.null, logical(1))]
        } else {
            list()
        }
    )
}

.phylogenetic_signal_statistic <- function(tree, values, method) {
    tree <- .ensure_tree_branch_lengths(tree)
    phylo_dist <- as.matrix(ape::cophenetic.phylo(tree))
    response_dist <- abs(outer(values, values, "-"))
    lower <- lower.tri(phylo_dist, diag = FALSE)
    x <- as.numeric(phylo_dist[lower])
    y <- as.numeric(response_dist[lower])

    if (length(x) == 0L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
        return(NA_real_)
    }

    stats::cor(x, y, method = method)
}

.permutation_p_value <- function(observed, permuted, alternative) {
    permuted <- as.numeric(permuted)
    permuted <- permuted[!is.na(permuted) & is.finite(permuted)]

    if (!is.finite(observed)) {
        return(NA_real_)
    }

    if (length(permuted) == 0L) {
        return(1)
    }

    if (identical(alternative, "greater")) {
        (sum(permuted >= observed) + 1) / (length(permuted) + 1)
    } else if (identical(alternative, "less")) {
        (sum(permuted <= observed) + 1) / (length(permuted) + 1)
    } else {
        (sum(abs(permuted) >= abs(observed)) + 1) / (length(permuted) + 1)
    }
}

.permuted_signal_statistics <- function(tree, values, method, n_perm) {
    if (n_perm == 0L) {
        return(numeric())
    }

    vapply(
        seq_len(n_perm),
        function(i) {
            .phylogenetic_signal_statistic(
                tree = tree,
                values = sample(values, replace = FALSE),
                method = method
            )
        },
        numeric(1)
    )
}

.clade_memberships <- function(tree, min_tips, max_tips, collapse) {
    if (!ape::is.rooted(tree)) {
        stop(
            "Clade-scanning workflows require a rooted genomeTree(object).",
            call. = FALSE
        )
    }

    n_tips <- ape::Ntip(tree)
    max_tips <- .normalize_tree_max_tips(max_tips, n_tips = n_tips)
    node_ids <- seq.int(n_tips + 1L, n_tips + tree$Nnode)

    memberships <- lapply(node_ids, function(node_id) {
        tips <- sort(ape::extract.clade(tree, node = node_id)$tip.label)
        if (length(tips) < min_tips || length(tips) >= n_tips) {
            return(NULL)
        }

        if (!is.null(max_tips) && length(tips) > max_tips) {
            return(NULL)
        }

        node_index <- node_id - n_tips
        node_label <- NA_character_
        if (!is.null(tree$node.label) && length(tree$node.label) >= node_index) {
            node_label <- as.character(tree$node.label[[node_index]])
        }

        list(
            cladeId = paste0("node_", node_id),
            cladeLabel = if (!is.na(node_label) && nzchar(node_label)) {
                node_label
            } else {
                paste0("node_", node_id)
            },
            nodeId = as.integer(node_id),
            memberGenomeIds = tips,
            memberGenomeLabel = paste(tips, collapse = collapse)
        )
    })

    Filter(Negate(is.null), memberships)
}

.clade_difference_statistic <- function(values, member_ids) {
    member_values <- values[names(values) %in% member_ids]
    complement_values <- values[!(names(values) %in% member_ids)]

    if (length(member_values) == 0L || length(complement_values) == 0L) {
        return(NA_real_)
    }

    mean(member_values) - mean(complement_values)
}

.permuted_clade_statistics <- function(values, member_ids, n_perm) {
    if (n_perm == 0L) {
        return(numeric())
    }

    vapply(
        seq_len(n_perm),
        function(i) {
            permuted <- sample(values, replace = FALSE)
            names(permuted) <- names(values)
            .clade_difference_statistic(permuted, member_ids = member_ids)
        },
        numeric(1)
    )
}

.signal_direction <- function(statistic) {
    if (is.na(statistic) || !is.finite(statistic)) {
        return(NA_character_)
    }

    if (statistic > 0) {
        "positive"
    } else if (statistic < 0) {
        "negative"
    } else {
        "zero"
    }
}

.phylo_signal_row <- function(
    analysis_id,
    feature_id_column,
    variable_info,
    statistic,
    p_value,
    n_available,
    n_tested,
    genome_ids,
    correlation_method,
    alternative,
    n_perm,
    status,
    error_message = NA_character_,
    warning_message = NA_character_
) {
    row <- S4Vectors::DataFrame(
        tested_term = variable_info$testedTerm,
        variable_type = variable_info$type,
        reference_level = variable_info$referenceLevel,
        contrast_level = variable_info$contrastLevel,
        effect_label = variable_info$effectLabel,
        signal_statistic = statistic,
        signal_direction = .signal_direction(statistic),
        p_value = p_value,
        q_value = NA_real_,
        n_genomes_available = as.integer(n_available),
        n_genomes_tested = as.integer(n_tested),
        tested_genome_ids = if (length(genome_ids) == 0L) NA_character_ else {
            paste(genome_ids, collapse = ";")
        },
        correlation_method = correlation_method,
        alternative = alternative,
        n_permutations = as.integer(n_perm),
        warning_message = warning_message,
        error_message = error_message,
        status = status,
        row.names = analysis_id
    )
    row[[feature_id_column]] <- analysis_id
    row[, c(feature_id_column, setdiff(names(row), feature_id_column))]
}

.clade_scan_row <- function(
    row_id,
    variable_info,
    estimate,
    p_value,
    clade_id,
    clade_label,
    node_id,
    feature_id_column,
    feature_id,
    member_ids,
    all_ids,
    n_perm,
    alternative,
    status,
    warning_message = NA_character_,
    error_message = NA_character_
) {
    member_values <- as.character(member_ids)
    complement_ids <- setdiff(as.character(all_ids), member_values)
    clade_id_value <- if (!is.na(clade_id) && nzchar(clade_id)) {
        clade_id
    } else {
        row_id
    }

    row <- S4Vectors::DataFrame(
        tested_term = variable_info$testedTerm,
        variable_type = variable_info$type,
        reference_level = variable_info$referenceLevel,
        contrast_level = variable_info$contrastLevel,
        effect_label = variable_info$effectLabel,
        clade_id = clade_id_value,
        clade_label = clade_label,
        node_id = as.integer(node_id),
        estimate = estimate,
        p_value = p_value,
        q_value = NA_real_,
        n_genomes_clade = as.integer(length(member_values)),
        n_genomes_complement = as.integer(length(complement_ids)),
        member_genome_ids = if (length(member_values) == 0L) NA_character_ else {
            paste(member_values, collapse = ";")
        },
        complement_genome_ids = if (length(complement_ids) == 0L) NA_character_ else {
            paste(complement_ids, collapse = ";")
        },
        alternative = alternative,
        n_permutations = as.integer(n_perm),
        warning_message = warning_message,
        error_message = error_message,
        status = status,
        row.names = row_id
    )

    if (!is.null(feature_id_column) && !is.null(feature_id)) {
        row[[feature_id_column]] <- feature_id
        row <- row[, c(
            feature_id_column,
            setdiff(names(row), feature_id_column)
        )]
    }

    row
}

#' Test for Overall Phylogenetic Signal in Genome-Level Responses
#'
#' `fitGenomePhylogeneticSignal()` follows [fitGenomeModel()] with a tree-based
#' summary of whether more closely related genomes tend to show more similar
#' fitted responses.
#'
#' It tests the association between pairwise phylogenetic distances and
#' pairwise absolute response differences across genomes. Positive statistics
#' indicate that more distantly related genomes tend to differ more strongly in
#' their fitted response estimates.
#'
#' @param x A genome-level `MTTKFit`, typically returned by [fitGenomeModel()].
#' @param object An `MTTKExperiment` supplying `genomeTree(object)`.
#' @param tree Optional `ape::phylo` object. When `NULL`, `genomeTree(object)`
#'   is used.
#' @param correlationMethod Correlation used to quantify the association
#'   between phylogenetic distance and response difference.
#' @param alternative Alternative hypothesis for the permutation p-value.
#' @param nPerm Number of tip-label permutations used to estimate the p-value.
#' @param status Optional genome-fit status filter applied before the signal
#'   test. The default keeps only genome rows with `status == "ok"`.
#'
#' @return An `MTTKFit` with one row summarizing the overall phylogenetic
#'   signal in the genome-level responses.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     genome_fit <- fitGenomeModel(x, variable = "condition")
#'     fitGenomePhylogeneticSignal(genome_fit, x, nPerm = 99)
#' }
#'
#' @export
fitGenomePhylogeneticSignal <- function(
    x,
    object,
    tree = NULL,
    correlationMethod = c("spearman", "pearson"),
    alternative = c("two.sided", "greater", "less"),
    nPerm = 999L,
    status = "ok"
) {
    genome_results <- .require_genome_fit(x)
    tree <- .require_genome_tree_object(object, tree = tree)
    correlation_method <- .normalize_phylo_correlation_method(correlationMethod)
    alternative <- .normalize_phylo_alternative(alternative)
    n_perm <- .normalize_phylo_permutations(nPerm)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    keep <- rep(TRUE, nrow(genome_results))
    if (!is.null(status_filter) && "status" %in% names(genome_results)) {
        keep <- as.character(genome_results$status) %in% status_filter
    }

    filtered <- genome_results[keep, , drop = FALSE]
    aligned <- .align_tree_values(
        tree = tree,
        genome_ids = filtered$genome_id,
        values = filtered$estimate
    )

    if (length(aligned$values) < 3L) {
        results <- .phylo_signal_row(
            analysis_id = "genome_response_signal",
            feature_id_column = "analysis_id",
            variable_info = variable_info,
            statistic = NA_real_,
            p_value = NA_real_,
            n_available = aligned$nAvailable,
            n_tested = length(aligned$values),
            genome_ids = aligned$genomeIds,
            correlation_method = correlation_method,
            alternative = alternative,
            n_perm = n_perm,
            status = "skipped",
            error_message = "At least three genomes with finite fitted responses are required."
        )
    } else {
        statistic <- .phylogenetic_signal_statistic(
            tree = aligned$tree,
            values = aligned$values,
            method = correlation_method
        )
        permuted <- .permuted_signal_statistics(
            tree = aligned$tree,
            values = aligned$values,
            method = correlation_method,
            n_perm = n_perm
        )
        results <- .phylo_signal_row(
            analysis_id = "genome_response_signal",
            feature_id_column = "analysis_id",
            variable_info = variable_info,
            statistic = statistic,
            p_value = .permutation_p_value(statistic, permuted, alternative),
            n_available = aligned$nAvailable,
            n_tested = length(aligned$values),
            genome_ids = aligned$genomeIds,
            correlation_method = correlation_method,
            alternative = alternative,
            n_perm = n_perm,
            status = if (is.finite(statistic)) "ok" else "skipped",
            error_message = if (is.finite(statistic)) NA_character_ else {
                "The fitted responses did not vary enough to compute a distance correlation."
            }
        )
    }

    MTTKFit(
        results = results,
        info = list(
            backend = "permutation_test",
            model = "genome_phylogenetic_signal",
            featureType = "analysis",
            featureIdColumn = "analysis_id",
            sourceFeatureType = "genome",
            sourceFeatureIdColumn = "genome_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            correlationMethod = correlation_method,
            alternative = alternative,
            nPermutations = n_perm,
            signalStatistic = "cor(phylogenetic_distance, abs(response_difference))",
            treeSource = "genomeTree(object)"
        ),
        models = list()
    )
}

#' Test for KO-Level Phylogenetic Signal in Genome-Specific Responses
#'
#' `fitKOPhylogeneticSignal()` follows [fitKORandomSlopeModel()] with a
#' tree-based summary for each KO. It asks whether genome-specific KO response
#' estimates tend to be more similar among closely related genomes.
#'
#' @param x A KO random-slope `MTTKFit`, from [fitKORandomSlopeModel()].
#' @param object An `MTTKExperiment` supplying `genomeTree(object)`.
#' @param tree Optional `ape::phylo` object. When `NULL`, `genomeTree(object)`
#'   is used.
#' @param correlationMethod Correlation used to quantify the association
#'   between phylogenetic distance and response difference.
#' @param alternative Alternative hypothesis for the permutation p-value.
#' @param nPerm Number of tip-label permutations used to estimate the p-value.
#' @param minGenomes Minimum number of genomes with finite KO genome effects
#'   required for a KO to receive an `"ok"` result.
#' @param status Optional KO-fit status filter applied before the signal test.
#'   The default keeps only KO rows with `status == "ok"`.
#'
#' @return An `MTTKFit` with one row per KO.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     ko_fit <- fitKORandomSlopeModel(x, variable = "condition")
#'     fitKOPhylogeneticSignal(ko_fit, x, nPerm = 99)
#' }
#'
#' @export
fitKOPhylogeneticSignal <- function(
    x,
    object,
    tree = NULL,
    correlationMethod = c("spearman", "pearson"),
    alternative = c("two.sided", "greater", "less"),
    nPerm = 999L,
    minGenomes = 3L,
    status = "ok"
) {
    parsed <- .require_ko_random_slope_fit(x)
    tree <- .require_genome_tree_object(object, tree = tree)
    correlation_method <- .normalize_phylo_correlation_method(correlationMethod)
    alternative <- .normalize_phylo_alternative(alternative)
    n_perm <- .normalize_phylo_permutations(nPerm)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    ko_results <- as.data.frame(parsed$koResults, stringsAsFactors = FALSE)
    ko_results$ko_id <- rownames(parsed$koResults)
    ko_status <- stats::setNames(
        if ("status" %in% names(ko_results)) as.character(ko_results$status) else {
            rep(NA_character_, nrow(ko_results))
        },
        ko_results$ko_id
    )

    group_effects <- as.data.frame(parsed$groupEffects, stringsAsFactors = FALSE)
    group_effects$ko_id <- as.character(group_effects$ko_id)
    group_effects$genome_id <- as.character(group_effects$genome_id)
    group_split <- split(group_effects, group_effects$ko_id)
    ko_ids <- rownames(parsed$koResults)

    results <- do.call(
        rbind,
        lapply(ko_ids, function(ko_id) {
            data <- group_split[[ko_id]]
            if (is.null(data)) {
                return(.phylo_signal_row(
                    analysis_id = ko_id,
                    feature_id_column = "ko_id",
                    variable_info = variable_info,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = 0L,
                    n_tested = 0L,
                    genome_ids = character(),
                    correlation_method = correlation_method,
                    alternative = alternative,
                    n_perm = n_perm,
                    status = "skipped",
                    error_message = "No genome-specific KO effects were available."
                ))
            }

            if (!is.null(status_filter) &&
                !(ko_status[[ko_id]] %in% status_filter)) {
                return(.phylo_signal_row(
                    analysis_id = ko_id,
                    feature_id_column = "ko_id",
                    variable_info = variable_info,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = nrow(data),
                    n_tested = 0L,
                    genome_ids = character(),
                    correlation_method = correlation_method,
                    alternative = alternative,
                    n_perm = n_perm,
                    status = "skipped",
                    error_message = "The KO did not pass the requested status filter."
                ))
            }

            aligned <- .align_tree_values(
                tree = tree,
                genome_ids = data$genome_id,
                values = data$conditional_effect_estimate
            )

            if (length(aligned$values) < min_genomes) {
                return(.phylo_signal_row(
                    analysis_id = ko_id,
                    feature_id_column = "ko_id",
                    variable_info = variable_info,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    n_available = aligned$nAvailable,
                    n_tested = length(aligned$values),
                    genome_ids = aligned$genomeIds,
                    correlation_method = correlation_method,
                    alternative = alternative,
                    n_perm = n_perm,
                    status = "skipped",
                    error_message = paste0(
                        "At least ",
                        min_genomes,
                        " genomes with finite KO effects are required."
                    )
                ))
            }

            statistic <- .phylogenetic_signal_statistic(
                tree = aligned$tree,
                values = aligned$values,
                method = correlation_method
            )
            permuted <- .permuted_signal_statistics(
                tree = aligned$tree,
                values = aligned$values,
                method = correlation_method,
                n_perm = n_perm
            )

            .phylo_signal_row(
                analysis_id = ko_id,
                feature_id_column = "ko_id",
                variable_info = variable_info,
                statistic = statistic,
                p_value = .permutation_p_value(statistic, permuted, alternative),
                n_available = aligned$nAvailable,
                n_tested = length(aligned$values),
                genome_ids = aligned$genomeIds,
                correlation_method = correlation_method,
                alternative = alternative,
                n_perm = n_perm,
                status = if (is.finite(statistic)) "ok" else "skipped",
                error_message = if (is.finite(statistic)) NA_character_ else {
                    "The genome-specific KO effects did not vary enough to compute a distance correlation."
                }
            )
        })
    )
    rownames(results) <- ko_ids

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = results,
        info = list(
            backend = "permutation_test",
            model = "ko_phylogenetic_signal",
            featureType = "KO",
            featureIdColumn = "ko_id",
            sourceFeatureType = "KO_genome_effect",
            sourceFeatureIdColumn = "ko_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            correlationMethod = correlation_method,
            alternative = alternative,
            nPermutations = n_perm,
            minGenomes = min_genomes,
            signalStatistic = "cor(phylogenetic_distance, abs(response_difference))",
            treeSource = "genomeTree(object)",
            effectSource = "conditional_effect_estimate"
        ),
        models = list()
    )
}

#' Test for Module-Level Phylogenetic Signal in Genome-Specific Responses
#'
#' `fitModulePhylogeneticSignal()` follows [fitModuleRandomSlopeModel()] with a
#' tree-based summary for each module. It asks whether genome-specific module
#' response estimates tend to be more similar among closely related genomes.
#'
#' @inheritParams fitKOPhylogeneticSignal
#' @param x A module random-slope `MTTKFit`, from [fitModuleRandomSlopeModel()].
#'
#' @return An `MTTKFit` with one row per module.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     module_fit <- fitModuleRandomSlopeModel(x, variable = "condition")
#'     fitModulePhylogeneticSignal(module_fit, x, nPerm = 99)
#' }
#'
#' @export
fitModulePhylogeneticSignal <- function(
    x,
    object,
    tree = NULL,
    correlationMethod = c("spearman", "pearson"),
    alternative = c("two.sided", "greater", "less"),
    nPerm = 999L,
    minGenomes = 3L,
    status = "ok"
) {
    parsed <- .require_feature_random_slope_fit(
        x = x,
        model_name = "module_random_slope_model",
        feature_id_column = "module_id",
        function_label = "fitModuleRandomSlopeModel"
    )
    tree <- .require_genome_tree_object(object, tree = tree)
    correlation_method <- .normalize_phylo_correlation_method(correlationMethod)
    alternative <- .normalize_phylo_alternative(alternative)
    n_perm <- .normalize_phylo_permutations(nPerm)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    results <- .fit_feature_phylogenetic_signal_results(
        parsed = parsed,
        tree = tree,
        feature_id_column = "module_id",
        correlation_method = correlation_method,
        alternative = alternative,
        n_perm = n_perm,
        min_genomes = min_genomes,
        status_filter = status_filter,
        variable_info = variable_info
    )

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = results,
        info = list(
            backend = "permutation_test",
            model = "module_phylogenetic_signal",
            featureType = "module",
            featureIdColumn = "module_id",
            sourceFeatureType = "module_genome_effect",
            sourceFeatureIdColumn = "module_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            correlationMethod = correlation_method,
            alternative = alternative,
            nPermutations = n_perm,
            minGenomes = min_genomes,
            signalStatistic = "cor(phylogenetic_distance, abs(response_difference))",
            treeSource = "genomeTree(object)",
            effectSource = "conditional_effect_estimate"
        ),
        models = list()
    )
}

#' Test for Pathway-Level Phylogenetic Signal in Genome-Specific Responses
#'
#' `fitPathwayPhylogeneticSignal()` follows [fitPathwayRandomSlopeModel()] with
#' a tree-based summary for each pathway. It asks whether genome-specific
#' pathway response estimates tend to be more similar among closely related
#' genomes.
#'
#' @inheritParams fitKOPhylogeneticSignal
#' @param x A pathway random-slope `MTTKFit`, from [fitPathwayRandomSlopeModel()].
#'
#' @return An `MTTKFit` with one row per pathway.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     pathway_fit <- fitPathwayRandomSlopeModel(x, variable = "condition")
#'     fitPathwayPhylogeneticSignal(pathway_fit, x, nPerm = 99)
#' }
#'
#' @export
fitPathwayPhylogeneticSignal <- function(
    x,
    object,
    tree = NULL,
    correlationMethod = c("spearman", "pearson"),
    alternative = c("two.sided", "greater", "less"),
    nPerm = 999L,
    minGenomes = 3L,
    status = "ok"
) {
    parsed <- .require_feature_random_slope_fit(
        x = x,
        model_name = "pathway_random_slope_model",
        feature_id_column = "pathway_id",
        function_label = "fitPathwayRandomSlopeModel"
    )
    tree <- .require_genome_tree_object(object, tree = tree)
    correlation_method <- .normalize_phylo_correlation_method(correlationMethod)
    alternative <- .normalize_phylo_alternative(alternative)
    n_perm <- .normalize_phylo_permutations(nPerm)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    results <- .fit_feature_phylogenetic_signal_results(
        parsed = parsed,
        tree = tree,
        feature_id_column = "pathway_id",
        correlation_method = correlation_method,
        alternative = alternative,
        n_perm = n_perm,
        min_genomes = min_genomes,
        status_filter = status_filter,
        variable_info = variable_info
    )

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = results,
        info = list(
            backend = "permutation_test",
            model = "pathway_phylogenetic_signal",
            featureType = "pathway",
            featureIdColumn = "pathway_id",
            sourceFeatureType = "pathway_genome_effect",
            sourceFeatureIdColumn = "pathway_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            correlationMethod = correlation_method,
            alternative = alternative,
            nPermutations = n_perm,
            minGenomes = min_genomes,
            signalStatistic = "cor(phylogenetic_distance, abs(response_difference))",
            treeSource = "genomeTree(object)",
            effectSource = "conditional_effect_estimate"
        ),
        models = list()
    )
}

#' Scan a Genome Phylogeny for Response-Shifted Clades
#'
#' `scanGenomeClades()` follows [fitGenomeModel()] with a subtree scan. Each
#' internal clade in the rooted genome phylogeny is tested for a shifted mean
#' fitted response relative to the complement of that clade.
#'
#' @param x A genome-level `MTTKFit`, typically returned by [fitGenomeModel()].
#' @param object An `MTTKExperiment` supplying `genomeTree(object)`.
#' @param tree Optional rooted `ape::phylo` object. When `NULL`,
#'   `genomeTree(object)` is used.
#' @param alternative Alternative hypothesis for the permutation p-value.
#' @param nPerm Number of tip-label permutations used to estimate the p-value.
#' @param minTips Minimum subtree size to test.
#' @param maxTips Optional maximum subtree size to test.
#' @param status Optional genome-fit status filter applied before the clade
#'   scan. The default keeps only genome rows with `status == "ok"`.
#'
#' @return An `MTTKFit` with one row per tested clade.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     genome_fit <- fitGenomeModel(x, variable = "condition")
#'     scanGenomeClades(genome_fit, x, nPerm = 99, minTips = 2)
#' }
#'
#' @export
scanGenomeClades <- function(
    x,
    object,
    tree = NULL,
    alternative = c("two.sided", "greater", "less"),
    nPerm = 999L,
    minTips = 2L,
    maxTips = NULL,
    status = "ok"
) {
    genome_results <- .require_genome_fit(x)
    tree <- .require_genome_tree_object(object, tree = tree)
    alternative <- .normalize_phylo_alternative(alternative)
    n_perm <- .normalize_phylo_permutations(nPerm)
    min_tips <- .normalize_tree_min_tips(minTips)
    collapse <- .normalize_tree_collapse(";")
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    keep <- rep(TRUE, nrow(genome_results))
    if (!is.null(status_filter) && "status" %in% names(genome_results)) {
        keep <- as.character(genome_results$status) %in% status_filter
    }

    filtered <- genome_results[keep, , drop = FALSE]
    aligned <- .align_tree_values(
        tree = tree,
        genome_ids = filtered$genome_id,
        values = filtered$estimate
    )

    if (length(aligned$values) < 3L) {
        results <- .clade_scan_row(
            row_id = "node_na",
            variable_info = variable_info,
            estimate = NA_real_,
            p_value = NA_real_,
            clade_id = NA_character_,
            clade_label = NA_character_,
            node_id = NA_integer_,
            feature_id_column = NULL,
            feature_id = NULL,
            member_ids = character(),
            all_ids = character(),
            n_perm = n_perm,
            alternative = alternative,
            status = "skipped",
            error_message = "At least three genomes with finite fitted responses are required."
        )
    } else {
        values <- stats::setNames(aligned$values, aligned$genomeIds)
        clades <- .clade_memberships(
            tree = aligned$tree,
            min_tips = min_tips,
            max_tips = maxTips,
            collapse = collapse
        )

        results <- if (length(clades) == 0L) {
            .clade_scan_row(
                row_id = "node_na",
                variable_info = variable_info,
                estimate = NA_real_,
                p_value = NA_real_,
                clade_id = "node_na",
                clade_label = NA_character_,
                node_id = NA_integer_,
                feature_id_column = NULL,
                feature_id = NULL,
                member_ids = character(),
                all_ids = aligned$genomeIds,
                n_perm = n_perm,
                alternative = alternative,
                status = "skipped",
                error_message = "No clades satisfied the requested subtree-size constraints."
            )
        } else {
            do.call(
                rbind,
                lapply(clades, function(clade) {
                    estimate <- .clade_difference_statistic(values, member_ids = clade$memberGenomeIds)
                    permuted <- .permuted_clade_statistics(
                        values = values,
                        member_ids = clade$memberGenomeIds,
                        n_perm = n_perm
                    )
                    .clade_scan_row(
                        row_id = clade$cladeId,
                        variable_info = variable_info,
                        estimate = estimate,
                        p_value = .permutation_p_value(estimate, permuted, alternative),
                        clade_id = clade$cladeId,
                        clade_label = clade$cladeLabel,
                        node_id = clade$nodeId,
                        feature_id_column = NULL,
                        feature_id = NULL,
                        member_ids = clade$memberGenomeIds,
                        all_ids = aligned$genomeIds,
                        n_perm = n_perm,
                        alternative = alternative,
                        status = if (is.finite(estimate)) "ok" else "skipped",
                        error_message = if (is.finite(estimate)) NA_character_ else {
                            "The clade statistic could not be computed."
                        }
                    )
                })
            )
        }
    }

    rownames(results) <- if ("clade_id" %in% names(results)) {
        ifelse(is.na(results$clade_id), rownames(results), as.character(results$clade_id))
    } else {
        rownames(results)
    }

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = results,
        info = list(
            backend = "permutation_test",
            model = "genome_clade_scan",
            featureType = "genome_clade",
            featureIdColumn = "clade_id",
            sourceFeatureType = "genome",
            sourceFeatureIdColumn = "genome_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            alternative = alternative,
            nPermutations = n_perm,
            minTips = min_tips,
            maxTips = maxTips,
            treeSource = "genomeTree(object)",
            effectStatistic = "mean(clade_response) - mean(complement_response)"
        ),
        models = list()
    )
}

#' Scan a Genome Phylogeny for KO-Specific Responding Clades
#'
#' `scanKOClades()` follows [fitKORandomSlopeModel()] with a subtree scan for
#' each KO. Each internal clade in the rooted genome phylogeny is tested for a
#' shifted mean KO genome effect relative to the complement of that clade.
#'
#' @param x A KO random-slope `MTTKFit`, from [fitKORandomSlopeModel()].
#' @param object An `MTTKExperiment` supplying `genomeTree(object)`.
#' @param tree Optional rooted `ape::phylo` object. When `NULL`,
#'   `genomeTree(object)` is used.
#' @param alternative Alternative hypothesis for the permutation p-value.
#' @param nPerm Number of tip-label permutations used to estimate the p-value.
#' @param minTips Minimum subtree size to test.
#' @param maxTips Optional maximum subtree size to test.
#' @param minGenomes Minimum number of genomes with finite KO effects required
#'   before clades are scanned for a KO.
#' @param status Optional KO-fit status filter applied before the clade scan.
#'   The default keeps only KO rows with `status == "ok"`.
#'
#' @return An `MTTKFit` with one row per KO/clade combination.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     ko_fit <- fitKORandomSlopeModel(x, variable = "condition")
#'     scanKOClades(ko_fit, x, nPerm = 99, minTips = 2)
#' }
#'
#' @export
scanKOClades <- function(
    x,
    object,
    tree = NULL,
    alternative = c("two.sided", "greater", "less"),
    nPerm = 999L,
    minTips = 2L,
    maxTips = NULL,
    minGenomes = 3L,
    status = "ok"
) {
    parsed <- .require_ko_random_slope_fit(x)
    tree <- .require_genome_tree_object(object, tree = tree)
    alternative <- .normalize_phylo_alternative(alternative)
    n_perm <- .normalize_phylo_permutations(nPerm)
    min_tips <- .normalize_tree_min_tips(minTips)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    collapse <- .normalize_tree_collapse(";")
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    ko_results <- as.data.frame(parsed$koResults, stringsAsFactors = FALSE)
    ko_results$ko_id <- rownames(parsed$koResults)
    ko_status <- stats::setNames(
        if ("status" %in% names(ko_results)) as.character(ko_results$status) else {
            rep(NA_character_, nrow(ko_results))
        },
        ko_results$ko_id
    )

    group_effects <- as.data.frame(parsed$groupEffects, stringsAsFactors = FALSE)
    group_effects$ko_id <- as.character(group_effects$ko_id)
    group_effects$genome_id <- as.character(group_effects$genome_id)
    group_split <- split(group_effects, group_effects$ko_id)
    ko_ids <- rownames(parsed$koResults)

    results <- do.call(
        rbind,
        lapply(ko_ids, function(ko_id) {
            data <- group_split[[ko_id]]

            if (is.null(data)) {
                return(.clade_scan_row(
                    row_id = paste0(ko_id, "::node_na"),
                    variable_info = variable_info,
                    estimate = NA_real_,
                    p_value = NA_real_,
                    clade_id = NA_character_,
                    clade_label = NA_character_,
                    node_id = NA_integer_,
                    feature_id_column = "ko_id",
                    feature_id = ko_id,
                    member_ids = character(),
                    all_ids = character(),
                    n_perm = n_perm,
                    alternative = alternative,
                    status = "skipped",
                    error_message = "No genome-specific KO effects were available."
                ))
            }

            if (!is.null(status_filter) &&
                !(ko_status[[ko_id]] %in% status_filter)) {
                return(.clade_scan_row(
                    row_id = paste0(ko_id, "::node_na"),
                    variable_info = variable_info,
                    estimate = NA_real_,
                    p_value = NA_real_,
                    clade_id = NA_character_,
                    clade_label = NA_character_,
                    node_id = NA_integer_,
                    feature_id_column = "ko_id",
                    feature_id = ko_id,
                    member_ids = character(),
                    all_ids = character(),
                    n_perm = n_perm,
                    alternative = alternative,
                    status = "skipped",
                    error_message = "The KO did not pass the requested status filter."
                ))
            }

            aligned <- .align_tree_values(
                tree = tree,
                genome_ids = data$genome_id,
                values = data$conditional_effect_estimate
            )

            if (length(aligned$values) < min_genomes) {
                return(.clade_scan_row(
                    row_id = paste0(ko_id, "::node_na"),
                    variable_info = variable_info,
                    estimate = NA_real_,
                    p_value = NA_real_,
                    clade_id = NA_character_,
                    clade_label = NA_character_,
                    node_id = NA_integer_,
                    feature_id_column = "ko_id",
                    feature_id = ko_id,
                    member_ids = character(),
                    all_ids = aligned$genomeIds,
                    n_perm = n_perm,
                    alternative = alternative,
                    status = "skipped",
                    error_message = paste0(
                        "At least ",
                        min_genomes,
                        " genomes with finite KO effects are required."
                    )
                ))
            }

            values <- stats::setNames(aligned$values, aligned$genomeIds)
            clades <- .clade_memberships(
                tree = aligned$tree,
                min_tips = min_tips,
                max_tips = maxTips,
                collapse = collapse
            )

            if (length(clades) == 0L) {
                return(.clade_scan_row(
                    row_id = paste0(ko_id, "::node_na"),
                    variable_info = variable_info,
                    estimate = NA_real_,
                    p_value = NA_real_,
                    clade_id = "node_na",
                    clade_label = NA_character_,
                    node_id = NA_integer_,
                    feature_id_column = "ko_id",
                    feature_id = ko_id,
                    member_ids = character(),
                    all_ids = aligned$genomeIds,
                    n_perm = n_perm,
                    alternative = alternative,
                    status = "skipped",
                    error_message = "No clades satisfied the requested subtree-size constraints."
                ))
            }

            do.call(
                rbind,
                lapply(clades, function(clade) {
                    estimate <- .clade_difference_statistic(values, member_ids = clade$memberGenomeIds)
                    permuted <- .permuted_clade_statistics(
                        values = values,
                        member_ids = clade$memberGenomeIds,
                        n_perm = n_perm
                    )

                    .clade_scan_row(
                        row_id = paste0(ko_id, "::", clade$cladeId),
                        variable_info = variable_info,
                        estimate = estimate,
                        p_value = .permutation_p_value(estimate, permuted, alternative),
                        clade_id = clade$cladeId,
                        clade_label = clade$cladeLabel,
                        node_id = clade$nodeId,
                        feature_id_column = "ko_id",
                        feature_id = ko_id,
                        member_ids = clade$memberGenomeIds,
                        all_ids = aligned$genomeIds,
                        n_perm = n_perm,
                        alternative = alternative,
                        status = if (is.finite(estimate)) "ok" else "skipped",
                        error_message = if (is.finite(estimate)) NA_character_ else {
                            "The clade statistic could not be computed."
                        }
                    )
                })
            )
        })
    )
    results$ko_clade_id <- rownames(results)
    results <- results[, c("ko_clade_id", setdiff(names(results), "ko_clade_id"))]

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = results,
        info = list(
            backend = "permutation_test",
            model = "ko_clade_scan",
            featureType = "KO_clade",
            featureIdColumn = "ko_clade_id",
            sourceFeatureType = "KO_genome_effect",
            sourceFeatureIdColumn = "ko_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            alternative = alternative,
            nPermutations = n_perm,
            minTips = min_tips,
            maxTips = maxTips,
            minGenomes = min_genomes,
            treeSource = "genomeTree(object)",
            effectSource = "conditional_effect_estimate",
            effectStatistic = "mean(clade_response) - mean(complement_response)"
        ),
        models = list()
    )
}

#' Scan a Genome Phylogeny for Module-Specific Responding Clades
#'
#' `scanModuleClades()` follows [fitModuleRandomSlopeModel()] with a subtree
#' scan for each module. Each internal clade in the rooted genome phylogeny is
#' tested for a shifted mean module genome effect relative to the complement of
#' that clade.
#'
#' @inheritParams scanKOClades
#' @param x A module random-slope `MTTKFit`, from [fitModuleRandomSlopeModel()].
#'
#' @return An `MTTKFit` with one row per module/clade combination.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     module_fit <- fitModuleRandomSlopeModel(x, variable = "condition")
#'     scanModuleClades(module_fit, x, nPerm = 99, minTips = 2)
#' }
#'
#' @export
scanModuleClades <- function(
    x,
    object,
    tree = NULL,
    alternative = c("two.sided", "greater", "less"),
    nPerm = 999L,
    minTips = 2L,
    maxTips = NULL,
    minGenomes = 3L,
    status = "ok"
) {
    parsed <- .require_feature_random_slope_fit(
        x = x,
        model_name = "module_random_slope_model",
        feature_id_column = "module_id",
        function_label = "fitModuleRandomSlopeModel"
    )
    tree <- .require_genome_tree_object(object, tree = tree)
    alternative <- .normalize_phylo_alternative(alternative)
    n_perm <- .normalize_phylo_permutations(nPerm)
    min_tips <- .normalize_tree_min_tips(minTips)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    results <- .scan_feature_clades_results(
        parsed = parsed,
        tree = tree,
        feature_id_column = "module_id",
        result_id_column = "module_clade_id",
        alternative = alternative,
        n_perm = n_perm,
        min_tips = min_tips,
        max_tips = maxTips,
        min_genomes = min_genomes,
        status_filter = status_filter,
        variable_info = variable_info
    )

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = results,
        info = list(
            backend = "permutation_test",
            model = "module_clade_scan",
            featureType = "module_clade",
            featureIdColumn = "module_clade_id",
            sourceFeatureType = "module_genome_effect",
            sourceFeatureIdColumn = "module_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            alternative = alternative,
            nPermutations = n_perm,
            minTips = min_tips,
            maxTips = maxTips,
            minGenomes = min_genomes,
            treeSource = "genomeTree(object)",
            effectSource = "conditional_effect_estimate",
            effectStatistic = "mean(clade_response) - mean(complement_response)"
        ),
        models = list()
    )
}

#' Scan a Genome Phylogeny for Pathway-Specific Responding Clades
#'
#' `scanPathwayClades()` follows [fitPathwayRandomSlopeModel()] with a subtree
#' scan for each pathway. Each internal clade in the rooted genome phylogeny is
#' tested for a shifted mean pathway genome effect relative to the complement of
#' that clade.
#'
#' @inheritParams scanKOClades
#' @param x A pathway random-slope `MTTKFit`, from [fitPathwayRandomSlopeModel()].
#'
#' @return An `MTTKFit` with one row per pathway/clade combination.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     pathway_fit <- fitPathwayRandomSlopeModel(x, variable = "condition")
#'     scanPathwayClades(pathway_fit, x, nPerm = 99, minTips = 2)
#' }
#'
#' @export
scanPathwayClades <- function(
    x,
    object,
    tree = NULL,
    alternative = c("two.sided", "greater", "less"),
    nPerm = 999L,
    minTips = 2L,
    maxTips = NULL,
    minGenomes = 3L,
    status = "ok"
) {
    parsed <- .require_feature_random_slope_fit(
        x = x,
        model_name = "pathway_random_slope_model",
        feature_id_column = "pathway_id",
        function_label = "fitPathwayRandomSlopeModel"
    )
    tree <- .require_genome_tree_object(object, tree = tree)
    alternative <- .normalize_phylo_alternative(alternative)
    n_perm <- .normalize_phylo_permutations(nPerm)
    min_tips <- .normalize_tree_min_tips(minTips)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    results <- .scan_feature_clades_results(
        parsed = parsed,
        tree = tree,
        feature_id_column = "pathway_id",
        result_id_column = "pathway_clade_id",
        alternative = alternative,
        n_perm = n_perm,
        min_tips = min_tips,
        max_tips = maxTips,
        min_genomes = min_genomes,
        status_filter = status_filter,
        variable_info = variable_info
    )

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = results,
        info = list(
            backend = "permutation_test",
            model = "pathway_clade_scan",
            featureType = "pathway_clade",
            featureIdColumn = "pathway_clade_id",
            sourceFeatureType = "pathway_genome_effect",
            sourceFeatureIdColumn = "pathway_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            alternative = alternative,
            nPermutations = n_perm,
            minTips = min_tips,
            maxTips = maxTips,
            minGenomes = min_genomes,
            treeSource = "genomeTree(object)",
            effectSource = "conditional_effect_estimate",
            effectStatistic = "mean(clade_response) - mean(complement_response)"
        ),
        models = list()
    )
}

#' Fit a Genome-Level Phylogenetic GLS Model
#'
#' `fitGenomePhylogeneticGLS()` fits an intercept-only phylogenetic GLS model to
#' the genome-level fitted response estimates returned by [fitGenomeModel()].
#' The covariance among genomes is defined by the rooted genome phylogeny under
#' a Brownian-motion model.
#'
#' This workflow is useful when the question is "what is the phylogenetically
#' corrected mean genome response?" rather than "is there overall phylogenetic
#' signal?" The response is the fitted genome-level effect estimate, not the raw
#' RNA count.
#'
#' @param x A genome-level `MTTKFit`, typically returned by [fitGenomeModel()].
#' @param object An `MTTKExperiment` supplying `genomeTree(object)`.
#' @param tree Optional `ape::phylo` object. When `NULL`, `genomeTree(object)`
#'   is used.
#' @param method Fitting method passed to `nlme::gls()`.
#' @param useWeights Logical; if `TRUE`, the phylogenetic GLS uses the fitted
#'   effect standard errors as inverse-variance weights.
#' @param status Optional genome-fit status filter applied before fitting. The
#'   default keeps only genome rows with `status == "ok"`.
#' @param keepFits Logical; should the backend `gls` fit object be stored in the
#'   returned `MTTKFit`?
#'
#' @return An `MTTKFit` with one row summarizing the phylogenetically corrected
#'   mean genome response.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     genome_fit <- fitGenomeModel(x, variable = "condition")
#'     fitGenomePhylogeneticGLS(genome_fit, x)
#' }
#'
#' @export
fitGenomePhylogeneticGLS <- function(
    x,
    object,
    tree = NULL,
    method = c("REML", "ML"),
    useWeights = TRUE,
    status = "ok",
    keepFits = FALSE
) {
    genome_results <- .require_genome_fit(x)
    tree <- .require_genome_tree_object(object, tree = tree)
    gls_method <- .normalize_phylo_gls_method(method)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    keep <- rep(TRUE, nrow(genome_results))
    if (!is.null(status_filter) && "status" %in% names(genome_results)) {
        keep <- as.character(genome_results$status) %in% status_filter
    }

    filtered <- as.data.frame(genome_results[keep, , drop = FALSE], stringsAsFactors = FALSE)
    aligned <- .align_tree_data_frame(
        tree = tree,
        data = data.frame(
            genome_id = filtered$genome_id,
            response = filtered$estimate,
            std_error = filtered$std_error,
            stringsAsFactors = FALSE
        ),
        genome_column = "genome_id",
        required_columns = c("response", "std_error")
    )

    models <- list()
    if (nrow(aligned$data) < 3L) {
        results <- .phylo_gls_row(
            analysis_id = "genome_phylogenetic_gls",
            feature_id_column = "analysis_id",
            variable_info = variable_info,
            estimate = NA_real_,
            std_error = NA_real_,
            statistic = NA_real_,
            p_value = NA_real_,
            n_available = aligned$nAvailable,
            n_tested = nrow(aligned$data),
            genome_ids = aligned$genomeIds,
            gls_method = gls_method,
            use_weights = useWeights,
            status = "skipped",
            error_message = "At least three genomes with finite fitted responses are required."
        )
    } else {
        fitted <- tryCatch(
            .fit_brownian_gls(
                aligned_data = aligned$data,
                tree = aligned$tree,
                method = gls_method,
                use_weights = useWeights
            ),
            error = function(e) e
        )

        if (inherits(fitted, "error")) {
            results <- .phylo_gls_row(
                analysis_id = "genome_phylogenetic_gls",
                feature_id_column = "analysis_id",
                variable_info = variable_info,
                estimate = NA_real_,
                std_error = NA_real_,
                statistic = NA_real_,
                p_value = NA_real_,
                n_available = aligned$nAvailable,
                n_tested = nrow(aligned$data),
                genome_ids = aligned$genomeIds,
                gls_method = gls_method,
                use_weights = useWeights,
                status = "error",
                error_message = conditionMessage(fitted)
            )
        } else {
            if (isTRUE(keepFits)) {
                models$genome_phylogenetic_gls <- fitted$fit
            }
            results <- .phylo_gls_row(
                analysis_id = "genome_phylogenetic_gls",
                feature_id_column = "analysis_id",
                variable_info = variable_info,
                estimate = fitted$estimate,
                std_error = fitted$std_error,
                statistic = fitted$statistic,
                p_value = fitted$p_value,
                n_available = aligned$nAvailable,
                n_tested = nrow(aligned$data),
                genome_ids = aligned$genomeIds,
                gls_method = gls_method,
                use_weights = useWeights,
                status = "ok"
            )
        }
    }

    ok_rows <- !is.na(results$p_value) & results$status == "ok"
    results$q_value <- NA_real_
    results$q_value[ok_rows] <- stats::p.adjust(results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = results,
        info = list(
            backend = "phylogenetic_gls",
            model = "genome_phylogenetic_gls",
            featureType = "analysis",
            featureIdColumn = "analysis_id",
            sourceFeatureType = "genome",
            sourceFeatureIdColumn = "genome_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            glsMethod = gls_method,
            correlationModel = "Brownian",
            useMeasurementErrorWeights = useWeights,
            treeSource = "genomeTree(object)",
            responseSource = "estimate"
        ),
        models = models
    )
}

#' Fit a KO-Level Phylogenetic GLS Model
#'
#' `fitKOPhylogeneticGLS()` fits one intercept-only phylogenetic GLS model per
#' KO using the genome-specific KO effects returned by [fitKORandomSlopeModel()].
#' The covariance among genomes is defined by the rooted genome phylogeny under
#' a Brownian-motion model.
#'
#' This workflow is useful when the question is "what is the phylogenetically
#' corrected mean KO response across genomes?" rather than whether there is
#' overall phylogenetic signal or a localized responding clade.
#'
#' @param x A KO random-slope `MTTKFit`, from [fitKORandomSlopeModel()].
#' @inheritParams fitGenomePhylogeneticGLS
#' @param minGenomes Minimum number of genomes with finite feature effects
#'   required for a KO to receive an `"ok"` result.
#'
#' @return An `MTTKFit` with one row per KO.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     ko_fit <- fitKORandomSlopeModel(x, variable = "condition")
#'     fitKOPhylogeneticGLS(ko_fit, x)
#' }
#'
#' @export
fitKOPhylogeneticGLS <- function(
    x,
    object,
    tree = NULL,
    method = c("REML", "ML"),
    useWeights = TRUE,
    minGenomes = 3L,
    status = "ok",
    keepFits = FALSE
) {
    parsed <- .require_feature_random_slope_fit(
        x = x,
        model_name = "ko_random_slope_model",
        feature_id_column = "ko_id",
        function_label = "fitKORandomSlopeModel"
    )
    tree <- .require_genome_tree_object(object, tree = tree)
    gls_method <- .normalize_phylo_gls_method(method)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    fitted <- .fit_feature_phylogenetic_gls_results(
        parsed = parsed,
        tree = tree,
        feature_id_column = "ko_id",
        gls_method = gls_method,
        use_weights = useWeights,
        min_genomes = min_genomes,
        status_filter = status_filter,
        variable_info = variable_info,
        keep_fits = keepFits
    )

    ok_rows <- !is.na(fitted$results$p_value) & fitted$results$status == "ok"
    fitted$results$q_value <- NA_real_
    fitted$results$q_value[ok_rows] <- stats::p.adjust(fitted$results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = fitted$results,
        info = list(
            backend = "phylogenetic_gls",
            model = "ko_phylogenetic_gls",
            featureType = "KO",
            featureIdColumn = "ko_id",
            sourceFeatureType = "KO_genome_effect",
            sourceFeatureIdColumn = "ko_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            glsMethod = gls_method,
            correlationModel = "Brownian",
            useMeasurementErrorWeights = useWeights,
            minGenomes = min_genomes,
            treeSource = "genomeTree(object)",
            responseSource = "conditional_effect_estimate"
        ),
        models = fitted$models
    )
}

#' Fit a Module-Level Phylogenetic GLS Model
#'
#' `fitModulePhylogeneticGLS()` fits one intercept-only phylogenetic GLS model
#' per module using the genome-specific module effects returned by
#' [fitModuleRandomSlopeModel()].
#'
#' @param x A module random-slope `MTTKFit`, from [fitModuleRandomSlopeModel()].
#' @inheritParams fitKOPhylogeneticGLS
#'
#' @return An `MTTKFit` with one row per module.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     module_fit <- fitModuleRandomSlopeModel(x, variable = "condition")
#'     fitModulePhylogeneticGLS(module_fit, x)
#' }
#'
#' @export
fitModulePhylogeneticGLS <- function(
    x,
    object,
    tree = NULL,
    method = c("REML", "ML"),
    useWeights = TRUE,
    minGenomes = 3L,
    status = "ok",
    keepFits = FALSE
) {
    parsed <- .require_feature_random_slope_fit(
        x = x,
        model_name = "module_random_slope_model",
        feature_id_column = "module_id",
        function_label = "fitModuleRandomSlopeModel"
    )
    tree <- .require_genome_tree_object(object, tree = tree)
    gls_method <- .normalize_phylo_gls_method(method)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    fitted <- .fit_feature_phylogenetic_gls_results(
        parsed = parsed,
        tree = tree,
        feature_id_column = "module_id",
        gls_method = gls_method,
        use_weights = useWeights,
        min_genomes = min_genomes,
        status_filter = status_filter,
        variable_info = variable_info,
        keep_fits = keepFits
    )

    ok_rows <- !is.na(fitted$results$p_value) & fitted$results$status == "ok"
    fitted$results$q_value <- NA_real_
    fitted$results$q_value[ok_rows] <- stats::p.adjust(fitted$results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = fitted$results,
        info = list(
            backend = "phylogenetic_gls",
            model = "module_phylogenetic_gls",
            featureType = "module",
            featureIdColumn = "module_id",
            sourceFeatureType = "module_genome_effect",
            sourceFeatureIdColumn = "module_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            glsMethod = gls_method,
            correlationModel = "Brownian",
            useMeasurementErrorWeights = useWeights,
            minGenomes = min_genomes,
            treeSource = "genomeTree(object)",
            responseSource = "conditional_effect_estimate"
        ),
        models = fitted$models
    )
}

#' Fit a Pathway-Level Phylogenetic GLS Model
#'
#' `fitPathwayPhylogeneticGLS()` fits one intercept-only phylogenetic GLS model
#' per pathway using the genome-specific pathway effects returned by
#' [fitPathwayRandomSlopeModel()].
#'
#' @param x A pathway random-slope `MTTKFit`, from [fitPathwayRandomSlopeModel()].
#' @inheritParams fitKOPhylogeneticGLS
#'
#' @return An `MTTKFit` with one row per pathway.
#'
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'     x <- makeShowcaseMTTKExperiment()
#'     pathway_fit <- fitPathwayRandomSlopeModel(x, variable = "condition")
#'     fitPathwayPhylogeneticGLS(pathway_fit, x)
#' }
#'
#' @export
fitPathwayPhylogeneticGLS <- function(
    x,
    object,
    tree = NULL,
    method = c("REML", "ML"),
    useWeights = TRUE,
    minGenomes = 3L,
    status = "ok",
    keepFits = FALSE
) {
    parsed <- .require_feature_random_slope_fit(
        x = x,
        model_name = "pathway_random_slope_model",
        feature_id_column = "pathway_id",
        function_label = "fitPathwayRandomSlopeModel"
    )
    tree <- .require_genome_tree_object(object, tree = tree)
    gls_method <- .normalize_phylo_gls_method(method)
    min_genomes <- .normalize_tree_min_genomes(minGenomes)
    status_filter <- .normalize_meta_status_filter(status)
    variable_info <- .meta_variable_info(x)

    fitted <- .fit_feature_phylogenetic_gls_results(
        parsed = parsed,
        tree = tree,
        feature_id_column = "pathway_id",
        gls_method = gls_method,
        use_weights = useWeights,
        min_genomes = min_genomes,
        status_filter = status_filter,
        variable_info = variable_info,
        keep_fits = keepFits
    )

    ok_rows <- !is.na(fitted$results$p_value) & fitted$results$status == "ok"
    fitted$results$q_value <- NA_real_
    fitted$results$q_value[ok_rows] <- stats::p.adjust(fitted$results$p_value[ok_rows], method = "BH")

    MTTKFit(
        results = fitted$results,
        info = list(
            backend = "phylogenetic_gls",
            model = "pathway_phylogenetic_gls",
            featureType = "pathway",
            featureIdColumn = "pathway_id",
            sourceFeatureType = "pathway_genome_effect",
            sourceFeatureIdColumn = "pathway_id",
            sourceFitModel = fitInfo(x)$model,
            sourceBackend = fitInfo(x)$backend,
            sourceRandomEffects = fitInfo(x)$randomEffects,
            sourceGroupEffectColumn = fitInfo(x)$groupEffectColumn,
            variable = variable_info$variable,
            variableType = variable_info$type,
            referenceLevel = variable_info$referenceLevel,
            contrastLevel = variable_info$contrastLevel,
            effectLabel = variable_info$effectLabel,
            glsMethod = gls_method,
            correlationModel = "Brownian",
            useMeasurementErrorWeights = useWeights,
            minGenomes = min_genomes,
            treeSource = "genomeTree(object)",
            responseSource = "conditional_effect_estimate"
        ),
        models = fitted$models
    )
}

#' @import ggplot2
NULL

# ── Shared helpers ────────────────────────────────────────────────────────────

.plot_fit_table <- function(fit, term, status) {
    tbl <- fitTable(fit, term = term, status = status)
    tbl <- as.data.frame(tbl)
    tbl$.feature_id <- rownames(tbl)
    tbl
}

.check_plot_columns <- function(tbl, required, fn) {
    missing_cols <- setdiff(required, names(tbl))
    if (length(missing_cols) > 0L) {
        stop(
            fn, "() requires column(s) not found in the fit: ",
            paste(missing_cols, collapse = ", "), ".",
            call. = FALSE
        )
    }
}

.classify_direction <- function(estimate, q_value, alpha) {
    sig       <- !is.na(q_value) & q_value <= alpha
    direction <- rep("ns", length(estimate))
    direction[sig & !is.na(estimate) & estimate > 0] <- "up"
    direction[sig & !is.na(estimate) & estimate < 0] <- "down"
    factor(direction, levels = c("up", "down", "ns"))
}

# ── Tree layout (ape only, no ggtree) ─────────────────────────────────────────

# Computes a rectangular phylogram layout for ggplot2 from an ape::phylo tree.
# Returns list(segs, tips) where segs is a data.frame of branch segments and
# tips is a data.frame of tip positions with columns label, x, y.
.tree_layout <- function(tree) {
    tree    <- ape::ladderize(tree)
    n_tip   <- length(tree$tip.label)
    n_node  <- tree$Nnode
    n_total <- n_tip + n_node

    # x-coordinates: distance from root via edge lengths (phylogram) or equal
    x_all <- if (!is.null(tree$edge.length)) {
        ape::node.depth.edgelength(tree)
    } else {
        d <- ape::node.depth(tree)
        as.numeric(max(d) - d)
    }

    # Normalize to [0, 1]
    x_max <- max(x_all)
    if (x_max > 0) x_all <- x_all / x_max

    # y-coordinates: sequential tip positions via preorder traversal
    y_all   <- numeric(n_total)
    tip_idx <- 0L
    .assign_y <- function(nd) {
        children <- tree$edge[tree$edge[, 1L] == nd, 2L]
        if (length(children) == 0L) {
            tip_idx   <<- tip_idx + 1L
            y_all[nd] <<- tip_idx
        } else {
            for (ch in children) .assign_y(ch)
            y_all[nd] <<- mean(y_all[children])
        }
    }
    .assign_y(n_tip + 1L)

    # Horizontal arm segments (one per edge: parent-x → child-x at child-y)
    h_segs <- data.frame(
        x    = x_all[tree$edge[, 1L]],
        xend = x_all[tree$edge[, 2L]],
        y    = y_all[tree$edge[, 2L]],
        yend = y_all[tree$edge[, 2L]],
        stringsAsFactors = FALSE
    )

    # Vertical bar segments (one per internal node, spanning its children's y range)
    int_nodes <- seq.int(n_tip + 1L, n_total)
    v_segs <- do.call(rbind, lapply(int_nodes, function(nd) {
        ch_y <- y_all[tree$edge[tree$edge[, 1L] == nd, 2L]]
        data.frame(
            x = x_all[nd], xend = x_all[nd],
            y = min(ch_y),  yend = max(ch_y),
            stringsAsFactors = FALSE
        )
    }))

    segs <- rbind(h_segs, v_segs)

    tips <- data.frame(
        label   = tree$tip.label,
        x       = x_all[seq_len(n_tip)],
        y       = y_all[seq_len(n_tip)],
        stringsAsFactors = FALSE
    )

    list(segs = segs, tips = tips)
}

# ── plotVolcano ───────────────────────────────────────────────────────────────

#' Volcano plot from an MTTKFit
#'
#' `plotVolcano()` draws a volcano plot with the effect estimate on the x-axis
#' and \eqn{-\log_{10}(q)} on the y-axis. Points are coloured by direction and
#' significance, and the top features by significance can be labelled.
#'
#' @param fit An `MTTKFit` object.
#' @param term Optional character scalar selecting a specific model term when
#'   the fit stores multiple terms.
#' @param alpha Significance threshold applied to `q_value`. Default `0.05`.
#' @param label Logical. Whether to label the top features. Default `TRUE`.
#' @param n Integer. Number of features to label, ranked by significance.
#'   Default `10L`.
#' @param status Optional character vector of fit statuses to retain before
#'   plotting. Pass `NULL` to keep all rows. Default `"ok"`.
#' @param colors Named character vector with elements `"up"`, `"down"`, and
#'   `"ns"` specifying point colours. A default palette is used when `NULL`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' fit <- fitKOMixedModel(x, variable = "condition")
#' plotVolcano(fit)
#' plotVolcano(fit, alpha = 0.1, n = 20)
#' }
#'
#' @export
plotVolcano <- function(
    fit,
    term   = NULL,
    alpha  = 0.05,
    label  = TRUE,
    n      = 10L,
    status = "ok",
    colors = NULL
) {
    tbl <- .plot_fit_table(fit, term = term, status = status)
    .check_plot_columns(tbl, c("estimate", "q_value"), "plotVolcano")

    tbl$.x   <- as.numeric(tbl$estimate)
    tbl$.y   <- -log10(as.numeric(tbl$q_value))
    tbl$.dir <- .classify_direction(tbl$.x, as.numeric(tbl$q_value), alpha)

    palette <- c(up = "#E64B35", down = "#4DBBD5", ns = "#BBBBBB")
    if (!is.null(colors)) {
        missing_keys <- setdiff(c("up", "down", "ns"), names(colors))
        if (length(missing_keys) > 0L) {
            stop("'colors' must have named elements: up, down, ns.", call. = FALSE)
        }
        palette <- colors
    }

    p <- ggplot(tbl, aes(x = .x, y = .y, color = .dir)) +
        geom_point(aes(size = .dir), alpha = 0.75, stroke = 0) +
        scale_color_manual(
            values = palette,
            labels = c(
                up   = paste0("Up (q ≤ ", alpha, ")"),
                down = paste0("Down (q ≤ ", alpha, ")"),
                ns   = "Not significant"
            ),
            name = NULL,
            drop = FALSE
        ) +
        scale_size_manual(
            values = c(up = 1.8, down = 1.8, ns = 1.2),
            guide  = "none"
        ) +
        geom_vline(
            xintercept = 0,
            linetype = "dashed", color = "grey60", linewidth = 0.4
        ) +
        geom_hline(
            yintercept = -log10(alpha),
            linetype = "dashed", color = "grey60", linewidth = 0.4
        ) +
        labs(
            x = "Effect estimate",
            y = expression(-log[10](italic(q)))
        ) +
        theme_bw(base_size = 11) +
        theme(legend.position = "top")

    if (isTRUE(label) && n > 0L) {
        finite_y  <- is.finite(tbl$.y)
        label_tbl <- tbl[finite_y, , drop = FALSE]
        label_tbl <- label_tbl[order(-label_tbl$.y), , drop = FALSE]
        label_tbl <- utils::head(label_tbl, n)

        p <- p + geom_text(
            data          = label_tbl,
            aes(label     = .feature_id),
            size          = 2.5,
            vjust         = -0.6,
            show.legend   = FALSE,
            check_overlap = TRUE
        )
    }

    p
}

# ── plotEffects ───────────────────────────────────────────────────────────────

#' Ranked effect plot from an MTTKFit
#'
#' `plotEffects()` draws a horizontal lollipop or bar chart of effect estimates,
#' ranked by estimate magnitude. Error bars (95 % CI) are added when
#' `std_error` is present in the fit.
#'
#' @param fit An `MTTKFit` object.
#' @param term Optional character scalar selecting a specific model term.
#' @param n Integer or `NULL`. Maximum number of features to display, selected
#'   by largest absolute estimate. Pass `NULL` to show all. Default `30L`.
#' @param alpha Significance threshold applied to `q_value` for colouring.
#'   Default `0.05`.
#' @param type `"lollipop"` (default) or `"bar"`.
#' @param status Optional character vector of fit statuses to retain before
#'   plotting. Pass `NULL` to keep all rows. Default `"ok"`.
#' @param colors Named character vector with elements `"significant"` and
#'   `"ns"`. A default palette is used when `NULL`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' fit <- fitKOMixedModel(x, variable = "condition")
#' plotEffects(fit)
#' plotEffects(fit, n = 50, type = "bar")
#' }
#'
#' @export
plotEffects <- function(
    fit,
    term   = NULL,
    n      = 30L,
    alpha  = 0.05,
    type   = c("lollipop", "bar"),
    status = "ok",
    colors = NULL
) {
    type <- match.arg(type)
    tbl  <- .plot_fit_table(fit, term = term, status = status)
    .check_plot_columns(tbl, "estimate", "plotEffects")

    tbl$.x <- as.numeric(tbl$estimate)
    has_se <- "std_error" %in% names(tbl)
    has_q  <- "q_value"   %in% names(tbl)

    if (has_se) tbl$.se <- as.numeric(tbl$std_error)

    if (has_q) {
        q <- as.numeric(tbl$q_value)
        tbl$.sig <- factor(
            ifelse(!is.na(q) & q <= alpha, "significant", "ns"),
            levels = c("significant", "ns")
        )
    } else {
        tbl$.sig <- factor(rep("ns", nrow(tbl)), levels = c("significant", "ns"))
    }

    if (!is.null(n) && n < nrow(tbl)) {
        ord <- order(-abs(tbl$.x), na.last = TRUE)
        tbl <- tbl[utils::head(ord, n), , drop = FALSE]
    }

    tbl <- tbl[order(tbl$.x, na.last = TRUE), , drop = FALSE]
    tbl$.feature_id <- factor(tbl$.feature_id, levels = tbl$.feature_id)

    palette <- c(significant = "#E64B35", ns = "#BBBBBB")
    if (!is.null(colors)) {
        missing_keys <- setdiff(c("significant", "ns"), names(colors))
        if (length(missing_keys) > 0L) {
            stop("'colors' must have named elements: significant, ns.", call. = FALSE)
        }
        palette <- colors
    }

    sig_labels <- c(
        significant = paste0("q ≤ ", alpha),
        ns          = "Not significant"
    )

    p <- ggplot(tbl, aes(x = .x, y = .feature_id, color = .sig))

    if (type == "lollipop") {
        p <- p +
            geom_segment(
                aes(x = 0, xend = .x, yend = .feature_id),
                color = "grey75", linewidth = 0.5
            ) +
            geom_point(size = 2.5)
    } else {
        p <- p +
            geom_col(aes(fill = .sig), color = NA, width = 0.7) +
            scale_fill_manual(
                values = palette,
                labels = sig_labels,
                name   = NULL,
                drop   = FALSE
            )
    }

    if (has_se) {
        p <- p + geom_errorbarh(
            aes(xmin = .x - 1.96 * .se, xmax = .x + 1.96 * .se),
            height = 0.3, linewidth = 0.4
        )
    }

    p <- p +
        scale_color_manual(
            values = palette,
            labels = sig_labels,
            name   = NULL,
            drop   = FALSE
        ) +
        geom_vline(
            xintercept = 0,
            linetype = "dashed", color = "grey60", linewidth = 0.4
        ) +
        labs(x = "Effect estimate", y = NULL) +
        theme_bw(base_size = 11) +
        theme(
            legend.position    = "top",
            panel.grid.major.y = element_blank()
        )

    p
}

# ── plotEffectHeatmap ─────────────────────────────────────────────────────────

#' Heatmap of conditional genome effects from a random-slope MTTKFit
#'
#' `plotEffectHeatmap()` draws a feature-by-genome heatmap of the
#' genome-specific conditional effects stored in a random-slope fit (e.g. from
#' [fitKORandomSlopeModel()], [fitModuleRandomSlopeModel()], or
#' [fitPathwayRandomSlopeModel()]). Rows are features (KO, module, or pathway)
#' and columns are genomes.
#'
#' @param fit An `MTTKFit` returned by a random-slope modelling function that
#'   stores per-genome conditional effects.
#' @param genomes Optional character vector. Subset to these genome IDs before
#'   plotting.
#' @param features Optional character vector. Subset to these feature IDs before
#'   plotting.
#' @param cluster Logical. Whether to cluster rows and columns by hierarchical
#'   clustering. Default `TRUE`.
#' @param midpoint Numeric. Midpoint for the diverging colour scale. Default
#'   `0`.
#' @param colors Optional character vector of length 2 giving the low and high
#'   colours for the diverging scale. Default blue–red.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' ko_fit <- fitKORandomSlopeModel(x, variable = "condition", keepFits = TRUE)
#' plotEffectHeatmap(ko_fit)
#' plotEffectHeatmap(ko_fit, cluster = FALSE)
#' }
#'
#' @export
plotEffectHeatmap <- function(
    fit,
    genomes  = NULL,
    features = NULL,
    cluster  = TRUE,
    midpoint = 0,
    colors   = NULL
) {
    group_eff <- .stored_group_effects(fit)
    group_eff <- as.data.frame(group_eff)

    info       <- fitInfo(fit)
    feat_col   <- if (!is.null(info$featureIdColumn))  as.character(info$featureIdColumn)  else "ko_id"
    genome_col <- if (!is.null(info$groupEffectColumn)) as.character(info$groupEffectColumn) else "genome_id"
    value_col  <- "conditional_effect_estimate"

    .check_plot_columns(group_eff, c(feat_col, genome_col, value_col), "plotEffectHeatmap")

    if (!is.null(genomes)) {
        group_eff <- group_eff[as.character(group_eff[[genome_col]]) %in% genomes, , drop = FALSE]
    }
    if (!is.null(features)) {
        group_eff <- group_eff[as.character(group_eff[[feat_col]]) %in% features, , drop = FALSE]
    }

    if (nrow(group_eff) == 0L) {
        stop(
            "No data remains after subsetting. Check 'genomes' and 'features'.",
            call. = FALSE
        )
    }

    feature_ids <- unique(as.character(group_eff[[feat_col]]))
    genome_ids  <- unique(as.character(group_eff[[genome_col]]))

    if (isTRUE(cluster) && length(feature_ids) > 1L && length(genome_ids) > 1L) {
        mat <- matrix(
            NA_real_,
            nrow     = length(feature_ids),
            ncol     = length(genome_ids),
            dimnames = list(feature_ids, genome_ids)
        )
        for (i in seq_len(nrow(group_eff))) {
            r       <- as.character(group_eff[[feat_col]][i])
            cc      <- as.character(group_eff[[genome_col]][i])
            mat[r, cc] <- as.numeric(group_eff[[value_col]][i])
        }
        mat_fill               <- mat
        mat_fill[is.na(mat_fill)] <- 0
        row_order              <- stats::hclust(stats::dist(mat_fill))$order
        col_order              <- stats::hclust(stats::dist(t(mat_fill)))$order
        feature_levels         <- feature_ids[row_order]
        genome_levels          <- genome_ids[col_order]
    } else {
        feature_levels <- feature_ids
        genome_levels  <- genome_ids
    }

    # Rename dynamic columns to fixed names for aes()
    plot_df            <- group_eff
    names(plot_df)[names(plot_df) == feat_col]   <- ".feature"
    names(plot_df)[names(plot_df) == genome_col] <- ".genome"
    names(plot_df)[names(plot_df) == value_col]  <- ".value"

    plot_df$.feature <- factor(as.character(plot_df$.feature), levels = feature_levels)
    plot_df$.genome  <- factor(as.character(plot_df$.genome),  levels = genome_levels)
    plot_df$.value   <- as.numeric(plot_df$.value)

    abs_max <- max(abs(plot_df$.value), na.rm = TRUE)

    low_col  <- "#4DBBD5"
    high_col <- "#E64B35"
    if (!is.null(colors) && length(colors) >= 2L) {
        low_col  <- colors[[1L]]
        high_col <- colors[[2L]]
    }

    ggplot(plot_df, aes(x = .genome, y = .feature, fill = .value)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_gradient2(
            low      = low_col,
            mid      = "white",
            high     = high_col,
            midpoint = midpoint,
            limits   = c(-abs_max, abs_max),
            name     = "Conditional\neffect",
            na.value = "grey90"
        ) +
        labs(x = NULL, y = NULL) +
        theme_bw(base_size = 10) +
        theme(
            axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
            axis.text.y  = element_text(size = 7),
            panel.border = element_blank(),
            panel.grid   = element_blank()
        )
}

# ── plotGenomeEffects ─────────────────────────────────────────────────────────

#' Tree plot of per-genome effects from a genome-level MTTKFit
#'
#' `plotGenomeEffects()` overlays genome-level effect estimates on a
#' phylogenetic tree. Tips are coloured by direction and scaled by effect
#' magnitude. Supply the same tree that was used for modelling (e.g.
#' `genomeTree(x)`).
#'
#' @param fit An `MTTKFit` returned by [fitGenomeModel()] or a similar
#'   genome-level modelling function. Must contain a `genome_id` column.
#' @param tree An `ape::phylo` object whose tip labels match the `genome_id`
#'   values in `fit`.
#' @param term Optional character scalar for multi-term fits.
#' @param alpha Significance threshold for `q_value` used for direction
#'   colouring. Default `0.05`.
#' @param status Fit status filter. Default `"ok"`.
#' @param tip_label Logical. Whether to draw tip labels. Default `TRUE`.
#' @param colors Named character vector with elements `"up"`, `"down"`, and
#'   `"ns"`. A default palette is used when `NULL`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' genome_fit <- fitGenomeModel(x, variable = "condition")
#' plotGenomeEffects(genome_fit, tree = genomeTree(x))
#' }
#'
#' @export
plotGenomeEffects <- function(
    fit,
    tree,
    term      = NULL,
    alpha     = 0.05,
    status    = "ok",
    tip_label = TRUE,
    colors    = NULL
) {
    tbl <- .plot_fit_table(fit, term = term, status = status)
    .check_plot_columns(tbl, c("estimate", "genome_id"), "plotGenomeEffects")

    layout <- .tree_layout(tree)
    tips   <- layout$tips

    # Join estimates to tree tips
    tbl$genome_id  <- as.character(tbl$genome_id)
    idx            <- match(tips$label, tbl$genome_id)
    tips$.estimate <- as.numeric(tbl$estimate)[idx]
    tips$.q_value  <- if ("q_value" %in% names(tbl)) as.numeric(tbl$q_value)[idx] else NA_real_
    tips$.dir      <- .classify_direction(tips$.estimate, tips$.q_value, alpha)
    tips$.size     <- ifelse(is.na(tips$.estimate), 1, 1 + 3 * (abs(tips$.estimate) / max(abs(tips$.estimate), na.rm = TRUE)))

    palette <- c(up = "#E64B35", down = "#4DBBD5", ns = "#BBBBBB")
    if (!is.null(colors)) {
        missing_keys <- setdiff(c("up", "down", "ns"), names(colors))
        if (length(missing_keys) > 0L) {
            stop("'colors' must have named elements: up, down, ns.", call. = FALSE)
        }
        palette <- colors
    }

    p <- ggplot() +
        geom_segment(
            data = layout$segs,
            aes(x = x, xend = xend, y = y, yend = yend),
            color = "grey50", linewidth = 0.4, na.rm = TRUE
        ) +
        geom_point(
            data = tips,
            aes(x = x, y = y, color = .dir, size = .size)
        ) +
        scale_color_manual(
            values = palette,
            labels = c(
                up   = paste0("Up (q ≤ ", alpha, ")"),
                down = paste0("Down (q ≤ ", alpha, ")"),
                ns   = "Not significant / missing"
            ),
            name = NULL,
            drop = FALSE
        ) +
        scale_size_identity()

    if (isTRUE(tip_label)) {
        p <- p + geom_text(
            data  = tips,
            aes(x = x + 0.02, y = y, label = label),
            hjust = 0, size = 2.2, color = "grey30"
        )
    }

    p <- p +
        scale_x_continuous(expand = expansion(add = c(0.02, if (isTRUE(tip_label)) 0.35 else 0.05))) +
        theme_minimal(base_size = 10) +
        theme(
            axis.title      = element_blank(),
            axis.text       = element_blank(),
            axis.ticks      = element_blank(),
            panel.grid      = element_blank(),
            legend.position = "bottom"
        )

    p
}

# ── plotCladeScan ─────────────────────────────────────────────────────────────

#' Tree plot of significant clades from a clade scan MTTKFit
#'
#' `plotCladeScan()` visualises the result of [scanGenomeClades()] or
#' [scanKOClades()] (filtered to a single KO) by colouring tree tips by the
#' most significant clade they belong to. Tips not in any significant clade are
#' drawn in a neutral colour.
#'
#' When the fit comes from [scanKOClades()], filter to a single KO before
#' passing to this function (e.g. `fit[fit$ko_id == "K00001", ]`).
#'
#' @param fit An `MTTKFit` returned by [scanGenomeClades()] or
#'   [scanKOClades()], optionally pre-filtered to a single feature.
#' @param tree An `ape::phylo` object whose tip labels match the genome IDs in
#'   the clade scan result.
#' @param alpha Significance threshold for `q_value`. Default `0.05`.
#' @param status Fit status filter. Default `"ok"`.
#' @param top_n Integer. Maximum number of significant clades to highlight,
#'   ranked by significance. Default `5L`.
#' @param tip_label Logical. Whether to draw tip labels. Default `TRUE`.
#' @param colors Named character vector with elements `"up"`, `"down"`, and
#'   `"background"`. A default palette is used when `NULL`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' genome_fit  <- fitGenomeModel(x, variable = "condition")
#' clade_fit   <- scanGenomeClades(genome_fit, x)
#' plotCladeScan(clade_fit, tree = genomeTree(x))
#' }
#'
#' @export
plotCladeScan <- function(
    fit,
    tree,
    alpha     = 0.05,
    status    = "ok",
    top_n     = 5L,
    tip_label = TRUE,
    colors    = NULL
) {
    tbl <- .plot_fit_table(fit, term = NULL, status = status)
    .check_plot_columns(tbl, c("estimate", "q_value", "member_genome_ids"), "plotCladeScan")

    # Select significant clades, ordered by q_value
    has_q   <- !is.na(as.numeric(tbl$q_value))
    is_sig  <- has_q & as.numeric(tbl$q_value) <= alpha
    sig     <- tbl[is_sig, , drop = FALSE]
    sig     <- sig[order(as.numeric(sig$q_value)), , drop = FALSE]
    if (!is.null(top_n) && nrow(sig) > top_n) sig <- utils::head(sig, top_n)

    layout <- .tree_layout(tree)
    tips   <- layout$tips

    # Assign each tip to its most significant clade (first match wins)
    tips$.clade_dir <- "background"
    tips$.clade_id  <- NA_character_

    if (nrow(sig) > 0L) {
        for (i in seq_len(nrow(sig))) {
            raw_members <- as.character(sig$member_genome_ids[i])
            members     <- strsplit(raw_members, ";", fixed = TRUE)[[1L]]
            members     <- members[nzchar(members) & !is.na(members)]
            unassigned  <- is.na(tips$.clade_id)
            in_clade    <- tips$label %in% members & unassigned
            if (!any(in_clade)) next
            est       <- as.numeric(sig$estimate[i])
            direction <- if (!is.na(est) && est > 0) "up" else "down"
            tips$.clade_dir[in_clade] <- direction
            tips$.clade_id[in_clade]  <- sig$.feature_id[i]
        }
    }

    tips$.clade_dir <- factor(tips$.clade_dir, levels = c("up", "down", "background"))

    palette <- c(up = "#E64B35", down = "#4DBBD5", background = "#DDDDDD")
    if (!is.null(colors)) {
        missing_keys <- setdiff(c("up", "down", "background"), names(colors))
        if (length(missing_keys) > 0L) {
            stop("'colors' must have named elements: up, down, background.", call. = FALSE)
        }
        palette <- colors
    }

    n_sig_total <- sum(is_sig)
    caption_txt <- if (n_sig_total == 0L) {
        paste0("No significant clades found (q ≤ ", alpha, ")")
    } else {
        paste0(
            n_sig_total, " significant clade(s) found (q ≤ ", alpha, ")",
            if (n_sig_total > top_n) paste0("; showing top ", top_n) else ""
        )
    }

    p <- ggplot() +
        geom_segment(
            data = layout$segs,
            aes(x = x, xend = xend, y = y, yend = yend),
            color = "grey60", linewidth = 0.4, na.rm = TRUE
        ) +
        geom_point(
            data = tips,
            aes(x = x, y = y, color = .clade_dir),
            size = 2.5
        ) +
        scale_color_manual(
            values = palette,
            labels = c(
                up         = paste0("Increased (q ≤ ", alpha, ")"),
                down       = paste0("Decreased (q ≤ ", alpha, ")"),
                background = "Not in significant clade"
            ),
            name = NULL,
            drop = FALSE
        )

    if (isTRUE(tip_label)) {
        p <- p + geom_text(
            data  = tips,
            aes(x = x + 0.02, y = y, label = label, color = .clade_dir),
            hjust = 0, size = 2.2, show.legend = FALSE
        )
    }

    p <- p +
        scale_x_continuous(expand = expansion(add = c(0.02, if (isTRUE(tip_label)) 0.35 else 0.05))) +
        labs(caption = caption_txt) +
        theme_minimal(base_size = 10) +
        theme(
            axis.title      = element_blank(),
            axis.text       = element_blank(),
            axis.ticks      = element_blank(),
            panel.grid      = element_blank(),
            legend.position = "bottom",
            plot.caption    = element_text(size = 8, color = "grey50")
        )

    p
}

# ── Internal helpers ───────────────────────────────────────────────────────────

.auto_sep <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext == "csv") "," else "\t"
}

.read_table_input <- function(x, sep = NULL, row_names = 1L) {
    if (is.character(x) && length(x) == 1L) {
        sep <- if (is.null(sep)) .auto_sep(x) else sep
        rn  <- if (is.null(row_names)) FALSE else row_names
        out <- utils::read.table(
            x,
            header           = TRUE,
            sep              = sep,
            row.names        = rn,
            check.names      = FALSE,
            stringsAsFactors = FALSE
        )
        return(out)
    }
    if (is.data.frame(x) || methods::is(x, "DataFrame")) return(as.data.frame(x))
    if (is.matrix(x)) return(as.data.frame(x))
    stop("Input must be a file path, data.frame, DataFrame, or matrix.", call. = FALSE)
}

.read_matrix_input <- function(x, sep = NULL) {
    df <- .read_table_input(x, sep = sep, row_names = 1L)
    m  <- as.matrix(df)
    storage.mode(m) <- "integer"
    m
}

.read_tree_input <- function(x) {
    if (is.null(x)) return(NULL)
    if (methods::is(x, "phylo")) return(x)
    if (is.character(x) && length(x) == 1L) {
        if (!file.exists(x)) stop("Tree file not found: ", x, call. = FALSE)
        return(ape::read.tree(x))
    }
    stop("'genomeTree' must be a file path (Newick) or an ape::phylo object.", call. = FALSE)
}

.align_samples <- function(count_mat, meta_df, label) {
    common <- intersect(colnames(count_mat), rownames(meta_df))
    if (length(common) == 0L) {
        stop(
            "No shared sample IDs between the count matrix and ", label, ". ",
            "Check that column names of the count matrix match row names of the metadata.",
            call. = FALSE
        )
    }
    if (length(common) < ncol(count_mat)) {
        warning(
            ncol(count_mat) - length(common), " sample(s) absent from ", label,
            " will be dropped.", call. = FALSE
        )
    }
    count_mat[, common, drop = FALSE]
}

# ── KO column detection ────────────────────────────────────────────────────────

.ko_pattern <- "^K[0-9]{5}$"

.known_ko_col_names <- c(
    "ko_id", "ko", "kegg", "kegg_id", "kegg_ko", "koid", "ko_ids",
    "kegg_orthology", "kegg_ortholog", "kegg_annotation",
    "k_number", "ortholog_id", "ko_number"
)

.is_ko_values <- function(values, min_fraction = 0.05) {
    values <- as.character(values)
    non_empty <- values[!is.na(values) & nzchar(values) & values != "NA"]
    if (length(non_empty) == 0L) return(FALSE)
    mean(grepl(.ko_pattern, non_empty)) >= min_fraction
}

#' Detect the KO identifier column in a data frame
#'
#' Detection proceeds in three steps:
#'   1. Exact match to `ko_id_col` (case-sensitive).
#'   2. Case-insensitive match against a list of known KO column names.
#'   3. Full scan: any column where ≥ 50 % of non-empty values match K[0-9]{5}.
#' Returns the matched column name, or NULL with an optional warning.
.detect_ko_column <- function(df, ko_id_col) {
    col_names   <- names(df)
    lower_names <- tolower(col_names)

    # Step 1: explicit column provided
    if (!is.null(ko_id_col) && !is.na(ko_id_col) && nzchar(ko_id_col)) {
        if (ko_id_col %in% col_names) {
            if (.is_ko_values(df[[ko_id_col]])) {
                return(ko_id_col)
            } else {
                warning(
                    "Column '", ko_id_col, "' was found but contains no values ",
                    "matching the KO pattern K[0-9]{5}. Skipping KO link construction.",
                    call. = FALSE
                )
                return(NULL)
            }
        }
        # Try case-insensitive match on the supplied name
        ci_idx <- which(lower_names == tolower(ko_id_col))
        if (length(ci_idx) > 0L) {
            candidate <- col_names[[ci_idx[[1L]]]]
            if (.is_ko_values(df[[candidate]])) {
                message("KO column matched case-insensitively: '", candidate, "'")
                return(candidate)
            }
        }
        warning(
            "Column '", ko_id_col, "' not found. Attempting auto-detection.",
            call. = FALSE
        )
    }

    # Step 2: known name variants (case-insensitive)
    for (known in .known_ko_col_names) {
        idx <- which(lower_names == known)
        if (length(idx) > 0L) {
            candidate <- col_names[[idx[[1L]]]]
            if (.is_ko_values(df[[candidate]])) {
                message("Auto-detected KO column by name: '", candidate, "'")
                return(candidate)
            }
        }
    }

    # Step 3: full value scan (stricter threshold to avoid false positives)
    for (col in col_names) {
        if (.is_ko_values(df[[col]], min_fraction = 0.5)) {
            message("Auto-detected KO column by value pattern: '", col, "'")
            return(col)
        }
    }

    NULL
}

# ── Link table construction ────────────────────────────────────────────────────

.build_links_from_annotations <- function(gene_df, gene_ids,
                                           genome_id_col, ko_id_col) {
    links <- list()

    # gene → genome
    if (!is.null(genome_id_col) &&
        genome_id_col %in% names(gene_df) &&
        !all(is.na(gene_df[[genome_id_col]]))) {
        links$gene_to_genome <- S4Vectors::DataFrame(
            gene_id   = gene_ids,
            genome_id = as.character(gene_df[[genome_id_col]])
        )
    }

    # gene → KO (flexible detection + value validation)
    detected_ko_col <- .detect_ko_column(gene_df, ko_id_col)
    if (!is.null(detected_ko_col)) {
        ko_vals <- as.character(gene_df[[detected_ko_col]])
        has_ko  <- !is.na(ko_vals) & nzchar(ko_vals) & grepl(.ko_pattern, ko_vals)
        n_invalid <- sum(!is.na(ko_vals) & nzchar(ko_vals) & !grepl(.ko_pattern, ko_vals))
        if (n_invalid > 0L) {
            warning(
                n_invalid, " value(s) in column '", detected_ko_col,
                "' do not match K[0-9]{5} and will be excluded from the gene_to_ko link.",
                call. = FALSE
            )
        }
        if (any(has_ko)) {
            links$gene_to_ko <- S4Vectors::DataFrame(
                gene_id = gene_ids[has_ko],
                ko_id   = ko_vals[has_ko]
            )
        }
    }

    S4Vectors::SimpleList(links)
}

# ── KEGG link fetching ─────────────────────────────────────────────────────────

.batch_kegg_link <- function(target, source_ids, batch_size = 10L, sleep = 0.3) {
    if (length(source_ids) == 0L) return(stats::setNames(character(0), character(0)))
    batches <- split(source_ids, ceiling(seq_along(source_ids) / batch_size))
    results <- vector("list", length(batches))
    for (i in seq_along(batches)) {
        results[[i]] <- tryCatch(
            KEGGREST::keggLink(target, batches[[i]]),
            error = function(e) {
                warning("KEGGREST error fetching '", target, "' links (batch ", i, "): ",
                        conditionMessage(e), call. = FALSE)
                stats::setNames(character(0), character(0))
            }
        )
        if (sleep > 0 && i < length(batches)) Sys.sleep(sleep)
    }
    unlist(results, use.names = TRUE)
}

#' Fetch KEGG Higher-Order Link Tables for a Set of KO Identifiers
#'
#' `buildKEGGLinks()` queries the KEGG REST API via `KEGGREST` to retrieve
#' `ko_to_module`, `ko_to_pathway`, and `module_to_pathway` mapping tables for
#' a supplied set of KO identifiers. The result is a named `SimpleList` that
#' can be passed directly to the `links` argument of [MTTKExperiment()], or
#' used alongside [readMTTKExperiment()] and [fromPhyloseq()].
#'
#' Only KEGG reference pathways (`map` prefix) are returned; organism-specific
#' pathway entries are excluded. Module-to-pathway links are derived from the
#' modules found for the supplied KOs, so only modules represented in your
#' dataset are included.
#'
#' Access to the KEGG REST API requires network connectivity and is subject to
#' KEGG's academic-use terms. The function batches requests (10 IDs per call)
#' with an optional inter-batch pause controlled by `sleep`.
#'
#' @param ko_ids Character vector of KO identifiers in the form `K[0-9]{5}`
#'   (e.g. `"K00001"`). Values that do not match this pattern are silently
#'   dropped.
#' @param sleep Seconds to wait between successive KEGG REST requests.
#'   Defaults to `0.3`. Set to `0` to disable.
#'
#' @return A `SimpleList` with three `DataFrame` entries: `ko_to_module`,
#'   `ko_to_pathway`, and `module_to_pathway`.
#'
#' @seealso [readMTTKExperiment()] which calls this automatically when
#'   `fetchKEGGLinks = TRUE`.
#'
#' @examples
#' \dontrun{
#' links <- buildKEGGLinks(c("K00001", "K00002", "K01006"))
#' names(links)
#' }
#'
#' @export
buildKEGGLinks <- function(ko_ids, sleep = 0.3) {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop(
            "Package 'KEGGREST' is required for buildKEGGLinks(). ",
            "Install it with: BiocManager::install('KEGGREST')",
            call. = FALSE
        )
    }
    ko_ids <- unique(as.character(ko_ids))
    ko_ids <- ko_ids[!is.na(ko_ids) & nzchar(ko_ids) & grepl(.ko_pattern, ko_ids)]
    if (length(ko_ids) == 0L) {
        stop("No valid KO identifiers (K[0-9]{5}) supplied.", call. = FALSE)
    }
    if (!is.numeric(sleep) || length(sleep) != 1L || is.na(sleep) || sleep < 0) {
        stop("'sleep' must be a single non-negative number.", call. = FALSE)
    }

    # KO → module
    message("Fetching KO-to-module links for ", length(ko_ids), " KOs ...")
    raw_ko_mod  <- .batch_kegg_link("module", paste0("ko:", ko_ids), sleep = sleep)
    ko_to_module <- if (length(raw_ko_mod) > 0L) {
        S4Vectors::DataFrame(
            ko_id     = sub("^ko:", "", names(raw_ko_mod)),
            module_id = sub("^md:", "", raw_ko_mod)
        )
    } else {
        message("  No KO-to-module links found.")
        S4Vectors::DataFrame(ko_id = character(0L), module_id = character(0L))
    }

    # KO → pathway (map/reference pathways only)
    message("Fetching KO-to-pathway links ...")
    raw_ko_path     <- .batch_kegg_link("pathway", paste0("ko:", ko_ids), sleep = sleep)
    raw_ko_path_map <- raw_ko_path[grepl("^path:map", raw_ko_path)]
    ko_to_pathway <- if (length(raw_ko_path_map) > 0L) {
        S4Vectors::DataFrame(
            ko_id      = sub("^ko:", "", names(raw_ko_path_map)),
            pathway_id = sub("^path:", "", raw_ko_path_map)
        )
    } else {
        message("  No KO-to-pathway links found.")
        S4Vectors::DataFrame(ko_id = character(0L), pathway_id = character(0L))
    }

    # module → pathway (derived from modules present in the dataset)
    module_ids <- unique(as.character(ko_to_module$module_id))
    module_to_pathway <- if (length(module_ids) > 0L) {
        message("Fetching module-to-pathway links for ", length(module_ids), " modules ...")
        raw_mod_path     <- .batch_kegg_link("pathway", paste0("md:", module_ids),
                                              sleep = sleep)
        raw_mod_path_map <- raw_mod_path[grepl("^path:map", raw_mod_path)]
        if (length(raw_mod_path_map) > 0L) {
            S4Vectors::DataFrame(
                module_id  = sub("^md:", "", names(raw_mod_path_map)),
                pathway_id = sub("^path:", "", raw_mod_path_map)
            )
        } else {
            S4Vectors::DataFrame(module_id = character(0L), pathway_id = character(0L))
        }
    } else {
        S4Vectors::DataFrame(module_id = character(0L), pathway_id = character(0L))
    }

    message("Done.")
    S4Vectors::SimpleList(
        ko_to_module      = ko_to_module,
        ko_to_pathway     = ko_to_pathway,
        module_to_pathway = module_to_pathway
    )
}

# ── readMTTKExperiment() ───────────────────────────────────────────────────────

#' Build an \code{MTTKExperiment} from Delimited Files
#'
#' `readMTTKExperiment()` constructs an `MTTKExperiment` from the flat files
#' produced by a genome-resolved metatranscriptomics pipeline: a gene-level RNA
#' count matrix, a gene annotation table, a sample metadata table, and
#' optionally genome-level DNA counts, genome metadata, and a genome phylogeny.
#'
#' **File format rules**
#'
#' - Count matrices (`geneCounts`, `genomeCounts`): rows are features, columns
#'   are samples. The **first column** must contain feature identifiers.
#' - Annotation/metadata tables (`geneAnnotations`, `sampleMetadata`,
#'   `genomeMetadata`): the **first column** must contain row identifiers.
#' - Tree file: plain Newick text (`.nwk`, `.tree`, or `.txt`).
#'
#' The separator is inferred from the file extension (`.csv` → comma, anything
#' else → tab). Pass `sep = ","` or `sep = "\t"` to override.
#'
#' **Link table construction**
#'
#' `gene_to_genome` and `gene_to_ko` are built automatically from
#' `geneAnnotations`. KO column detection is flexible: the function first
#' looks for a column named by `koIdCol`, then tries common name variants
#' (`ko`, `KO`, `KEGG`, `kegg_id`, etc.) case-insensitively, and finally
#' scans all columns for values matching the pattern `K[0-9]{5}`. Only values
#' that match the pattern are included in the link table.
#'
#' When `fetchKEGGLinks = TRUE` and KO identifiers are present, the function
#' automatically fetches `ko_to_module`, `ko_to_pathway`, and
#' `module_to_pathway` from the KEGG REST API via `KEGGREST`. This requires
#' network access and the `KEGGREST` package.
#'
#' @param geneCounts Required. File path to a gene-level RNA count matrix, or
#'   an already-loaded `matrix` or `data.frame`. Rows are genes, columns are
#'   samples. The first column must contain gene identifiers.
#' @param geneAnnotations Optional. File path or `data.frame` with one row per
#'   gene. The first column must match the row names of `geneCounts`. Columns
#'   such as `genome_id`, `ko_id`, and `gene_name` are recognised automatically.
#' @param sampleMetadata Optional. File path or `data.frame`. The first column
#'   must contain sample identifiers matching the column names of `geneCounts`.
#' @param genomeCounts Optional. File path or matrix with genome-level counts.
#'   Rows are genomes, first column contains genome identifiers.
#' @param genomeMetadata Optional. File path or `data.frame` with genome
#'   metadata. First column must contain genome identifiers.
#' @param genomeTree Optional. File path (Newick) or `ape::phylo` object.
#' @param fetchKEGGLinks Logical; if `TRUE` and KO identifiers are detected,
#'   automatically fetch `ko_to_module`, `ko_to_pathway`, and
#'   `module_to_pathway` from the KEGG REST API. Requires `KEGGREST`.
#'   Defaults to `FALSE`.
#' @param geneAssayName Name for the gene-level count assay. Default
#'   `"rna_gene_counts"`.
#' @param genomeAssayName Name for the genome-level count assay. Default
#'   `"dna_genome_counts"`.
#' @param genomeIdCol Column name in `geneAnnotations` mapping genes to genomes.
#'   Defaults to `"genome_id"`.
#' @param koIdCol Hint for the KO annotation column name in `geneAnnotations`.
#'   The function accepts many common variants and also searches by value
#'   pattern when the name is not found. Set to `NULL` to skip KO detection.
#'   Defaults to `"ko_id"`.
#' @param sleep Seconds between KEGG REST requests when `fetchKEGGLinks = TRUE`.
#'   Defaults to `0.3`.
#' @param sep Column separator. `NULL` auto-detects from file extension.
#' @param metadata Named list of experiment-level metadata.
#'
#' @return A valid `MTTKExperiment`.
#'
#' @seealso [fromPhyloseq()], [buildKEGGLinks()], [MTTKExperiment()].
#'
#' @examples
#' \dontrun{
#' x <- readMTTKExperiment(
#'     geneCounts      = "gene_counts.tsv",
#'     geneAnnotations = "gene_annotations.tsv",
#'     sampleMetadata  = "sample_metadata.tsv",
#'     genomeCounts    = "genome_counts.tsv",
#'     genomeMetadata  = "genome_metadata.tsv",
#'     genomeTree      = "genome_tree.nwk",
#'     fetchKEGGLinks  = TRUE
#' )
#' }
#'
#' @export
readMTTKExperiment <- function(
    geneCounts,
    geneAnnotations = NULL,
    sampleMetadata  = NULL,
    genomeCounts    = NULL,
    genomeMetadata  = NULL,
    genomeTree      = NULL,
    fetchKEGGLinks  = FALSE,
    geneAssayName   = "rna_gene_counts",
    genomeAssayName = "dna_genome_counts",
    genomeIdCol     = "genome_id",
    koIdCol         = "ko_id",
    sleep           = 0.3,
    sep             = NULL,
    metadata        = list()
) {
    # ── gene counts ────────────────────────────────────────────────────────────
    gene_mat   <- .read_matrix_input(geneCounts, sep = sep)
    gene_ids   <- rownames(gene_mat)
    sample_ids <- colnames(gene_mat)

    if (is.null(gene_ids) || any(is.na(gene_ids)) || any(!nzchar(gene_ids))) {
        stop(
            "Gene count matrix must have non-empty row names. ",
            "Ensure the first column of the file contains gene identifiers.",
            call. = FALSE
        )
    }

    # ── gene annotations → rowData + gene-level links ──────────────────────────
    row_data  <- NULL
    link_list <- list()

    if (!is.null(geneAnnotations)) {
        anno_df <- .read_table_input(geneAnnotations, sep = sep)
        common  <- intersect(gene_ids, rownames(anno_df))

        if (length(common) == 0L) {
            stop(
                "No gene IDs in 'geneAnnotations' match the row names of 'geneCounts'. ",
                "The first column of the annotation file must contain gene identifiers.",
                call. = FALSE
            )
        }
        if (length(common) < length(gene_ids)) {
            warning(
                length(gene_ids) - length(common),
                " gene(s) in 'geneCounts' have no annotation and will have NA metadata.",
                call. = FALSE
            )
        }

        anno_aligned           <- anno_df[gene_ids, , drop = FALSE]
        rownames(anno_aligned) <- gene_ids
        row_data               <- S4Vectors::DataFrame(anno_aligned, check.names = FALSE)

        auto_links <- .build_links_from_annotations(
            gene_df       = anno_aligned,
            gene_ids      = gene_ids,
            genome_id_col = genomeIdCol,
            ko_id_col     = koIdCol
        )
        for (nm in names(auto_links)) link_list[[nm]] <- auto_links[[nm]]
    }

    # ── sample metadata → colData ──────────────────────────────────────────────
    col_data <- NULL

    if (!is.null(sampleMetadata)) {
        meta_df    <- .read_table_input(sampleMetadata, sep = sep)
        gene_mat   <- .align_samples(gene_mat, meta_df, "'sampleMetadata'")
        sample_ids <- colnames(gene_mat)
        col_data   <- S4Vectors::DataFrame(
            meta_df[sample_ids, , drop = FALSE], check.names = FALSE
        )
    }

    # ── KEGG higher-order links ────────────────────────────────────────────────
    if (isTRUE(fetchKEGGLinks)) {
        ko_ids <- if (!is.null(link_list$gene_to_ko)) {
            unique(as.character(link_list$gene_to_ko$ko_id))
        } else {
            character(0L)
        }
        if (length(ko_ids) == 0L) {
            warning(
                "fetchKEGGLinks = TRUE but no KO identifiers were detected. ",
                "Skipping KEGG link retrieval.",
                call. = FALSE
            )
        } else {
            kegg_links <- buildKEGGLinks(ko_ids, sleep = sleep)
            for (nm in names(kegg_links)) link_list[[nm]] <- kegg_links[[nm]]
        }
    }

    # ── genome-level data ──────────────────────────────────────────────────────
    genome_assays_list <- NULL
    genome_data_df     <- NULL

    if (!is.null(genomeCounts)) {
        genome_mat     <- .read_matrix_input(genomeCounts, sep = sep)
        common_samples <- intersect(colnames(genome_mat), sample_ids)
        if (length(common_samples) == 0L) {
            stop(
                "No shared sample IDs between 'genomeCounts' and 'geneCounts'.",
                call. = FALSE
            )
        }
        genome_mat         <- genome_mat[, common_samples, drop = FALSE]
        genome_assays_list <- stats::setNames(list(genome_mat), genomeAssayName)
    }

    if (!is.null(genomeMetadata)) {
        gm_df          <- .read_table_input(genomeMetadata, sep = sep)
        genome_data_df <- S4Vectors::DataFrame(gm_df, check.names = FALSE)
    }

    tree      <- .read_tree_input(genomeTree)
    links_out <- if (length(link_list) > 0L) S4Vectors::SimpleList(link_list) else NULL

    MTTKExperiment(
        assays       = stats::setNames(list(gene_mat), geneAssayName),
        rowData      = row_data,
        colData      = col_data,
        genomeAssays = genome_assays_list,
        genomeData   = genome_data_df,
        genomeTree   = tree,
        links        = links_out,
        metadata     = metadata
    )
}

# ── fromPhyloseq() ─────────────────────────────────────────────────────────────

#' Build an \code{MTTKExperiment} from phyloseq Objects
#'
#' `fromPhyloseq()` converts one or two `phyloseq` objects into an
#' `MTTKExperiment`. The primary object `x` is treated as the **gene-level**
#' layer (RNA counts). An optional `genomePhyloseq` provides the
#' **genome-level** layer (DNA counts, genome metadata, and genome phylogeny).
#'
#' **Mapping from phyloseq slots**
#'
#' For `x` (gene level):
#' - `otu_table` → gene RNA count assay.
#' - `sam_data` → sample metadata (`colData`).
#' - `tax_table` → gene metadata (`rowData`). KO column detection is flexible:
#'   the function searches by common name variants and by value pattern
#'   (`K[0-9]{5}`) when the column name is ambiguous.
#'
#' For `genomePhyloseq` (genome level):
#' - `otu_table` → genome count assay.
#' - `tax_table` → genome metadata.
#' - `phy_tree` → genome phylogeny.
#'
#' @param x A `phyloseq` object representing gene-level RNA data.
#' @param genomePhyloseq Optional `phyloseq` object with genome-level DNA
#'   counts. Sample names must overlap with those of `x`.
#' @param fetchKEGGLinks Logical; if `TRUE` and KO identifiers are detected,
#'   automatically fetch higher-order KEGG links via `buildKEGGLinks()`.
#'   Defaults to `FALSE`.
#' @param geneAssayName Name for the gene-level count assay. Default
#'   `"rna_gene_counts"`.
#' @param genomeAssayName Name for the genome-level count assay. Default
#'   `"dna_genome_counts"`.
#' @param genomeIdCol Column name in `tax_table(x)` mapping genes to genomes.
#'   Defaults to `"genome_id"`.
#' @param koIdCol Hint for the KO annotation column. Flexible detection is
#'   applied when the exact name is not found. `NULL` to skip. Default
#'   `"ko_id"`.
#' @param sleep Seconds between KEGG REST requests when `fetchKEGGLinks = TRUE`.
#'   Defaults to `0.3`.
#' @param metadata Named list of experiment-level metadata.
#'
#' @return A valid `MTTKExperiment`.
#'
#' @seealso [readMTTKExperiment()], [buildKEGGLinks()], [MTTKExperiment()].
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' x <- fromPhyloseq(gene_physeq, genomePhyloseq = genome_physeq,
#'                   fetchKEGGLinks = TRUE)
#' }
#'
#' @export
fromPhyloseq <- function(
    x,
    genomePhyloseq  = NULL,
    fetchKEGGLinks  = FALSE,
    geneAssayName   = "rna_gene_counts",
    genomeAssayName = "dna_genome_counts",
    genomeIdCol     = "genome_id",
    koIdCol         = "ko_id",
    sleep           = 0.3,
    metadata        = list()
) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
        stop(
            "The 'phyloseq' package must be installed to use fromPhyloseq(). ",
            "Install it with: BiocManager::install('phyloseq')",
            call. = FALSE
        )
    }
    if (!methods::is(x, "phyloseq")) stop("'x' must be a phyloseq object.", call. = FALSE)

    # ── gene count matrix ──────────────────────────────────────────────────────
    otu      <- phyloseq::otu_table(x)
    gene_mat <- if (phyloseq::taxa_are_rows(otu)) as.matrix(otu) else t(as.matrix(otu))
    storage.mode(gene_mat) <- "integer"
    gene_ids   <- rownames(gene_mat)
    sample_ids <- colnames(gene_mat)

    # ── sample metadata ────────────────────────────────────────────────────────
    col_data <- NULL
    sd       <- tryCatch(phyloseq::sample_data(x), error = function(e) NULL)
    if (!is.null(sd)) {
        sd_df  <- as.data.frame(sd)
        common <- intersect(sample_ids, rownames(sd_df))
        if (length(common) > 0L) {
            col_data <- S4Vectors::DataFrame(
                sd_df[sample_ids[sample_ids %in% common], , drop = FALSE],
                check.names = FALSE
            )
        }
    }

    # ── gene metadata → rowData + links ───────────────────────────────────────
    row_data  <- NULL
    link_list <- list()

    tt <- tryCatch(phyloseq::tax_table(x), error = function(e) NULL)
    if (!is.null(tt)) {
        tt_df <- as.data.frame(as.matrix(tt), stringsAsFactors = FALSE)
        tt_df <- tt_df[gene_ids[gene_ids %in% rownames(tt_df)], , drop = FALSE]
        if (nrow(tt_df) > 0L) {
            aligned           <- tt_df[gene_ids, , drop = FALSE]
            rownames(aligned) <- gene_ids
            row_data          <- S4Vectors::DataFrame(aligned, check.names = FALSE)
            auto_links <- .build_links_from_annotations(
                gene_df       = aligned,
                gene_ids      = gene_ids,
                genome_id_col = genomeIdCol,
                ko_id_col     = koIdCol
            )
            for (nm in names(auto_links)) link_list[[nm]] <- auto_links[[nm]]
        }
    }

    # ── KEGG higher-order links ────────────────────────────────────────────────
    if (isTRUE(fetchKEGGLinks)) {
        ko_ids <- if (!is.null(link_list$gene_to_ko)) {
            unique(as.character(link_list$gene_to_ko$ko_id))
        } else {
            character(0L)
        }
        if (length(ko_ids) == 0L) {
            warning("fetchKEGGLinks = TRUE but no KO identifiers detected. Skipping.",
                    call. = FALSE)
        } else {
            kegg_links <- buildKEGGLinks(ko_ids, sleep = sleep)
            for (nm in names(kegg_links)) link_list[[nm]] <- kegg_links[[nm]]
        }
    }

    # ── genome-level phyloseq ──────────────────────────────────────────────────
    genome_assays_list <- NULL
    genome_data_df     <- NULL
    genome_tree        <- tryCatch(phyloseq::phy_tree(x), error = function(e) NULL)

    if (!is.null(genomePhyloseq)) {
        if (!methods::is(genomePhyloseq, "phyloseq")) {
            stop("'genomePhyloseq' must be a phyloseq object.", call. = FALSE)
        }
        g_otu      <- phyloseq::otu_table(genomePhyloseq)
        genome_mat <- if (phyloseq::taxa_are_rows(g_otu)) as.matrix(g_otu) else t(as.matrix(g_otu))
        storage.mode(genome_mat) <- "integer"

        common_samples <- intersect(colnames(genome_mat), sample_ids)
        if (length(common_samples) == 0L) {
            stop("No shared sample IDs between 'x' and 'genomePhyloseq'.", call. = FALSE)
        }
        genome_mat         <- genome_mat[, common_samples, drop = FALSE]
        genome_assays_list <- stats::setNames(list(genome_mat), genomeAssayName)

        g_tt <- tryCatch(phyloseq::tax_table(genomePhyloseq), error = function(e) NULL)
        if (!is.null(g_tt)) {
            genome_data_df <- S4Vectors::DataFrame(
                as.data.frame(as.matrix(g_tt), stringsAsFactors = FALSE),
                check.names = FALSE
            )
        }
        g_tree <- tryCatch(phyloseq::phy_tree(genomePhyloseq), error = function(e) NULL)
        if (!is.null(g_tree)) genome_tree <- g_tree
    }

    links_out <- if (length(link_list) > 0L) S4Vectors::SimpleList(link_list) else NULL

    MTTKExperiment(
        assays       = stats::setNames(list(gene_mat), geneAssayName),
        rowData      = row_data,
        colData      = col_data,
        genomeAssays = genome_assays_list,
        genomeData   = genome_data_df,
        genomeTree   = genome_tree,
        links        = links_out,
        metadata     = metadata
    )
}

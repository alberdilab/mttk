# ── Internal helpers ──────────────────────────────────────────────────────────

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
    if (is.data.frame(x) || methods::is(x, "DataFrame")) {
        return(as.data.frame(x))
    }
    if (is.matrix(x)) {
        return(as.data.frame(x))
    }
    stop("Input must be a file path, data.frame, DataFrame, or matrix.", call. = FALSE)
}

.read_matrix_input <- function(x, sep = NULL) {
    df <- .read_table_input(x, sep = sep, row_names = 1L)
    m  <- as.matrix(df)
    storage.mode(m) <- "integer"
    m
}

.build_links_from_annotations <- function(gene_df, gene_ids,
                                           genome_id_col, ko_id_col) {
    links <- list()

    if (!is.null(genome_id_col) &&
        genome_id_col %in% names(gene_df) &&
        !all(is.na(gene_df[[genome_id_col]]))) {
        links$gene_to_genome <- S4Vectors::DataFrame(
            gene_id   = gene_ids,
            genome_id = as.character(gene_df[[genome_id_col]])
        )
    }

    if (!is.null(ko_id_col) &&
        ko_id_col %in% names(gene_df) &&
        !all(is.na(gene_df[[ko_id_col]]))) {
        ko_vals    <- as.character(gene_df[[ko_id_col]])
        has_ko     <- !is.na(ko_vals) & ko_vals != ""
        links$gene_to_ko <- S4Vectors::DataFrame(
            gene_id = gene_ids[has_ko],
            ko_id   = ko_vals[has_ko]
        )
    }

    S4Vectors::SimpleList(links)
}

.read_tree_input <- function(x) {
    if (is.null(x)) return(NULL)
    if (methods::is(x, "phylo")) return(x)
    if (is.character(x) && length(x) == 1L) {
        if (!file.exists(x)) stop("Tree file not found: ", x, call. = FALSE)
        return(ape::read.tree(x))
    }
    stop("'genomeTree' must be a file path (Newick) or an ape::phylo object.",
         call. = FALSE)
}

.align_samples <- function(count_mat, meta_df, label) {
    common <- intersect(colnames(count_mat), rownames(meta_df))
    if (length(common) == 0L) {
        stop(
            "No shared sample IDs between the count matrix and ", label, ". ",
            "Check that column names of the count matrix match row names of ",
            "the metadata table.",
            call. = FALSE
        )
    }
    if (length(common) < ncol(count_mat)) {
        warning(
            ncol(count_mat) - length(common), " sample(s) in the count matrix ",
            "are absent from ", label, " and will be dropped.",
            call. = FALSE
        )
    }
    count_mat[, common, drop = FALSE]
}

# ── readMTTKExperiment() ──────────────────────────────────────────────────────

#' Build an \code{MTTKExperiment} from Delimited Files
#'
#' `readMTTKExperiment()` constructs an `MTTKExperiment` from the flat files
#' that are typically produced by a genome-resolved metatranscriptomics
#' pipeline: a gene-level RNA count matrix, a gene annotation table, a sample
#' metadata table, and optionally genome-level DNA counts, genome metadata, a
#' genome phylogeny, and functional link tables.
#'
#' **File format rules**
#'
#' - Count matrices (`geneCounts`, `genomeCounts`): rows are features, columns
#'   are samples. The first column must contain the feature identifiers; all
#'   remaining columns are treated as samples and must be numeric.
#' - Annotation and metadata tables (`geneAnnotations`, `sampleMetadata`,
#'   `genomeMetadata`): the first column must contain the row identifiers.
#'   Additional columns become metadata.
#' - Link tables (`koToModule`, `koToPathway`, `moduleToPathway`): the first
#'   column is the source identifier and the second column is the target
#'   identifier.
#' - Tree file: plain Newick text (`.nwk`, `.tree`, or `.txt`).
#'
#' The separator is inferred from the file extension (`.csv` → comma, anything
#' else → tab). Pass `sep = ","` or `sep = "\t"` to override.
#'
#' **Link table construction**
#'
#' When `geneAnnotations` is supplied and contains a column matching
#' `genomeIdCol`, a `gene_to_genome` link table is built automatically.
#' If a column matching `koIdCol` is also present, a `gene_to_ko` table is
#' built. Higher-order links (`koToModule`, `koToPathway`, `moduleToPathway`)
#' can be provided as separate files.
#'
#' @param geneCounts Required. File path to a gene-level RNA count matrix, or
#'   an already-loaded `matrix` or `data.frame`. Rows are genes, columns are
#'   samples. The first column must contain gene identifiers.
#' @param geneAnnotations Optional. File path or `data.frame` with one row per
#'   gene. The first column must contain gene identifiers matching those in
#'   `geneCounts`. Additional columns (e.g. `genome_id`, `ko_id`, `gene_name`)
#'   become gene-level metadata and are used to build link tables automatically.
#' @param sampleMetadata Optional. File path or `data.frame` with one row per
#'   sample. The first column must contain sample identifiers matching the
#'   column names of `geneCounts`.
#' @param genomeCounts Optional. File path or matrix with genome-level DNA (or
#'   RNA) counts. Rows are genomes, columns are samples. The first column must
#'   contain genome identifiers.
#' @param genomeMetadata Optional. File path or `data.frame` with one row per
#'   genome. The first column must contain genome identifiers.
#' @param genomeTree Optional. File path to a Newick tree or an `ape::phylo`
#'   object. Tip labels must match the genome identifiers in `genomeCounts` or
#'   `geneAnnotations`.
#' @param koToModule Optional. File path or `data.frame` with two columns:
#'   `ko_id` and `module_id`.
#' @param koToPathway Optional. File path or `data.frame` with two columns:
#'   `ko_id` and `pathway_id`.
#' @param moduleToPathway Optional. File path or `data.frame` with two columns:
#'   `module_id` and `pathway_id`.
#' @param geneAssayName Name to assign to the gene-level count assay.
#'   Defaults to `"rna_gene_counts"`.
#' @param genomeAssayName Name to assign to the genome-level count assay.
#'   Defaults to `"dna_genome_counts"`.
#' @param genomeIdCol Column name in `geneAnnotations` that maps each gene to
#'   its parent genome. Used to build the `gene_to_genome` link table and
#'   `rowData(x)$genome_id`. Defaults to `"genome_id"`.
#' @param koIdCol Column name in `geneAnnotations` that gives the KO annotation
#'   of each gene. Used to build the `gene_to_ko` link table. Defaults to
#'   `"ko_id"`. Set to `NULL` to skip KO link construction.
#' @param sep Column separator for file reading. `NULL` (default) auto-detects
#'   from file extension: `.csv` → `","`, everything else → `"\t"`.
#' @param metadata Named list of experiment-level metadata stored in
#'   `S4Vectors::metadata(x)`.
#'
#' @return A valid `MTTKExperiment`.
#'
#' @seealso [fromPhyloseq()] for importing from phyloseq objects,
#'   [MTTKExperiment()] for manual construction.
#'
#' @examples
#' \dontrun{
#' x <- readMTTKExperiment(
#'     geneCounts      = "gene_counts.tsv",
#'     geneAnnotations = "gene_annotations.tsv",
#'     sampleMetadata  = "sample_metadata.tsv",
#'     genomeCounts    = "genome_counts.tsv",
#'     genomeMetadata  = "genome_metadata.tsv",
#'     genomeTree      = "genome_tree.nwk"
#' )
#' }
#'
#' @export
readMTTKExperiment <- function(
    geneCounts,
    geneAnnotations  = NULL,
    sampleMetadata   = NULL,
    genomeCounts     = NULL,
    genomeMetadata   = NULL,
    genomeTree       = NULL,
    koToModule       = NULL,
    koToPathway      = NULL,
    moduleToPathway  = NULL,
    geneAssayName    = "rna_gene_counts",
    genomeAssayName  = "dna_genome_counts",
    genomeIdCol      = "genome_id",
    koIdCol          = "ko_id",
    sep              = NULL,
    metadata         = list()
) {
    # ── gene counts ────────────────────────────────────────────────────────────
    gene_mat <- .read_matrix_input(geneCounts, sep = sep)
    gene_ids   <- rownames(gene_mat)
    sample_ids <- colnames(gene_mat)

    if (is.null(gene_ids) || any(is.na(gene_ids)) || any(gene_ids == "")) {
        stop(
            "Gene count matrix must have non-empty row names. ",
            "Ensure the first column of the file contains gene identifiers.",
            call. = FALSE
        )
    }

    # ── gene annotations → rowData + links ────────────────────────────────────
    row_data <- NULL
    link_list <- list()

    if (!is.null(geneAnnotations)) {
        anno_df  <- .read_table_input(geneAnnotations, sep = sep)
        gene_ids_common <- intersect(gene_ids, rownames(anno_df))

        if (length(gene_ids_common) == 0L) {
            stop(
                "No gene IDs in 'geneAnnotations' match the row names of 'geneCounts'. ",
                "The first column of the annotation file must contain gene identifiers.",
                call. = FALSE
            )
        }
        if (length(gene_ids_common) < length(gene_ids)) {
            warning(
                length(gene_ids) - length(gene_ids_common),
                " gene(s) in 'geneCounts' have no annotation and will have NA metadata.",
                call. = FALSE
            )
        }

        anno_aligned             <- anno_df[gene_ids, , drop = FALSE]
        rownames(anno_aligned)   <- gene_ids
        row_data                 <- S4Vectors::DataFrame(anno_aligned, check.names = FALSE)

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
        meta_df      <- .read_table_input(sampleMetadata, sep = sep)
        gene_mat     <- .align_samples(gene_mat, meta_df, "'sampleMetadata'")
        sample_ids   <- colnames(gene_mat)
        meta_aligned <- meta_df[sample_ids, , drop = FALSE]
        col_data     <- S4Vectors::DataFrame(meta_aligned, check.names = FALSE)
    }

    # ── additional link tables ─────────────────────────────────────────────────
    if (!is.null(koToModule)) {
        df <- .read_table_input(koToModule, sep = sep, row_names = NULL)
        names(df)[seq_len(min(2L, ncol(df)))] <- c("ko_id", "module_id")[seq_len(min(2L, ncol(df)))]
        link_list$ko_to_module <- S4Vectors::DataFrame(df)
    }
    if (!is.null(koToPathway)) {
        df <- .read_table_input(koToPathway, sep = sep, row_names = NULL)
        names(df)[seq_len(min(2L, ncol(df)))] <- c("ko_id", "pathway_id")[seq_len(min(2L, ncol(df)))]
        link_list$ko_to_pathway <- S4Vectors::DataFrame(df)
    }
    if (!is.null(moduleToPathway)) {
        df <- .read_table_input(moduleToPathway, sep = sep, row_names = NULL)
        names(df)[seq_len(min(2L, ncol(df)))] <- c("module_id", "pathway_id")[seq_len(min(2L, ncol(df)))]
        link_list$module_to_pathway <- S4Vectors::DataFrame(df)
    }

    # ── genome-level data ──────────────────────────────────────────────────────
    genome_assays_list <- NULL
    genome_data_df     <- NULL

    if (!is.null(genomeCounts)) {
        genome_mat <- .read_matrix_input(genomeCounts, sep = sep)
        common_samples <- intersect(colnames(genome_mat), sample_ids)
        if (length(common_samples) == 0L) {
            stop(
                "No shared sample IDs between 'genomeCounts' and 'geneCounts'. ",
                "Column names must match.",
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

    tree <- .read_tree_input(genomeTree)

    links_out <- if (length(link_list) > 0L) S4Vectors::SimpleList(link_list) else NULL

    MTTKExperiment(
        assays         = stats::setNames(list(gene_mat), geneAssayName),
        rowData        = row_data,
        colData        = col_data,
        genomeAssays   = genome_assays_list,
        genomeData     = genome_data_df,
        genomeTree     = tree,
        links          = links_out,
        metadata       = metadata
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
#' For the gene-level object `x`:
#' - `otu_table` → gene RNA count assay (name set by `geneAssayName`).
#' - `sam_data` → sample metadata (`colData`).
#' - `tax_table` → gene metadata (`rowData`). Columns named `genomeIdCol` and
#'   `koIdCol` are used to build `gene_to_genome` and `gene_to_ko` link tables
#'   automatically.
#'
#' For the genome-level object `genomePhyloseq`:
#' - `otu_table` → genome count assay (name set by `genomeAssayName`).
#' - `tax_table` → genome metadata (`genomeData`).
#' - `phy_tree` → genome phylogeny (`genomeTree`).
#' - `sam_data` → used for alignment only; sample names must overlap with `x`.
#'
#' If `genomePhyloseq` is omitted but `x` contains a `phy_tree`, that tree is
#' stored as the genome tree when `genomeIdCol` is set (assuming tips are genome
#' identifiers). If `genomePhyloseq` is omitted and there is no `phy_tree`,
#' only gene-level data is stored.
#'
#' @param x A `phyloseq` object representing gene-level RNA data.
#' @param genomePhyloseq Optional `phyloseq` object representing genome-level
#'   DNA abundance data. Sample names must overlap with those of `x`.
#' @param geneAssayName Name to assign to the gene-level count assay.
#'   Defaults to `"rna_gene_counts"`.
#' @param genomeAssayName Name to assign to the genome-level count assay.
#'   Defaults to `"dna_genome_counts"`.
#' @param genomeIdCol Column name in `tax_table(x)` that maps each gene to its
#'   parent genome. Used to build the `gene_to_genome` link and
#'   `rowData$genome_id`. `NULL` to skip. Defaults to `"genome_id"`.
#' @param koIdCol Column name in `tax_table(x)` that gives the KO annotation.
#'   Used to build the `gene_to_ko` link. `NULL` to skip. Defaults to `"ko_id"`.
#' @param metadata Named list of experiment-level metadata.
#'
#' @return A valid `MTTKExperiment`.
#'
#' @seealso [readMTTKExperiment()] for file-based import,
#'   [MTTKExperiment()] for manual construction.
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' # gene-level phyloseq with genome_id and ko_id in tax_table
#' x <- fromPhyloseq(gene_physeq, genomePhyloseq = genome_physeq)
#' }
#'
#' @export
fromPhyloseq <- function(
    x,
    genomePhyloseq  = NULL,
    geneAssayName   = "rna_gene_counts",
    genomeAssayName = "dna_genome_counts",
    genomeIdCol     = "genome_id",
    koIdCol         = "ko_id",
    metadata        = list()
) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
        stop(
            "The 'phyloseq' package must be installed to use fromPhyloseq(). ",
            "Install it with: BiocManager::install('phyloseq')",
            call. = FALSE
        )
    }
    if (!methods::is(x, "phyloseq")) {
        stop("'x' must be a phyloseq object.", call. = FALSE)
    }

    # ── gene count matrix ──────────────────────────────────────────────────────
    otu <- phyloseq::otu_table(x)
    gene_mat <- if (phyloseq::taxa_are_rows(otu)) {
        as.matrix(otu)
    } else {
        t(as.matrix(otu))
    }
    storage.mode(gene_mat) <- "integer"

    gene_ids   <- rownames(gene_mat)
    sample_ids <- colnames(gene_mat)

    # ── sample metadata ────────────────────────────────────────────────────────
    col_data <- NULL
    sd <- tryCatch(phyloseq::sample_data(x), error = function(e) NULL)
    if (!is.null(sd)) {
        sd_df    <- as.data.frame(sd)
        common   <- intersect(sample_ids, rownames(sd_df))
        if (length(common) > 0L) {
            col_data <- S4Vectors::DataFrame(sd_df[sample_ids[sample_ids %in% common], , drop = FALSE],
                                             check.names = FALSE)
        }
    }

    # ── gene metadata (tax_table) → rowData + links ────────────────────────────
    row_data  <- NULL
    link_list <- list()

    tt <- tryCatch(phyloseq::tax_table(x), error = function(e) NULL)
    if (!is.null(tt)) {
        tt_df <- as.data.frame(as.matrix(tt), stringsAsFactors = FALSE)
        tt_df <- tt_df[gene_ids[gene_ids %in% rownames(tt_df)], , drop = FALSE]
        if (nrow(tt_df) > 0L) {
            aligned              <- tt_df[gene_ids, , drop = FALSE]
            rownames(aligned)    <- gene_ids
            row_data             <- S4Vectors::DataFrame(aligned, check.names = FALSE)
            auto_links <- .build_links_from_annotations(
                gene_df       = aligned,
                gene_ids      = gene_ids,
                genome_id_col = genomeIdCol,
                ko_id_col     = koIdCol
            )
            for (nm in names(auto_links)) link_list[[nm]] <- auto_links[[nm]]
        }
    }

    # ── genome phylogeny from gene-level phyloseq ──────────────────────────────
    gene_tree <- tryCatch(phyloseq::phy_tree(x), error = function(e) NULL)

    # ── genome-level phyloseq ──────────────────────────────────────────────────
    genome_assays_list <- NULL
    genome_data_df     <- NULL
    genome_tree        <- gene_tree  # fall back to gene tree if no genome phyloseq

    if (!is.null(genomePhyloseq)) {
        if (!methods::is(genomePhyloseq, "phyloseq")) {
            stop("'genomePhyloseq' must be a phyloseq object.", call. = FALSE)
        }

        g_otu <- phyloseq::otu_table(genomePhyloseq)
        genome_mat <- if (phyloseq::taxa_are_rows(g_otu)) {
            as.matrix(g_otu)
        } else {
            t(as.matrix(g_otu))
        }
        storage.mode(genome_mat) <- "integer"

        common_samples <- intersect(colnames(genome_mat), sample_ids)
        if (length(common_samples) == 0L) {
            stop(
                "No shared sample IDs between 'x' and 'genomePhyloseq'. ",
                "Sample names must overlap.",
                call. = FALSE
            )
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

.mttk_example_components <- function() {
    gene_ids <- paste0("gene_", seq_len(6))
    sample_ids <- paste0("sample_", seq_len(4))
    genome_ids <- paste0("genome_", seq_len(3))

    rna_counts <- matrix(
        c(
            120L, 110L, 210L, 230L,
            35L, 30L, 70L, 80L,
            90L, 95L, 130L, 140L,
            12L, 10L, 18L, 22L,
            60L, 58L, 40L, 38L,
            18L, 20L, 12L, 10L
        ),
        nrow = length(gene_ids),
        byrow = TRUE,
        dimnames = list(gene_ids, sample_ids)
    )

    dna_counts <- matrix(
        c(
            50L, 48L, 52L, 54L,
            50L, 48L, 52L, 54L,
            30L, 32L, 31L, 29L,
            30L, 32L, 31L, 29L,
            20L, 22L, 19L, 18L,
            20L, 22L, 19L, 18L
        ),
        nrow = length(gene_ids),
        byrow = TRUE,
        dimnames = list(gene_ids, sample_ids)
    )

    row_data <- S4Vectors::DataFrame(
        gene_id = gene_ids,
        gene_name = c(
            "rpoB_like",
            "atpA_like",
            "glnA_like",
            "nirK_like",
            "acsA_like",
            "fdhF_like"
        ),
        genome_id = c(
            "genome_1",
            "genome_1",
            "genome_2",
            "genome_2",
            "genome_3",
            "genome_3"
        ),
        ko_id = c("K03043", "K02111", "K01915", "K00368", "K00198", "K00123")
    )
    rownames(row_data) <- gene_ids

    col_data <- S4Vectors::DataFrame(
        sample_id = sample_ids,
        condition = c("control", "control", "treated", "treated"),
        site = c("reef", "reef", "estuary", "estuary"),
        pH = c(8.1, 8.0, 7.4, 7.3)
    )
    rownames(col_data) <- sample_ids

    genome_data <- S4Vectors::DataFrame(
        genome_name = c("Genome A", "Genome B", "Genome C"),
        clade = c("Clade Alpha", "Clade Beta", "Clade Gamma"),
        taxonomy = c(
            "Bacteria; Proteobacteria",
            "Bacteria; Bacteroidota",
            "Archaea; Thermoproteota"
        )
    )
    rownames(genome_data) <- genome_ids

    link_tables <- S4Vectors::SimpleList(
        gene_to_genome = S4Vectors::DataFrame(
            gene_id = gene_ids,
            genome_id = row_data$genome_id
        ),
        gene_to_ko = S4Vectors::DataFrame(
            gene_id = gene_ids,
            ko_id = row_data$ko_id
        ),
        ko_to_module = S4Vectors::DataFrame(
            ko_id = c("K03043", "K02111", "K01915", "K00368", "K00198", "K00123"),
            module_id = c("M001", "M001", "M002", "M003", "M004", "M004")
        ),
        module_to_pathway = S4Vectors::DataFrame(
            module_id = c("M001", "M002", "M003", "M004"),
            pathway_id = c("map00970", "map00910", "map00920", "map00680")
        )
    )

    list(
        assays = list(
            rna_counts = rna_counts,
            dna_counts = dna_counts
        ),
        rowData = row_data,
        colData = col_data,
        genomeData = genome_data,
        links = link_tables,
        activeHierarchies = c("biological_origin", "functional_annotation"),
        metadata = list(
            study_name = "Synthetic genome-resolved metatranscriptomics example",
            study_design = "Two control and two treated community samples"
        )
    )
}

.make_mttk_example <- function() {
    components <- .mttk_example_components()

    MTTKExperiment(
        assays = components$assays,
        rowData = components$rowData,
        colData = components$colData,
        genomeData = components$genomeData,
        links = components$links,
        activeHierarchies = components$activeHierarchies,
        metadata = components$metadata
    )
}

#' Create the Packaged Example `MTTKExperiment`
#'
#' `makeExampleMTTKExperiment()` returns a small `MTTKExperiment` that can be
#' used in examples, tests, and interactive exploration of the package.
#'
#' The returned object contains:
#'
#' - six genes measured across four samples,
#' - paired `rna_counts` and `dna_counts` assays,
#' - feature, sample, and genome metadata,
#' - explicit biological and functional mapping tables.
#'
#' @return A valid `MTTKExperiment`.
#'
#' @examples
#' x <- makeExampleMTTKExperiment()
#' x
#' names(links(x))
#'
#' @export
makeExampleMTTKExperiment <- function() {
    .make_mttk_example()
}

#' Example `MTTKExperiment`
#'
#' `MTTKExample` is a packaged `MTTKExperiment` object for examples,
#' documentation, and quick interactive use.
#'
#' The object contains two assays (`rna_counts` and `dna_counts`), six genes,
#' four samples, three genomes, and a small set of biological and functional
#' link tables.
#'
#' @docType data
#' @format An `MTTKExperiment` object with 6 rows and 4 columns.
#'
#' @source Generated from `data-raw/make-example-data.R`.
#'
#' @keywords datasets
#'
#' @examples
#' data("MTTKExample")
#' MTTKExample
#' genomeData(MTTKExample)
#'
"MTTKExample"

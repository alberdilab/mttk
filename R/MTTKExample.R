.mttk_example_components <- function() {
    gene_ids <- paste0("gene_", seq_len(6))
    sample_ids <- paste0("sample_", seq_len(4))
    genome_ids <- paste0("genome_", seq_len(3))

    rna_gene_counts <- matrix(
        c(
            120L, 110L, 230L, 250L,
            40L, 35L, 42L, 38L,
            90L, 95L, 170L, 180L,
            70L, 75L, 40L, 35L,
            55L, 60L, 58L, 56L,
            30L, 28L, 15L, 12L
        ),
        nrow = length(gene_ids),
        byrow = TRUE,
        dimnames = list(gene_ids, sample_ids)
    )

    dna_genome_counts <- matrix(
        c(
            50L, 48L, 52L, 54L,
            30L, 32L, 31L, 29L,
            20L, 22L, 19L, 18L
        ),
        nrow = length(genome_ids),
        byrow = TRUE,
        dimnames = list(genome_ids, sample_ids)
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
        ko_id = c("K03043", "K02111", "K03043", "K00368", "K02111", "K00368")
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
            ko_id = c("K03043", "K02111", "K00368"),
            module_id = c("M001", "M002", "M003")
        ),
        module_to_pathway = S4Vectors::DataFrame(
            module_id = c("M001", "M002", "M003"),
            pathway_id = c("map00970", "map00910", "map00920")
        )
    )

    list(
        assays = list(
            rna_gene_counts = rna_gene_counts
        ),
        rowData = row_data,
        colData = col_data,
        genomeAssays = list(
            dna_genome_counts = dna_genome_counts
        ),
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
        genomeAssays = components$genomeAssays,
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
#' - a gene-level `rna_gene_counts` assay,
#' - a genome-level `dna_genome_counts` assay stored in `genomeExperiment(x)`,
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
#' The object contains one gene-level assay (`rna_gene_counts`), one
#' genome-level assay (`dna_genome_counts`), six genes, four samples, three
#' genomes, and a small set of biological and functional link tables.
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
#' dnaGenomeCounts(MTTKExample)
#'
"MTTKExample"

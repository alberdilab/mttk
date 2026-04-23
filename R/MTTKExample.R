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
        station_id = c("station_1", "station_2", "station_1", "station_2"),
        site = c("reef", "reef", "estuary", "estuary"),
        pH = c(8.1, 8.0, 7.4, 7.3)
    )
    rownames(col_data) <- sample_ids

    genome_data <- S4Vectors::DataFrame(
        genome_name = c("Genome A", "Genome B", "Genome C"),
        domain = c("Bacteria", "Bacteria", "Archaea"),
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
            ko_id = c(
                "K03043", "K03043",
                "K02111", "K02111",
                "K00368", "K00368"
            ),
            module_id = c(
                "M00001", "M00002",
                "M00009", "M00115",
                "M00010", "M00011"
            )
        ),
        ko_to_pathway = S4Vectors::DataFrame(
            ko_id = c(
                "K03043", "K03043",
                "K02111", "K02111",
                "K00368", "K00368"
            ),
            pathway_id = c(
                "map00010", "map00020",
                "map00020", "map00190",
                "map00910", "map00680"
            )
        ),
        module_to_pathway = S4Vectors::DataFrame(
            module_id = c(
                "M00001", "M00002",
                "M00009", "M00115",
                "M00010", "M00011"
            ),
            pathway_id = c(
                "map00010", "map00020",
                "map00020", "map00190",
                "map00910", "map00680"
            )
        )
    )

    genome_tree <- ape::read.tree(
        text = "((genome_1:0.1,genome_2:0.1)Bacteria_clade:0.2,genome_3:0.3)root;"
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
        genomeTree = genome_tree,
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
        genomeTree = components$genomeTree,
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
#' - a rooted genome phylogeny stored in `genomeTree(x)`,
#' - feature, sample, and genome metadata,
#' - explicit biological and functional mapping tables, including direct
#'   many-to-many `ko_to_module` and `ko_to_pathway` links.
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
#' genome-level assay (`dna_genome_counts`), a rooted genome phylogeny, six
#' genes, four samples, three genomes, and a small set of biological and
#' functional link tables.
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
#' genomeTree(MTTKExample)
#'
"MTTKExample"

.mttk_showcase_components <- function() {
    gene_ids <- paste0("gene_", seq_len(16))
    sample_ids <- paste0("sample_", seq_len(6))
    genome_ids <- paste0("genome_", seq_len(4))

    rna_gene_counts <- matrix(
        c(
            120L, 118L, 122L, 220L, 230L, 225L,
            80L, 78L, 82L, 140L, 145L, 142L,
            70L, 72L, 68L, 35L, 34L, 36L,
            45L, 44L, 46L, 50L, 52L, 51L,
            110L, 108L, 112L, 205L, 210L, 208L,
            75L, 73L, 77L, 130L, 128L, 132L,
            62L, 64L, 60L, 30L, 31L, 29L,
            40L, 39L, 41L, 42L, 43L, 41L,
            85L, 87L, 84L, 160L, 165L, 162L,
            65L, 66L, 64L, 40L, 42L, 41L,
            58L, 57L, 59L, 28L, 27L, 29L,
            35L, 36L, 34L, 39L, 38L, 40L,
            95L, 94L, 96L, 180L, 182L, 185L,
            68L, 70L, 67L, 44L, 45L, 43L,
            55L, 54L, 56L, 26L, 25L, 27L,
            38L, 37L, 39L, 70L, 72L, 71L
        ),
        nrow = length(gene_ids),
        byrow = TRUE,
        dimnames = list(gene_ids, sample_ids)
    )

    dna_genome_counts <- matrix(
        c(
            100L, 102L, 98L, 110L, 111L, 109L,
            90L, 88L, 92L, 84L, 82L, 83L,
            60L, 62L, 59L, 65L, 67L, 64L,
            70L, 72L, 68L, 75L, 74L, 76L
        ),
        nrow = length(genome_ids),
        byrow = TRUE,
        dimnames = list(genome_ids, sample_ids)
    )

    row_data <- S4Vectors::DataFrame(
        gene_id = gene_ids,
        gene_name = c(
            "rbcL_A", "uspA_A", "narG_A", "cyoB_A",
            "rbcL_B", "uspA_B", "narG_B", "ectB_B",
            "rbcL_C", "uspA_C", "narG_C", "cyoB_C",
            "rbcL_D", "uspA_D", "narG_D", "ectB_D"
        ),
        genome_id = rep(genome_ids, each = 4L),
        ko_id = c(
            "K10001", "K10002", "K10003", "K10004",
            "K10001", "K10002", "K10003", "K10005",
            "K10001", "K10002", "K10003", "K10004",
            "K10001", "K10002", "K10003", "K10005"
        ),
        row.names = gene_ids
    )

    col_data <- S4Vectors::DataFrame(
        sample_id = sample_ids,
        condition = factor(
            c("control", "control", "control", "oxygen_pulse", "oxygen_pulse", "oxygen_pulse"),
            levels = c("control", "oxygen_pulse")
        ),
        station_id = c("station_1", "station_2", "station_3", "station_1", "station_2", "station_3"),
        salinity = c(32.1, 31.8, 32.4, 24.3, 23.9, 24.1),
        oxygen = c(7.8, 7.6, 7.9, 3.1, 2.9, 3.2),
        row.names = sample_ids
    )

    genome_data <- S4Vectors::DataFrame(
        genome_name = c(
            "Pelagic Alpha",
            "Pelagic Beta",
            "Methano Gamma",
            "Halo Delta"
        ),
        domain = c("Bacteria", "Bacteria", "Archaea", "Archaea"),
        clade = c(
            "Surface bloom responders",
            "Particle-associated heterotrophs",
            "Methanogenic archaea",
            "Halophilic archaea"
        ),
        taxonomy = c(
            "Bacteria; Proteobacteria",
            "Bacteria; Bacteroidota",
            "Archaea; Methanobacteriota",
            "Archaea; Halobacteriota"
        ),
        row.names = genome_ids
    )

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
            ko_id = c(
                "K10001", "K10001",
                "K10002", "K10002",
                "K10003", "K10003",
                "K10004",
                "K10005"
            ),
            module_id = c(
                "M00001", "M00002",
                "M00001", "M00003",
                "M00002", "M00004",
                "M00004",
                "M00003"
            )
        ),
        ko_to_pathway = S4Vectors::DataFrame(
            ko_id = c(
                "K10001", "K10001",
                "K10002", "K10002",
                "K10003", "K10003",
                "K10004",
                "K10005", "K10005"
            ),
            pathway_id = c(
                "map00010", "map00190",
                "map00010", "map02020",
                "map00190", "map00910",
                "map00910",
                "map02020", "map00500"
            )
        ),
        module_to_pathway = S4Vectors::DataFrame(
            module_id = c("M00001", "M00002", "M00003", "M00004"),
            pathway_id = c("map00010", "map00190", "map02020", "map00910")
        )
    )

    genome_tree <- ape::read.tree(
        text = paste0(
            "((genome_1:0.1,genome_2:0.1)Bacteria_clade:0.35,",
            "(genome_3:0.1,genome_4:0.1)Archaea_clade:0.35)root;"
        )
    )

    rna_genome_counts <- base::rowsum(
        rna_gene_counts,
        group = as.character(row_data$genome_id),
        reorder = FALSE
    )
    rownames(rna_genome_counts) <- genome_ids

    list(
        assays = list(
            rna_gene_counts = rna_gene_counts
        ),
        rowData = row_data,
        colData = col_data,
        genomeAssays = list(
            rna_genome_counts = rna_genome_counts,
            dna_genome_counts = dna_genome_counts
        ),
        genomeData = genome_data,
        genomeTree = genome_tree,
        links = link_tables,
        activeHierarchies = c("biological_origin", "functional_annotation"),
        metadata = list(
            study_name = "Synthetic coastal oxygen-pulse metatranscriptomics showcase",
            study_design = paste(
                "Four genomes followed across six incubations spanning",
                "control and oxygen-pulse conditions."
            )
        )
    )
}

.make_mttk_showcase <- function() {
    components <- .mttk_showcase_components()

    MTTKExperiment(
        assays = components$assays,
        rowData = components$rowData,
        colData = components$colData,
        genomeAssays = components$genomeAssays,
        genomeData = components$genomeData,
        genomeTree = components$genomeTree,
        links = components$links,
        activeHierarchies = components$activeHierarchies,
        metadata = components$metadata
    )
}

#' Create the Packaged Showcase `MTTKExperiment`
#'
#' `makeShowcaseMTTKExperiment()` returns a richer `MTTKExperiment` designed
#' for the package vignette and for trying the main modeling workflows on one
#' compact but more realistic object.
#'
#' The returned object contains:
#'
#' - sixteen genes measured across six samples,
#' - both `rna_gene_counts` and genome-level `rna_genome_counts`,
#' - a genome-level `dna_genome_counts` assay stored in `genomeExperiment(x)`,
#' - a rooted genome phylogeny stored in `genomeTree(x)`,
#' - four genomes spanning two domains,
#' - sample metadata for a synthetic oxygen-pulse experiment, including a
#'   repeated-measures `station_id` blocking variable,
#' - direct many-to-many links from KO to module and pathway.
#'
#' @return A valid `MTTKExperiment`.
#'
#' @examples
#' x <- makeShowcaseMTTKExperiment()
#' x
#' names(genomeAssays(x))
#'
#' @export
makeShowcaseMTTKExperiment <- function() {
    .make_mttk_showcase()
}

#' Showcase `MTTKExperiment`
#'
#' `MTTKShowcase` is a packaged `MTTKExperiment` intended for vignette-style
#' walkthroughs of the main mttk workflows.
#'
#' The object contains one gene-level assay (`rna_gene_counts`), two
#' genome-level assays (`rna_genome_counts` and `dna_genome_counts`), a rooted
#' genome phylogeny, sixteen genes, six samples, four genomes, and a small set
#' of overlapping functional mappings chosen to support KO-, module-,
#' pathway-, genome-, and group-comparison examples.
#'
#' @docType data
#' @format An `MTTKExperiment` object with 16 rows and 6 columns.
#'
#' @source Generated from `data-raw/make-example-data.R`.
#'
#' @keywords datasets
#'
#' @examples
#' data("MTTKShowcase")
#' MTTKShowcase
#' genomeData(MTTKShowcase)
#' genomeTree(MTTKShowcase)
#'
"MTTKShowcase"

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(MTTK)
x <- makeExampleMTTKExperiment()

x

## -----------------------------------------------------------------------------
data("MTTKExample")
MTTKExample

## -----------------------------------------------------------------------------
genomeData(x)

## -----------------------------------------------------------------------------
names(links(x))

## -----------------------------------------------------------------------------
activeHierarchies(x)

## -----------------------------------------------------------------------------
rnaCounts(x)

## -----------------------------------------------------------------------------
genome_counts <- aggregateToGenome(x)
genome_counts

## -----------------------------------------------------------------------------
module_counts <- aggregateByLink(
    x,
    path = c("gene_to_ko", "ko_to_module"),
    assays = "rna_counts"
)

module_counts

## -----------------------------------------------------------------------------
feature_activity <- summarizeActivity(x, by = "feature")
genome_activity <- summarizeActivity(x, by = "genome")

feature_activity
genome_activity


# MTTK

Metatranscriptomics toolkit for nested gene, genome, and community analysis.

## What MTTK is for

MTTK is designed for genome-resolved metatranscriptomics, where the data are
not only a gene-by-sample count matrix.

- genes are nested within genomes,
- genome abundance may be measured separately from gene RNA counts,
- genes may also belong to functional hierarchies such as KO, module, and
  pathway.

The package therefore focuses on:

- an explicit container that stores gene-level and genome-level measurements
  together,
- transparent link tables between genes, genomes, and annotations,
- aggregation and modeling workflows that can use those nested relationships.

## How MTTK Differs From edgeR, DESeq2, and limma-voom

`edgeR`, `DESeq2`, and `limma-voom` are established tools for differential
expression on feature-by-sample expression matrices. They are strong choices
when the main task is gene-wise or transcript-wise differential expression in a
mostly flat design.

MTTK is different in two main ways.

First, MTTK is built around a different data structure. It keeps:

- gene-level RNA data,
- genome-level data such as DNA abundance,
- explicit gene-to-genome nesting,
- optional functional hierarchies such as gene-to-KO, KO-to-module, and
  KO-to-pathway links, including many-to-many mappings

in one Bioconductor-friendly object.

Second, MTTK is aiming at questions that are awkward to express in a flat
RNA-seq workflow, for example:

- comparing activity while accounting for the parent genome of each gene,
- aggregating genes to KO or higher annotation levels before modeling,
- combining gene-level RNA with genome-level abundance offsets,
- moving between gene-, genome-, and function-level analyses in one framework.

The difference is also in the modeling target, not just in the container.

- `edgeR` typically fits one negative-binomial GLM per feature in a flat
  feature-by-sample design, with empirical-Bayes shrinkage for dispersion.
- `DESeq2` typically fits one negative-binomial GLM per feature in a flat
  feature-by-sample design, with shrinkage for dispersion and fold changes.
- `limma-voom` transforms counts to log-expression with precision weights, then
  fits weighted linear models per feature.
- MTTK's first implemented workflow instead fits one `glmmTMB`
  negative-binomial mixed model per KO after aggregating genes to
  KO-within-genome counts, with `genome_id` as a random intercept and optional
  genome-abundance offsets.

So the current MTTK workflow answers a different question:

- standard RNA-seq workflows mostly ask whether each feature changes across
  samples,
- the first MTTK modeling workflow asks whether a function, represented as a
  KO, is associated with the condition across genomes while accounting for
  between-genome heterogeneity.

So the current position of MTTK is:

- not a replacement for `edgeR`, `DESeq2`, or `limma-voom`,
- complementary infrastructure for nested metatranscriptomic data,
- a place to implement hierarchy-aware analyses that those packages were not
  primarily designed around.

## Current scope

The current version already supports:

- `MTTKExperiment` for explicit gene/genome data storage,
- aggregation to genomes and annotation-linked groups,
- KO-level mixed models with `glmmTMB`, where KO counts are modeled across
  genomes with genome random effects,
- gene-level negative-binomial models with optional parent-genome abundance
  offsets,
- optional KEGG module annotation retrieval through `KEGGREST` for joining
  module names and classes onto KO-level results.

Gene-wise differential expression with `edgeR`, `DESeq2`, or `limma-voom`
could still be used alongside MTTK, with MTTK handling the container,
alignment, and hierarchical summaries around those analyses.

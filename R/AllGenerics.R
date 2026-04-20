#' @rdname geneExperiment
#' @export
methods::setGeneric("geneExperiment", function(x) {
    standardGeneric("geneExperiment")
})

#' @rdname geneExperiment
#' @export
methods::setGeneric("geneExperiment<-", function(x, value) {
    standardGeneric("geneExperiment<-")
})

#' @rdname genomeExperiment
#' @export
methods::setGeneric("genomeExperiment", function(x) {
    standardGeneric("genomeExperiment")
})

#' @rdname genomeExperiment
#' @export
methods::setGeneric("genomeExperiment<-", function(x, value) {
    standardGeneric("genomeExperiment<-")
})

#' @rdname geneAssays
#' @export
methods::setGeneric("geneAssays", function(x, withDimnames = TRUE) {
    standardGeneric("geneAssays")
})

#' @rdname geneAssays
#' @export
methods::setGeneric("geneAssays<-", function(x, value) {
    standardGeneric("geneAssays<-")
})

#' @rdname genomeAssays
#' @export
methods::setGeneric("genomeAssays", function(x, withDimnames = TRUE) {
    standardGeneric("genomeAssays")
})

#' @rdname genomeAssays
#' @export
methods::setGeneric("genomeAssays<-", function(x, value) {
    standardGeneric("genomeAssays<-")
})

#' @rdname genomeData
#' @export
methods::setGeneric("genomeData", function(x) {
    standardGeneric("genomeData")
})

#' @rdname genomeData
#' @export
methods::setGeneric("genomeData<-", function(x, value) {
    standardGeneric("genomeData<-")
})

#' @rdname links
#' @export
methods::setGeneric("links", function(x) {
    standardGeneric("links")
})

#' @rdname links
#' @export
methods::setGeneric("links<-", function(x, value) {
    standardGeneric("links<-")
})

#' @rdname activeHierarchies
#' @export
methods::setGeneric("activeHierarchies", function(x) {
    standardGeneric("activeHierarchies")
})

#' @rdname activeHierarchies
#' @export
methods::setGeneric("activeHierarchies<-", function(x, value) {
    standardGeneric("activeHierarchies<-")
})

#' @rdname rnaGeneCounts
#' @export
methods::setGeneric("rnaGeneCounts", function(x, withDimnames = TRUE) {
    standardGeneric("rnaGeneCounts")
})

#' @rdname rnaGeneCounts
#' @export
methods::setGeneric("rnaGeneCounts<-", function(x, value) {
    standardGeneric("rnaGeneCounts<-")
})

#' @rdname rnaGenomeCounts
#' @export
methods::setGeneric("rnaGenomeCounts", function(x, withDimnames = TRUE) {
    standardGeneric("rnaGenomeCounts")
})

#' @rdname rnaGenomeCounts
#' @export
methods::setGeneric("rnaGenomeCounts<-", function(x, value) {
    standardGeneric("rnaGenomeCounts<-")
})

#' @rdname dnaGenomeCounts
#' @export
methods::setGeneric("dnaGenomeCounts", function(x, withDimnames = TRUE) {
    standardGeneric("dnaGenomeCounts")
})

#' @rdname dnaGenomeCounts
#' @export
methods::setGeneric("dnaGenomeCounts<-", function(x, value) {
    standardGeneric("dnaGenomeCounts<-")
})

#' @rdname fitInfo
#' @export
methods::setGeneric("fitInfo", function(x) {
    standardGeneric("fitInfo")
})

#' @rdname modelObjects
#' @export
methods::setGeneric("modelObjects", function(x) {
    standardGeneric("modelObjects")
})

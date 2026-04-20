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

#' @rdname rnaCounts
#' @export
methods::setGeneric("rnaCounts", function(x, withDimnames = TRUE) {
    standardGeneric("rnaCounts")
})

#' @rdname rnaCounts
#' @export
methods::setGeneric("rnaCounts<-", function(x, value) {
    standardGeneric("rnaCounts<-")
})

#' @rdname dnaCounts
#' @export
methods::setGeneric("dnaCounts", function(x, withDimnames = TRUE) {
    standardGeneric("dnaCounts")
})

#' @rdname dnaCounts
#' @export
methods::setGeneric("dnaCounts<-", function(x, value) {
    standardGeneric("dnaCounts<-")
})

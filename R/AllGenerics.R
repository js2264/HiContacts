#' @rdname plotMatrix
#' @param x Can be a `HiCExperiment`, a GInteractions or an array/matrix object. 
#' @param ... Extra arguments passed to the corresponding method.

setGeneric("plotMatrix", function(x, ...) {standardGeneric("plotMatrix")})

#' @rdname plotMatrix

setGeneric("montage", function(x, ...) {standardGeneric("montage")})

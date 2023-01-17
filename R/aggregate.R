#' @rdname arithmetics
#' @export

setMethod("aggregate", signature(x = "HiCExperiment"), function(x, ...) {
    params <- list(...)
    if (!'targets' %in% names(params)) 
        stop("Please provide a `targets` argument (`GRanges` or `GInteractions`)")
    targets <- params[['targets']]
    if ('flanking_bins' %in% names(params)) {flanking_bins <- params[['flanking_bins']]} 
    else {flanking_bins <- 51}
    if ('BPPARAM' %in% names(params)) {BPPARAM <- params[['BPPARAM']]} 
    else {BPPARAM <- BiocParallel::bpparam()}
    HiCExperiment::AggrHiCExperiment(
        file = fileName(x), 
        resolution = resolution(x), 
        targets = targets,  
        flanking_bins = flanking_bins, 
        metadata = S4Vectors::metadata(x), 
        topologicalFeatures = topologicalFeatures(x), 
        pairsFile = pairsFile(x), 
        BPPARAM = BPPARAM,
        bed = metadata(x)$bed
    )
})

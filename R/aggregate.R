#' HiContacts arithmetics functionalities
#' @name arithmetics
#' @rdname arithmetics
#' @export

setMethod("aggregate", signature(x = "HiCExperiment"), function(
    x, 
    targets, 
    flankingBins = 51, 
    maxDistance = NULL,
    BPPARAM = BiocParallel::bpparam()
) {
    HiCExperiment::AggrHiCExperiment(
        file = fileName(x), 
        resolution = resolution(x), 
        targets = targets,  
        flankingBins = flankingBins, 
        metadata = S4Vectors::metadata(x), 
        topologicalFeatures = topologicalFeatures(x), 
        pairsFile = pairsFile(x), 
        maxDistance = maxDistance, 
        BPPARAM = BPPARAM,
        bed = metadata(x)$regions
    )
})

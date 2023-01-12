#' @name arithmetics
#' @aliases aggregate,HiCExperiment-method
#' @importFrom S4Vectors aggregate
#' @importFrom BiocParallel bpparam
#' @return a `AggrHiCExperiment` object
#' @export
#' @examples 
#' #### -----
#' #### Aggregate a contact matrix over centromeres, at different scales
#' #### -----
#' 
#' contacts <- full_contacts_yeast() |> zoom(resolution = 1000)
#' centros <- topologicalFeatures(contacts, 'centromeres')
#' aggr <- aggregate(contacts, targets = centros, flanking_bins = 50)
#' plotMatrix(aggr, 'detrended', scale = 'linear', limits = c(-1, 1))
#' 
#' contacts <- full_contacts_yeast() |> zoom(resolution = 8000)
#' centros <- topologicalFeatures(contacts, 'centromeres')
#' aggr <- aggregate(contacts, targets = centros, flanking_bins = 20)
#' plotMatrix(aggr, 'detrended', scale = 'linear', limits = c(-1, 1))
#' 

setMethod("aggregate", signature(x = "HiCExperiment"), function(x, ...) {
    params <- list(...)
    if (!'targets' %in% names(params)) 
        stop("Please provide a `targets` argument (`GRanges` or `GInteractions`)")
    targets <- params[['targets']]
    if ('flanking_bins' %in% names(params)) {flanking_bins <- params[['flanking_bins']]} 
    else {flanking_bins <- 50}
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

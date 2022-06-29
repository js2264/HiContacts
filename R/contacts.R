################################################################################
#                                                                              #
#                               CONSTRUCTORS                                   #
#                                                                              #
################################################################################

setClassUnion("GRangesOrcharacterOrNULL", members = c("GRanges", "character", "NULL"))

methods::setClass("contacts", 
    contains = c("Annotated"), 
    slots = c(
        focus = "GRangesOrcharacterOrNULL", 
        metadata = "list", 
        seqinfo = "Seqinfo", 
        resolutions = "numeric", 
        current_resolution = "numeric", 
        bins = "GRanges",
        interactions = "GInteractions",
        assays = "SimpleList"
    )
)

contacts <- function(cool_path, resolution, focus = NULL, metadata = NULL) {
    
    ## -- Check that provided file is valid 
    check_cool_file(cool_path)

    ## -- Read resolutions
    resolutions <- lsCoolResolutions(cool_path, full.list = FALSE, silent = TRUE)
    first_res <- resolutions[1]
    if (!missing(resolution)) {
        current_res <- resolution
    }
    else {
        current_res <- first_res
    }

    ## -- Read interactions
    gis <- cool2gi(cool_path, res = current_res, coords = focus)

    ## -- Read seqinfo
    si <- cool2seqinfo(cool_path, first_res)

    ## -- Tile the genome
    bins <- GenomicRanges::tileGenome(
        seqlengths = GenomeInfoDb::seqlengths(si), 
        tilewidth = current_res, 
        cut.last.tile.in.chrom = TRUE
    )
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- NULL

    ## -- Create contact object
    x <- methods::new("contacts", 
        focus = focus, 
        metadata = purrr::flatten(c(list(path = cool_path), metadata)), 
        seqinfo = si, 
        resolutions = resolutions, 
        current_resolution = current_res, 
        bins = bins, 
        interactions = gis, 
        assays = S4Vectors::SimpleList(
            'raw' = tibble(score = mcols$count), 
            'balanced' = tibble(score = mcols$score)
        )
    )
    validObject(x)
    return(x)
}

setValidity("contacts",
    function(object) {
        if (!is(object@focus, "GRangesOrcharacterOrNULL"))
            return("'focus' slot must be a GRanges")
        if (!is(resolutions(object), "numeric"))
            return("'resolutions' slot must be an numeric vector")
        if (!is(bins(object), "GRanges"))
            return("'bins' slot must be a GRanges")
        if (!is(assays(object), "SimpleList"))
            return("'assays' slot must be a SimpleList")
        TRUE
    }
)

################################################################################
#                                                                              #
#                                 FUNCTIONS                                    #
#                                                                              #
################################################################################

focus <- function(x, focus) {
    y <- contacts(
        path(x), 
        resolution = currentResolution(x), 
        focus = focus, 
        metadata = S4Vectors::metadata(x))
    validObject(y)
    y
}
switchResolution <- function(x, resolution) {
    check_resolution(x, resolution)
    y <- contacts(
        path(x), 
        resolution = resolution, 
        focus = x@focus, 
        metadata = S4Vectors::metadata(x))
    validObject(y)
    y
}

################################################################################
#                                                                              #
#                                 METHODS                                      #
#                                                                              #
################################################################################

# setGeneric("filter", function(x, locus) {standardGeneric("filter")})
# setMethod("filter", signature("contacts", locus = "GRanges"), function(x, locus) {
#     anchors <- InteractionSet::anchors(x@interactions)
#     left_anchor_in_locus <- seq_along(anchors[[1]]) %in% S4Vectors::queryHits(findOverlaps(anchors[[1]], locus))
#     right_anchor_in_locus <- seq_along(anchors[[2]]) %in% S4Vectors::queryHits(findOverlaps(anchors[[2]], locus))
#     subset_interactions <- x@interactions[left_anchor_in_locus & right_anchor_in_locus]
#     subset_interactions <- reduceRegions(subset_interactions)
#     subset_assays <- lapply(x@assays, function(assay) {
#         subset_assay <- assay[left_anchor_in_locus & right_anchor_in_locus,]
#     }) %>% S4Vectors::SimpleList()
#     x@interactions <- subset_interactions
#     x@assays <- subset_assays
#     return(x)
# })

################################################################################
#                                                                              #
#                                 SETTERS                                      #
#                                                                              #
################################################################################

################################################################################
#                                                                              #
#                                 GETTERS                                      #
#                                                                              #
################################################################################

# setMethod("length", signature(x = "contacts"), definition = function(x) length(x@regions))
setMethod("length", "contacts", function(x) length(regions(x)))

setGeneric("path", function(x) {standardGeneric("path")})
setMethod("path", "contacts", function(x) S4Vectors::metadata(x)$path)
setGeneric("seqinfo", function(x) {standardGeneric("seqinfo")})
setMethod("seqinfo", "contacts", function(x) x@seqinfo)
setGeneric("resolutions", function(x) {standardGeneric("resolutions")})
setMethod("resolutions", "contacts", function(x) x@resolutions)
setGeneric("currentResolution", function(x) {standardGeneric("currentResolution")})
setMethod("currentResolution", "contacts", function(x) x@current_resolution)
setGeneric("bins", function(x) {standardGeneric("bins")})
setMethod("bins", "contacts", function(x) x@bins)
setMethod("interactions", "contacts", function(x) x@interactions)
setGeneric("assays", function(x) {standardGeneric("assays")})
setMethod("assays", "contacts", function(x) x@assays)
setGeneric("assay", function(x, name) {standardGeneric("assay")})
setMethod("assay", signature(x = "contacts", name = "numeric"), function(x, name) {
    gis <- x@interactions
    gis$score <- x@assays[[name]]
    return(gis)
})
setMethod("assay", signature(x = "contacts", name = "character"), function(x, name) {
    gis <- x@interactions
    gis$score <- x@assays[[name]]
    return(gis)
})

setMethod("anchors", "contacts", function(x) anchors(assay(x, 1)))
setMethod("regions", "contacts", function(x) regions(assay(x, 1)))

setMethod("show", signature("contacts"), function(object) {

    cat(glue::glue('`contacts` object with {format(length(interactions(object)), big.mark = ",")} interactions over {format(length(object), big.mark = ",")} regions'), '\n')
    cat(glue::glue('file: {path(object)}'), '\n')
    cat(glue::glue('focus: {ifelse(is.null(object@focus), "whole genome", formatCoords(object@focus))}'), '\n')
    cat('------------\n')

    ## Metadata
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))

    ## Resolutions
    S4Vectors::coolcat("resolutions(%d): %s\n", resolutions(object))
    cat(glue::glue('current resolution({object@current_resolution}): {length(interactions(object))} interactions'), '\n')

    ## Bins
    cat(glue::glue('bins: {length(bins(object))}'), '\n')

    ## Regions
    cat(glue::glue('regions: {length(regions(object))}'), '\n')

    ## Assays
    cat(glue::glue('assays({length(assays(object))}): {paste(names(assays(object)), collapse = " ")}'), '\n')

})


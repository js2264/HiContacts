################################################################################
#                                                                              #
#                               CONSTRUCTORS                                   #
#                                                                              #
################################################################################

setClassUnion("GRangesOrGInteractions", members = c("GRanges", "GInteractions"))
setClassUnion("GRangesOrcharacterOrNULL", members = c("GRanges", "character", "NULL"))
setClassUnion("numericOrcharacter", members = c("numeric", "character"))

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
        assays = "SimpleList", 
        features = "SimpleList"
    )
)

contacts <- function(cool_path, resolution = NULL, focus = NULL, metadata = NULL, features = NULL) {
    
    ## -- Check that provided file is valid 
    check_cool_file(cool_path)

    ## -- Read resolutions
    resolutions <- lsCoolResolutions(cool_path, full.list = FALSE, silent = TRUE)
    res <- resolutions[length(resolutions)]
    if (!is.null(resolution)) {
        current_res <- resolution
    }
    else {
        current_res <- res
    }

    ## -- Read interactions
    gis <- cool2gi(cool_path, res = current_res, coords = focus)

    ## -- Read seqinfo
    si <- cool2seqinfo(cool_path, res)

    ## -- Import features
    if (is.null(features)) {
        features <- S4Vectors::SimpleList(
            'loops' = GenomicRanges::GRanges(), 
            'borders' = GenomicRanges::GRanges(), 
            'compartments' = GenomicRanges::GRanges(), 
            'viewpoints' = GenomicRanges::GRanges()
        )
    }

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
        metadata = purrr::flatten(c(list(path = cool_path), metadata[names(metadata)!='path'])), 
        seqinfo = si, 
        resolutions = resolutions, 
        current_resolution = current_res, 
        bins = bins, 
        interactions = gis, 
        assays = S4Vectors::SimpleList(
            'raw' = mcols$count, 
            'balanced' = mcols$score
        ), 
        features = features
    )
    validObject(x)
    return(x)
} 

setValidity("contacts",
    function(object) {
        if (!is(focus(object), "GRangesOrcharacterOrNULL"))
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

zoom <- function(x, focus = NULL, resolution = NULL) {
    res <- ifelse(!is.null(resolution), resolution, resolution(x))
    foc <- case_when(
        is.null(focus) & !is.null(focus(x)) ~ list(focus(x)), 
        is.null(focus) & is.null(focus(x)) ~ list(NULL), 
        !is.null(focus) ~ list(as.character(focus))
    )[[1]]
    y <- contacts(
        path(x), 
        resolution = res, 
        focus = foc, 
        metadata = S4Vectors::metadata(x)
    )
    y@features <- features(x)
    validObject(y)
    y
}
setGeneric("addFeature", function(x, feature, feature.name) {standardGeneric("addFeature")})
setMethod("addFeature", 
    signature("contacts", "GRangesOrGInteractions", "character"), 
    function(x, feature, feature.name) {
        x@features[[feature.name]] <- feature
        x
    }
)

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
setMethod("dim", "contacts", function(x) dim(gi2cm(assay(x, 1))))
setGeneric("type", function(x) {standardGeneric("type")})
setMethod("type", "contacts", function(x) ifelse(length(x)^2 != length(interactions(x)), 'sparse', 'full'))

setGeneric("path", function(x) {standardGeneric("path")})
setMethod("path", "contacts", function(x) S4Vectors::metadata(x)$path)
setGeneric("seqinfo", function(x) {standardGeneric("seqinfo")})
setMethod("seqinfo", "contacts", function(x) x@seqinfo)
setGeneric("resolutions", function(x) {standardGeneric("resolutions")})
setMethod("resolutions", "contacts", function(x) x@resolutions)
setGeneric("resolution", function(x) {standardGeneric("resolution")})
setMethod("resolution", "contacts", function(x) x@current_resolution)
setGeneric("bins", function(x) {standardGeneric("bins")})
setMethod("bins", "contacts", function(x) x@bins)
setGeneric("focus", function(x) {standardGeneric("focus")})
setMethod("focus", "contacts", function(x) x@focus)
setMethod("interactions", "contacts", function(x) x@interactions)
setGeneric("assays", function(x) {standardGeneric("assays")})
setMethod("assays", "contacts", function(x) x@assays)
setGeneric("assay", function(x, name) {standardGeneric("assay")})
setMethod("assay", signature(x = "contacts", name = "numericOrcharacter"), function(x, name) {
    gis <- x@interactions
    gis$score <- x@assays[[name]]
    return(gis)
})
setMethod("assay", signature(x = "contacts", name = "missing"), function(x, name) {
    gis <- x@interactions
    gis$score <- x@assays[[1]]
    return(gis)
})
setGeneric("features", function(x) {standardGeneric("features")})
setMethod("features", "contacts", function(x) {
    S4Vectors::SimpleList(as.list(x@features))
})
setGeneric("feature", function(x, name) {standardGeneric("feature")})
setMethod("feature", signature(x = "contacts", name = "numericOrcharacter"), function(x, name) {
    features(x)[[name]]
})
setMethod("feature", signature(x = "contacts", name = "missing"), function(x, name) {
    features(x)[[1]]
})

setMethod("anchors", "contacts", function(x) anchors(assay(x, 1)))
setMethod("regions", "contacts", function(x) regions(assay(x, 1)))

setMethod("show", signature("contacts"), function(object) {

    if (is.null(focus(object))) {
        focus_str <- "whole genome"
    } 
    else if (is(focus(object), 'GRanges') & length(focus(object)) == 1) {
        focus_str <- formatCoords(focus(object))
    } 
    else if (is(focus(object), 'GRanges') & length(focus(object)) > 1) {
        focus_str <- glue::glue('multiple GRanges({length(focus(object))})')
    }
    else if (is(focus(object), 'character')) {
        focus_str <- formatCoords(focus(object))
    }

    cat(glue::glue('{type(object)} `contacts` object with {format(length(interactions(object)), big.mark = ",")} interactions over {format(length(object), big.mark = ",")} regions'), '\n')
    cat(glue::glue('file: {path(object)}'), '\n')
    cat(glue::glue('focus: {focus_str}'), '\n')
    cat('------------\n')

    ## Metadata
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))

    ## Resolutions
    S4Vectors::coolcat("resolutions(%d): %s\n", resolutions(object))
    cat(glue::glue('current resolution({object@current_resolution}): {length(interactions(object))} interactions in memory'), '\n')

    ## Bins
    cat(glue::glue('bins: {length(bins(object))}'), '\n')

    ## Regions
    cat(glue::glue('regions: {length(regions(object))}'), '\n')

    ## Assays
    cat(glue::glue('assays({length(assays(object))}): {paste(names(assays(object)), collapse = " ")}'), '\n')

    ## Features
    cat(glue::glue('features: {paste(paste0(names(features(object)), "(", lengths(features(object)), ")"), collapse = " ")}'), '\n')

})


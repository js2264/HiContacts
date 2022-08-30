################################################################################
#                                                                              #
#                               CONSTRUCTORS                                   #
#                                                                              #
################################################################################

setClassUnion("GRangesOrGInteractions", members = c("GRanges", "GInteractions"))
setClassUnion("GRangesOrPairsOrcharacterOrNULL", members = c("GRanges", "Pairs", "character", "NULL"))
setClassUnion("numericOrcharacter", members = c("numeric", "character"))
setClassUnion("characterOrNULL", members = c("character", "NULL"))

#' contacts S4 class 
#' 
#' An S4 class to represent a file-stored (as .(m)cool) contact matrix 
#' imported in R.
#'
#' @slot focus Chr. coordinates for which interaction counts are extracted 
#'   from the .(m)cool file.
#' @slot metadata <etadata associated with the .(m)cool file.
#' @slot seqinfo Seqinfo, deduced from the largest resolution available 
#'   in the .(m)cool file.
#' @slot resolutions Resolutions available in the .(m)cool file.
#' @slot current_resolution Current resolution
#' @slot bins Genomic bins (computed on-the-fly) of seqinfo at the current 
#'   resolution
#' @slot interactions Genomic Interactions extracted from the .(m)cool object
#' @slot assays Available interaction scores. 
#' @slot features Genomic features associated with the dataset (e.g. 
#'   loops, borders, etc...)
#' @slot pairsFile path for the .pairs file associated with the .(m)cool file
#' @slot type Optional. Type of contacts matrix (sparse, full, aggr, ratio, ...)

methods::setClass("contacts", 
    contains = c("Annotated"), 
    slots = c(
        focus = "GRangesOrPairsOrcharacterOrNULL", 
        metadata = "list", 
        seqinfo = "Seqinfo", 
        resolutions = "numeric", 
        current_resolution = "numeric", 
        bins = "GRanges",
        interactions = "GInteractions",
        assays = "SimpleList", 
        features = "SimpleList",
        pairsFile = "characterOrNULL", 
        type = "characterOrNULL"
    )
)

#' contacts
#' 
#' @import methods
#' @export

contacts <- function(cool_path, resolution = NULL, focus = NULL, metadata = NULL, features = NULL, pairs = NULL) {
    
    ## -- Check that provided file is valid 
    check_cool_file(cool_path)
    check_cool_format(cool_path, resolution)

    ## -- Read resolutions
    resolutions <- lsCoolResolutions(cool_path, silent = TRUE)
    if (is_mcool(cool_path)) {
        res <- resolutions[length(resolutions)]
        if (!is.null(resolution)) {
            current_res <- resolution
        }
        else {
            current_res <- res
        }
    } 
    else {
        res <- resolutions[length(resolutions)]
        current_res <- NULL
    }

    ## -- Read seqinfo
    if (is_mcool(cool_path)) {
        si <- cool2seqinfo(cool_path, res)
    }
    else {
        si <- cool2seqinfo(cool_path)
    }
    
    ## -- Tile the genome
    if (is_mcool(cool_path)) {
        bins <- GenomicRanges::tileGenome(
            seqlengths = GenomeInfoDb::seqlengths(si), 
            tilewidth = current_res, 
            cut.last.tile.in.chrom = TRUE
        )    
    }
    else {
        bins <- GenomicRanges::tileGenome(
            seqlengths = GenomeInfoDb::seqlengths(si), 
            tilewidth = res, 
            cut.last.tile.in.chrom = TRUE
        )
    }

    ## -- Read interactions
    gis <- cool2gi(cool_path, res = current_res, coords = focus)
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- NULL

    ## -- Import features
    if (is.null(features)) {
        features <- S4Vectors::SimpleList(
            'loops' = GenomicRanges::GRanges(), 
            'borders' = GenomicRanges::GRanges(), 
            'compartments' = GenomicRanges::GRanges(), 
            'viewpoints' = GenomicRanges::GRanges()
        )
    }

    ## -- Check pairs file
    if (!is.null(pairs)) {
        if (!file.exists(pairs)) {
            stop("Provided pairs file does not exist. Aborting now.")
        }
        pairsFile <- pairs
    }
    else {
        pairsFile <- NULL
    }

    ## -- Create contact object
    if (methods::is(focus, 'Pairs')) {
        foc <- paste0(
            as.character(S4Vectors::first(focus)), 
            ' x ', 
            as.character(S4Vectors::second(focus))
        )
    }
    else {
        foc <- focus
    }
    x <- methods::new("contacts", 
        focus = foc, 
        metadata = purrr::flatten(c(
            list(path = cool_path), 
            metadata[names(metadata)!='path'])
        ), 
        seqinfo = si, 
        resolutions = resolutions, 
        current_resolution = ifelse(is_mcool(cool_path), current_res, res), 
        bins = bins, 
        interactions = gis, 
        assays = S4Vectors::SimpleList(
            'raw' = mcols$count, 
            'balanced' = mcols$score
        ), 
        features = features, 
        pairsFile = pairsFile, 
        type = "sparse"
    )
    methods::validObject(x)
    return(x)
} 

setValidity("contacts",
    function(object) {
        if (!is(focus(object), "GRangesOrPairsOrcharacterOrNULL"))
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

#' zoom
#' 
#' @import methods
#' @importFrom dplyr case_when
#' @importFrom S4Vectors metadata
#' @export

zoom <- function(x, focus = NULL, resolution = NULL) {
    res <- ifelse(!is.null(resolution), resolution, resolution(x))
    foc <- dplyr::case_when(
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

#' addFeature
#'
#' @param x a contacts object
#' @param feature a GRanges or GInteractions object
#' @param feature.name a name to give to the new features
#' @return The input contacts object, with a newly added feature
#' @export
#' @docType methods
#' @rdname addFeature-methods

setGeneric("addFeature", function(x, feature, feature.name) {
    standardGeneric("addFeature")
})

#' @rdname addFeature-methods
#' @aliases addFeature,contacts,GRangesOrGInteractions-method

setMethod("addFeature", 
    signature("contacts", "GRangesOrGInteractions", "character"), 
    function(x, feature, feature.name) {
        x@features[[feature.name]] <- feature
        x
    }
)

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

#' length method for objects of class \code{contacts}.
#'
#' @docType methods
#' @name length-contacts
#' @rdname length-contacts
#' @aliases length-contacts length,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export

setMethod("length", "contacts", function(x) length(regions(x)))

#' dim method for objects of class \code{contacts}.
#'
#' @docType methods
#' @name dim-contacts
#' @rdname dim-contacts
#' @aliases dim-contacts dim,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export

setMethod("dim", "contacts", function(x) dim(gi2cm(assay(x, 1))))

setMethod('[', signature(x="contacts"), function(x, i) {
    x@interactions <- interactions(x)[i]
    for (K in seq_along(assays(x))) {
        x@assays[[K]] <- x@assays[[K]][i]
    }
    return(x)
})

setGeneric("type", function(x) {standardGeneric("type")})
setMethod("type", "contacts", function(x) x@type)
setGeneric("path", function(x) {standardGeneric("path")})
setMethod("path", "contacts", function(x) S4Vectors::metadata(x)$path)
setMethod("seqinfo", "contacts", function(x) x@seqinfo)
setGeneric("resolutions", function(x) {standardGeneric("resolutions")})
setMethod("resolutions", "contacts", function(x) x@resolutions)
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

setGeneric("assay<-", function(x, name, value) {standardGeneric("assay<-")})
setMethod("assay<-", c(x = "contacts", name = "numericOrcharacter", value = "numeric"), function(x, name, value) {
    x@assays[[name]] <- value
    return(x)
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
setGeneric("feature<-", function(x, name, value) {standardGeneric("feature<-")})
setMethod("feature<-", signature(x = "contacts", name = "numericOrcharacter", value = "GRangesOrGInteractions"), function(x, name, value) {
    x@features[[name]] <- value
    return(x)
})

setGeneric("pairsFile", function(x, name) {standardGeneric("pairsFile")})
setMethod("pairsFile", "contacts", function(x) {
    x@pairsFile
})
setGeneric("pairsFile<-", function(x, value) {standardGeneric("pairsFile<-")})
setMethod("pairsFile<-", signature(x = "contacts", value = "character"), function(x, value) {
    if (!file.exists(value)) {
        stop("Provided pairs file does not exist. Aborting now.")
    }
    x@pairsFile <- value
    x
})

setMethod("anchors", "contacts", function(x) anchors(assay(x, 1)))
setMethod("regions", "contacts", function(x) regions(assay(x, 1)))

#' Show method for objects of class \code{contacts}.
#'
#' @docType methods
#' @name show-contacts
#' @rdname show-contacts
#' @aliases show-contacts show,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export

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

    ## Pairs
    cat(glue::glue('pairs: {ifelse(is.null(object@pairsFile), "N/A", object@pairsFile)}'), '\n')

})


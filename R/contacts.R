################################################################################
#                                                                              #
#                               CONSTRUCTORS                                   #
#                                                                              #
################################################################################

setClassUnion("GRangesOrGInteractions", members = c("GRanges", "GInteractions"))
setClassUnion("GRangesOrPairsOrcharacterOrNULL", members = c("GRanges", "Pairs", "character", "NULL"))
setClassUnion("numericOrcharacter", members = c("numeric", "character"))
setClassUnion("characterOrNULL", members = c("character", "NULL"))

#' @title `contacts` S4 class and methods
#'
#' @slot focus Chr. coordinates for which interaction counts are extracted 
#'   from the .(m)cool file.
#' @slot metadata metadata associated with the .(m)cool file.
#' @slot seqinfo Seqinfo, deduced from the largest resolution available 
#'   in the .(m)cool file.
#' @slot resolutions Resolutions available in the .(m)cool file.
#' @slot current_resolution Current resolution
#' @slot bins Genomic bins (computed on-the-fly) of seqinfo at the current 
#'   resolution
#' @slot interactions Genomic Interactions extracted from the .(m)cool object
#' @slot scores Available interaction scores. 
#' @slot features Genomic features associated with the dataset (e.g. 
#'   loops, borders, etc...)
#' @slot pairsFile path for the .pairs file associated with the .(m)cool file
#' @slot matrixType Type of contacts matrix (sparse, full, aggr, ratio, ...)
#' @slot coolPath Path of (m)cool file
#' 
#' @importClassesFrom S4Vectors Pairs
#' @importClassesFrom S4Vectors Annotated
#' @importFrom methods setClass
#' @export
#' @rdname contacts 

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
        scores = "SimpleList", 
        features = "SimpleList",
        pairsFile = "characterOrNULL", 
        matrixType = "characterOrNULL",
        coolPath = "character"
    )
)

#' contacts constructor
#' 
#' @param cool_path Path of a (m)cool file
#' @param resolution Resolution to use with mcool file
#' @param focus focus Chr. coordinates for which 
#'   interaction counts are extracted from the .(m)cool file.
#'   Can be provided as a character string or as a GRanges object
#' @param metadata list of metadata
#' @param features features provided as a named SimpleList
#' @param pairs Path to an associated .pairs file
#' @return a new `contacts` object.
#' 
#' @import methods
#' @rdname contacts
#' @export
#' @examples 
#' library(HiContacts)
#' data(contacts_yeast)
#' contacts_yeast

contacts <- function(
    cool_path, 
    resolution = NULL, 
    focus = NULL, 
    metadata = list(), 
    features = S4Vectors::SimpleList(
        'loops' = GenomicRanges::GRanges(), 
        'borders' = GenomicRanges::GRanges(), 
        'compartments' = GenomicRanges::GRanges(), 
        'viewpoints' = GenomicRanges::GRanges()
    ), 
    pairs = NULL
) {
    
    ## -- Check that provided file is valid 
    check_cool_file(cool_path)
    check_cool_format(cool_path, resolution)

    ## -- Read resolutions
    resolutions <- lsCoolResolutions(cool_path)
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
    gis <- cool2gi(cool_path, resolution = current_res, coords = focus)
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- NULL

    ## -- Check pairs file
    if (!is.null(pairs)) {
        if (!file.exists(pairs)) {
            stop("Provided pairs file does not exist. Aborting now.")
        }
    }
    pairsFile <- pairs

    ## -- Create contact object
    x <- methods::new("contacts", 
        focus = focus, 
        metadata = metadata, 
        seqinfo = si, 
        resolutions = resolutions, 
        current_resolution = ifelse(is_mcool(cool_path), current_res, res), 
        bins = bins, 
        interactions = gis, 
        scores = S4Vectors::SimpleList(
            'raw' = mcols$count, 
            'balanced' = mcols$score
        ), 
        features = features, 
        pairsFile = pairsFile, 
        matrixType = "sparse", 
        coolPath = as.character(cool_path)
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
        if (!is(scores(object), "SimpleList"))
            return("'scores' slot must be a SimpleList")
        TRUE
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
#' @name length
#' @docType methods
#' @rdname contacts
#' @aliases length,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' length(contacts_yeast)

setMethod("length", "contacts", function(x) length(interactions(x)))

#' [ method for objects of class \code{contacts}.
#'
#' @name [
#' @docType methods
#' @rdname contacts
#' @aliases [,contacts,ANY,ANY,ANY-method
#'
#' @param x A \code{contacts} object.
#' @param i a range or boolean vector.
#'
#' @export
#' @examples 
#' contacts_yeast[1:10]

setMethod("[", signature("contacts"), function(x, i) {
    x@interactions <- interactions(x)[i]
    for (K in seq_along(scores(x))) {
        x@scores[[K]] <- x@scores[[K]][i]
    }
    return(x)
})

#' matrixType method for objects of class \code{contacts}.
#'
#' @name matrixType
#' @docType methods
#' @rdname contacts
#' @aliases matrixType,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' matrixType(contacts_yeast)

setGeneric("matrixType", function(x) {standardGeneric("matrixType")})
setMethod("matrixType", "contacts", function(x) x@matrixType)

#' coolPath method for objects of class \code{contacts}.
#'
#' @name coolPath
#' @docType methods
#' @rdname contacts
#' @aliases coolPath,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' coolPath(contacts_yeast)

setGeneric("coolPath", function(x) {standardGeneric("coolPath")})
setMethod("coolPath", "contacts", function(x) {x@coolPath})

#' seqinfo method for objects of class \code{contacts}.
#'
#' @name seqinfo,contacts-method
#' @docType methods
#' @rdname contacts
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' seqinfo(contacts_yeast)

setMethod("seqinfo", "contacts", function(x) x@seqinfo)

#' resolutions method for objects of class \code{contacts}.
#'
#' @name resolutions
#' @docType methods
#' @rdname contacts
#' @aliases resolutions,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' resolutions(contacts_yeast)

setGeneric("resolutions", function(x) {standardGeneric("resolutions")})
setMethod("resolutions", "contacts", function(x) x@resolutions)

#' resolution method for objects of class \code{contacts}.
#'
#' @name resolution
#' @docType methods
#' @rdname contacts
#' @aliases resolution,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' resolution(contacts_yeast)

setMethod("resolution", "contacts", function(x) x@current_resolution)

#' bins method for objects of class \code{contacts}.
#'
#' @name bins
#' @docType methods
#' @rdname contacts
#' @aliases bins,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' bins(contacts_yeast)

setGeneric("bins", function(x) {standardGeneric("bins")})
setMethod("bins", "contacts", function(x) x@bins)

#' focus method for objects of class \code{contacts}.
#'
#' @name focus
#' @docType methods
#' @rdname contacts
#' @aliases focus,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' focus(contacts_yeast)

setGeneric("focus", function(x) {standardGeneric("focus")})
setMethod("focus", "contacts", function(x) x@focus)

#' interactions method for objects of class \code{contacts}.
#'
#' @name interactions
#' @docType methods
#' @rdname contacts
#' @aliases interactions,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' interactions(contacts_yeast)

setMethod("interactions", "contacts", function(x) x@interactions)

#' scores method for objects of class \code{contacts}.
#'
#' @name scores
#' @docType methods
#' @rdname contacts
#' @aliases scores,contacts,missing-method
#' @aliases scores,contacts,character-method
#' @aliases scores,contacts,numeric-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' scores(contacts_yeast)
#' scores(contacts_yeast, 1)
#' scores(contacts_yeast, 'balanced')

setGeneric("scores", function(x, name) {standardGeneric("scores")})
setMethod("scores", signature(x = "contacts", name = "missing"), function(x) x@scores)
setMethod("scores", signature(x = "contacts", name = "character"), function(x, name) {
    gis <- x@interactions
    if (!name %in% names(scores(x))) {
        stop(paste0(name, ' not in scores.'))
    }
    gis$score <- x@scores[[name]]
    return(gis)
})
setMethod("scores", signature(x = "contacts", name = "numeric"), function(x, name) {
    gis <- x@interactions
    if (name > length(scores(x))) {
        stop(paste0('Only ', length(scores(x)), ' scores in x.'))
    }
    gis$score <- x@scores[[name]]
    return(gis)
})

#' `scores<-` method for objects of class \code{contacts}.
#'
#' @name scores<-
#' @docType methods
#' @rdname contacts
#' @aliases scores<-,contacts,character,numeric-method
#'
#' @param x A \code{contacts} object.
#' @param name name
#' @param value value
#'
#' @export
#' @examples 
#' scores(contacts_yeast, 'test') <- runif(length(contacts_yeast))
#' scores(contacts_yeast, 'test')

setGeneric("scores<-", function(x, name, value) {standardGeneric("scores<-")})
setMethod("scores<-", c(x = "contacts", name = "character", value = "numeric"), function(x, name, value) {
    x@scores[[name]] <- value
    return(x)
})

#' features method for objects of class \code{contacts}.
#'
#' @name features
#' @docType methods
#' @rdname contacts
#' @aliases features,contacts,missing-method
#' @aliases features,contacts,character-method
#' @aliases features,contacts,numeric-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' data(full_contacts_yeast)
#' features(full_contacts_yeast)
#' features(full_contacts_yeast, 1)
#' features(full_contacts_yeast, 'centromeres')

setGeneric("features", function(x, name) {standardGeneric("features")})
setMethod("features", signature(x = "contacts", name = "missing"), function(x) {
    S4Vectors::SimpleList(as.list(x@features))
})
setMethod("features", signature(x = "contacts", name = "character"), function(x, name) {
    if (!name %in% names(features(x))) {
        stop(paste0(name, ' not in features.'))
    }
    x@features[[name]]
})
setMethod("features", signature(x = "contacts", name = "numeric"), function(x, name) {
    if (name > length(features(x))) {
        stop(paste0('Only ', length(features(x)), ' features in x.'))
    }
    x@features[[name]]
})

#' `features<-` method for objects of class \code{contacts}.
#'
#' @name features<-
#' @docType methods
#' @rdname contacts
#' @aliases features<-,contacts,character,GRangesOrGInteractions-method
#'
#' @param x A \code{contacts} object.
#' @param name name
#' @param value value
#'
#' @export
#' @examples 
#' data(centros_yeast)
#' features(contacts_yeast, 'centromeres') <- centros_yeast
#' features(contacts_yeast, 'centromeres')

setGeneric("features<-", function(x, name, value) {standardGeneric("features<-")})
setMethod("features<-", signature(x = "contacts", name = "character", value = "GRangesOrGInteractions"), function(x, name, value) {
    x@features[[name]] <- value
    return(x)
})

#' pairsFile method for objects of class \code{contacts}.
#'
#' @name pairsFile
#' @docType methods
#' @rdname contacts
#' @aliases pairsFile,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' pairsFile(full_contacts_yeast)

setGeneric("pairsFile", function(x, name) {standardGeneric("pairsFile")})
setMethod("pairsFile", "contacts", function(x) {
    x@pairsFile
})

#' `pairsFile<-` method for objects of class \code{contacts}.
#'
#' @name pairsFile<-
#' @docType methods
#' @rdname contacts
#' @aliases pairsFile<-,contacts,character-method
#'
#' @param x A \code{contacts} object.
#' @param name name
#' @param value value
#'
#' @export

setGeneric("pairsFile<-", function(x, value) {standardGeneric("pairsFile<-")})
setMethod("pairsFile<-", signature(x = "contacts", value = "character"), function(x, value) {
    if (!file.exists(value)) {
        stop("Provided pairs file does not exist. Aborting now.")
    }
    x@pairsFile <- value
    x
})

#' anchors method for objects of class \code{contacts}.
#'
#' @name anchors
#' @docType methods
#' @rdname contacts
#' @aliases anchors,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' anchors(contacts_yeast)

setMethod("anchors", "contacts", function(x) anchors(scores(x, 1)))

#' regions method for objects of class \code{contacts}.
#'
#' @name regions
#' @docType methods
#' @rdname contacts
#' @aliases regions,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' regions(contacts_yeast)

setMethod("regions", "contacts", function(x) regions(scores(x, 1)))

#' summary method for objects of class \code{contacts}.
#'
#' @name summary
#' @docType methods
#' @rdname contacts
#' @aliases summary,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' summary(contacts_yeast)

setMethod("summary", "contacts", function(object) {
    cat(glue::glue(
        '{matrixType(object)} `contacts` object with {format(length(interactions(object)), big.mark = ",")} interactions over {format(length(regions(object)), big.mark = ",")} regions'
    ), '\n')
})


#' Show method for objects of class \code{contacts}.
#'
#' @name show
#' @docType methods
#' @rdname contacts
#' @aliases show,contacts-method
#'
#' @param object A \code{contacts} object.
#'
#' @export
#' @examples 
#' show(contacts_yeast)

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
    else if (is(focus(object), 'Pairs')) {
        focus_str <- paste0(
            as.character(S4Vectors::first(focus(object))), 
            ' x ', 
            as.character(S4Vectors::second(focus(object)))
        )
    }
    else if (is(focus(object), 'character')) {
        focus_str <- formatCoords(focus(object))
    }

    cat(summary(object))
    cat(glue::glue('coolPath: {coolPath(object)}'), '\n')
    cat(glue::glue('focus: {focus_str}'), '\n')
    cat('------------\n')

    ## Metadata
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))

    ## Resolutions
    S4Vectors::coolcat("resolutions(%d): %s\n", resolutions(object))
    cat(glue::glue('current resolution({resolution(object)}): {length(interactions(object))} interactions in memory'), '\n')

    ## Bins
    cat(glue::glue('bins: {length(bins(object))}'), '\n')

    ## Regions
    cat(glue::glue('regions: {length(regions(object))}'), '\n')

    ## Scores
    cat(glue::glue('scores({length(scores(object))}): {paste(names(scores(object)), collapse = " ")}'), '\n')

    ## Features
    cat(glue::glue('features: {paste(paste0(names(features(object)), "(", lengths(features(object)), ")"), collapse = " ")}'), '\n')

    ## Pairs
    cat(glue::glue('pairs: {ifelse(is.null(pairsFile(object)), "N/A", pairsFile(object))}'), '\n')

})


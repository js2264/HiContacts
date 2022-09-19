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
#' @import methods
#' @importClassesFrom S4Vectors Pairs
#' @importClassesFrom S4Vectors Annotated
#' @importFrom methods setClass
#' @export
#' @rdname contacts 

methods::setClass("contacts", 
    contains = c("Annotated"), 
    slots = c(
        focus = "characterOrNULL", 
        metadata = "list", 
        seqinfo = "Seqinfo", 
        resolutions = "numeric", 
        current_resolution = "numeric", 
        bins = "GRanges",
        interactions = "GInteractions",
        scores = "SimpleList", 
        features = "SimpleList",
        pairsFile = "characterOrNULL", 
        matrixType = "character",
        coolPath = "character"
    )
)

#' @rdname contacts
#' 
#' @param path Path of a (m)cool file
#' @param resolution Resolution to use with mcool file
#' @param focus focus Chr. coordinates for which 
#'   interaction counts are extracted from the .(m)cool file, provided
#'   as a character string (e.g. "II:4000-5000"). If not provided, 
#'   the entire (m)cool file will be imported. 
#' @param metadata list of metadata
#' @param features features provided as a named SimpleList
#' @param pairs Path to an associated .pairs file
#' @return a new `contacts` object.
#' 
#' @export
#' @examples 
#' library(HiContacts)
#' data(contacts_yeast)
#' contacts_yeast

contacts <- function(
    path, 
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
    check_cool_file(path)
    check_cool_format(path, resolution)

    ## -- Read resolutions
    resolutions <- lsCoolResolutions(path)
    if (is_mcool(path)) {
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
    if (is_mcool(path)) {
        si <- cool2seqinfo(path, res)
    }
    else {
        si <- cool2seqinfo(path)
    }
    
    ## -- Tile the genome
    if (is_mcool(path)) {
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
    gis <- cool2gi(path, resolution = current_res, coords = focus)
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
        current_resolution = ifelse(is_mcool(path), current_res, res), 
        bins = bins, 
        interactions = gis, 
        scores = S4Vectors::SimpleList(
            'raw' = mcols$count, 
            'balanced' = mcols$score
        ), 
        features = features, 
        pairsFile = pairsFile, 
        matrixType = "sparse", 
        coolPath = as.character(path)
    )
    methods::validObject(x)
    return(x)
} 

setValidity("contacts",
    function(object) {
        if (!is(focus(object), "characterOrNULL"))
            return("'focus' slot must be a characterOrNULL")
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
#                                 ACCESSORS                                    #
#                                                                              #
################################################################################

#' @rdname contacts
#'
#' @name length
#' @docType methods
#' @aliases length,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' length(contacts_yeast)

setMethod("length", "contacts", function(x) length(interactions(x)))

#' @rdname contacts
#'
#' @name [
#' @docType methods
#' @aliases [,contacts,numeric,ANY,ANY-method
#' @aliases [,contacts,logical,ANY,ANY-method
#' @aliases [,contacts,character,ANY,ANY-method
#'
#' @param x A \code{contacts} object.
#' @param i a range or boolean vector.
#'
#' @export
#' @examples 
#' contacts_yeast[1:10]

setMethod("[", signature("contacts", "numeric"), function(x, i) {
    interactions(x) <- interactions(x)[i]
    for (K in seq_along(scores(x))) {
        x@scores[[K]] <- x@scores[[K]][i]
    }
    matrixType(x) <- "subset"
    return(x)
})
setMethod("[", signature("contacts", "logical"), function(x, i) {
    interactions(x) <- interactions(x)[i]
    for (K in seq_along(scores(x))) {
        x@scores[[K]] <- x@scores[[K]][i]
    }
    matrixType(x) <- "subset"
    return(x)
})
setMethod("[", signature("contacts", "character"), function(x, i) {
    `%over%` <- IRanges::`%over%`
    i_ <- char2pair(i)
    sub <- anchors(x)[['first']] %over% S4Vectors::first(i_) & 
        anchors(x)[['second']] %over% S4Vectors::second(i_)
    interactions(x) <- interactions(x)[sub]
    for (K in seq_along(scores(x))) {
        x@scores[[K]] <- x@scores[[K]][sub]
    }
    matrixType(x) <- "subset"
    return(x)
})

#' @rdname contacts
#'
#' @name matrixType
#' @docType methods
#' @aliases matrixType,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' matrixType(contacts_yeast)

setGeneric("matrixType", function(x) {standardGeneric("matrixType")})
setMethod("matrixType", "contacts", function(x) x@matrixType)

#' @rdname contacts
#'
#' @name matrixType<-
#' @docType methods
#' @aliases matrixType<-,contacts,character-method
#'
#' @param x A \code{contacts} object.
#' @param name name
#' @param value value
#'
#' @export
#' @examples 
#' matrixType(contacts_yeast) <- "custom"
#' matrixType(contacts_yeast)

setGeneric("matrixType<-", function(x, value) {standardGeneric("matrixType<-")})
setMethod("matrixType<-", c(x = "contacts", value = "character"), function(x, value) {
    x@matrixType <- value
    return(x)
})

#' @rdname contacts
#'
#' @name coolPath
#' @docType methods
#' @aliases coolPath,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' coolPath(contacts_yeast)

setGeneric("coolPath", function(x) {standardGeneric("coolPath")})
setMethod("coolPath", "contacts", function(x) {x@coolPath})

#' @rdname contacts
#'
#' @name seqinfo,contacts-method
#' @docType methods
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' seqinfo(contacts_yeast)

setMethod("seqinfo", "contacts", function(x) x@seqinfo)

#' @rdname contacts
#'
#' @name resolutions
#' @docType methods
#' @aliases resolutions,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' resolutions(contacts_yeast)

setGeneric("resolutions", function(x) {standardGeneric("resolutions")})
setMethod("resolutions", "contacts", function(x) x@resolutions)

#' @rdname contacts
#'
#' @name resolution
#' @docType methods
#' @aliases resolution,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' resolution(contacts_yeast)

setMethod("resolution", "contacts", function(x) x@current_resolution)

#' @rdname contacts
#'
#' @name bins
#' @docType methods
#' @aliases bins,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' bins(contacts_yeast)

setGeneric("bins", function(x) {standardGeneric("bins")})
setMethod("bins", "contacts", function(x) x@bins)

#' @rdname contacts
#'
#' @name focus
#' @docType methods
#' @aliases focus,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' focus(contacts_yeast)

setGeneric("focus", function(x) {standardGeneric("focus")})
setMethod("focus", "contacts", function(x) x@focus)

#' @rdname contacts
#'
#' @name interactions
#' @docType methods
#' @aliases interactions,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' interactions(contacts_yeast)

setMethod("interactions", "contacts", function(x) x@interactions)

#' @rdname contacts
#'
#' @name interactions<-
#' @docType methods
#' @aliases interactions<-,contacts,GInteractions-method
#'
#' @param x A \code{contacts} object.
#' @param name name
#' @param value value
#'
#' @export

setMethod("interactions<-", signature(x = "contacts", value = "GInteractions"), function(x, value) {
    x@interactions <- value
    x
})

#' @rdname contacts
#'
#' @name scores
#' @docType methods
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

#' @rdname contacts
#'
#' @name scores<-
#' @docType methods
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

#' @rdname contacts
#' 
#' @name features
#' @docType methods
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

#' @rdname contacts
#' 
#' @name features<-
#' @docType methods
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

#' @rdname contacts
#' 
#' @name pairsFile
#' @docType methods
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

#' @rdname contacts
#' 
#' @name pairsFile<-
#' @docType methods
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

#' @rdname contacts
#' 
#' @name anchors
#' @docType methods
#' @aliases anchors,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' anchors(contacts_yeast)

setMethod("anchors", "contacts", function(x) anchors(scores(x, 1)))

#' @rdname contacts
#' 
#' @name regions
#' @docType methods
#' @aliases regions,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' regions(contacts_yeast)

setMethod("regions", "contacts", function(x) regions(scores(x, 1)))

#' @rdname contacts
#' 
#' @name summary
#' @docType methods
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


#' @rdname contacts
#' 
#' @name show
#' @docType methods
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
    else {
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


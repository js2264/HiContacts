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
#' @slot fileName Path of (m)cool file
#' @slot focus Chr. coordinates for which interaction counts are extracted 
#'   from the .(m)cool file.
#' @slot resolutions Resolutions available in the .(m)cool file.
#' @slot resolution Current resolution
#' @slot interactions Genomic Interactions extracted from the .(m)cool object
#' @slot scores Available interaction scores. 
#' @slot topologicalFeatures Topological features associated with the dataset 
#'   (e.g. loops (\<Pairs\>), borders (\<GRanges\>), 
#'   viewpoints (\<GRanges\>), etc...)
#' @slot pairsFile Path to the .pairs file associated with the .(m)cool file
#' @slot metadata metadata associated with the .(m)cool file.
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
        fileName = "character",
        focus = "characterOrNULL", 
        resolutions = "numeric", 
        resolution = "numeric", 
        interactions = "GInteractions",
        scores = "SimpleList", 
        topologicalFeatures = "SimpleList",
        pairsFile = "characterOrNULL",
        metadata = "list"
    )
)

#' @rdname contacts
#' 
#' @param file Path to a (m)cool file
#' @param resolution Resolution to use with mcool file
#' @param focus focus Chr. coordinates for which 
#'   interaction counts are extracted from the .(m)cool file, provided
#'   as a character string (e.g. "II:4000-5000"). If not provided, 
#'   the entire (m)cool file will be imported. 
#' @param metadata list of metadata
#' @param topologicalFeatures topologicalFeatures provided as a named SimpleList
#' @param pairsFile Path to an associated .pairs file
#' @return a new `contacts` object.
#' 
#' @export
#' @examples 
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' contacts_yeast

contacts <- function(
    file, 
    resolution = NULL, 
    focus = NULL, 
    metadata = list(), 
    topologicalFeatures = S4Vectors::SimpleList(
        'loops' = S4Vectors::Pairs(
            GenomicRanges::GRanges(), 
            GenomicRanges::GRanges()
        ), 
        'borders' = GenomicRanges::GRanges(), 
        'compartments' = GenomicRanges::GRanges(), 
        'viewpoints' = GenomicRanges::GRanges()
    ), 
    pairsFile = NULL
) {
    
    ## -- Check that provided file is valid 
    check_cool_file(file)
    check_cool_format(file, resolution)

    ## -- Read resolutions
    resolutions <- lsCoolResolutions(file)
    if (is_mcool(file)) {
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

    ## -- Read interactions
    gis <- cool2gi(file, resolution = current_res, coords = focus)
    mcols <- GenomicRanges::mcols(gis)
    GenomicRanges::mcols(gis) <- mcols[, c('bin_id1', 'bin_id2')]

    ## -- Check pairs file
    if (!is.null(pairsFile)) {
        if (!file.exists(pairsFile)) {
            stop("Provided pairsFile does not exist. Aborting now.")
        }
    }

    ## -- Create contact object
    x <- methods::new("contacts", 
        fileName = as.character(file),
        focus = focus, 
        resolutions = resolutions, 
        resolution = ifelse(is_mcool(file), current_res, res), 
        interactions = gis, 
        scores = S4Vectors::SimpleList(
            'raw' = as.numeric(mcols$count),
            'balanced' = as.numeric(mcols$score)
        ), 
        topologicalFeatures = topologicalFeatures, 
        pairsFile = pairsFile, 
        metadata = metadata
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
#' @name fileName
#' @docType methods
#' @aliases fileName,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @examples 
#' fileName(contacts_yeast)

setMethod("fileName", "contacts", function(object) object@fileName)

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

setMethod("resolution", "contacts", function(x) x@resolution)

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
#' tail(scores(contacts_yeast, 1))
#' tail(scores(contacts_yeast, 'balanced'))

setGeneric("scores", function(x, name) {standardGeneric("scores")})
setMethod("scores", signature(x = "contacts", name = "missing"), function(x) x@scores)
setMethod("scores", signature(x = "contacts", name = "character"), function(x, name) {
    if (!name %in% names(scores(x))) {
        stop(paste0(name, ' not in scores.'))
    }
    return(x@scores[[name]])
})
setMethod("scores", signature(x = "contacts", name = "numeric"), function(x, name) {
    if (name > length(scores(x))) {
        stop(paste0('Only ', length(scores(x)), ' scores in x.'))
    }
    return(x@scores[[name]])
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
#' tail(scores(contacts_yeast, 'test'))

setGeneric("scores<-", function(x, name, value) {standardGeneric("scores<-")})
setMethod("scores<-", c(x = "contacts", name = "character", value = "numeric"), function(x, name, value) {
    x@scores[[name]] <- value
    return(x)
})

#' @rdname contacts
#' 
#' @name topologicalFeatures
#' @docType methods
#' @aliases topologicalFeatures,contacts,missing-method
#' @aliases topologicalFeatures,contacts,character-method
#' @aliases topologicalFeatures,contacts,numeric-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' full_contacts_yeast <- full_contacts_yeast()
#' topologicalFeatures(full_contacts_yeast)
#' topologicalFeatures(full_contacts_yeast, 1)
#' topologicalFeatures(full_contacts_yeast, 'centromeres')

setGeneric("topologicalFeatures", function(x, name) {standardGeneric("topologicalFeatures")})
setMethod("topologicalFeatures", signature(x = "contacts", name = "missing"), function(x) {
    S4Vectors::SimpleList(as.list(x@topologicalFeatures))
})
setMethod("topologicalFeatures", signature(x = "contacts", name = "character"), function(x, name) {
    if (!name %in% names(topologicalFeatures(x))) {
        stop(paste0(name, ' not in topologicalFeatures.'))
    }
    x@topologicalFeatures[[name]]
})
setMethod("topologicalFeatures", signature(x = "contacts", name = "numeric"), function(x, name) {
    if (name > length(topologicalFeatures(x))) {
        stop(paste0('Only ', length(topologicalFeatures(x)), ' topologicalFeatures in x.'))
    }
    x@topologicalFeatures[[name]]
})

#' @rdname contacts
#' 
#' @name topologicalFeatures<-
#' @docType methods
#' @aliases topologicalFeatures<-,contacts,character,GRangesOrGInteractions-method
#'
#' @param x A \code{contacts} object.
#' @param name name
#' @param value value
#'
#' @export
#' @examples 
#' data(centros_yeast)
#' topologicalFeatures(contacts_yeast, 'centromeres') <- centros_yeast
#' topologicalFeatures(contacts_yeast, 'centromeres')

setGeneric("topologicalFeatures<-", function(x, name, value) {standardGeneric("topologicalFeatures<-")})
setMethod("topologicalFeatures<-", signature(x = "contacts", name = "character", value = "GRangesOrGInteractions"), function(x, name, value) {
    x@topologicalFeatures[[name]] <- value
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
        stop("Provided pairsFile does not exist. Aborting now.")
    }
    x@pairsFile <- value
    x
})

#' @rdname contacts
#' 
#' @name metadata<-
#' @docType methods
#' @aliases metadata<-,contacts,list-method
#'
#' @param x A \code{contacts} object.
#' @param name name
#' @param value value
#'
#' @export

setGeneric("metadata<-", function(x, value) {standardGeneric("metadata<-")})
setMethod("metadata<-", signature(x = "contacts", value = "list"), function(x, value) {
    x@metadata <- value
    x
})

################################################################################
#                                                                              #
#                                 OTHER METHODS                                #
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
#' @importFrom InteractionSet reduceRegions
#' @export
#' @examples 
#' contacts_yeast[seq_len(10)]

setMethod("[", signature("contacts", "numeric"), function(x, i) {
    interactions(x) <- InteractionSet::reduceRegions(
        interactions(x)[i]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[i]
    }
    return(x)
})
setMethod("[", signature("contacts", "logical"), function(x, i) {
    interactions(x) <- InteractionSet::reduceRegions(
        interactions(x)[i]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[i]
    }
    return(x)
})
setMethod("[", signature("contacts", "character"), function(x, i) {
    i_ <- char2pair(i)
    bi_ <- bins(x)
    ints_ <- interactions(x)
    valid_bins_first <- subsetByOverlaps(
        bi_, S4Vectors::first(i_), type = 'within'
    )$bin_id
    valid_bins_second <- subsetByOverlaps(
        bi_, S4Vectors::second(i_), type = 'within'
    )$bin_id
    sub <- ints_$bin_id1 %in% valid_bins_first & ints_$bin_id2 %in% valid_bins_second
    interactions(x) <- InteractionSet::reduceRegions(
        ints_[sub]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[sub]
    }
    return(x)
})

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

setMethod("seqinfo", "contacts", function(x) {
    if (is_mcool(fileName(x))) {
        si <- cool2seqinfo(fileName(x), resolution(x))
    }
    else {
        si <- cool2seqinfo(fileName(x))
    }
    return(si)
})

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
setMethod("bins", "contacts", function(x) {
    bins <- GenomicRanges::tileGenome(
        seqlengths = GenomeInfoDb::seqlengths(seqinfo(x)), 
        tilewidth = resolution(x), 
        cut.last.tile.in.chrom = TRUE
    )
    seqinfo(bins) <- seqinfo(x)
    bins$bin_id <- seq_along(bins)
    # bins <- bins[GenomicRanges::width(bins) == resolution(x)]
    return(bins)
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

setMethod("anchors", "contacts", function(x) anchors(interactions(x)))

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

setMethod("regions", "contacts", function(x) regions(interactions(x)))

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
        '`contacts` object with {format(length(interactions(object)), big.mark = ",")} interactions over {format(length(regions(object)), big.mark = ",")} regions'
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
    cat('-------\n')
    cat(glue::glue('fileName: "{fileName(object)}"'), '\n')
    cat(glue::glue('focus: "{focus_str}"'), '\n')

    ## Resolutions
    S4Vectors::coolcat("resolutions(%d): %s\n", resolutions(object))
    cat(glue::glue('current resolution: {resolution(object)}'), '\n')

    ## Interactions
    cat(glue::glue('interactions: {length(interactions(object))}'), '\n')

    ## Scores
    cat(glue::glue('scores({length(scores(object))}): {paste(names(scores(object)), collapse = " ")}'), '\n')

    ## topologicalFeatures
    cat(glue::glue('topologicalFeatures: {paste(paste0(names(topologicalFeatures(object)), "(", lengths(topologicalFeatures(object)), ")"), collapse = " ")}'), '\n')

    ## Pairs
    cat(glue::glue('pairsFile: {ifelse(is.null(pairsFile(object)), "N/A", pairsFile(object))}'), '\n')

    ## Metadata
    S4Vectors::coolcat("metadata(%d): %s\n", names(S4Vectors::metadata(object)))

})

################################################################################
#                                                                              #
#                                 COERCING                                    -#
#                                                                              #
################################################################################

#' @rdname contacts
#' 
#' @name setAs
#' @docType methods
#' @aliases setAs,contacts-method
#'
#' @param x A \code{contacts} object.
#'
#' @export
#' @examples 
#' as(contacts_yeast, 'GInteractions')
#' as(contacts_yeast, 'ContactMatrix')
#' as(contacts_yeast, 'matrix')[seq_len(10), seq_len(10)]

setAs("contacts", "GInteractions", function(from) interactions(from))
setAs("contacts", "ContactMatrix", function(from) {
    if ('balanced' %in% names(scores(from))) {
        x <- interactions(from)
        x$score <- scores(from, 'balanced')
        gi2cm(x)
    } 
    else {
        x <- interactions(from)
        x$score <- scores(from, 1)
        gi2cm(x)
    }
})
setAs("contacts", "matrix", function(from) {
    as(from, "ContactMatrix") |> cm2matrix()
})

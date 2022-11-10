################################################################################
################################################################################
###############                                                  ###############
###############          METHODS FOR NEW GENERICS                ###############
###############                                                  ###############
################################################################################
################################################################################

#' @rdname Contacts
#'
#' @name resolutions
#' @docType methods
#' @aliases resolutions,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' resolutions(contacts_yeast)

setMethod("resolutions", "Contacts", function(x) x@resolutions)

#' @rdname Contacts
#'
#' @name resolution
#' @docType methods
#' @aliases resolution,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' resolution(contacts_yeast)

setMethod("resolution", "Contacts", function(x) x@resolution)

#' @rdname Contacts
#'
#' @name focus
#' @docType methods
#' @aliases focus,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' focus(contacts_yeast)

setMethod("focus", "Contacts", function(x) x@focus)

#' @rdname Contacts
#'
#' @name focus<-
#' @docType methods
#' @aliases focus<-,Contacts,character-method
#'
#' @param x A \code{Contacts} object.
#' @param name name
#' @param value value
#'
#' @export

setMethod("focus<-", signature(x = "Contacts", value = "character"), function(x, value) {
    x@focus <- value
    x
})

#' @rdname Contacts
#'
#' @name scores
#' @docType methods
#' @aliases scores,Contacts,missing-method
#' @aliases scores,Contacts,character-method
#' @aliases scores,Contacts,numeric-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' scores(contacts_yeast)
#' tail(scores(contacts_yeast, 1))
#' tail(scores(contacts_yeast, 'balanced'))

setMethod("scores", signature(x = "Contacts", name = "missing"), function(x) x@scores)
setMethod("scores", signature(x = "Contacts", name = "character"), function(x, name) {
    if (!name %in% names(scores(x))) {
        stop(paste0(name, ' not in scores.'))
    }
    return(x@scores[[name]])
})
setMethod("scores", signature(x = "Contacts", name = "numeric"), function(x, name) {
    if (name > length(scores(x))) {
        stop(paste0('Only ', length(scores(x)), ' scores in x.'))
    }
    return(x@scores[[name]])
})

#' @rdname Contacts
#'
#' @name scores<-
#' @docType methods
#' @aliases scores<-,Contacts,character,numeric-method
#'
#' @param x A \code{Contacts} object.
#' @param name name
#' @param value value
#'
#' @export
#' @examples 
#' scores(contacts_yeast, 'test') <- runif(length(contacts_yeast))
#' tail(scores(contacts_yeast, 'test'))

setMethod("scores<-", c(x = "Contacts", name = "character", value = "numeric"), function(x, name, value) {
    x@scores[[name]] <- value
    return(x)
})

#' @rdname Contacts
#' 
#' @name topologicalFeatures
#' @docType methods
#' @aliases topologicalFeatures,Contacts,missing-method
#' @aliases topologicalFeatures,Contacts,character-method
#' @aliases topologicalFeatures,Contacts,numeric-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' full_contacts_yeast <- full_contacts_yeast()
#' topologicalFeatures(full_contacts_yeast)
#' topologicalFeatures(full_contacts_yeast, 1)
#' topologicalFeatures(full_contacts_yeast, 'centromeres')

setMethod("topologicalFeatures", signature(x = "Contacts", name = "missing"), function(x) {
    S4Vectors::SimpleList(as.list(x@topologicalFeatures))
})
setMethod("topologicalFeatures", signature(x = "Contacts", name = "character"), function(x, name) {
    if (!name %in% names(topologicalFeatures(x))) {
        stop(paste0(name, ' not in topologicalFeatures.'))
    }
    x@topologicalFeatures[[name]]
})
setMethod("topologicalFeatures", signature(x = "Contacts", name = "numeric"), function(x, name) {
    if (name > length(topologicalFeatures(x))) {
        stop(paste0('Only ', length(topologicalFeatures(x)), ' topologicalFeatures in x.'))
    }
    x@topologicalFeatures[[name]]
})

#' @rdname Contacts
#' 
#' @name topologicalFeatures<-
#' @docType methods
#' @aliases topologicalFeatures<-,Contacts,character,GRangesOrGInteractions-method
#'
#' @param x A \code{Contacts} object.
#' @param name name
#' @param value value
#'
#' @export
#' @examples 
#' data(centros_yeast)
#' topologicalFeatures(contacts_yeast, 'centromeres') <- centros_yeast
#' topologicalFeatures(contacts_yeast, 'centromeres')

setMethod("topologicalFeatures<-", signature(x = "Contacts", name = "character", value = "GRangesOrGInteractions"), function(x, name, value) {
    x@topologicalFeatures[[name]] <- value
    return(x)
})

#' @rdname Contacts
#' 
#' @name pairsFile
#' @docType methods
#' @aliases pairsFile,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' pairsFile(full_contacts_yeast)

setMethod("pairsFile", "Contacts", function(x) {
    x@pairsFile
})

#' @rdname Contacts
#' 
#' @name pairsFile<-
#' @docType methods
#' @aliases pairsFile<-,Contacts,character-method
#'
#' @param x A \code{Contacts} object.
#' @param name name
#' @param value value
#'
#' @export

setMethod("pairsFile<-", signature(x = "Contacts", value = "character"), function(x, value) {
    if (!file.exists(value)) {
        stop("Provided pairsFile does not exist. Aborting now.")
    }
    x@pairsFile <- value
    x
})

#' @rdname Contacts
#' 
#' @name metadata<-
#' @docType methods
#' @aliases metadata<-,Contacts,list-method
#'
#' @param x A \code{Contacts} object.
#' @param name name
#' @param value value
#'
#' @export

setMethod("metadata<-", signature(x = "Contacts", value = "list"), function(x, value) {
    x@metadata <- value
    x
})

################################################################################
################################################################################
###############                                                  ###############
###############          METHODS FOR EXISTING GENERICS           ###############
###############                                                  ###############
################################################################################
################################################################################

#' @rdname Contacts
#'
#' @name fileName
#' @docType methods
#' @aliases fileName,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @importMethodsFrom BiocGenerics fileName
#' @export
#' @examples 
#' fileName(contacts_yeast)

setMethod("fileName", "Contacts", function(object) object@fileName)

#' @rdname Contacts
#'
#' @name interactions
#' @docType methods
#' @aliases interactions,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' interactions(contacts_yeast)

setMethod("interactions", "Contacts", function(x) x@interactions)

#' @rdname Contacts
#'
#' @name interactions<-
#' @docType methods
#' @aliases interactions<-,Contacts,GInteractions-method
#'
#' @param x A \code{Contacts} object.
#' @param name name
#' @param value value
#'
#' @export

setMethod("interactions<-", signature(x = "Contacts", value = "GInteractions"), function(x, value) {
    x@interactions <- value
    x
})

#' @rdname Contacts
#'
#' @name length
#' @docType methods
#' @aliases length,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' length(contacts_yeast)

setMethod("length", "Contacts", function(x) length(interactions(x)))

#' @rdname Contacts
#'
#' @name [
#' @docType methods
#' @aliases [,Contacts,numeric,ANY,ANY-method
#' @aliases [,Contacts,logical,ANY,ANY-method
#' @aliases [,Contacts,character,ANY,ANY-method
#'
#' @param x A \code{Contacts} object.
#' @param i a range or boolean vector.
#'
#' @importFrom InteractionSet reduceRegions
#' @export
#' @examples 
#' contacts_yeast[seq_len(10)]

setMethod("[", signature("Contacts", "numeric"), function(x, i) {
    interactions(x) <- InteractionSet::reduceRegions(
        interactions(x)[i]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[i]
    }
    return(x)
})
setMethod("[", signature("Contacts", "logical"), function(x, i) {
    interactions(x) <- InteractionSet::reduceRegions(
        interactions(x)[i]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[i]
    }
    return(x)
})
setMethod("[", signature("Contacts", "character"), function(x, i) {
    re_ <- regions(x)
    ints_ <- interactions(x)
    if (length(i) == 1) {
        if (grepl(
            '[A-Za-z0-9]*:[0-9]*-[0-9]* [xX/-;] [A-Za-z0-9]*:[0-9]*-[0-9]*$', i
        )) {
            i_ <- char2coords(i)
            valid_regions_first <- subsetByOverlaps(
                re_, S4Vectors::first(i_), type = 'within'
            )$bin_id
            valid_regions_second <- subsetByOverlaps(
                re_, S4Vectors::second(i_), type = 'within'
            )$bin_id
        }
        else if (grepl(
            '[A-Za-z0-9]*:[0-9]*-[0-9]*$', i
        )) {
            i_ <- char2coords(i)
            valid_regions_first <- subsetByOverlaps(
                re_, S4Vectors::first(i_), type = 'within'
            )$bin_id
            valid_regions_second <- valid_regions_first
        }
        else if (
            i %in% seqnames(seqinfo(x))
        ){
            valid_regions_first <- re_$bin_id[as.vector(seqnames(re_)) %in% i]
            valid_regions_second <- valid_regions_first
        }
        else {
            stop("Failed to coerce i into a Pairs/GRanges/chr.")
        }
    }
    else {
        if (
            all(i %in% seqnames(seqinfo(x)))
        ){
            valid_regions_first <- re_$bin_id[as.vector(seqnames(re_)) %in% i]
            valid_regions_second <- valid_regions_first
        }
        else {
            stop("Failed to coerce i into a valid Pairs/GRanges/chr.")
        }
    }
    sub <- ints_$bin_id1 %in% valid_regions_first & ints_$bin_id2 %in% valid_regions_second
    interactions(x) <- InteractionSet::reduceRegions(
        ints_[sub]
    )
    for (n in names(scores(x))) {
        scores(x, n) <- scores(x, n)[sub]
    }
    focus(x) <- i
    return(x)
})

#' @rdname Contacts
#'
#' @name seqinfo,Contacts-method
#' @docType methods
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' seqinfo(contacts_yeast)

setMethod("seqinfo", "Contacts", function(x) {
    if (is_mcool(fileName(x))) {
        si <- cool2seqinfo(fileName(x), resolution(x))
    }
    else {
        si <- cool2seqinfo(fileName(x))
    }
    return(si)
})

#' @rdname Contacts
#' 
#' @name anchors
#' @docType methods
#' @aliases anchors,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' anchors(contacts_yeast)

setMethod("anchors", "Contacts", function(x) anchors(interactions(x)))

#' @rdname Contacts
#' 
#' @name regions
#' @docType methods
#' @aliases regions,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' regions(contacts_yeast)

setMethod("regions", "Contacts", function(x) regions(interactions(x)))

#' @rdname Contacts
#' 
#' @name summary
#' @docType methods
#' @aliases summary,Contacts-method
#'
#' @param x A \code{Contacts} object.
#'
#' @export
#' @examples 
#' summary(contacts_yeast)

setMethod("summary", "Contacts", function(object) {
    cat(glue::glue(
        '`Contacts` object with {format(length(interactions(object)), big.mark = ",")} interactions over {format(length(regions(object)), big.mark = ",")} regions'
    ), '\n')
})

#' @rdname Contacts
#' 
#' @name show
#' @docType methods
#' @aliases show,Contacts-method
#'
#' @param object A \code{Contacts} object.
#'
#' @export
#' @examples 
#' show(contacts_yeast)

setMethod("show", signature("Contacts"), function(object) {

    if (is.null(focus(object))) {
        focus_str <- "whole genome"
    } 
    else {
        focus_str <- coords2char(focus(object))
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

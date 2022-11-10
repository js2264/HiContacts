################################################################################
#                                                                              #
#                               CONSTRUCTORS                                   #
#                                                                              #
################################################################################

#' @rdname Contacts
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
#' @return a new `Contacts` object.
#' 
#' @export
#' @examples 
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' contacts_yeast

Contacts <- function(
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
    x <- methods::new("Contacts", 
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

setValidity("Contacts",
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


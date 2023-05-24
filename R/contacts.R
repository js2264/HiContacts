################################################################################
#                                                                              #
#                               CONSTRUCTORS                                   #
#                                                                              #
################################################################################

#' @title Contacts
#' 
#' @rdname Contacts
#' 
#' @description 
#' 
#' This function has been deprecated in favor of the generic `HiCExperiment()`
#' constructor (from HiCExperiment package).
#' 
#' @param file Path to a (m)cool file
#' @param resolution Resolution to use with mcool file
#' @param focus focus Chr. coordinates for which 
#'   interaction counts are extracted from the .(m)cool file, provided
#'   as a character string (e.g. "II:4001-5000"). If not provided, 
#'   the entire (m)cool file will be imported. 
#' @param metadata list of metadata
#' @param topologicalFeatures topologicalFeatures provided as a named SimpleList
#' @param pairsFile Path to an associated .pairs file
#' @return a new `HiCExperiment` object.
#' 
#' @import HiCExperiment
#' @export
#' @examples 
#' library(HiContacts)
#' library(HiContactsData)
#' mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#' Contacts(mcool_path, resolution = 1000)

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
    .Deprecated(
            "HiCExperiment::HiCExperiment", 
            msg = paste(
                "`HiContacts::Contacts` is deprecated;", 
                "see '?HiCExperiment::HiCExperiment' constructor instead."
            )
    )
    HiCExperiment::HiCExperiment(
        file = file, 
        resolution = resolution, 
        focus = focus, 
        metadata = metadata, 
        topologicalFeatures = topologicalFeatures, 
        pairsFile = pairsFile
    )
} 

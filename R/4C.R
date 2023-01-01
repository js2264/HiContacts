#' Computing virtual 4C profiles
#' 
#' From a (m)cool pre-imported in memory, computes a 4C profile 
#' using a user-specified `viewpoint`. 
#'
#' @param x a `HiCExperiment` object
#' @param viewpoint viewpoint, defined as a GRanges
#' @param use.scores use.scores
#' @return A tibble with the contact frequency of the viewpoint, per bin 
#'   along the imported genomic range.
#' 
#' @import tibble
#' @importFrom scales unit_format
#' @importFrom S4Vectors queryHits
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps
#' @export
#' @examples 
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' v4C <- virtual4C(contacts_yeast, GenomicRanges::GRanges('II:490000-510000'))
#' v4C

virtual4C <- function(x, viewpoint, use.scores = 'balanced') {
    gis <- InteractionSet::interactions(x)
    gis$score <- HiCExperiment::scores(x, use.scores)
    cm <- cm2matrix(gi2cm(gis))
    regions <- InteractionSet::regions(gis)
    regions_in_viewpoint <- seq_along(regions) %in% S4Vectors::queryHits(
        GenomicRanges::findOverlaps(regions, viewpoint)
    )
    score <- rowSums(cm[, regions_in_viewpoint], na.rm = TRUE)
    GenomicRanges::GRanges(
        seqnames = as.vector(GenomicRanges::seqnames(regions)),
        IRanges::IRanges(
            GenomicRanges::start(regions), 
            GenomicRanges::end(regions)
        ), 
        score = score,
        viewpoint = as.character(viewpoint), 
        center = GenomicRanges::start(regions) + 
            (GenomicRanges::end(regions) - GenomicRanges::start(regions))/2, 
        in_viewpoint = regions_in_viewpoint
    )
}


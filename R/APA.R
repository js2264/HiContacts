#' Aggregated Plot Analysis
#'
#' @param x x
#' @param coords coords
#' @param bins bins
#' @param use.scores use.scores
#' @return an aggregated contact matrix
#'
#' @import InteractionSet
#' @importFrom ggrastr geom_tile_rast
#' @import ggplot2
#' @import tibble
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr arrange
#' @importFrom dplyr left_join
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges resize
#' @importFrom GenomicRanges nearest
#' @importFrom GenomicRanges start
#' @importFrom methods as
#' @importFrom S4Vectors SimpleList

APA <- function(x, coords, nbins = 50, use.scores = 'balanced') {
    `%within%` <- IRanges::`%within%`

    ## -- Resize targets
    expand <- resize(coords, fix = 'center', width = resolution(x) * nbins)
    GenomeInfoDb::seqlevels(expand) <- GenomeInfoDb::seqlevels(x)
    seqinfo(expand) <- seqinfo(x)
    expand <- trim(expand)
    expand <- expand[width(expand) == resolution(x) * nbins]

    ## -- Create aggregated GInteractions obj
    # ints <- fullContactInteractions(
    #     chr = 'aggr', start = -expand/2, end = {expand/2}, 
    #     binning = resolution(x)
    # )

    lapply(seq_along(expand), function(K) {
        gr <- expand[K]
        x[as.character(gr)]
    })
    
    return(x)
}

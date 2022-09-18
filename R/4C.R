#' virtual4C
#'
#' @param x a `contacts` object
#' @param viewpoint viewpoint defined as a GRanges
#' @param use.scores use.scores
#' @return A tibble with the contact frequency of the viewpoint, per bin 
#'   along the imported genomic range.
#' 
#' @import ggplot2
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
#' data(contacts_yeast)
#' virtual4C(contacts_yeast, GenomicRanges::GRanges('II:490000-510000'))

virtual4C <- function(x, viewpoint, use.scores = 'balanced') {
    gis <- scores(x, use.scores)
    cm <- cm2matrix(gi2cm(gis))
    regions <- regions(gis)
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

#' plot4C
#'
#' @param x GRanges, generally the output of `virtual4C()`
#' @param mapping aes to pass on to ggplot2
#' @return ggplot
#'
#' @import ggplot2
#' @import tibble
#' @importFrom scales unit_format
#' @export
#' @examples 
#' library(HiContacts)
#' data(contacts_yeast)
#' v4C <- virtual4C(contacts_yeast, GenomicRanges::GRanges('II:490000-510000'))
#' plot4C(v4C, ggplot2::aes(x = center, y = score))

plot4C <- function(x, mapping) {
    x <- tibble::as_tibble(x)
    p <- ggplot2::ggplot(x, mapping) + 
        ggplot2::geom_line() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::labs(x = "Position", y = "Contacts with viewpoint") + 
        ggplot2::scale_x_continuous(
            expand = c(0, 0), 
            labels = scales::unit_format(unit = "M", scale = 1e-6)
        )
    p
}

#' virtual4C
#'
#' @param x x
#' @param viewpoint viewpoint
#' @param use.assay use.assay
#'
#' @import ggplot2
#' @import tibble
#' @importFrom scales unit_format
#' @importFrom S4Vectors queryHits
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps
#' @export

virtual4C <- function(x, viewpoint, use.assay = 'balanced') {
    gis <- assay(x, use.assay)
    cm <- cm2matrix(gi2cm(gis))
    regions <- regions(gis)
    regions_in_viewpoint <- seq_along(regions) %in% S4Vectors::queryHits(findOverlaps(regions, viewpoint))
    score <- rowSums(cm[, regions_in_viewpoint], na.rm = TRUE)
    tibble::tibble(
        viewpoint = as.character(viewpoint), 
        chr = as.vector(seqnames(regions)), 
        center = GenomicRanges::start(regions) + 
            (GenomicRanges::end(regions) - GenomicRanges::start(regions))/2, 
        score = score
    )
}

#' gg4C
#'
#' @param ... ...
#'
#' @import ggplot2
#' @import tibble
#' @importFrom scales unit_format
#' @export

gg4C <- function(x, mapping) {
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

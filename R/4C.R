virtual4C <- function(x, viewpoint, use.assay = 'balanced') {
    gis <- assay(x, use.assay)
    cm <- cm2matrix(gi2cm(gis))
    regions <- regions(gis)
    regions_in_viewpoint <- seq_along(regions) %in% S4Vectors::queryHits(findOverlaps(regions, viewpoint))
    score <- rowSums(cm[, regions_in_viewpoint], na.rm = TRUE)
    tibble::tibble(
        viewpoint = as.character(viewpoint), 
        chr = as.vector(seqnames(regions)), 
        center = start(regions) + (end(regions) - start(regions))/2, 
        score = score
    )
}

#' gg4C
#'
#' @param ... ...
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @export

gg4C <- function(x, aes) {
    p <- ggplot(x, aes) + 
        geom_line() + 
        theme_minimal() + 
        theme(panel.border = element_rect(fill = NA)) + 
        theme(panel.grid.minor = element_blank()) +
        labs(x = "Position", y = "Contacts with viewpoint")
    p
}

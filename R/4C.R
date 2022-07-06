virtual4C <- function(x, viewpoint, use.assay = 'balanced') {
    gis <- assay(x, use.assay)
    cm <- cm2matrix(gi2cm(gis))
    regions <- regions(gis)
    regions_in_viewpoint <- seq_along(regions) %in% S4Vectors::queryHits(findOverlaps(regions, viewpoint))
    score <- rowSums(cm[, regions_in_viewpoint], na.rm = TRUE)
    regions$viewpoint <- as.character(viewpoint)
    regions$score <- score
    return(regions)
}


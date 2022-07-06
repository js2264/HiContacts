#' APA
#'
#' @param x x
#' @param coords coords
#' @param res res
#' @param limits limits
#' @param dpi dpi
#' @param rasterize rasterize
#' @param symmetrical symmetrical
#' @param BPPARAM BPPARAM
#' @param scale scale
#' @param cmap cmap
#'
#' @import InteractionSet
#' @import ggrastr
#' @import ggplot2
#' @import BiocParallel
#' @import zeallot
#' @import tibble
#' @import dplyr
#' @importFrom GenomicRanges width
#' @export

APA <- function(x, coords, bins = 50, use.assay = 'balanced', BPPARAM = BiocParallel::bpparam()) {

    ## -- Resize coordinates
    expand <- bins*unique(width(coords))
    expanded_coords <- GenomicRanges::resize(coords, fix = 'center', width = expand)
    centered_coords <- GenomicRanges::resize(coords, fix = 'center', width = 1)

    ## -- Create aggregated GInteractions obj
    apa <- fullContactInteractions(chr = 'aggr', start = -expand/2, end = {expand/2}, binning = resolution(x))

    ## -- Get scores ~ x/y, per chromosome
    dists <- BiocParallel::bplapply(BPPARAM = BPPARAM, unique(seqnames(coords)), function(chr) {
        
        ## -- Pre-load the cool file for subcoords (in chromosome `chr`)
        subcoords <- expanded_coords[GenomeInfoDb::seqnames(expanded_coords) == chr]
        subcoords <- IRanges::subsetByOverlaps(subcoords, as(seqinfo(x), 'GRanges'), type = 'within')
        subcoords_centers <- GenomicRanges::resize(subcoords, fix = 'center', width = 1)
        xx <- contacts(path(x), focus = subcoords, resolution = resolution(x))
        gis <- assay(xx, use.assay)
        reg <- regions(gis)
        ans <- anchors(gis)
        
        # -- For each `row`/`column` anchor, get the distance to centered_coords
        nearest_to_row <- GenomicRanges::nearest(ans[['first']], subcoords_centers)
        nearest_to_column <- GenomicRanges::nearest(ans[['second']], subcoords_centers)
        df <- tibble::tibble(
            chr = chr,
            row_center = start(ans[['first']]) + resolution(x)/2 - 1,
            dist_row = row_center - start(subcoords_centers[nearest_to_row]),
            column_center = start(ans[['second']]) + resolution(x)/2 - 1,
            dist_column = column_center - start(subcoords_centers[nearest_to_column]),
            score = gis$score
        ) %>% 
        dplyr::filter(
            abs(dist_row) <= expand/2, 
            abs(dist_column) <= expand/2, 
            dist_row < expand/2, 
            dist_column < expand/2
        ) 
        df$dist_column[df$dist_column < df$dist_row & abs(df$dist_column) > df$dist_row] <- -{df$dist_column[df$dist_column < df$dist_row & abs(df$dist_column) > df$dist_row]}
        df$dist_row[df$dist_column < df$dist_row & df$dist_column < abs(df$dist_row)] <- -{df$dist_row[df$dist_column < df$dist_row & df$dist_column < abs(df$dist_row)]}
        # df$dist_row[df$dist_column < df$dist_row & df$dist_column <= abs(df$dist_row)] <- df$dist_column[df$dist_column < df$dist_row & df$dist_column <= abs(df$dist_row)] * sample(c(-1, 1), 1)
        # df <- df[df$dist_column != df$dist_row, ]
        
        # df <- rbind(df, df %>% dplyr::mutate(x2 = dist_column, dist_column = dist_row, dist_row = x2) %>% dplyr::select(-x2))
        message(glue::glue('Finishing {chr}...'))
        return(df)

    }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::group_by(dist_row, dist_column) %>% 
        dplyr::summarize(score = mean(score, na.rm = TRUE), .groups = 'drop') %>% 
        dplyr::arrange(dist_column) 
    x@interactions <- apa
    x@assays <- S4Vectors::SimpleList(APA = left_join(as_tibble(apa), dists, by = c(start1 = 'dist_row', start2 = 'dist_column'))$score)
    x@features <- c(features(x), S4Vectors::SimpleList(APA = coords))
    return(x)
}

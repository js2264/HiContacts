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
#' @export

APA <- function(x, coords, bins = 50, use.assay = 'balanced', BPPARAM = BiocParallel::bpparam()) {
    `%within%` <- IRanges::`%within%`

    ## -- Resize targets
    expand <- bins*unique(GenomicRanges::width(coords))
    expanded_coords <- GenomicRanges::resize(coords, fix = 'center', width = expand)
    centered_coords <- GenomicRanges::resize(coords, fix = 'center', width = 1)

    ## -- Filter targets
    sub <- expanded_coords %within% GenomicRanges::reduce(regions(x))
    coords <- coords[sub]
    expanded_coords <- expanded_coords[sub]
    centered_coords <- centered_coords[sub]

    ## -- Create aggregated GInteractions obj
    ints <- fullContactInteractions(
        chr = 'aggr', start = -expand/2, end = {expand/2}, 
        binning = resolution(x)
    )

    ## -- Get scores ~ x/y, per chromosome
    dists <- BiocParallel::bplapply(
        BPPARAM = BPPARAM, 
        unique(GenomeInfoDb::seqnames(coords)), 
        function(chr) {
        
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
                row_center = GenomicRanges::start(ans[['first']]) + resolution(x)/2 - 1,
                dist_row = row_center - GenomicRanges::start(subcoords_centers[nearest_to_row]),
                column_center = GenomicRanges::start(ans[['second']]) + resolution(x)/2 - 1,
                dist_column = column_center - GenomicRanges::start(subcoords_centers[nearest_to_column]),
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
            # message(glue::glue('Finishing {chr}...'))
            return(df)
        }
    ) %>% 
        dplyr::bind_rows() %>% 
        dplyr::group_by(dist_row, dist_column) %>% 
        dplyr::summarize(score = mean(score, na.rm = TRUE), .groups = 'drop') %>% 
        dplyr::arrange(dist_column) 
    x@interactions <- ints
    x@assays <- S4Vectors::SimpleList(
        APA = dplyr::left_join(
            tibble::as_tibble(ints), 
            dists, 
            by = c(start1 = 'dist_row', start2 = 'dist_column')
        )$score
    )
    
    x@features <- c(features(x), S4Vectors::SimpleList(APA = coords))
    x@type <- 'aggr.'
    return(x)
}

APA_ <- function(x, coords, bins = 50, use.assay = 'balanced', BPPARAM = BiocParallel::bpparam()) {
    `%within%` <- IRanges::`%within%`
    `%over%` <- IRanges::`%over%`

    ## -- Resize targets
    expand <- bins*unique(GenomicRanges::width(coords))
    expanded_coords <- GenomicRanges::resize(coords, fix = 'center', width = expand)
    centered_coords <- GenomicRanges::resize(coords, fix = 'center', width = 1)

    ## -- Filter targets
    sub <- expanded_coords %within% GenomicRanges::reduce(bins(x))
    coords <- coords[sub]
    expanded_coords <- expanded_coords[sub]
    centered_coords <- centered_coords[sub]

    ## -- Get scores ~ x/y, per chromosome
    gis_aggr <- BiocParallel::bplapply(
        BPPARAM = BPPARAM, 
        unique(GenomeInfoDb::seqnames(coords)), 
        function(chr) {
            ## -- Pre-load the cool file for interactions in coords in chromosome `chr`
            subcoords <- expanded_coords[GenomeInfoDb::seqnames(expanded_coords) == chr]
            subcoords_centers <- GenomicRanges::resize(subcoords, fix = 'center', width = 1)
            xx <- contacts(path(x), focus = subcoords, resolution = resolution(x))
            gis <- assay(xx, use.assay)
            reg <- regions(gis)
            ans <- anchors(gis)
            
            # -- For each subcoord, fetch overlapping interactions
            lapply(seq_along(subcoords), function(K) {
                gr <- subcoords[K]
                mid <- start(subcoords_centers[K])
                sub_first <- ans[['first']] %within% gr
                sub_second <- ans[['second']] %within% gr
                as_tibble(gis[sub_first & sub_second]) %>% 
                    mutate(
                        seqnames1 = 'aggr', 
                        seqnames2 = 'aggr', 
                        start1 = start1 - mid, 
                        end1 = end1 - mid, 
                        start2 = start2 - mid, 
                        end2 = end2 - mid
                    ) %>%
                    as_GInteractions()
            }) %>% 
                do.call(c, .) %>%
                as_tibble() %>% 
                group_by(seqnames1, seqnames2, start1, end1, start2, end2) %>% 
                summarize(score = sum(score, na.rm = TRUE)) %>% 
                as_GInteractions()
        }
    ) %>% do.call(c, .)
    x@interactions <- gis_aggr
    x@assays <- S4Vectors::SimpleList(APA = gis_aggr$score) 
    x@features <- c(features(x), S4Vectors::SimpleList(APA = coords))
    x@type <- 'aggr.'
    return(x)
}


APA2 <- function(x, coords, bins = 50, use.assay = 'balanced', BPPARAM = BiocParallel::bpparam()) {
    `%within%` <- IRanges::`%within%`

    ## -- Resize targets
    expand <- bins*unique(GenomicRanges::width(coords))
    expanded_coords <- GenomicRanges::resize(coords, fix = 'center', width = expand)
    centered_coords <- GenomicRanges::resize(coords, fix = 'center', width = 1)

    ## -- Filter targets
    sub <- expanded_coords %within% GenomicRanges::reduce(regions(x))
    coords <- coords[sub]
    expanded_coords <- expanded_coords[sub]
    centered_coords <- centered_coords[sub]

    ## -- Iterate over each target




    x@interactions <- ints
    x@assays <- S4Vectors::SimpleList(
        APA = dplyr::left_join(
            tibble::as_tibble(ints), 
            dists, 
            by = c(start1 = 'dist_row', start2 = 'dist_column')
        )$score
    )
    
    x@features <- c(features(x), S4Vectors::SimpleList(APA = coords))
    x@type <- 'aggr.'
    return(x)
}

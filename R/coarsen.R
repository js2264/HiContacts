#' @rdname arithmetics
#' @export

coarsen <- function(x, bin.size) {

    res <- resolution(x)
    if (bin.size <= res) stop("Provided `bin.size` <= resolution(x). Please coarsen only to a coarser resolution.")
    gis <- HiCExperiment::interactions(x)
    an <- anchors(gis)
    re <- regions(gis)
    rebinned_re <- GenomicRanges::tileGenome(
        GenomeInfoDb::seqinfo(re), tilewidth = bin.size, 
        cut.last.tile.in.chrom = TRUE
    )
    rebinned_re$chr <- as.vector(GenomicRanges::seqnames(rebinned_re))
    rebinned_re$center <- floor(as.vector(GenomicRanges::start(rebinned_re) + {GenomicRanges::end(rebinned_re) - GenomicRanges::start(rebinned_re)}/2))
    rebinned_re$bin_id <- seq_along(rebinned_re) - 1
    names(rebinned_re) <- paste(
        as.vector(GenomicRanges::seqnames(rebinned_re)), 
        as.vector(GenomicRanges::start(rebinned_re)), 
        as.vector(GenomicRanges::end(rebinned_re)), 
        sep = '_'
    )
    regions(gis)$bin_id <- rebinned_re[GenomicRanges::nearest(
        GenomicRanges::resize(regions(gis), fix = 'center', width = 1), rebinned_re
    )]$bin_id

    # - Coarsen GIs
    score_names <- names(scores(x))
    gis$bin_id1 <- rebinned_re[GenomicRanges::nearest(
        GenomicRanges::resize(an[[1]], fix = 'center', width = 1), 
        rebinned_re
    )]$bin_id
    gis$bin_id2 <- rebinned_re[GenomicRanges::nearest(
        GenomicRanges::resize(an[[2]], fix = 'center', width = 1), 
        rebinned_re
    )]$bin_id
    df <- as.data.frame(S4Vectors::mcols(gis))
    gis2 <- dplyr::select(df, bin_id1, bin_id2) |> 
        dplyr::distinct() |>
        dplyr::left_join(as.data.frame(rebinned_re), by = c('bin_id1' = 'bin_id')) |> 
        dplyr::select(seqnames, start, end, bin_id1, bin_id2) |> 
        dplyr::rename(seqnames1 = seqnames, start1 = start, end1 = end) |> 
        dplyr::left_join(as.data.frame(rebinned_re), by = c('bin_id2' = 'bin_id')) |> 
        dplyr::select(seqnames1, start1, end1, seqnames, start, end, bin_id1, bin_id2) |> 
        dplyr::rename(seqnames2 = seqnames, start2 = start, end2 = end) |> 
        HiCExperiment::df2gi()
    replaceRegions(gis2) <- rebinned_re

    # - Coarsen scores
    scores <- dplyr::group_by(df, bin_id1, bin_id2) |> 
        dplyr::summarise(dplyr::across(all_of(score_names), function(x)sum(x)), .groups = "drop")
    scores <- left_join(as.data.frame(S4Vectors::mcols(gis2)), scores, by = c('bin_id1', 'bin_id2'))

    # - Replace interactions, resolution and scores
    x@resolution <- bin.size
    x@interactions <- gis2
    x@scores <- S4Vectors::SimpleList(as.list(scores[, score_names]))

    return(x)
}
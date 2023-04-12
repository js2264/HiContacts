#' @title Aligning tracks with HiCExperiment objects
#' @name tracks
#' @rdname tracks
#' @description
#' 
#' Aligning tracks with HiCExperiment objects
#'
#' @param x A `HiCExperiment` object over a full genome
#' @param use.pairs logical. Whether to use pairsFile to compute Hi-C coverage
#' @param bin.size if `use.pairs == TRUE`, to which resolution 
#' @return A `HiCExperiment` object with 2 added columns in `regions(x)`
#' 
#' @examples 
#' mcool_file <- HiContactsData::HiContactsData('yeast_wt', format = 'mcool')
#' hic <- import(mcool_file, format = 'mcool', resolution = 1000)
#' coverage(hic)
NULL

#' @param rle a bigwig track imported as 'Rle'
#' @param track.name Name to give to the track in regions metadata
#' @noRd

alignTrack <- function(x, rle, track.name = 'track') {
    re <- regions(interactions(x))
    mc <- S4Vectors::mcols(re)
    mc[[track.name]] <- rle[re]
    mc[[paste0(track.name, '.mean')]] <- BiocGenerics::mean(
        mc[[track.name]], na.rm = TRUE
    )
    S4Vectors::mcols(regions(interactions(x))) <- mc
    return(x)
}

#' @rdname tracks
#' @importFrom GenomicRanges tileGenome
#' @importFrom GenomicRanges countOverlaps
#' @importFrom IRanges coverage
#' @export

setMethod("coverage", signature(x = "HiCExperiment"), function(x, use.pairs = FALSE, bin.size = resolution(x)) {
    pairsFile <- HiCExperiment::pairsFile(x)
    
    ## -- Define tiles on which to count coverage
    if (bin.size > resolution(x)) {
        tiles <- GenomicRanges::tileGenome(
            seqinfo(x), tilewidth = bin.size, cut.last.tile.in.chrom = TRUE
        )
    } 
    else {
        tiles <- bins(x) |> as("GRanges")
    }

    if (use.pairs & !is.null(pairsFile)) { ## Pairs provided: easy countOverlaps
        message("Importing pairs file ", pairsFile, " in memory. This may take a while...")
        pairs <- BiocIO::import(pairsFile, format = 'pairs') 
        an <- anchors(pairs)
        reads <- c(an[[1]], an[[2]])
        tiles$cnt <- GenomicRanges::countOverlaps(tiles, reads)
        tiles$CPM <- tiles$cnt / sum(tiles$cnt) * 1e6
        IRanges::coverage(tiles, weight = 'CPM')
    } 
    else { ## Pairs not provided: recover each ends from binned ints, findOverlaps
        ints <- interactions(x)
        an <- anchors(ints)
        an[[1]]$cnt <- ints$count
        an[[2]]$cnt <- ints$count
        reads <- c(an[[1]], an[[2]])
        df <- dplyr::left_join(
            as.data.frame(reads), 
            as.data.frame(tiles),
            by = dplyr::join_by(seqnames, start, end)
        ) |> dplyr::group_by(seqnames, start, end) |> 
            dplyr::summarize(cnt = sum(cnt), .groups = 'drop') 
        tiles <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
        tiles$CPM <- tiles$cnt / sum(tiles$cnt) * 1e6
        IRanges::coverage(tiles, weight = 'CPM')
    }

})

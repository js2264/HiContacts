#' getAnchors
#'
#' This was adapted from `dovetail-genomics/coolR`
#'
#' @param file file
#' @param res res
#' @param balanced balanced
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom IRanges IRanges
#' @export

getAnchors <- function(file, res = NULL, balanced = "cooler") {
    bins <- fetchCool(file, "bins", res)
    anchors <- GenomicRanges::GRanges(
        bins$chr,
        IRanges::IRanges(bins$start + 1, bins$end),
        seqinfo = cool2seqinfo(file, res)
    )
    names(anchors) <- paste(GenomicRanges::seqnames(anchors), GenomicRanges::start(anchors), GenomicRanges::end(anchors), sep = "_")
    if ("weight" %in% names(bins) & {balanced == "cooler" | balanced == TRUE}) {
        anchors$weight <- bins$weight
    }
    else {
        weight <- 1
    }
    return(anchors)
}

#' getCounts
#'
#' @param file file
#' @param coords coords
#' @param anchors anchors
#' @param coords2 coords2
#' @param res res
#'
#' @import zeallot
#' @import tidyr
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom glue glue
#' @importFrom S4Vectors subjectHits
#' @export

getCounts <- function(file,
                      coords,
                      anchors,
                      coords2 = NULL,
                      res = NULL) {
    ## Process coordinates
    `%<-%` <- zeallot::`%<-%`
    c(coords_chr, coords_start, coords_end) %<-% splitCoords(coords)
    if (any(is.na(coords_start) & !is.na(coords_chr))) {
        coords_start <- rep(1, length(coords_start))
        coords_end <- GenomeInfoDb::seqlengths(anchors)[coords_chr]
    }

    ## Check that queried chr. exists
    if (any(!coords_chr %in% as.vector(GenomicRanges::seqnames(anchors)) & !is.na(coords_chr))) {
        sn <- paste0(unique(as.vector(GenomicRanges::seqnames(anchors))), collapse = ", ")
        stop(glue::glue("{coords_chr} not in file. Available seqnames: {sn}"))
    }

    ## Get chunks to parse ------------------------- This was adapted from `dovetail-genomics/coolR`
    if (is.na(coords_chr)) {
        chunks <- NULL
    } else {
        start <- GenomicRanges::GRanges(coords_chr, IRanges::IRanges(coords_start + 1, width = 1))
        end <- GenomicRanges::GRanges(coords_chr, IRanges::IRanges(coords_end, width = 1))
        start_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(start, anchors))
        end_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(end, anchors))
        bin1_idx <- fetchCool(file, path = "indexes/bin1_offset", res, idx = seq(start_idx, end_idx))
        slice <- sum(bin1_idx[-1] - bin1_idx[-length(bin1_idx)])
        chunks <- seq(
            bin1_idx[1] + 1,
            bin1_idx[1] + 1 + slice
        )
    }

    ## Reading the chunks from the cool file
    df <- tidyr::tibble(
        bin1_id = fetchCool(file, path = "pixels/bin1_id", res, idx = chunks) + 1,
        bin2_id = fetchCool(file, path = "pixels/bin2_id", res, idx = chunks) + 1,
        count = fetchCool(file, path = "pixels/count", res, idx = chunks)
    )

    ## Filter ranges if there is a coords or coords2
    if (is.null(coords2)) {
        coords2 <- coords
    }
    c(coords_chr2, coords_start2, coords_end2) %<-% splitCoords(coords2)
    if (is.na(coords_start2) & !is.na(coords_chr2)) {
        coords_start2 <- 1
        coords_end2 <- GenomeInfoDb::seqlengths(anchors)[coords_chr2]
    }
    if (is.na(coords_chr2)) {
        chunks <- NULL
    } else {
        start <- GenomicRanges::GRanges(coords_chr2, IRanges::IRanges(coords_start2 + 1, width = 1))
        end <- GenomicRanges::GRanges(coords_chr2, IRanges::IRanges(coords_end2, width = 1))
        start_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(start, anchors))
        end_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(end, anchors))
        bin2_idx <- fetchCool(file, path = "indexes/bin1_offset", res, idx = seq(start_idx, end_idx))
        slice <- sum(bin2_idx[-1] - bin2_idx[-length(bin2_idx)])
        chunks <- seq(
            bin2_idx[1] + 1,
            bin2_idx[1] + 1 + slice
        )
        df2 <- tidyr::tibble(
            bin1_id = fetchCool(file, path = "pixels/bin1_id", res, idx = chunks),
            bin2_id = fetchCool(file, path = "pixels/bin2_id", res, idx = chunks),
            count = fetchCool(file, path = "pixels/count", res, idx = chunks)
        )
        df <- df[df$bin2_id %in% df2$bin1_id, ]
    }

    ## Return df
    return(df)
}

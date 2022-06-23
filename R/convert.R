#' cool2seqinfo
#'
#' @param file file
#' @param res res
#'
#' @importFrom GenomeInfoDb Seqinfo
#' @export

cool2seqinfo <- function(file, res = NULL) {
    chroms <- fetchCool(file, "chroms", res)
    seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames = as.vector(chroms$name),
        seqlengths = as.vector(chroms$length)
    )
    return(seqinfo)
}

#' cool2gi
#'
#' @param file file
#' @param coords coords
#' @param coords2 coords2
#' @param res res
#'
#' @import InteractionSet
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges resize
#' @export

cool2gi <- function(file, coords = NULL, coords2 = NULL, res = NULL) {
    anchors <- getAnchors(file, res)
    cnts <- getCounts(file, coords = coords, anchors = anchors, coords2 = coords2, res = res)
    gi <- InteractionSet::GInteractions(
        anchors[cnts$bin1_id],
        anchors[cnts$bin2_id],
        count = cnts$count
    )
    gi$bin1 <- cnts$bin1_id
    gi$bin2 <- cnts$bin2_id

    if (!is.null(gi$anchor1.weight) & !is.null(gi$anchor2.weight)) {
        gi$score <- gi$count * gi$anchor1.weight * gi$anchor2.weight
    } else {
        gi$score <- gi$count
    }
    InteractionSet::regions(gi)$chr <- GenomicRanges::seqnames(InteractionSet::regions(gi))
    InteractionSet::regions(gi)$center <- GenomicRanges::start(GenomicRanges::resize(InteractionSet::regions(gi), fix = "center", width = 1))
    return(gi)
}

#' gi2cm
#'
#' @param gi gi
#'
#' @import InteractionSet
#' @importFrom GenomicRanges mcols
#' @export

gi2cm <- function(gi, fill = "count") {
    InteractionSet::inflate(
        gi,
        rows = 1:length(InteractionSet::regions(gi)),
        columns = 1:length(InteractionSet::regions(gi)),
        fill = GenomicRanges::mcols(gi)[[fill]]
    )
}

cm2matrix <- function(cm, replace_NA = NA) {
    m <- Matrix::as.matrix(cm)
    m[is.na(m)] <- replace_NA
    m
}

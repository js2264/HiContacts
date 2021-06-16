cool2seqinfo <- function(file, res = NULL) {
    chroms <- fetchCool(file, 'chroms', res)
    seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames = as.vector(chroms$name),
        seqlengths = as.vector(chroms$length)
    )
    return(seqinfo)
}

cool2gi <- function(file, balanced = 'cooler', coords = NULL, coords2 = NULL, res = NULL) {
    anchors <- getAnchors(file, res, balanced = balanced)
    cnts <- getCounts(file, coords = coords, anchors = anchors, coords2 = coords2, res = res)
    gi <- InteractionSet::GInteractions(
        anchors[cnts$bin1_id],
        anchors[cnts$bin2_id],
        count = cnts$count
    )
    gi$bin1 <- cnts$bin1_id
    gi$bin2 <- cnts$bin2_id

    if ({!is.null(gi$anchor1.weight) & !is.null(gi$anchor2.weight)} & balanced == 'cooler') {
        gi$score <- gi$count * gi$anchor1.weight * gi$anchor2.weight
    } 
    else if (balanced == 'ICE') {
        gi <- iceGis(gi)
    } 
    else {
        gi$score <- gi$count
    }
    InteractionSet::regions(gi)$chr <- GenomicRanges::seqnames(InteractionSet::regions(gi))
    InteractionSet::regions(gi)$center <- GenomicRanges::start(GenomicRanges::resize(InteractionSet::regions(gi), fix = 'center', width = 1))
    return(gi)
}

gi2cm <- function(gi) {
    InteractionSet::inflate(
        gi, 
        rows = 1:length(InteractionSet::regions(gi)), 
        columns = 1:length(InteractionSet::regions(gi)), 
        fill = ifelse('norm_count' %in% names(GenomicRanges::mcols(gi)), gi$norm_count, gi$count)
    )
}

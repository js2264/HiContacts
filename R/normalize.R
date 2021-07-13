#' @import reticulate
#' @import InteractionSet
#' @import SummarizedExperiment
#' @import dplyr
#' @import tidyr
#' @export

iceGis <- function(gi) {
    icedn <- reticulate::import("iced.normalization")
    cm <- gi2cm(gi)
    m <- cm@matrix
    m_toremove <- order(rowSums(is.na(m)), decreasing = TRUE)[1:{
        nrow(m) * 0.04
    }]
    # m_toremove <- NA
    m_norm <- icedn$ICE_normalization(m[!{
        1:nrow(m) %in% m_toremove
    }, !{
        1:nrow(m) %in% m_toremove
    }])
    new_cm <- cm[!{
        1:nrow(m) %in% m_toremove
    }, !{
        1:nrow(m) %in% m_toremove
    }]
    InteractionSet::as.matrix(new_cm) <- m_norm
    is <- InteractionSet::deflate(new_cm)
    new_gis <- InteractionSet::GInteractions(
        InteractionSet::anchors(is)[[1]],
        InteractionSet::anchors(is)[[2]],
        score = SummarizedExperiment::assay(is)[, 1]
    )
    new_gis$bin1 <- tidyr::as_tibble(new_gis) %>%
        dplyr::left_join(
            tidyr::as_tibble(gi),
            by = c("seqnames1", "start1", "end1", "width1", "strand1", "seqnames2", "start2", "end2", "width2", "strand2")
        ) %>%
        dplyr::pull(bin1)
    new_gis$bin2 <- tidyr::as_tibble(new_gis) %>%
        dplyr::left_join(
            tidyr::as_tibble(gi),
            by = c("seqnames1", "start1", "end1", "width1", "strand1", "seqnames2", "start2", "end2", "width2", "strand2")
        ) %>%
        dplyr::pull(bin2)
    return(new_gis)
}

#' @import dplyr
#' @export

normalizeOverExpected <- function(mat) {

    # Get the binning resolution
    binsize <- sort(unique(mat$y - mat$x))[2]

    # Get the number of bins in the matrix
    nbins <-
        {
            range(mat$y - mat$x) / binsize
        }[2] + 1

    expected <- mat %>%
        dplyr::mutate(
            diag = (mat$y - mat$x) / binsize,
            pixelsPerDiag = nbins - diag
        ) %>%
        dplyr::group_by(diag) %>%
        dplyr::summarize(
            expected = sum(score, na.rm = TRUE) / pixelsPerDiag,
            .groups = "keep"
        ) %>%
        dplyr::distinct()

    mat <- mat %>%
        dplyr::mutate(diag = (mat$y - mat$x) / binsize) %>%
        dplyr::left_join(expected, by = "diag") %>%
        dplyr::mutate(scoreOverExpected = -log2(score / expected)) ## `-` because both `score` and `expected` have to be < 0
}

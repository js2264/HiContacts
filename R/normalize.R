#' iceGis
#'
#' @param gi gi
#'
#' @import reticulate
#' @import InteractionSet
#' @import dplyr
#' @import tidyr
#' @importFrom SummarizedExperiment assay
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

#' normalizeOverExpected
#'
#' @param mat mat
#'
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

#' correlateMatrix
#'
#' @param mat mat
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import corrr
#' @export

correlateMatrix <- function(mat) {
    x <- mat %>%
        dplyr::select(x, y, score) %>%
        distinct() %>%
        tidyr::pivot_wider(names_from = y, values_from = score) %>%
        tibble::column_to_rownames("x")
    x <- x[rownames(x), rownames(x)]
    # x <- x[colSums(!is.na(x)) > 0, colSums(!is.na(x)) > 0]
    # x[lower.tri(x)] <- t(x)[lower.tri(x)]

    co <- corrr::correlate(x, diagonal = 0, method = "pearson", quiet = TRUE)
    mat2 <- co %>%
        tidyr::pivot_longer(-term, names_to = "y", values_to = "corr") %>%
        dplyr::rename("x" = "term") %>%
        mutate(
            x = as.numeric(x),
            y = as.numeric(y)
        )
    mat %>%
        left_join(mat2) %>%
        mutate(score = corr)
}

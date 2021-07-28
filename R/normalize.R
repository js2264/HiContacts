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
#' @importFrom scales rescale
#' @importFrom tibble tibble
#' @import dplyr
#' @export

normalizeOverExpected <- function(mat) {

    # Get the binning resolution
    binsize <- sort(unique(mat$y - mat$x))[2]

    # Get the number of bins in the matrix
    mBin <- min(c(mat$bin2, mat$bin1))
    MBin <- max(c(mat$bin2, mat$bin1))
    bins <- seq(mBin, MBin)
    nbins <- length(unique(abs(bins)))
    bins_df <- tibble::tibble(
        bin1 = rep(bins, each = length(bins)),
        bin2 = rep(bins, length(bins))
    ) %>% dplyr::filter(bin2 >= bin1)

    df <- left_join(bins_df, mat) %>%
        dplyr::mutate(
            seqnames1 = unique(mat$seqnames1),
            start1 = bin1 * binsize - (mat$bin1[1] * binsize - mat$start1[1]),
            end1 = bin1 * binsize - (mat$bin1[1] * binsize - mat$end1[1]),
            width1 = binsize,
            strand1 = "*",
            chr1 = unique(mat$chr1),
            center1 = end1 - (end1 - start1) / 2,
            seqnames2 = unique(mat$seqnames2),
            start2 = bin2 * binsize - (mat$bin2[1] * binsize - mat$start2[1]),
            end2 = bin2 * binsize - (mat$bin2[1] * binsize - mat$end2[1]),
            width2 = binsize,
            strand2 = "*",
            chr2 = unique(mat$chr2),
            center2 = end2 - (end2 - start2) / 2,
            count = 0,
            anchor1.weight = 0,
            anchor2.weight = 0,
            x = center1,
            y = center2,
            diag = (y - x) / binsize,
            pixelsPerDiag = nbins - diag
        )

    expected <- df %>%
        dplyr::group_by(diag) %>%
        dplyr::summarize(
            expected = sum(score, na.rm = TRUE) / pixelsPerDiag,
            .groups = "keep"
        ) %>%
        dplyr::distinct()
    # expected$expected <- scales::rescale(expected$expected, to = range(mat$score, na.rm = TRUE))

    res <- df %>%
        dplyr::left_join(expected, by = "diag") %>%
        dplyr::mutate(scoreOverExpected = score - expected) ## `-` because both `score` and `expected` have to be < 0

    return(res)
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

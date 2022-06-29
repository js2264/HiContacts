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

#' detrend
#'
#' @param gis gis
#'
#' @importFrom scales rescale
#' @importFrom tibble tibble
#' @import dplyr
#' @export

detrend <- function(gis, update_scores = FALSE) {
    gis$diag <- abs(gis$bin2 - gis$bin1)
    expected <- as_tibble(gis) %>% 
        group_by(diag) %>% 
        summarize(average_interaction_per_diag = mean(score, na.rm = TRUE)) %>% 
        mutate(average_interaction_per_diag = average_interaction_per_diag / 2)
    gis$expected <- as_tibble(gis) %>% left_join(expected, by = 'diag') %>% pull(average_interaction_per_diag)
    gis$score_over_expected <- log2(gis$score / gis$expected)
    if (update_scores) {
        gis$score <- gis$score_over_expected
    }
    return(gis)
}

smoothen <- function(gis, use_serpentine_trend = TRUE, serpentine_niter = 100L, serpentine_ncores = 10L) {
 
    options(reticulate.repl.quiet = TRUE)
    reticulate::use_condaenv('tm')
    sp <- reticulate::import('serpentine')

    ## Check that only 1 chromosome is present in the gis object
    seqnames <- unique(as.vector(GenomicRanges::seqnames(InteractionSet::anchors(gis)[[1]])))
    if (length(seqnames) > 1) {
        stop('Smoothing maps across multiple chromosomes is not supported. Aborting now.')
    }
    binsize <- width(regions(gis)[1])[1]

    ## Run Serpentine
    B <- cm2matrix(gi2cm(gis, fill = 'score'), replace_NA = 0)
    A <- matrix(data = 1, nrow = nrow(B), ncol = ncol(B))
    c(sm1, sm2, sK) %<-% sp$serpentin_binning(A, B, verbose = FALSE, iterations = serpentine_niter, parallel = serpentine_ncores)
    
    ## Re-center smoothened matrix 
    if (use_serpentine_trend) {
        c(trend, threshold) %<-% sp$MDbefore(m1, m2, show = FALSE)
        sK <- sK - trend
    }
    else {
        sK <- sK - mean(sK, na.rm = TRUE)
    }

    ## Make a full-featured interactions (storing smoothed scores in `score`)
    gis_smoothened <- sK %>%
        as_tibble() %>% 
        setNames(start(anchors(gi2cm(gis, fill = 'score'))$row)) %>%
        mutate(start1 = start(anchors(gi2cm(gis, fill = 'score'))$row)) %>% 
        pivot_longer(-start1, names_to = 'start2', values_to = 'score') %>% 
        mutate(start2 = as.numeric(start2))
    an1 <- GenomicRanges::GRanges(
        seqnames = seqnames, 
        IRanges::IRanges(gis_smoothened$start1, width = binsize)
    )
    an2 <- GenomicRanges::GRanges(
        seqnames = seqnames, 
        IRanges::IRanges(gis_smoothened$start2, width = binsize)
    )
    reg <- unique(c(an1, an2))
    bins <- as_tibble(gis) %>% 
        select(seqnames1, start1, end1, bin1) %>% 
        distinct() %>% 
        GenomicRanges::makeGRangesFromDataFrame(seqnames.field = 'seqnames1', start.field = 'start1', end.field = 'end1', keep.extra.columns = TRUE)
    gi <- InteractionSet::GInteractions(
        anchor1 = an1, 
        anchor2 = an2, 
        regions = reg, 
        count = NA, 
        anchor1.weight = NA, 
        anchor2.weight = NA
    )
    gi$bin1 <- plyranges::join_overlap_left(anchors(gi)[[1]], bins)$bin1
    gi$bin2 <- plyranges::join_overlap_left(anchors(gi)[[2]], bins)$bin1
    gi$score <- gis_smoothened$score
    return(gi)
}

autocorrelate <- function(gis) {

}

divide <- function(x, by, serpentine = FALSE, ...) {
    `%>%` <- tidyr::`%>%`
    `%<-%` <- zeallot::`%<-%`
    
    ## -- If regions are different, manually merge them 
    InteractionSet::replaceRegions(x) <- unique(c(InteractionSet::regions(x), InteractionSet::regions(by)))
    InteractionSet::replaceRegions(by) <- unique(c(InteractionSet::regions(x), InteractionSet::regions(by)))

    ## -- Convert to matrices 
    m1 <- cm2matrix(gi2cm(x, fill = 'count'), replace_NA = 0)
    m2 <- cm2matrix(gi2cm(by, fill = 'count'), replace_NA = 0)
    binsize <- width(regions(x)[1])[1]

    ## Compute ratio
    if (serpentine) {
        ## -- Run serpentine
        options(reticulate.repl.quiet = TRUE)
        reticulate::use_condaenv('tm')
        sp <- reticulate::import('serpentine')
        c(trend, threshold) %<-% sp$MDbefore(m1, m2, show = FALSE)
        c(sm1, sm2, sK) %<-% sp$serpentin_binning(m1, m2, threshold = threshold, minthreshold = threshold/5, verbose = FALSE, iterations = serpentine_niter, parallel = serpentine_ncores)
        sK <- sK - trend
    }
    else {
        sK <- m2/m1
    }

    ## Make a full-featured interactions (storing divided scores in `score`)
    seqnames <- unique(seqnames(regions(x)))
    mat <- sK %>%
        as_tibble() %>% 
        setNames(start(anchors(gi2cm(x, fill = 'score'))$row)) %>%
        mutate(start2 = start(anchors(gi2cm(by, fill = 'score'))$row)) %>% 
        pivot_longer(-start2, names_to = 'start1', values_to = 'score') %>% 
        mutate(start1 = as.numeric(start1))
    an1 <- GenomicRanges::GRanges(
        seqnames = seqnames, 
        IRanges::IRanges(mat$start1, width = binsize)
    )
    an2 <- GenomicRanges::GRanges(
        seqnames = seqnames, 
        IRanges::IRanges(mat$start2, width = binsize)
    )
    reg <- unique(c(an1, an2))
    bins <- as_tibble(c(x, by)) %>% 
        select(seqnames1, start1, end1, bin1) %>% 
        distinct() %>% 
        GenomicRanges::makeGRangesFromDataFrame(seqnames.field = 'seqnames1', start.field = 'start1', end.field = 'end1', keep.extra.columns = TRUE)
    mat <- sK %>%
        as_tibble() %>% 
        setNames(start(anchors(gi2cm(x, fill = 'score'))$row)) %>%
        mutate(start1 = start(anchors(gi2cm(by, fill = 'score'))$row)) %>% 
        pivot_longer(-start1, names_to = 'start2', values_to = 'score') %>% 
        mutate(start2 = as.numeric(start2)) %>%
        dplyr::mutate(
            end1 = start1 + binsize, 
            end2 = start2 + binsize
        )
    gi <- InteractionSet::GInteractions(
        anchor1 = an1, 
        anchor2 = an2, 
        regions = reg, 
        count = NA, 
        anchor1.weight = NA, 
        anchor2.weight = NA
    )
    gi$bin1 <- plyranges::join_overlap_left(anchors(gi)[[1]], bins)$bin1
    gi$bin2 <- plyranges::join_overlap_left(anchors(gi)[[2]], bins)$bin1
    gi$gis1_v_gis2 <- mat$score

    ## -- Filter ratio
    gi <- gi[!is.na(gi$gis1_v_gis2) & is.finite(gi$gis1_v_gis2)]
    return(gi)

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

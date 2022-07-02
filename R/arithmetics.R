#' detrend
#'
#' @param gis gis
#'
#' @importFrom scales rescale
#' @importFrom tibble tibble
#' @import dplyr
#' @export

detrend <- function(x, use.assay = 'balanced') {
    gis <- assay(x, use.assay)
    gis$diag <- pairdist(gis) / resolution(x)
    expected <- as_tibble(gis) %>% 
        group_by(diag) %>% 
        summarize(average_interaction_per_diag = mean(score, na.rm = TRUE)) %>% 
        mutate(average_interaction_per_diag = average_interaction_per_diag / 2)
    gis$expected <- as_tibble(gis) %>% left_join(expected, by = 'diag') %>% pull(average_interaction_per_diag)
    gis$score_over_expected <- log2(gis$score / gis$expected)
    x@assays[['expected']] <- gis$expected
    x@assays[['detrended']] <- gis$score_over_expected
    return(x)
}

smoothen <- function(x, use.assay = 'balanced', use_serpentine_trend = TRUE, serpentine_niter = 10L, serpentine_ncores = 16L) {
 
    # options(reticulate.repl.quiet = TRUE)
    # reticulate::use_condaenv('tm')
    sp <- reticulate::import('serpentine')
    gis <- assay(x, use.assay)

    ## Check that only 1 chromosome is present in the gis object
    seqnames <- unique(as.vector(GenomicRanges::seqnames(InteractionSet::anchors(gis)[[1]])))
    if (length(seqnames) > 1) {
        stop('Smoothing maps across multiple chromosomes is not supported. Aborting now.')
    }
    binsize <- resolution(x)

    ## Run Serpentine
    B <- cm2matrix(gi2cm(gis), replace_NA = 0)
    A <- matrix(data = 1, nrow = nrow(B), ncol = ncol(B))
    c(sm1, sm2, sK) %<-% sp$serpentin_binning(A, B, verbose = FALSE, iterations = serpentine_niter, parallel = serpentine_ncores)
    
    ## Re-center smoothened matrix 
    if (use_serpentine_trend) {
        c(trend, threshold) %<-% sp$MDbefore(B, A, show = FALSE)
        if (is.na(trend)) trend <- mean(sK, na.rm = TRUE)
        sK <- sK - trend
    }
    else {
        sK <- sK - mean(sK, na.rm = TRUE)
    }

    ## Make a full-featured interactions (storing smoothed scores in `score`)
    gis_smoothened <- sK %>%
        as_tibble() %>% 
        setNames(start(anchors(gi2cm(gis))$row)) %>%
        mutate(start1 = start(anchors(gi2cm(gis))$row)) %>% 
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
    gi <- InteractionSet::GInteractions(
        anchor1 = an1, 
        anchor2 = an2, 
        regions = reg
    )
    x@interactions <- gi
    x@assays <- S4Vectors::SimpleList(smoothen = gis_smoothened$score)
    return(x)
}

autocorrelate <- function(x, use.assay = 'balanced') {
    gis <- assay(x, use.assay)
    reg <- regions(gis)
    mat <- cm2matrix(gi2cm(gis))
    co <- corrr::correlate(log10(mat), diagonal = 0, method = "pearson", quiet = TRUE)
    colnames(co) <- c('term', names(reg))
    co$term <- names(reg)
    mat2 <- co %>%
        tidyr::pivot_longer(-term, names_to = "y", values_to = "corr") %>%
        dplyr::rename("x" = "term")
    gis2 <- InteractionSet::GInteractions(
        anchor1 = stringr::str_replace(mat2$x, '_', ':') %>% stringr::str_replace('_', '-') %>% as('GRanges'), 
        anchor2 = stringr::str_replace(mat2$y, '_', ':') %>% stringr::str_replace('_', '-') %>% as('GRanges'), 
        score = as.array(mat2$corr)
    )
    x@interactions <- gis2
    x@assays <- S4Vectors::SimpleList(autocorrelation = gis2$score)
    return(x)
}

divide <- function(x, by) {
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

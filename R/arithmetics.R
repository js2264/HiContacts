#' detrend
#'
#' @param gis gis
#'
#' @importFrom scales rescale
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom InteractionSet pairdist
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @import GenomicRanges
#' @export

detrend <- function(x, use.assay = 'balanced') {
    gis <- assay(x, use.assay)
    gis$diag <- InteractionSet::pairdist(gis) / resolution(x)
    expected <- tibble::as_tibble(gis) %>% 
        dplyr::group_by(diag) %>% 
        dplyr::summarize(average_interaction_per_diag = mean(score, na.rm = TRUE)) %>% 
        dplyr::mutate(average_interaction_per_diag = average_interaction_per_diag / 2)
    gis$expected <- tibble::as_tibble(gis) %>% 
        dplyr::left_join(expected, by = 'diag') %>% 
        dplyr::pull(average_interaction_per_diag)
    gis$score_over_expected <- log2(gis$score / gis$expected)
    x@assays[['expected']] <- gis$expected
    x@assays[['detrended']] <- gis$score_over_expected
    return(x)
}

#' smoothen
#'
#' @import reticulate
#' @import GenomicRanges
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom InteractionSet anchors
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors SimpleList
#' @export

smoothen <- function(x, use.assay = 'balanced', use_serpentine_trend = TRUE, serpentine_niter = 10L, serpentine_ncores = 16L) {

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
        tibble::as_tibble() %>% 
        stats::setNames(GenomicRanges::start(anchors(gi2cm(gis))$row)) %>%
        dplyr::mutate(start1 = GenomicRanges::start(anchors(gi2cm(gis))$row)) %>% 
        tidyr::pivot_longer(-start1, names_to = 'start2', values_to = 'score') %>% 
        dplyr::mutate(start2 = as.numeric(start2))
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
    x@type <- 'smoothed'
    return(x)
}

#' autocorrelate
#'
#' @import InteractionSet
#' @import stringr
#' @importFrom tidyr pivot_longer
#' @importFrom corrr correlate
#' @importFrom dplyr rename
#' @importFrom S4Vectors SimpleList
#' @export

autocorrelate <- function(x, use.assay = 'balanced', ignore_ndiags = 3) {
    gis <- assay(x, use.assay)
    reg <- regions(gis)
    mat <- cm2matrix(gi2cm(gis))
    sdiag(mat, 0) <- NA
    for (K in seq(-ignore_ndiags, ignore_ndiags, by = 1)) {
        sdiag(mat, K) <- NA
    }
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
    x@type <- 'autocorr.'
    return(x)
}

#' divide
#'
#' @import tidyr
#' @import zeallot
#' @import reticulate
#' @import plyranges
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom GenomicRanges seqnames
#' @importFrom InteractionSet regions
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors SimpleList
#' @export

divide <- function(x, by, use.assay = 'balanced') {
    `%>%` <- tidyr::`%>%`
    `%<-%` <- zeallot::`%<-%`
    
    ## -- Check that all objects are comparable (bins, regions, resolution, seqinfo)
    is_comparable(x, by)

    x_gis <- assay(x, use.assay)
    by_gis <- assay(by, use.assay)

    ## -- If regions are different, manually merge them 
    InteractionSet::replaceRegions(x_gis) <- unique(
        c(InteractionSet::regions(x_gis), InteractionSet::regions(by_gis))
    )
    InteractionSet::replaceRegions(by_gis) <- unique(
        c(InteractionSet::regions(x_gis), InteractionSet::regions(by_gis))
    )

    ## -- Convert to matrices 
    m1 <- cm2matrix(gi2cm(x_gis), replace_NA = 0)
    m2 <- cm2matrix(gi2cm(by_gis), replace_NA = 0)
    binsize <- resolution(x)

    serpentine <- FALSE
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
        sK <- m1/m2
    }

    ## Make a full-featured interactions (storing divided scores in `score`)
    seqnames <- unique(GenomicRanges::seqnames(regions(x)))
    mat <- sK %>%
        tibble::as_tibble() %>% 
        setNames(GenomicRanges::start(anchors(gi2cm(x_gis))$row)) %>%
        dplyr::mutate(start2 = GenomicRanges::start(anchors(gi2cm(by_gis))$row)) %>% 
        tidyr::pivot_longer(-start2, names_to = 'start1', values_to = 'score') %>% 
        dplyr::mutate(start1 = as.numeric(start1)) %>%
        dplyr::mutate(
            end1 = start1 + binsize, 
            end2 = start2 + binsize
        )
    an1 <- GenomicRanges::GRanges(
        seqnames = seqnames, 
        IRanges::IRanges(mat$start1, width = binsize)
    )
    an2 <- GenomicRanges::GRanges(
        seqnames = seqnames, 
        IRanges::IRanges(mat$start2, width = binsize)
    )
    reg <- unique(c(an1, an2))
    bins <- c(bins(x), bins(by)) %>% 
        unique() 
    gi <- InteractionSet::GInteractions(
        anchor1 = an1, 
        anchor2 = an2, 
        regions = reg, 
        count = NA, 
        anchor1.weight = NA, 
        anchor2.weight = NA
    )
    gi$gis1_v_gis2 <- mat$score

    ## -- Filter ratio
    gi <- gi[!is.na(mat$score) & is.finite(mat$score)]

    ## -- Create 'in silico' contacts
    path <- paste0(
        basename(metadata(x)$path), ' / ', basename(metadata(by)$path)
    )
    res <- methods::new("contacts", 
        focus = focus(x), 
        metadata = list(
            path = path, x_path = metadata(x)$path, by_path = metadata(by)$path
        ), 
        seqinfo = seqinfo(x), 
        resolutions = binsize, 
        current_resolution = binsize, 
        bins = bins, 
        interactions = gi, 
        assays = S4Vectors::SimpleList(
            'ratio' = mat$score[!is.na(mat$score) & is.finite(mat$score)]
        ), 
        features = S4Vectors::SimpleList(), 
        pairsFile = NULL, 
        type = 'ratio'
    )
    return(res)

}

#' merge
#'
#' @import tidyr
#' @import zeallot
#' @import reticulate
#' @import plyranges
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom GenomicRanges seqnames
#' @importFrom InteractionSet regions
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors SimpleList
#' @export

merge <- function(..., use.assay = 'balanced') {
    `%>%` <- tidyr::`%>%`
    `%<-%` <- zeallot::`%<-%`
    contacts_list <- list(...)
    
    ## -- Check that all objects are comparable (bins, regions, resolution, seqinfo)
    is_comparable(...)

    # Unify all the interactions
    ints <- do.call(
        c, lapply(contacts_list, FUN = interactions) 
    ) %>% 
        unique() %>% 
        sort()
    
    # Set all scores for each assay to 0
    asss <- lapply(names(assays(contacts_list[[1]])), function(name) {
        rep(0, length(ints))
    })
    names(asss) <- names(assays(contacts_list[[1]]))
    asss <- S4Vectors::SimpleList(asss)

    ## -- Iterate over each contacts in `contacts_list`
    for (idx in seq_along(contacts_list)) {
        sub <- S4Vectors::subjectHits(
            GenomicRanges::findOverlaps(
                interactions(contacts_list[[idx]]),
                ints
            )
        )
        sub <- seq_along(ints) %in% sub
        for (K in seq_along(asss)) {
            vals <- contacts_list[[idx]]@assays[[K]]
            vals[is.na(vals)] <- 0
            asss[[K]][sub] <- asss[[K]][sub] + 
                vals
        }
    }

    ## -- Create 'in silico' contacts
    files <- paste0(
        basename(unlist(lapply(contacts_list, path))), 
        collapse = ' + '
    )
    res <- methods::new("contacts", 
        focus = focus(x), 
        metadata = list(
            path = '',
            merging = files, 
            operation = 'sum'
        ), 
        seqinfo = seqinfo(contacts_list[[1]]), 
        resolutions = resolutions(contacts_list[[1]]), 
        current_resolution = resolution(contacts_list[[1]]), 
        bins = bins(contacts_list[[1]]), 
        interactions = ints, 
        assays = asss, 
        features = S4Vectors::SimpleList(), 
        pairsFile = NULL, 
        type = 'merged'
    )
    return(res)

}

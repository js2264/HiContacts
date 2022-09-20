#' Arithmetics with (m)cool file(s)
#' 
#' Different operations can be performed:  
#'  - Detrending a contact matrix, i.e. removing the distance-dependent 
#' contact trend;
#'  - Autocorrelate a contact matrix: this is is typically done to highlight 
#' large-scale compartments;
#'  - Divide one contact matrix by another; 
#'  - Merge multiple contact matrices;
#'  - Serpentinify, or smooth a contact matrix out. This requires `serpentine` 
#' python package to be installed.
#' 
#' @rdname arithmetics
#'
#' @param x a `contacts` object
#' @param use.scores use.scores
#' @return a `contacts` object with two additional scoress: `expected` and
#'   `detrended`
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
#' @examples 
#' #### -----
#' #### Detrending a contact matrix
#' #### -----
#' 
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' contacts_yeast <- detrend(contacts_yeast)
#' scores(contacts_yeast)

detrend <- function(x, use.scores = 'balanced') {
    gis <- scores(x, use.scores)
    gis$diag <- InteractionSet::pairdist(gis) / resolution(x)
    expected <- tibble::as_tibble(gis) |> 
        dplyr::group_by(diag) |> 
        dplyr::summarize(average_interaction_per_diag = mean(score, na.rm = TRUE)) |> 
        dplyr::mutate(average_interaction_per_diag = average_interaction_per_diag / 2)
    gis$expected <- tibble::as_tibble(gis) |> 
        dplyr::left_join(expected, by = 'diag') |> 
        dplyr::pull(average_interaction_per_diag)
    gis$score_over_expected <- log2(gis$score / gis$expected)
    scores(x, "expected") <- gis$expected
    scores(x, "detrended") <- gis$score_over_expected
    return(x)
}

#' @rdname arithmetics
#'
#' @param x a `contacts` object
#' @param use.scores use.scores
#' @param ignore_ndiags ignore N diagonals when calculating correlations
#' @return a `contacts` object with a single `autocorrelation` scores
#' 
#' @import InteractionSet
#' @import stringr
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom dplyr rename
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom S4Vectors SimpleList
#' @export
#' @examples 
#' #### -----
#' #### Auto-correlate a contact matrix
#' #### -----
#' 
#' contacts_yeast <- autocorrelate(contacts_yeast)
#' scores(contacts_yeast)
#' plotMatrix(contacts_yeast, scale = 'linear', limits = c(-1, 1), cmap = bwrColors())

autocorrelate <- function(x, use.scores = 'balanced', ignore_ndiags = 3) {
    gis <- scores(x, use.scores)
    reg <- regions(gis)
    mat <- cm2matrix(gi2cm(gis))
    sdiag(mat, 0) <- NA
    for (K in seq(-ignore_ndiags, ignore_ndiags, by = 1)) {
        sdiag(mat, K) <- NA
    }
    # co <- corrr::correlate(log10(mat), diagonal = 0, method = "pearson", quiet = TRUE)
    co <- stats::cor(log10(mat), use = 'pairwise.complete.obs', method = "pearson")
    co <- tibble::as_tibble(co) |> 
        dplyr::mutate(term = paste0('V', seq_len(dplyr::n()))) |> 
        dplyr::relocate(term)
    colnames(co) <- c('term', names(reg))
    co$term <- names(reg)
    mat2 <- co |>
        tidyr::pivot_longer(-term, names_to = "y", values_to = "corr") |>
        dplyr::rename("x" = "term")
    gis2 <- InteractionSet::GInteractions(
        anchor1 = stringr::str_replace(mat2$x, '_', ':') |> stringr::str_replace('_', '-') |> as('GRanges'), 
        anchor2 = stringr::str_replace(mat2$y, '_', ':') |> stringr::str_replace('_', '-') |> as('GRanges'), 
        score = as.array(mat2$corr)
    )

    res <- methods::new("contacts", 
        focus = focus(x), 
        metadata = metadata(x),
        resolutions = resolutions(x), 
        resolution = resolution(x), 
        interactions = gis2, 
        scores = S4Vectors::SimpleList(autocorrelation = gis2$score),
        topologicalFeatures = topologicalFeatures(x), 
        pairsFile = pairsFile(x), 
        fileName = fileName(x)
    )

    return(res)
}

#' @rdname arithmetics
#'
#' @param x a `contacts` object
#' @param by a `contacts` object
#' @param use.scores use.scores
#' @return a `contacts` object with a single `ratio` scores
#'
#' @import tidyr
#' @import reticulate
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom GenomicRanges seqnames
#' @importFrom InteractionSet regions
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors SimpleList
#' @export
#' @examples 
#' #### -----
#' #### Divide 2 contact matrices
#' #### -----
#' 
#' contacts_yeast <- contacts_yeast()
#' contacts_yeast_eco1 <- contacts_yeast_eco1()
#' div_contacts <- divide(contacts_yeast_eco1, by = contacts_yeast)
#' div_contacts
#' plotMatrix(div_contacts, scale = 'log2', limits = c(-2, 2), cmap = bwrColors())

divide <- function(x, by, use.scores = 'balanced') {
    
    ## -- Check that all objects are comparable (bins, regions, resolution, seqinfo)
    is_comparable(x, by)

    x_gis <- scores(x, use.scores)
    by_gis <- scores(by, use.scores)

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
        sp <- reticulate::import('serpentine')
        .v <- sp$serpentin_binning(m1, m2, threshold = threshold, minthreshold = threshold/5, verbose = FALSE, iterations = serpentine_niter, parallel = serpentine_ncores)
        sm1 <- .v[[1]]
        sm2 <- .v[[2]]
        sK <- .v[[3]]
        .v <- sp$MDbefore(m1, m2, show = FALSE)
        trend <- .v[[1]]
        threshold <- .v[[2]]
        sK <- sK - trend
    }
    else {
        sK <- m1/m2
    }

    ## Make a full-featured interactions (storing divided scores in `score`)
    seqnames <- unique(GenomicRanges::seqnames(regions(x)))
    mat <- sK |>
        tibble::as_tibble() |> 
        setNames(GenomicRanges::start(anchors(gi2cm(x_gis))$row)) |>
        dplyr::mutate(start2 = GenomicRanges::start(anchors(gi2cm(by_gis))$row)) |> 
        tidyr::pivot_longer(-start2, names_to = 'start1', values_to = 'score') |> 
        dplyr::mutate(start1 = as.numeric(start1)) |>
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
    bins <- c(bins(x), bins(by)) |> 
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
    res <- methods::new("contacts", 
        focus = focus(x), 
        metadata = list(),
        resolutions = binsize, 
        resolution = binsize, 
        interactions = gi, 
        scores = S4Vectors::SimpleList(
            'ratio' = mat$score[!is.na(mat$score) & is.finite(mat$score)]
        ), 
        topologicalFeatures = S4Vectors::SimpleList(), 
        pairsFile = NULL, 
        fileName = paste0(basename(fileName(x)), ' / ', basename(fileName(by)))
    )
    return(res)

}

#' @rdname arithmetics
#'
#' @param ... `contacts` objects
#' @param use.scores use.scores
#' @return a `contacts` object. Each returned scores is the sum of the
#'   corresponding scores from input `contacts`.
#'
#' @import tidyr
#' @import reticulate
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom GenomicRanges seqnames
#' @importFrom InteractionSet regions
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors SimpleList
#' @export
#' @examples 
#' #### -----
#' #### Merge 2 contact matrices
#' #### -----
#' 
#' merged_contacts <- merge(contacts_yeast_eco1, contacts_yeast)
#' merged_contacts

merge <- function(..., use.scores = 'balanced') {
    contacts_list <- list(...)
    
    ## -- Check that at least 2 `contacts` objects are passed to `merge()`
    if (length(contacts_list) < 2) {
        stop("Please provide at least 2 `contacts` objects.")
    } 
    else if (!all(unlist(lapply(contacts_list, is, 'contacts')))) {
        stop("Please only provide `contacts` objects.")
    }

    ## -- Check that all objects are comparable (bins, regions, resolution, seqinfo)
    is_comparable(...)

    # Unify all the interactions
    ints <- do.call(
        c, lapply(contacts_list, FUN = interactions) 
    ) |> 
        unique() |> 
        sort()
    
    # Set all scores for each scores to 0
    asss <- lapply(names(scores(contacts_list[[1]])), function(name) {
        rep(0, length(ints))
    })
    names(asss) <- names(scores(contacts_list[[1]]))
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
            vals <- scores(contacts_list[[idx]])[[K]]
            vals[is.na(vals)] <- 0
            asss[[K]][sub] <- asss[[K]][sub] + 
                vals
        }
    }

    ## -- Create 'in silico' contacts
    files <- paste0(
        basename(unlist(lapply(contacts_list, fileName))), 
        collapse = ', '
    )
    res <- methods::new("contacts", 
        focus = paste0(
            basename(unlist(lapply(contacts_list, focus))), 
            collapse = ', '
        ), 
        metadata = list(
            merging = files, 
            operation = 'sum'
        ), 
        resolutions = resolutions(contacts_list[[1]]), 
        resolution = resolution(contacts_list[[1]]), 
        interactions = ints, 
        scores = asss, 
        topologicalFeatures = S4Vectors::SimpleList(), 
        pairsFile = NULL, 
        fileName = ""
    )
    return(res)

}

#' @rdname arithmetics
#'
#' @param x a `contacts` object
#' @param use.scores use.scores
#' @param use_serpentine_trend whether to use the trend estimated with 
#'   serpentine (this requires `reticulate` and the python package `serpentine`)
#' @param serpentine_niter number of iterations to use for serpentine
#' @param serpentine_ncores number of CPUs to use for serpentine
#' @return a `contacts` object with a single `smoothen` scores
#' 
#' @import reticulate
#' @import GenomicRanges
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom InteractionSet anchors
#' @importFrom InteractionSet GInteractions
#' @importFrom S4Vectors SimpleList

serpentinify <- function(x, use.scores = 'balanced', 
    use_serpentine_trend = TRUE, serpentine_niter = 10L, serpentine_ncores = 16L
) {

    sp <- reticulate::import('serpentine')
    gis <- scores(x, use.scores)

    ## Check that only 1 chromosome is present in the gis object
    seqnames <- unique(as.vector(GenomicRanges::seqnames(InteractionSet::anchors(gis)[['first']])))
    if (length(seqnames) > 1) {
        stop('Smoothing maps across multiple chromosomes is not supported. Aborting now.')
    }
    binsize <- resolution(x)

    ## Run Serpentine
    B <- cm2matrix(gi2cm(gis), replace_NA = 0)
    A <- matrix(data = 1, nrow = nrow(B), ncol = ncol(B))
    .v <- sp$serpentin_binning(A, B, verbose = FALSE, iterations = serpentine_niter, parallel = serpentine_ncores)
    sm1 <- .v[[1]]
    sm2 <- .v[[2]]
    sK <- .v[[3]]
    
    ## Re-center smoothened matrix 
    if (use_serpentine_trend) {
        .v <- sp$MDbefore(B, A, show = FALSE)
        trend <- .v[[1]]
        threshold <- .v[[2]]
        if (is.na(trend)) trend <- mean(sK, na.rm = TRUE)
        sK <- sK - trend
    }
    else {
        sK <- sK - mean(sK, na.rm = TRUE)
    }

    ## Make a full-featured interactions (storing smoothed scores in `score`)
    gis_smoothened <- sK |>
        tibble::as_tibble() |> 
        stats::setNames(GenomicRanges::start(anchors(gi2cm(gis))$row)) |>
        dplyr::mutate(start1 = GenomicRanges::start(anchors(gi2cm(gis))$row)) |> 
        tidyr::pivot_longer(-start1, names_to = 'start2', values_to = 'score') |> 
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

    res <- methods::new("contacts", 
        focus = focus(x), 
        metadata = metadata(x),
        resolutions = resolutions(x), 
        resolution = resolution(x), 
        interactions = gi, 
        scores = S4Vectors::SimpleList(smoothen = gis_smoothened$score),
        topologicalFeatures = topologicalFeatures(x), 
        pairsFile = pairsFile(x), 
        fileName = fileName(x)
    )

    return(res)
}


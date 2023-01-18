#' HiContacts arithmetics functionalities
#' 
#' @name arithmetics
#' @aliases detrend
#' @aliases autocorrelate
#' @aliases divide
#' @aliases merge
#' @aliases despeckle
#' @aliases aggregate,HiCExperiment-method
#' @aliases boost
#' @aliases subsample
#' @aliases normalize,HiCExperiment-method
#' 
#' @description 
#' Different arithmetic operations can be performed on Hi-C contact matrices:  
#'  - `normalize` a contact matrix using iterative correction;
#'  - `detrend` a contact matrix, i.e. remove the distance-dependent 
#' contact trend;
#'  - `autocorrelate` a contact matrix: this is is typically done to highlight 
#' large-scale compartments;
#'  - `divide` one contact matrix by another; 
#'  - `merge` multiple contact matrices;
#'  - `despeckle` (i.e. smooth out) a contact matrix out;
#'  - `aggregate` (average) a contact matrices over a set of genomic loci of 
#' interest;
#'  - `boost` Hi-C signal by enhancing long-range interactions while preserving short-
#' range interactions (this is adapted from Boost-HiC);
#'  - `subsample` interactions using a proportion or a fixed number of final 
#' interactions.
#' 
#' @param x,object a `HiCExperiment` object
#' @param use.scores Which scores to use to perform operations
#' @param ... `HiCExperiment` objects. For `aggregate`, `targets` (a set of 
#' GRanges or GInteractions).
#' @param detrend Detrend matrix before performing autocorrelation
#' @param ignore_ndiags ignore N diagonals when calculating correlations
#' @param by a `HiCExperiment` object
#' @param focal.size Size of the smoothing rectangle
#' @param alpha Power law scaling factor. As indicated in Boost-HiC documentation, 
#' the alpha parameter influences the weighting of contacts: if alpha < 1 
#' long-range interactions are prioritized; if alpha >> 1 short-range 
#' interactions have more weight when computing the distance matrix.
#' @param full.replace Whether to replace the entire set of contacts, 
#' rather than only filling the missing interactions (count=0) (Default: FALSE)
#' @param prop Float between 0 and 1, or integer corresponding to the # of 
#' @param niters Number of iterations for ICE matrix balancing
#' @param min.nnz Filter bins with less than `min.nnz` non-zero elements when 
#' performing ICE matrix balancing
#' @param mad.max Filter out bins whose log coverage is less than `mad.max` 
#' median absolute deviations below the median bin log coverage.
#' 
#' @return a `HiCExperiment` object with extra scores
#' 
#' @import InteractionSet
#' @import stringr
#' @import tidyr
#' @importMethodsFrom BiocGenerics normalize
#' @importFrom scales rescale
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr pull
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr rename
#' @importFrom dplyr n
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors aggregate
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocParallel bpparam
#' 
#' @examples 
#' #### -----
#' #### Normalize a contact matrix
#' #### -----
#' 
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' normalize(contacts_yeast)
#'
#' #### -----
#' #### Detrending a contact matrix
#' #### -----
#' 
#' detrend(contacts_yeast)
#'
#' #### -----
#' #### Auto-correlate a contact matrix
#' #### -----
#' 
#' autocorrelate(contacts_yeast)
#' 
#' #### -----
#' #### Divide 2 contact matrices
#' #### -----
#' 
#' contacts_yeast <- refocus(contacts_yeast, 'II')
#' contacts_yeast_eco1 <- contacts_yeast_eco1() |> refocus('II')
#' divide(contacts_yeast_eco1, by = contacts_yeast)
#' 
#' #### -----
#' #### Merge 2 contact matrices
#' #### -----
#' 
#' merge(contacts_yeast_eco1, contacts_yeast)
#' 
#' #### -----
#' #### Despeckle (smoothen) a contact map
#' #### -----
#' 
#' despeckle(contacts_yeast)
#' 
#' #### -----
#' #### Aggregate a contact matrix over centromeres, at different scales
#' #### -----
#' 
#' contacts <- contacts_yeast() |> zoom(resolution = 1000)
#' centros <- topologicalFeatures(contacts, 'centromeres')
#' aggregate(contacts, targets = centros, flanking_bins = 50)
#' 
#' #### -----
#' #### Enhance long-range interaction signal
#' #### -----
#' 
#' contacts <- contacts_yeast() |> zoom(resolution = 1000) |> refocus('II')
#' boost(contacts, alpha = 1)
#' 
#' #### -----
#' #### Subsample interactions 
#' #### -----
#' 
#' subsample(contacts, prop = 100000)
NULL

#' @rdname arithmetics
#' @export

detrend <- function(x, use.scores = 'balanced') {
    gis <- interactions(x)
    gis$score <- scores(x, use.scores)
    gis$diag <- gis$bin_id2 - gis$bin_id1
    expected <- tibble::as_tibble(gis) |> 
        dplyr::group_by(diag) |> 
        dplyr::summarize(average_interaction_per_diag = mean(score, na.rm = TRUE)) |> 
        dplyr::mutate(average_interaction_per_diag = average_interaction_per_diag / 2)
    gis$expected <- tibble::as_tibble(gis) |> 
        dplyr::left_join(expected, by = 'diag') |> 
        dplyr::pull(average_interaction_per_diag)
    gis$score_over_expected <- log2( {gis$score/sum(gis$score, na.rm = TRUE)} / {gis$expected/sum(gis$expected, na.rm = TRUE)} )
    scores(x, "expected") <- gis$expected
    scores(x, "detrended") <- gis$score_over_expected
    return(x)
}

#' @rdname arithmetics
#' @export

autocorrelate <- function(x, use.scores = 'balanced', detrend = TRUE, ignore_ndiags = 3) {
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
        message("Install WGCNA package to perform Boost-HiC.")
        message("install.packages('WGCNA')")
    }
    if (detrend) {
        if (!{'detrend' %in% names(scores(x))}) {
            x <- detrend(x, use.scores = use.scores)
        }
        use.scores <- 'detrended'
    }
    gis <- interactions(x)
    gis$score <- scores(x, use.scores)
    reg <- regions(gis)
    mat <- cm2matrix(gi2cm(gis))
    for (K in seq(-ignore_ndiags, ignore_ndiags, by = 1)) {
        sdiag(mat, K) <- NA
    }
    co <- WGCNA::cor(mat, use = 'pairwise.complete.obs', method = "pearson")
    colnames(co) <- names(regions(gis))
    rownames(co) <- names(regions(gis))
    mat2 <- tibble::as_tibble(co, rownames = 'region', .name_repair = "minimal") |> 
        tidyr::pivot_longer(-region, names_to = "y", values_to = "corr") |>
        dplyr::rename("x" = "region")
    mat2$bin_id1 <- dplyr::left_join(
        mat2, 
        data.frame(region = names(regions(gis)), id = regions(gis)$bin_id), 
        by = c(x = 'region')
    ) |> pull(id)
    mat2$bin_id2 <- dplyr::left_join(
        mat2, 
        data.frame(region = names(regions(gis)), id = regions(gis)$bin_id), 
        by = c(y = 'region')
    ) |> pull(id)
    scores(x, 'autocorrelated') <- left_join(
        as_tibble(gis), mat2, by = c("bin_id1", "bin_id2")
    ) |> pull(corr)
    return(x)
}

#' @rdname arithmetics
#' @export

divide <- function(x, by, use.scores = 'balanced') {
    
    ## -- Check that all objects are comparable (bins, regions, resolution, seqinfo)
    is_comparable(x, by)

    x_gis <- interactions(x)
    x_gis$score <- scores(x, use.scores)
    by_gis <- interactions(by)
    by_gis$score <- scores(by, use.scores)

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
    colnames(sK) <- paste0('x', GenomicRanges::start(anchors(gi2cm(x_gis))$row))
    rownames(sK) <- paste0('x', GenomicRanges::start(anchors(gi2cm(x_gis))$row))

    ## Make a full-featured interactions (storing divided scores in `score`)
    seqnames <- unique(GenomicRanges::seqnames(regions(x)))
    mat <- sK |>
        base::as.matrix() |> 
        tibble::as_tibble(.name_repair = "universal") |> 
        dplyr::mutate(start2 = GenomicRanges::start(anchors(gi2cm(by_gis))$row)) |> 
        tidyr::pivot_longer(-start2, names_to = 'start1', values_to = 'score') |> 
        dplyr::mutate(start1 = gsub('x', '', start1) |> as.numeric()) |>
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
    gi <- InteractionSet::GInteractions(
        anchor1 = an1, 
        anchor2 = an2, 
        regions = reg, 
        count = NA, 
        anchor1.weight = NA, 
        anchor2.weight = NA
    )
    gi$gis1_v_gis2 <- mat$score

    ## -- Create HiCExperiment
    res <- methods::new("HiCExperiment", 
        focus = HiCExperiment::focus(x), 
        metadata = S4Vectors::metadata(x),
        resolutions = binsize, 
        resolution = binsize, 
        interactions = gi[!is.na(mat$score) & is.finite(mat$score)], 
        scores = S4Vectors::SimpleList(
            'ratio' = mat$score[!is.na(mat$score) & is.finite(mat$score)]
        ), 
        topologicalFeatures = HiCExperiment::topologicalFeatures(x), 
        pairsFile = NULL, 
        fileName = paste0(basename(fileName(x)), ' / ', basename(fileName(by)))
    )
    return(res)

}

#' @rdname arithmetics
#' @export

merge <- function(..., use.scores = 'balanced') {
    contacts_list <- list(...)
    
    ## -- Check that at least 2 `HiCExperiment` objects are passed to `merge()`
    are_HiCExperiment(...)
    if (length(contacts_list) < 2) {
        stop("Please provide at least 2 `HiCExperiment` objects.")
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

    ## -- Iterate over each Contacts in `contacts_list`
    for (idx in seq_along(contacts_list)) {
        sub <- S4Vectors::subjectHits(
            GenomicRanges::findOverlaps(
                interactions(contacts_list[[idx]]),
                ints
            )
        )
        sub <- seq_along(ints) %in% sub
        for (K in seq_along(asss)) {
            vals <- scores(contacts_list[[idx]], K)
            vals[is.na(vals)] <- 0
            asss[[K]][sub] <- asss[[K]][sub] + 
                vals
        }
    }

    ## -- Create 'in silico' Contacts
    files <- paste0(
        basename(unlist(lapply(contacts_list, fileName))), 
        collapse = ', '
    )
    res <- methods::new("HiCExperiment", 
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
#' @export

despeckle <- function(x, use.scores = 'balanced', focal.size = 1) {
    if (!requireNamespace("terra", quietly = TRUE)) {
        message("Install terra package to despeckle contact matrix.")
        message("install.packages('terra')")
    }

    gis <- HiCExperiment::interactions(x)
    gis$score <- HiCExperiment::scores(x, use.scores)
    cm <- HiCExperiment::gi2cm(gis)
    an <- InteractionSet::anchors(cm)
    mat <- HiCExperiment::cm2matrix(cm)
    r <- terra::rast(mat)
    gauss <- matrix(data = 1/{(focal.size*2+1)*(focal.size*2+1)-1}, ncol = focal.size*2+1, nrow = focal.size*2+1)
    gauss[focal.size+1, focal.size+1] <- 0
    mat_despeckled <- terra::focal(
        x = r, 
        w = gauss, 
        fun = sum,
        na.rm = TRUE
    )#/{(focal.size*2+1)*(focal.size*2+1)-1}
    mat_despeckled <- terra::as.array(mat_despeckled)[, , 1]
    cm_despeckled <- InteractionSet::ContactMatrix(
        mat_despeckled, 
        anchor1 = an[[1]], 
        anchor2 = an[[2]], 
        regions = InteractionSet::regions(cm)
    )
    is_despeckled <- InteractionSet::deflate(cm_despeckled)
    gis_despeckled <- InteractionSet::interactions(is_despeckled)
    gis_despeckled$bin_id1 <- HiCExperiment::anchors(gis_despeckled)[[1]]$bin_id
    gis_despeckled$bin_id2 <- HiCExperiment::anchors(gis_despeckled)[[2]]$bin_id
    HiCExperiment::interactions(x) <- gis_despeckled
    m <- dplyr::left_join(
        as.data.frame(mcols(gis_despeckled)), 
        as.data.frame(mcols(gis)), 
        by = c('bin_id1', 'bin_id2')
    ) |> dplyr::select(-bin_id1, -bin_id2, -score)
    m$despeckled <- SummarizedExperiment::assay(is_despeckled, 1)[, 1]
    l <- as.list(m) |> S4Vectors::SimpleList()
    x@scores <- l
    return(x)
}


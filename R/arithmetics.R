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
#' @param targets Set of chromosome coordinates for which 
#'   interaction counts are extracted from the Hi-C contact file, provided
#'   as a GRanges object (for diagnoal-centered loci) or as a GInteractions
#'   object (for off-diagonal coordinates).
#' @param flankingBins Number of bins on each flank of the bins containing 
#'   input targets.
#' @param maxDistance Maximum distance to use when compiling distance decay
#' @param BPPARAM BiocParallel parameters
#' 
#' @return a `HiCExperiment` object with extra scores
#' 
#' @import InteractionSet
#' @import stringr
#' @import tidyr
#' @importFrom BiocGenerics normalize
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
#' aggregate(contacts, targets = centros, flankingBins = 51)
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
    if (is.null(metadata(x)[['detrending_model']])) {
        expected <- tibble::as_tibble(gis) |> 
            dplyr::group_by(diag) |> 
            dplyr::summarize(average_interaction_per_diag = mean(score, na.rm = TRUE)) |> 
            dplyr::mutate(average_interaction_per_diag = average_interaction_per_diag / 2)
    }
    else {
        expected <- metadata(x)[['detrending_model']]
        expected$average_interaction_per_diag <- expected$score
        expected$distance <- NULL
        expected$score <- NULL
    }
    gis$expected <- tibble::as_tibble(gis) |> 
        dplyr::left_join(expected, by = 'diag') |> 
        dplyr::pull(average_interaction_per_diag)
    gis$score_over_expected <- log2( {gis$score/sum(gis$score, na.rm = TRUE)} / {gis$expected/sum(gis$expected, na.rm = TRUE)} )
    # gis$score_over_expected <- log2( {scale(gis$score)} / {scale(gis$expected)} )[,1]
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
    colnames(co) <- names(reg)
    rownames(co) <- names(reg)
    mat2 <- tibble::as_tibble(co, rownames = 'region', .name_repair = "minimal") |> 
        tidyr::pivot_longer(-region, names_to = "y", values_to = "corr") |>
        dplyr::rename("x" = "region")
    mat2$bin_id1 <- dplyr::left_join(
        mat2, 
        data.frame(region = names(reg), id = reg$bin_id), 
        by = c(x = 'region')
    ) |> pull(id)
    mat2$bin_id2 <- dplyr::left_join(
        mat2, 
        data.frame(region = names(reg), id = reg$bin_id), 
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

    x_gis <- interactions(x)
    by_gis <- interactions(by)

    ## -- If regions are different, manually merge them 
    re <- unique(
        c(InteractionSet::regions(x_gis), InteractionSet::regions(by_gis))
    )
    InteractionSet::replaceRegions(x_gis) <- re
    InteractionSet::replaceRegions(by_gis) <- re
    
    ## -- Check that all objects are comparable (bins, resolution, seqinfo)
    is_same_seqinfo(x, by)
    is_same_resolution(x, by)
    is_same_bins(x, by)

    ## -- Convert to matrices 
    m1 <- cm2matrix(gi2cm(x_gis, use.scores), replace_NA = 0)
    m2 <- cm2matrix(gi2cm(by_gis, use.scores), replace_NA = 0)
    binsize <- resolution(x)

    # - Compute ratio matrix
    sK <- m1/m2

    # - Make a full-featured interactions (storing divided scores in `score`)
    gis <- dplyr::full_join(as_tibble(x_gis), as_tibble(by_gis), 
        by = c("seqnames1", "start1", "end1", "width1", "strand1", "chr1", 
        "center1", "bin_id1", "seqnames2", "start2", "end2", "width2",
        "strand2", "chr2", "center2", "bin_id2")
    ) |> dplyr::select(!ends_with(c('1.1.x', '1.1.y', '2.1.x', '2.1.y')))
    mat_ratio <- sK
    cm_ratio <- InteractionSet::ContactMatrix(
        mat_ratio, 
        anchor1 = re, 
        anchor2 = re, 
        regions = re
    )
    is_ratio <- InteractionSet::deflate(cm_ratio)
    gis_ratio <- InteractionSet::interactions(is_ratio)
    gis_ratio$bin_id1 <- HiCExperiment::anchors(gis_ratio)[[1]]$bin_id
    gis_ratio$bin_id2 <- HiCExperiment::anchors(gis_ratio)[[2]]$bin_id
    m <- dplyr::left_join(
        as.data.frame(mcols(gis_ratio)), 
        gis, 
        by = c('bin_id1', 'bin_id2')
    ) |> dplyr::select(ends_with(c('.x', '.y', 'ratio')))
    m <- dplyr::rename(m, count = 'count.x', balanced = 'balanced.x')
    m$ratio <- SummarizedExperiment::assay(is_ratio, 1)[, 1]
    scores <- as.list(m) |> S4Vectors::SimpleList()

    ## -- Create HiCExperiment
    res <- methods::new("HiCExperiment", 
        focus = HiCExperiment::focus(x), 
        metadata = S4Vectors::metadata(x),
        resolutions = binsize, 
        resolution = binsize, 
        interactions = gis_ratio, 
        scores = scores, 
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
    is_same_seqinfo(...)
    is_same_resolution(...)
    is_same_bins(...)

    # Unify all the interactions
    score_names <- names(scores(contacts_list[[1]]))
    ints <- do.call(
        c, lapply(contacts_list, FUN = interactions) 
    ) |> sort()
    ints_df <- as.data.frame(ints)
    merged_ints <- dplyr::select(ints_df, !any_of(score_names)) |> 
        dplyr::distinct() |> 
        HiCExperiment:::asGInteractions()

    # Group by bin_id1/bin_id2 and merge scores
    FUN_mean <- function(x) mean(x, na.rm = TRUE)
    FUN_sum <- function(x) sum(x, na.rm = TRUE)
    asss <- dplyr::select(ints_df, dplyr::any_of(c('bin_id1', 'bin_id2', score_names))) |> 
        dplyr::group_by(bin_id1, bin_id2) |>
        dplyr::summarise(dplyr::across(score_names, FUN_mean), .groups = "drop") |> 
        dplyr::select(-bin_id1, -bin_id2) |> 
        as("SimpleList")

    ## -- Create a new HiCExperiment object
    files <- paste0(
        basename(unlist(lapply(contacts_list, fileName))), 
        collapse = ', '
    )
    res <- methods::new("HiCExperiment", 
        focus = paste0(
            unlist(lapply(contacts_list, focus)), 
            collapse = ', '
        ), 
        metadata = list(
            merging = files, 
            operation = 'sum'
        ), 
        resolutions = resolutions(contacts_list[[1]]), 
        resolution = resolution(contacts_list[[1]]), 
        interactions = merged_ints, 
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


#' @importFrom stats rbinom
#' @rdname arithmetics
#' @export

subsample <- function(x, prop) {

    gis <- HiCExperiment::interactions(x)
    gis$score <- HiCExperiment::scores(x, 'count')
    gis$scaling <- gis$balanced/gis$score
    cm <- HiCExperiment::gi2cm(gis)
    an <- InteractionSet::anchors(cm)
    mat <- HiCExperiment::cm2matrix(cm)
    scores_names <- names(scores(x))
    if (prop > 1) prop <- prop/sum(mat, na.rm = TRUE)

    subsample_counts <- apply(mat, 2, function(x) stats::rbinom(nrow(mat), x, prop))
    subsample_counts[subsample_counts == 0] <- NA
    subsample_counts[lower.tri(subsample_counts)] <- NA
    # subsample_counts <- 2 * subsample_counts
    cm_subsampled <- InteractionSet::ContactMatrix(
        subsample_counts, 
        anchor1 = an[[1]], 
        anchor2 = an[[2]], 
        regions = InteractionSet::regions(cm)
    )
    is_subsampled <- InteractionSet::deflate(cm_subsampled)
    gis_subsampled <- InteractionSet::interactions(is_subsampled)
    gis_subsampled$bin_id1 <- HiCExperiment::anchors(gis_subsampled)[[1]]$bin_id
    gis_subsampled$bin_id2 <- HiCExperiment::anchors(gis_subsampled)[[2]]$bin_id
    HiCExperiment::interactions(x) <- gis_subsampled
    m <- dplyr::left_join(
        as.data.frame(GenomicRanges::mcols(gis_subsampled)), 
        as.data.frame(GenomicRanges::mcols(gis)), 
        by = c('bin_id1', 'bin_id2')
    ) |> dplyr::select(-bin_id1, -bin_id2, -score)
    n <- data.frame(
        count = SummarizedExperiment::assay(is_subsampled, 1)[, 1]
    )
    n$balanced <- n$count * m$scaling
    n$balanced <- scales::rescale(n$balanced, to = c(
        min(gis$balanced, na.rm = TRUE),
        max(gis$balanced, na.rm = TRUE)
    ))
    l <- as.list(n) |> S4Vectors::SimpleList()
    x@scores <- l

    return(x)

}

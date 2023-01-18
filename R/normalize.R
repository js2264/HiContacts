#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom stats median
#' @rdname arithmetics
#' @export

setMethod("normalize", signature(object = "HiCExperiment"), function(
    object, use.scores = 'count', niters = 200, min.nnz = 10, mad.max = 3
) {

    x <- object
    gis <- HiCExperiment::interactions(x)
    gis$score <- HiCExperiment::scores(x, use.scores)
    cm <- HiCExperiment::gi2cm(gis)
    an <- InteractionSet::anchors(cm)
    re <- InteractionSet::regions(cm)
    mat0 <- HiCExperiment::cm2matrix(cm, sparse = TRUE)
    mat0 <- Matrix::triu(mat0)
    mat <- mat0

    # - Filter bins 
    bins <- Matrix::colSums(mat)
    bins[bins == 0] <- 1
    lbins <- log10(bins)
    med <- stats::median(lbins)
    sig <- 1.4826 * stats::mad(lbins)
    s_min <- med - mad.max * sig
    filtered_bins <- {lbins > s_min} & 
        {Matrix::colSums(mat > 0, na.rm = TRUE) > min.nnz}
    invalid_bins <- re[which(!filtered_bins)]$bin_id

    # - Iterate over max_iter
    mat_iced <- mat[filtered_bins, filtered_bins]
    nonzeros <- data.frame(col = mat_iced@i+1, row = mat_iced@j+1) |> as.matrix()
    pb <- utils::txtProgressBar(min = 0, max = niters, style = 3, width = 50, char = "-") 
    for (i in seq_len(niters)) {
        utils::setTxtProgressBar(pb, i)
        bin_sums <- Matrix::colSums(mat_iced)
        bin_sums <- bin_sums / stats::median(bin_sums)
        mat_iced@x <- mat_iced@x / {bin_sums[nonzeros[, 'row']] * bin_sums[nonzeros[, 'col']]}
    }

    # - Rescale ice-d scores
    bin_sums <- Matrix::colSums(mat_iced)
    mat_iced@x <- mat_iced@x / stats::median(bin_sums)

    # - Make upper-tri
    mat_iced <- Matrix::triu(mat_iced)

    # - Recover scores 
    cm_iced <- InteractionSet::ContactMatrix(
        mat_iced, 
        anchor1 = an[[1]][filtered_bins],
        anchor2 = an[[2]][filtered_bins],
        regions = re[filtered_bins]
    )
    is_iced <- InteractionSet::deflate(cm_iced)
    gis_iced <- InteractionSet::interactions(is_iced)
    gis_iced$bin_id1 <- HiCExperiment::anchors(gis_iced)[[1]]$bin_id
    gis_iced$bin_id2 <- HiCExperiment::anchors(gis_iced)[[2]]$bin_id
    gis_iced$ICE <- SummarizedExperiment::assay(is_iced, 1)[, 1]
    m <- dplyr::left_join(
        as.data.frame(mcols(gis)), 
        as.data.frame(mcols(gis_iced)), 
        by = c('bin_id1', 'bin_id2')
    )
    m$ICE[{m$bin_id1 %in% invalid_bins} | {m$bin_id2 %in% invalid_bins}] <- NA
    m <- dplyr::select(m, -bin_id1, -bin_id2, -score)
    # - Save new scores
    x@scores <- as.list(m) |> S4Vectors::SimpleList()

    return(x)
})
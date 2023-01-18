#' @rdname arithmetics
#' @export

boost <- function(x, use.scores = 'balanced', alpha = 1, full.replace = FALSE) {

    if (!requireNamespace("Rfast", quietly = TRUE)) {
        message("Install Rfast package to perform Boost-HiC.")
        message("install.packages('Rfast')")
    }
    gis <- HiCExperiment::interactions(x)
    gis$score <- HiCExperiment::scores(x, use.scores)
    cm <- HiCExperiment::gi2cm(gis)
    an <- InteractionSet::anchors(cm)
    mat <- HiCExperiment::cm2matrix(cm)
    
    d0 <- (1/(mat^(1/alpha)))
    distances <- Rfast::floyd(d0)
    F <- (1/(distances^(1/alpha)))

    if (!full.replace) {
        mat_boosted <- mat
        mat_boosted[is.na(mat) | mat == 0] <- F[is.na(mat) | mat == 0]
    }
    else {
        mat_boosted <- F
    }

    cm_boosted <- InteractionSet::ContactMatrix(
        mat_boosted, 
        anchor1 = an[[1]], 
        anchor2 = an[[2]], 
        regions = InteractionSet::regions(cm)
    )
    is_boosted <- InteractionSet::deflate(cm_boosted)
    gis_boosted <- InteractionSet::interactions(is_boosted)
    gis_boosted$bin_id1 <- HiCExperiment::anchors(gis_boosted)[[1]]$bin_id
    gis_boosted$bin_id2 <- HiCExperiment::anchors(gis_boosted)[[2]]$bin_id
    HiCExperiment::interactions(x) <- gis_boosted
    m <- dplyr::left_join(
        as.data.frame(mcols(gis_boosted)), 
        as.data.frame(mcols(gis)), 
        by = c('bin_id1', 'bin_id2')
    ) |> dplyr::select(-bin_id1, -bin_id2, -score)
    m$boosted <- SummarizedExperiment::assay(is_boosted, 1)[, 1]
    l <- as.list(m) |> S4Vectors::SimpleList()
    x@scores <- l
    metadata(x)[['alpha']] <- alpha

    return(x)
}



getCompartments <- function(x, genome, neigens = 3, autocorrelation = FALSE, BPPARAM = BiocParallel::bpparam()) {
    
    BPPARAM <- BiocParallel::SerialParam()
    
    ## -- Compute eigens on each chr. separately
    compts <- BiocParallel::bplapply(
        BPPARAM = BPPARAM, 
        seqnames(seqinfo(x)), 
        function(chr) {
            x_chr <- x[chr]
            if (length(interactions(x_chr)) <= neigens) {
                gr <- regions(x_chr)
                if (length(gr)) {
                    for (k in seq_len(neigens)) {
                        GenomicRanges::mcols(gr)[[paste0('E', k)]] <- 0
                    }
                    gr$eigen <- 0
                    gr$GC <- 0
                }
                else {
                    for (k in seq_len(neigens)) {
                        GenomicRanges::mcols(gr)[[paste0('E', k)]] <- NULL
                    }
                    gr$eigen <- NULL
                    gr$GC <- NULL
                }
                return(gr)
            } else {
                .eigChr(
                    x_chr, 
                    genome, 
                    neigens, 
                    autocorrelation = autocorrelation,
                    ignore_diags = 2, 
                    clip_percentile = 0.999
                )
            }
        }
    ) |> GRangesList() |> unlist()

    ## -- Return initial HiCExperiment object with eigens in metadata
    A <- GenomicRanges::reduce(compts[compts$eigen > 0])
    A$compartment <- 'A'
    B <- GenomicRanges::reduce(compts[compts$eigen < 0])
    B$compartment <- 'B'
    cpts <- sort(c(A, B))
    topologicalFeatures(x, 'compartments') <- cpts
    metadata(x)$eigens <- compts
    return(x)

}

#' @importFrom RSpectra eigs_sym

.eigChr <- function(x_chr, genome, neigens, autocorrelation, clip_percentile, ignore_diags) {

    ## -- Get dense detrended matrix (with non-log2-detrended scores)
    if (!autocorrelation) {
        dx <- detrend(x_chr)
        gis <- InteractionSet::interactions(dx)
        gis$score <- HiCExperiment::scores(dx, 'detrended')
        gr <- regions(gis)
        m <- cm2matrix(gi2cm(gis))
        m <- m - colMeans(m, na.rm = TRUE)
    }
    else {
        dx <- autocorrelate(x_chr, ignore_ndiags = ignore_diags)
        gis <- InteractionSet::interactions(dx)
        gis$score <- HiCExperiment::scores(dx, 'autocorrelated')
        gr <- regions(gis)
        m <- cm2matrix(gi2cm(gis))
    }

    ## -- Remove ignored diagonals and NAs
    for (K in seq(-ignore_diags, ignore_diags, by = 1)) {
        sdiag(m, K) <- 0
    }
    m[is.na(m)] <- 0

    ## -- Mask white lines from matrix
    unmasked_bins <- rowSums(m) != 0
    m_no0 <- m[unmasked_bins, unmasked_bins]

    ## -- Clamp values at extremes
    m_no0 <- scales::oob_squish(
        m_no0, 
        quantile(m_no0, probs = c(
            1-clip_percentile, clip_percentile
        ), na.rm = TRUE)
    )

    ## -- Get eigen values and vectors 
    eigen <- RSpectra::eigs_sym(m_no0, k = neigens)
    
    ## -- Normalize eigenvector
    eigs <- eigen$vectors
    eigs <- eigs / apply(eigs, 2, function(x) as.numeric(sqrt(sum(x * x))))
    eigs <- lapply(seq_along(eigen$values), function(K) {
        eigs[, K] * as.numeric(sqrt(abs(eigen$values[K])))
    }) 
    names(eigs) <- paste0('E', seq_along(eigen$values))
    eigs <- dplyr::bind_cols(eigs)
    
    ## -- Recover full length eigenvector
    eigs_final <- as.data.frame(
        matrix(nrow = nrow(m), ncol = neigens, data = 0)
    )
    for (k in seq_len(neigens)) {
        eigs_final[unmasked_bins, k] <- eigs[, k]
        GenomicRanges::mcols(gr)[, paste0('E', k)] <- eigs_final[, k]
    }

    ## -- Get GC content
    gr$GC <- as.numeric(Biostrings::letterFrequency(
        genome[gr], letters = 'GC') / GenomicRanges::width(gr))

    ## -- Get matching eigenvector and phase it with GC
    gr <- .eigGCPhasing(gr, neigens)

    return(gr)

}

.eigGCPhasing <- function(gr, neigens) {
    cors <- lapply(seq_len(neigens), function(K) {
        cor(gr$GC, GenomicRanges::mcols(gr)[, paste0('E', K)])
    }) |> unlist()
    best_cor <- which.max(abs(cors))
    if (cors[best_cor] < 0) {
        message(glue::glue(
            "seqnames `{seqnames(gr)[1]}`: eigen #{best_cor} selected (flipped)"
        ))
        gr$eigen <- -mcols(gr)[, paste0('E', best_cor)]
    }
    else {
        message(glue::glue(
            "seqnames `{seqnames(gr)[1]}`: eigen #{best_cor} selected"
        ))
        gr$eigen <- mcols(gr)[, paste0('E', best_cor)]
    }
    return(gr)
}

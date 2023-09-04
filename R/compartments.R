#' @title Contact map compartments
#' @name getCompartments
#' @rdname getCompartments
#' @description
#' 
#' Computes eigen vectors for each chromosome using cis contacts and extract 
#' chromosome compartments. 
#'
#' @param x A `HiCExperiment` object over a full genome
#' @param resolution Which resolution to use to compute eigen vectors
#' @param genome a BSgenome of DNAStringSet object associated with the Hi-C 
#'   contact matrix. 
#' @param chromosomes character or integer vector indicating which 
#' @param neigens Numver of eigen vectors to extract
#' @param sort_eigens Can be FALSE or one of c('Spearman', 'Pearson')
#' @param BPPARAM BiocParallel parallelization settings
#' @return A `HiCExperiment` object with additional `eigens` metadata containing the
#' normalized eigenvectors and a new "compartments" topologicalFeatures 
#' storing A and B compartments as a GRanges object. 
#'
#' @importFrom BiocParallel bplapply
#' @importFrom stats cor
#' @export
#' @examples 
#' library(HiContacts)
#' full_contacts_yeast <- contacts_yeast(full = TRUE)
#' comps <- getCompartments(full_contacts_yeast)
#' metadata(comps)$eigens

getCompartments <- function(
    x, 
    resolution = NULL, 
    genome = NULL, 
    chromosomes = NULL, 
    neigens = 3, 
    sort_eigens = FALSE, 
    BPPARAM = BiocParallel::bpparam()
) {
    message( "Going through preflight checklist..." )
    # - Check resolutions and chromosome subset
    if (is.null(resolution)) resolution <- resolution(x)
    chrs <- GenomicRanges::seqnames(GenomeInfoDb::seqinfo(x))
    names(chrs) <- chrs
    if (!is.null(chromosomes)) chrs <- chrs[chromosomes]

    # - If HiCExperiment comes from HiC-Pro object, make sure that the 
    #   imported matrix is genome-wide, and that `balanced` scores are available
    if ('regions' %in% names(metadata(x))) {
        if (!is.null(focus(x))) {
            message('Re-importing genome-wide HiCExperiment from HiC-Pro data files.')
            x <- HiCExperiment::HiCExperiment(
                fileName(x), bed = metadata(x)$regions
            )
        }
        if (!"balanced" %in% names(x)) {
            message('Balancing genome-wide matrix.')
            x <- normalize(x)
            scores(x, 'balanced') <- scores(x, 'ICE') 
        }
    }

    ## -- Parse contact matrix for each chromosome 
    message( "Parsing intra-chromosomal contacts for each chromosome..." )
    if (interactive()) {
        bpparam <- BiocParallel::SerialParam(progressbar = TRUE)
    }
    else {
        bpparam <- BiocParallel::SerialParam(progressbar = FALSE)
    }
    l_subs <- BiocParallel::bplapply(
        BPPARAM = bpparam, 
        chrs, 
        function(chr) {
            if ('regions' %in% names(metadata(x))) {
                x[chr]
            }
            else {
                HiCExperiment::HiCExperiment(
                    fileName(x), resolution = resolution, focus = chr
                )
            }
        }
    )
    names(l_subs) <- chrs

    ## -- Compute eigens on each chr. separately
    message( "Computing eigenvectors for each chromosome..." )
    compts <- BiocParallel::bplapply(
        BPPARAM = BPPARAM, 
        l_subs, 
        function(x_chr) {
            if (length(HiCExperiment::regions(x_chr)) <= neigens) {
                gr <- HiCExperiment::regions(x_chr)
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
                    sort_eigens = sort_eigens,
                    autocorrelation = FALSE,
                    clip_percentile = 0.999, 
                    ignore_diags = 2
                )
            }
        }
    ) |> GenomicRanges::GRangesList() |> unlist()

    ## -- Return initial HiCExperiment object with eigens in metadata
    A <- GenomicRanges::reduce(compts[compts$eigen > 0])
    A$compartment <- 'A'
    B <- GenomicRanges::reduce(compts[compts$eigen < 0])
    B$compartment <- 'B'
    cpts <- sort(c(A, B))
    HiCExperiment::topologicalFeatures(x, 'compartments') <- cpts
    metadata(x)$eigens <- compts
    return(x)

}

#' @importFrom RSpectra eigs_sym

.eigChr <- function(
    x_chr, 
    genome = NULL, 
    neigens = 3, 
    sort_eigens = FALSE,
    autocorrelation = FALSE, 
    clip_percentile = 0.99, 
    ignore_diags = 2
) {

    ## -- Get dense detrended matrix (with non-log2-detrended scores)
    if (!autocorrelation) { ## THIS IS THE DEFAULT (NO CORRELATION MATRIX)
        dx <- detrend(x_chr)
        gis <- InteractionSet::interactions(dx)
        gis$score <- HiCExperiment::scores(dx, 'detrended')
        gr <- HiCExperiment::regions(gis)
        m <- HiCExperiment::cm2matrix(HiCExperiment::gi2cm(gis))
        m <- 2^m ## Return to linear scale
        m <- m - 1 ## -- Subtract 1
    }
    else {
        dx <- autocorrelate(x_chr, ignore_ndiags = ignore_diags)
        gis <- InteractionSet::interactions(dx)
        gis$score <- HiCExperiment::scores(dx, 'autocorrelated')
        gr <- HiCExperiment::regions(gis)
        m <- HiCExperiment::cm2matrix(HiCExperiment::gi2cm(gis))
    }
    
    ## -- Remove ignored diagonals and NAs
    for (K in seq(-ignore_diags, ignore_diags, by = 1)) {
        sdiag(m, K) <- 0
    }
    m[is.na(m)] <- 0
    if (all(m == 0)) stop("Autocorrelation matrix is not complex enough. Try with a finer resolution.")

    ## -- Mask white lines from matrix
    m <- as.matrix(m)
    unmasked_bins <- rowSums(m) != 0
    m_no0 <- m[unmasked_bins, unmasked_bins]

    ## -- Clamp values at extremes
    m_no0 <- scales::oob_squish(
        m_no0, 
        stats::quantile(m_no0, probs = c(
            1-clip_percentile, clip_percentile
        ), na.rm = TRUE)
    )

    ## -- Get eigen values and vectors
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
        message("RSpectra package not found. Falling back to base eigenvectors computation.")
        message("Install RSpectra package to perform faster eigenvectors computation on larger matrices.")
        message("install.packages('RSpectra')")
        eigen <- base::eigen(m_no0, symmetric = TRUE)
    } 
    else {
        eigen <- RSpectra::eigs_sym(m_no0, k = neigens)
    }
    
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

    if (!is.null(genome)) {
        if (is(genome, "BSgenome") | is(genome, "DNAStringSet")) {
            ## -- Get GC content
            gr$phasing <- as.numeric(Biostrings::letterFrequency(
                Biostrings::getSeq(genome, gr), letters = 'GC') / 
                GenomicRanges::width(gr))
        }
        else if (is(genome, "TxDb")) {
            ## -- Get gene density
            if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
                message("Install GenomicFeatures package to phase eigenvector using a `TxDb` object.")
                message("install.packages('GenomicFeatures')")
                gr$phasing <- 1
            } 
            else {
                cov <- GenomicFeatures::genes(genome) |> 
                    GenomicRanges::coverage()
                gr$phasing <- BiocGenerics::mean(cov[gr])
            }
        }
        else if (is(genome, "RleList")) {
            ## -- Get average coverage per bin
            gr$phasing <- BiocGenerics::mean(cov[gr])
        }
        else {
            stop("`genome` argument not recognized.")
        }
        ## -- Get matching eigenvector and phase it with GC
        gr <- .eigGCPhasing(gr, neigens, unmasked_bins, sort_eigens)
    }
    else {
        message("Caution! No genome is provided. The first eigenvector has ", 
        "been selected and no phasing has been performed.")
        gr$eigen <- gr$E1
    }

    return(gr)

}

.eigGCPhasing <- function(gr, neigens, unmasked_bins, sort_eigens) {
    if (isFALSE(sort_eigens)) {
        method <- 'spearman'
    }
    else if (sort_eigens == 'Spearman') {
        method <- 'spearman'
    }
    else if (sort_eigens == 'Pearson') {
        method <- 'pearson'
    }
    else {
        stop("`sort_eigens` can be FALSE or one of 'Spearman' or 'Pearson'")
    }
    cors <- lapply(seq_len(neigens), function(K) {
        vec <- GenomicRanges::mcols(gr)[, paste0('E', K)]
        stats::cor(
            gr$phasing[unmasked_bins], 
            vec[unmasked_bins], 
            method = method
        )
    }) |> unlist()

    # -- Make all eigens positively correlate with phasing
    for (i in seq_len(neigens)) {
        if (cors[i] < 0) {
            GenomicRanges::mcols(gr)[, paste0('E', i)] <- -GenomicRanges::mcols(gr)[, paste0('E', i)]
            cors[i] <- -cors[i]
        }
    }

    # -- Sort eigens by highest corr.
    if (!isFALSE(sort_eigens)) {
        gr0 <- gr
        cors0 <- cors
        o <- order(cors, decreasing = TRUE)
        for (i in seq_len(neigens)) {
            GenomicRanges::mcols(gr0)[, paste0('E', i)] <- GenomicRanges::mcols(gr)[, paste0('E', o[i])]
            cors0[i] <- cors[o[i]]
        }
        gr <- gr0
        cors <- cors0
        gr$eigen <- gr$E1
    }
    else {
        gr$eigen <- gr$E1
    }
    
    # -- Print correlations
    fcors <- paste(round(cors, 2), collapse = ' / ')
    message(paste0(
        "seqnames `", 
        GenomicRanges::seqnames(gr)[1], "` : ", fcors
    ))

    return(gr)
}

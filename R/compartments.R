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
    BPPARAM = BiocParallel::bpparam()
) {
    message( "Going through preflight checklist..." )
    # - Check resolutions and chromosome subset
    if (is.null(resolution)) resolution <- resolution(x)
    chrs <- GenomicRanges::seqnames(GenomeInfoDb::seqinfo(x))
    names(chrs) <- chrs
    if (!is.null(chromosomes)) chrs <- chrs[chromosomes]

    ## -- Parse contact matrix for each chromosome 
    message( "Parsing intra-chromosomal contacts for each chromosome..." )
    l_subs <- BiocParallel::bplapply(
        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE), 
        chrs, 
        function(chr) {
            HiCExperiment::HiCExperiment(
                fileName(x), resolution = resolution, focus = chr
            )
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
    genome, 
    neigens = 3, 
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
        ## -- Get GC content
        gr$GC <- as.numeric(Biostrings::letterFrequency(
            Biostrings::getSeq(genome, gr), letters = 'GC') / 
            GenomicRanges::width(gr))

        ## -- Get matching eigenvector and phase it with GC
        gr <- .eigGCPhasing(gr, neigens)
    }
    else {
        message("Caution! No genome is provided. The first eigenvector has ", 
        "been selected and no phasing has been performed.")
        gr$eigen <- gr$E1
    }

    return(gr)

}

.eigGCPhasing <- function(gr, neigens) {
    cors <- lapply(seq_len(neigens), function(K) {
        stats::cor(gr$GC, GenomicRanges::mcols(gr)[, paste0('E', K)], method = 'spearman')
    }) |> unlist()
    best_cor <- which.max(abs(cors))
    fcors <- paste(round(cors, 2), collapse = ' / ')
    if (cors[best_cor] < 0) {
        message(paste0(
            "seqnames `", 
            GenomicRanges::seqnames(gr)[1], "` | correlations : ", fcors, " | eigen #", 
            best_cor,
            " selected (flipped)"
        ))
        gr$eigen <- -GenomicRanges::mcols(gr)[, paste0('E', best_cor)]
    }
    else {
        message(paste0(
            "seqnames `", 
            GenomicRanges::seqnames(gr)[1], "` | correlations : ", fcors, "| eigen #", 
            best_cor,
            " selected"
        ))
        gr$eigen <- GenomicRanges::mcols(gr)[, paste0('E', best_cor)]
    }
    return(gr)
}

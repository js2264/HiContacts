getDiamondInsulation <- function(x, resolution, window_size = 10 * resolution, BPPARAM, ...) {

    message( "Going through preflight checklist..." )
    # - Get si and bins from contact matrix
    file <- fileName(x)
    si <- seqinfo(x)
    if (is_cool(file) | is_mcool(file)) {
        all_bins <- .getCoolAnchors(file, resolution = resolution)
        si <- .cool2seqinfo(file, resolution)
    }
    else if (is_hic(file)) {
        all_bins <- .getHicAnchors(file, resolution = resolution)
        si <- .hic2seqinfo(file)
    }
    else if (is_hicpro_matrix(file) & is_hicpro_regions(...)) {
        all_bins <- .getHicproAnchors(...)
        si <- .hicpro2seqinfo(file)
    }
    fact <- window_size / resolution
    # - Filter non-relevant bins
    .bins <- lapply(seqnames(si), function(chr) {
        all_bins[seqnames(all_bins) == chr] |> 
            head(-(fact+1)) |> tail(-(fact+1))
    })
    bins <- do.call(c, .bins)
    bins <- bins[1:400]

    # - !!! HEAVY LOAD !!! Parse ALL pixels and convert to sparse matrix
    # - Full parsing has to be done since parallelized access to HDF5 is not supported 
    # - Once full parsing is done, parallelization is trivial
    message( "Parsing the entire contact matrice as a sparse matrix..." )
    if (is_cool(file) | is_mcool(file)) {
        l <- .dumpCool(file, resolution = resolution)
    }
    else if (is_hic(file)) {
        l <- .dumpHic(file, resolution = resolution)
    }
    else if (is_hic_pro(file)) {
        l <- .dumpHicpro(file, ...)
    }
    l <- Matrix::sparseMatrix(
        i= l[['pixels']]$bin1_id + 1,
        j= l[['pixels']]$bin2_id + 1,
        x= l[['pixels']]$count
    )

    # - Scan each bin individually
    message( "Scan each window and compute diamond insulation score..." )
    scores <- BiocParallel::bplapply(BPPARAM = BPPARAM, 
        seq_along(bins), function(K) {
            .getDiamondScore(l, bins[K])
        }
    ) |> unlist()
    bins$insulation <- log2( scores / mean(scores) )
    fbins <- bins[!is.na(bins$insulation)]

    # - Annotate boundaries
    mins <- which(diff(sign(diff(fbins$insulation)))==+2)+1
    maxs <- which(diff(sign(diff(fbins$insulation)))==-2)+1
    if (length(maxs) < length(mins)) maxs[[length(maxs)+1]] <- length(fbins)
    fbins$min <- seq_along(fbins) %in% mins
    fbins$prominence <- NA
    fbins$prominence[fbins$min] <- -{fbins$insulation[fbins$min] - fbins$insulation[seq_along(fbins) %in% maxs]}
    
    ## -- Return initial HiCExperiment object with insulation in metadata
    metadata(x)$insulation <- fbins
    return(x)

}

getBorders <- function(x, weak_threshold = 0.2, strong_threshold = 0.5) {

    fbins <- metadata(x)$insulation
    fbins$border <- fbins$prominence >= weak_threshold & fbins$min
    fbins$type <- ifelse(fbins$prominence < weak_threshold, NA, ifelse(fbins$prominence < strong_threshold, 'weak', 'strong'))
    gr <- fbins[which(fbins$border)]
    gr$score <- gr$prominence
    names(gr) <- gr$type
    S4Vectors::mcols(gr) <- data.frame(score = gr$score)
    topologicalFeatures(x, 'borders') <- gr
    return(x)

}

.getDiamondScore <- function(l, bin) {
    ids <- seq(bin$bin_id - fact, bin$bin_id + fact)
    weights <- all_bins[ids+1]$weight
    counts <- l[ids+1, ids+1]
    balanced <- t(apply(
        apply(counts, 1, `*`, weights), 2, `*`, weights
    ))
    s <- sum(
        balanced[seq(1, fact) + 1 - 1, seq(fact, fact + fact - 1) + 1 + 1],
        na.rm = TRUE
    )
    return(s)
}

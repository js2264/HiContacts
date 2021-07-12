plotMatrix <- function(gis, limits = NULL, dpi = 500, rasterize = TRUE, symmetrical = TRUE) {

    `%>%` <- magrittr::`%>%`

    ## -- Matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    } 
    else {
        M <- max(gis$score[abs(gis$bin1 - gis$bin2) >= 0], na.rm = TRUE)
        m <- min(gis$score[abs(gis$bin1 - gis$bin2) >= 0], na.rm = TRUE)
        limits <- c(m, M)
    }

    # -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(raster.dpi = dpi)
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    # -- Check number of chromosomes that were extracted
    nseqnames <- length(unique(as.vector(GenomicRanges::seqnames(InteractionSet::anchors(gis)[[1]]))))
    
    if (nseqnames == 1) { ## Single chromosome coordinates to plot (easy scenario)

        ## -- Convert gis to table and extract x/y
        mat <- gis %>% 
            tibble::as_tibble() %>%
            dplyr::mutate(
                x = floor(end1 - (end1 - start1)/2), 
                y = floor(end2 - (end2 - start2)/2)
            )

        ## -- Clamp scores to limits
        mat <- dplyr::mutate(mat, score = ifelse(score > M, M, ifelse(score < m, m, score)))

        ## -- Add lower triangular matrix scores (if symetrical)
        if (symmetrical) {
            mat <- rbind(mat, mat %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% dplyr::select(-x2))   
        }

        ## -- Plot matrix
        p <- ggmatrix(mat, cols = afmhotr_colors) + 
            plotFun + 
            ggplot2::labs(
                x = unique(mat$seqnames1),
                y = unique(mat$seqnames2)
            )

    }

    else {

        chroms <- gis %>% 
            tidyr::as_tibble() %>% 
            dplyr::group_by(seqnames2) %>% 
            dplyr::slice_max(order_by = end2, n = 1) %>% 
            dplyr::select(seqnames2, end2) %>% 
            dplyr::distinct()
        chroms$cumlength <- cumsum(c(0, chroms$end2)[1:nrow(chroms)])
        chroms$end2 <- NULL
        mat <- gis %>% 
            tibble::as_tibble() %>%
            dplyr::left_join(chroms, by = c(seqnames1 = 'seqnames2')) %>% dplyr::rename(cumlength_x = cumlength) %>% 
            dplyr::left_join(chroms, by = c(seqnames2 = 'seqnames2')) %>% dplyr::rename(cumlength_y = cumlength) %>% 
            dplyr::mutate(
                x = floor(end1 - (end1 - start1)/2) + cumlength_x, 
                y = floor(end2 - (end2 - start2)/2) + cumlength_y
            ) 
        mat <- dplyr::mutate(mat, score = ifelse(score > M, M, ifelse(score < m, m, score)))
        mat <- rbind(mat, mat %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% dplyr::select(-x2))

        ## -- Plot matrix
        p <- ggmatrix(mat, cols = afmhotr_colors) + 
            plotFun + 
            ggplot2::labs(
                x = 'Genome coordinates',
                y = 'Genome coordinates'
            ) + 
            ggplot2::geom_hline(yintercept = chroms$cumlength[-1], colour = 'black', alpha = 0.75, size = 0.15) + 
            ggplot2::geom_vline(xintercept = chroms$cumlength[-1], colour = 'black', alpha = 0.75, size = 0.15) 

    }

    p

}

plotOverExpected <- function(
    gis, 
    limits = c(-1, 1), 
    dpi = 500, rasterize = TRUE, 
    return_expected = FALSE, 
    return_all = FALSE
) {

    `%>%` <- magrittr::`%>%`

    # -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(raster.dpi = dpi)
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    ## -- Convert gis to table and extract x/y
    mat <- gis %>% 
        tibble::as_tibble() %>%
        dplyr::mutate(
            x = floor(end1 - (end1 - start1)/2), 
            y = floor(end2 - (end2 - start2)/2)
        ) 
    
    ## -- Compute expected and overExpected score
    mat <- normalizeOverExpected(mat)

    if (return_expected) {
        mat$score <- mat$expected
    } 
    # else if (return_all = TRUE) {
    #     mat$score <- mat$scoreOverExpected
    # }
    else {
        mat$score <- mat$scoreOverExpected
    }

    ## -- Matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
        limits <- c(m, M)
    } else {
        M <- max(mat$score, na.rm = TRUE, is.finite = TRUE)
        m <- min(mat$score, na.rm = TRUE, is.finite = TRUE)
        mm <- max(abs(c(m, M)))
        limits <- c(-abs(mm), abs(mm))
    }

    ## -- Clamp scores to limits
    mat <- dplyr::mutate(mat, overExpected = ifelse(score > M, M, ifelse(score < m, m, score)))

    ## -- Add lower triangular matrix scores
    mat <- rbind(mat, mat %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% dplyr::select(-x2))

    ## -- Plot matrix
    p <- ggmatrix(mat) + 
        plotFun + 
        ggplot2::labs(
            x = unique(mat$seqnames1),
            y = unique(mat$seqnames1)
        ) + ggplot2::scale_fill_gradientn(
            colors = bwr_colors, 
            na.value = '#FFFFFF', 
            limits = limits
        )

    ## -- Add bluring to avoid speckles 
    # p <- ggfx::with_blur(p, sigma = 10)
    
    p

}

plotTriangularMatrix <- function(gis, limits = NULL, truncate_tip = 0.2) {

    `%>%` <- magrittr::`%>%`

    ## -- Convert gis to table and extract x/y
    mat <- gis %>% 
        tibble::as_tibble() %>%
        dplyr::mutate(
            x = floor(end1 - (end1 - start1)/2), 
            y = floor(end2 - (end2 - start2)/2)
        ) 
    binsize <- GenomicRanges::width(InteractionSet::regions(gis)[1])
    max_bin_dist <- (max(mat$bin2 - mat$bin1) + 1) * truncate_tip
    mat_ <- gi2mat(gis, limits, truncate_tip)

    ## -- Compute triangular truncated border
    M_y <- max(mat_$y)
    uptri <- tibble::tibble(
        x = c(min(mat_$x), min(mat_$x) + max_bin_dist*binsize/2, max(mat_$x) - max_bin_dist*binsize/2, max(mat_$x)), 
        y = c(0, M_y, M_y, 0)
    )

    ## -- Plot matrix
    p <- ggtiltedmatrix(mat_, cols = afmhotr_colors, truncate_tip = truncate_tip) + 
        ggplot2::geom_polygon() + 
        ggplot2::labs(
            x = unique(mat$seqnames1),
            y = ''
        ) + 
        ggplot2::geom_polygon(
            data = uptri, 
            ggplot2::aes(x, y), 
            col = 'black', 
            fill = NA, 
            size = 0.1, 
            inherit.aes = FALSE
        )

    p

}

plotMatrixList <- function(ls, truncate_tip = 0.2, ...) {

    `%>%` <- magrittr::`%>%`

    ## -- Convert gis to table and extract x/y, for each element in `ls`
    mat <- lapply(ls, function(df) {
        df %>% 
        tibble::as_tibble() %>%
        dplyr::mutate(
            x = floor(end1 - (end1 - start1)/2), 
            y = floor(end2 - (end2 - start2)/2)
        )
    }) %>% dplyr::bind_rows(.id = 'class')
    max_bin_dist <- (max(mat$bin2 - mat$bin1) + 1) * truncate_tip
    mat_ <- lapply(names(ls), function(n) {
        gis <- ls[[n]]
        m <- gi2mat(gis, ...) %>% dplyr::mutate(item = n)
    }) %>% dplyr::bind_rows() %>% 
        dplyr::mutate(item = factor(item, names(ls)))

    ## -- Plot matrix
    p <- ggtiltedmatrix(mat_, cols = afmhotr_colors, truncate_tip = truncate_tip, nmatrices = length(ls)) + 
        ggplot2::geom_polygon() + 
        ggplot2::labs(
            x = unique(mat$seqnames1),
            y = ''
        ) + 
        ggplot2::facet_grid(item~.)
    
}

plotAggregatedMatrix <- function(file, coords, res = NULL, limits = NULL, dpi = 500, rasterize = TRUE, symmetrical = TRUE, BPPARAM = BiocParallel::bpparam()) {

    `%>%` <- magrittr::`%>%`

    ## -- Coerce loops into coords and coords2 
    if (class(coords) == 'GInteractions') {
        isLoops <- TRUE
        loops <- coords
        coords <- InteractionSet::anchors(loops)[[1]]
        coords2 <- InteractionSet::anchors(loops)[[2]]
    }
    else {
        isLoops <- FALSE
        coords2 <- NULL
    }

    ## -- Check that coords are all the same width
    if (length(unique(GenomicRanges::width(coords))) != 1) {
        stop("Please provide GRanges that are all the same width. Aborting now.")
    }
    else {
        wi <- unique(GenomicRanges::width(coords))
        if (is.null(res)) {
            res0 <- unique(GenomicRanges::width(cool2gi(file, res = res, coords = coords[1], coords2 = coords[1]))[[1]])
        }
        else {
            res0 <- res
        }
        breaks <- seq({-wi/2}-res0, {wi/2}+res0, length.out = wi/res0 + 3)
    }

    ## -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(raster.dpi = dpi)
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    ## -- Calculate scores over each coords (+coords2)
    mats <- BiocParallel::bplapply(BPPARAM = BPPARAM, seq_along(coords), function(K) {
    # mats <- lapply(seq_along(coords), function(K) {
        range <- coords[K]
        # print(range)
        `%<-%` <- zeallot::`%<-%`
        c(coords_chr, coords_start, coords_end) %<-% splitCoords(range)
        midpoint <- coords_start + (coords_end - coords_start) / 2
        if (isLoops) {
            range2 <- coords2[K]
            `%<-%` <- zeallot::`%<-%`
            c(coords_chr2, coords_start2, coords_end2) %<-% splitCoords(range2)
            midpoint2 <- coords_start2 + (coords_end2 - coords_start2) / 2
        }
        else {
            range2 <- range
            midpoint2 <- midpoint
        }

        ## -- Filter to interactions which are within the range
        gis_sub <- cool2gi(file, res = res, coords = range, coords2 = range2)
        
        ## -- Convert gis_sub to table and extract x/y and adjust center to 0
        mat <- gis_sub %>% 
            tibble::as_tibble() %>%
            dplyr::mutate(
                x = center1, 
                y = center2, 
                ID = K
            ) %>% 
            dplyr::mutate(
                x = breaks[as.numeric(cut(x - midpoint, breaks, include.lowest = TRUE)) + 1], 
                y = breaks[as.numeric(cut(y - midpoint2, breaks, include.lowest = TRUE)) + 1], 
                score = scale(score)
            )
        
        # ggmatrix(mat, cols = afmhotr_colors) + 
        #     plotFun + 
        #     ggplot2::labs(
        #         x = unique(mat$seqnames1),
        #         y = unique(mat$seqnames1)
        #     )

        mat

    }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::group_by(x, y) %>% 
        dplyr::filter(x >= min(breaks[2:{length(breaks)-1}]), x <= max(breaks[2:{length(breaks)-1}]), y >= min(breaks[2:{length(breaks)-1}]), y <= max(breaks[2:{length(breaks)-1}]))
    
    ## -- Average scores for each bin
    mats <- mats %>% 
        dplyr::summarize(score = mean(score, na.rm = TRUE)) %>% 
        dplyr::mutate(seqnames1 = 'Aggr. coordinates') 
    
    ## -- Matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    } 
    else {
        M <- max(mats$score, na.rm = TRUE)
        m <- min(mats$score, na.rm = TRUE)
        limits <- c(m, M)
    }

    ## -- Clamp scores to limits
    mats <- mats %>% 
        dplyr::mutate(score = ifelse(score > M, M, ifelse(score < m, m, score)))

    ## -- Add lower triangular matrix scores 
    if (symmetrical) {
        mats <- rbind(mats, mats %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% 
            dplyr::select(-x2))
    }
    
    ## -- Plot matrix
    p <- ggmatrix(mats, cols = afmhotr_colors) + 
        plotFun + 
        ggplot2::labs(
            x = unique(mats$seqnames1),
            y = unique(mats$seqnames1)
        )

    p 

}
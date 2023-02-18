#' HiContacts plotting functionalities
#' 
#' @name HiContacts-plots
#' 
#' @description 
#' Several plots can be generated in HiContacts: 
#' - Hi-C contact matrices
#' - Distance-dependant interaction frequency decay (a.k.a. "Distance law" or "P(s)")
#' - Virtual 4C profiles
#' - Scalograms
#' - Saddle plots
NULL 

################################################################################
#                                                                              #
#                                 Matrix                                       #
#                                                                              #
################################################################################

#' Plotting a contact matrix
#'
#' @name plotMatrix
#' @aliases plotMatrix,HiCExperiment-method
#' @aliases plotMatrix,HiCExperiment-method
#' @aliases plotMatrix,GInteractions-method
#' @aliases plotMatrix,matrix-method
#' 
#' @param x A HiCExperiment object
#' @param compare.to Compare to a second HiC matrix in the lower left corner
#' @param use.scores Which scores to use in the heatmap
#' @param scale Any of 'log10', 'log2', 'linear', 'exp0.2' (Default: 'log10')
#' @param limits color map limits
#' @param maxDistance maximum distance. If provided, the heatmap is plotted 
#'   horizontally
#' @param loops Loops to plot on top of the heatmap, provided as `GInteractions`
#' @param borders Borders to plot on top of the heatmap, provided as `GRanges`
#' @param tracks Named list of bigwig tracks imported as `Rle`
#' @param dpi DPI to create the plot (Default: 500)
#' @param rasterize Whether the generated heatmap is rasterized or vectorized 
#' (Default: TRUE)
#' @param symmetrical Whether to enforce a symetrical heatmap (Default: TRUE)
#' @param chrom_lines Whether to display separating lines between chromosomes, 
#' should any be necessary (Default: TRUE)
#' @param show_grid Whether to display an underlying grid (Default: FALSE)
#' @param caption Whether to display a caption (Default: TRUE)
#' @param cmap Color scale to use. (Default: bgrColors() if limits are c(-1, 1)
#' and coolerColors() otherwise)
#' @return ggplot object
#'
#' @importFrom ggrastr geom_tile_rast
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr slice_max
#' @importFrom dplyr distinct
#' @importFrom dplyr left_join
#' @importFrom dplyr rename
#' @importFrom dplyr left_join
#' @importFrom dplyr rename
#' @importFrom GenomicRanges seqnames
#' @importFrom scales oob_squish
#' @importFrom scales unit_format
#' @examples 
#' contacts_yeast <- contacts_yeast()
#' plotMatrix(
#'     contacts_yeast, 
#'     use.scores = 'balanced', 
#'     scale = 'log10', 
#'     limits = c(-4, -1)
#' )
NULL 

#' @rdname plotMatrix
#' @export

setMethod("plotMatrix", "HiCExperiment", function(
    x, 
    compare.to = NULL, 
    use.scores = 'balanced', 
    scale = 'log10', 
    maxDistance = NULL, 
    loops = NULL, 
    borders = NULL, 
    tracks = NULL, 
    limits = NULL, 
    dpi = 500, 
    rasterize = TRUE, 
    symmetrical = TRUE, 
    chrom_lines = TRUE, 
    show_grid = FALSE, 
    cmap = NULL, 
    caption = TRUE  
) {
    if (!use.scores %in% names(scores(x))) 
        stop(paste0("Queried scores not found."))

    if (!is.null(compare.to)) {
        gis_x <- interactions(x)
        gis_y <- interactions(compare.to)
        gis <- c(
            as(gis_y, 'GInteractions'), 
            InteractionSet::swapAnchors(
                as(gis_x, 'GInteractions'), 
                mode = 'reverse'
            )
        )
        gis <- gis[pairdist(gis) != 0]
        p <- plotMatrix(
            gis, 
            use.scores = use.scores, 
            scale = scale, 
            maxDistance = maxDistance, 
            loops = loops, 
            borders = borders, 
            tracks = tracks, 
            limits = limits, 
            dpi = dpi, 
            rasterize = rasterize, 
            symmetrical = FALSE, 
            chrom_lines = chrom_lines, 
            show_grid = show_grid, 
            cmap = cmap  
        )
    }
    else {
        p <- plotMatrix(
            interactions(x), 
            use.scores = use.scores, 
            scale = scale, 
            maxDistance = maxDistance, 
            loops = loops, 
            borders = borders, 
            tracks = tracks, 
            limits = limits, 
            dpi = dpi, 
            rasterize = rasterize, 
            symmetrical = symmetrical, 
            chrom_lines = chrom_lines, 
            show_grid = show_grid, 
            cmap = cmap  
        )
    }
    if (caption & is.null(tracks)) {
        p <- p + ggplot2::labs(
            caption = paste(
                sep = '\n',
                paste0('file: ', fileName(x)), 
                paste0('resolution: ', resolution(x)), 
                paste0('focus: ', focus(x)),
                paste0('scores: ', use.scores),
                paste0('scale: ', scale)
            )
        )
    }
    return(p)
})

#' @rdname plotMatrix
#' @export

setMethod("plotMatrix", "GInteractions", function(
    x, 
    use.scores = NULL, 
    scale = 'log10', 
    maxDistance = NULL, 
    loops = NULL, 
    borders = NULL, 
    tracks = NULL, 
    limits = NULL, 
    dpi = 500, 
    rasterize = TRUE, 
    symmetrical = TRUE, 
    chrom_lines = TRUE, 
    show_grid = FALSE, 
    cmap = NULL  
) {
    `%over%` <- IRanges::`%over%`

    ## -- Extract scores
    gis <- x
    if (!is.null(use.scores)) {
        gis$score <- S4Vectors::mcols(gis)[, use.scores]
    }
    else {
        if ("balanced" %in% colnames(S4Vectors::mcols(gis))) {
            gis$score <- S4Vectors::mcols(gis)[, 'balanced']
        } 
        else {
            gis$score <- S4Vectors::mcols(gis)[, 'count']
        }
    }

    ## -- Put metric to plot in `score` column
    has_negative_scores <- any(gis$score < 0, na.rm = TRUE)

    ## -- Scale score
    if (scale == 'log10') {
        gis$score <- log10(gis$score)
    } 
    else if (scale == 'log2') {
        gis$score <- log2(gis$score)
    }
    else if (scale == 'exp0.2') {
        gis$score <- gis$score^0.2
    }
    else if (scale == 'linear') {
        gis$score <- gis$score
    }

    ## -- Set matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    }
    else {
        .scores <- gis$score[pairdist(x) != 0 & !is.na(pairdist(x) != 0)]
        .scores <- .scores[!is.na(.scores)]
        .scores <- .scores[!is.infinite(.scores)]
        M <- max(.scores)
        m <- min(.scores)
        limits <- c(m, M)
    }

    ## -- Choose color map 
    if (is.null(cmap)) {
        if (has_negative_scores | {m == -1 & M == 1}) {
            cmap <- bgrColors()
        }
        else {
            cmap <- coolerColors()
        }
    }
    
    # -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(
            raster.dpi = dpi,
            width = GenomicRanges::width(gis)[[1]][1], 
            height = GenomicRanges::width(gis)[[1]][1]
        )
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    ## -- If loops are provided, filter them and add
    if (!is.null(loops) & is.null(maxDistance)) {
        filtered_loops <- tibble::as_tibble(
            loops[anchors(loops)[['first']] %over% gis & anchors(loops)[['second']] %over% gis]
        )
        p_loops <- ggplot2::geom_point(
            data = filtered_loops, 
            mapping = ggplot2::aes(
                x = start1+(end1-start1)/2, y = start2+(end2-start2)/2
            ), 
            inherit.aes = FALSE, 
            shape = 21, 
            fill = "#76eaff00", 
            color = "#000000", 
            size = 3
        )
    }
    else {
        p_loops <- NULL
    }

    ## -- If borders are provided, filter them and add
    if (!is.null(borders) & is.null(maxDistance)) {
        filtered_borders <- tibble::as_tibble(
            borders[borders %over% gis]
        )
        p_borders <- ggplot2::geom_point(
            data = filtered_borders, 
            mapping = ggplot2::aes(x = start+(end-start)/2, y = start+(end-start)/2),
            inherit.aes = FALSE, 
            shape = 23, 
            fill = "#76eaff", 
            color = "#000000", 
            size = 1
        )
    }
    else {
        p_borders <- NULL
    }

    ## -- If tracks are provided, filter them and add
    if (!is.null(tracks) & is.null(maxDistance)) {
        re <- reduce(regions(gis))
        # if (length(re) > 1) {
        #     warning("Input track cannot be aligned to interactions with non-contiguous regions.")
        #     p_tracks <- ggplot2::ggplot() + ggplot2::theme_void()
        # }
        # else {
            ntracks <- length(tracks)
            step <- GenomicRanges::width(InteractionSet::anchors(gis[1])[[1]])/5
            filtered_tracks <- purrr::imap(tracks, ~ .x[re] |> 
                tibble::as_tibble() |> 
                dplyr::group_by(group=(0:(n()-1))%/%{step}) |> 
                dplyr::summarize(value = mean(value, na.rm = TRUE)) |> 
                dplyr::mutate(
                    # pos = seq(IRanges::start(re), IRanges::end(re), by = step) + step/2, 
                    track = .y
                )
            ) |> dplyr::bind_rows(.id = 'track')
            p_tracks <- ggplot2::ggplot() + 
                ggplot2::geom_area(
                    data = filtered_tracks, 
                    mapping = ggplot2::aes(x = group, y = value),
                    fill = "#393939", 
                    color = NA, 
                ) + 
                ggplot2::scale_x_continuous(expand = c(0, 0)) +
                ggplot2::scale_y_reverse(expand = c(0, 0)) +
                ggplot2::facet_grid(track~., scales = 'free') + 
                ggplot2::theme_void() + 
                ggplot2::theme(
                    strip.background = ggplot2::element_blank(), 
                    strip.text = ggplot2::element_blank()
                )
        # }
    }
    else {
        p_tracks <- ggplot2::ggplot() + ggplot2::theme_void()
    }

    # -- Check number of chromosomes that were extracted
    nseqnames <- length(unique(as.vector(
        GenomicRanges::seqnames(InteractionSet::regions(gis))
    )))

    if (nseqnames == 1) { ## SINGLE CHROMOSOME MAP
        
        if (is.null(maxDistance)) { ##### REGULAR SQUARE MATRIX
            ## -- Convert gis to table and extract x/y
            mat <- gis |>
                tibble::as_tibble() |>
                dplyr::mutate(
                    x = floor(end1 - (end1 - start1) / 2),
                    y = floor(end2 - (end2 - start2) / 2)
                ) |> 
                tidyr::drop_na(score)

            ## -- Clamp scores to limits
            mat$score = scales::oob_squish(mat$score, c(m, M))

            ## -- Add lower triangular matrix scores (if symetrical)
            if (symmetrical) {
                mat <- rbind(
                    mat, 
                    mat |> 
                        dplyr::mutate(x2 = y, y = x, x = x2) |> 
                        dplyr::select(-x2)
                )
                bb <- InteractionSet::boundingBox(x)
                bbOne <- InteractionSet::anchors(bb)[[1]]
                bbTwo <- InteractionSet::anchors(bb)[[2]]
                fraction_ov <- GenomicRanges::width(GenomicRanges::pintersect(bbOne, bbTwo)) /
                    {{GenomicRanges::width(bbOne) + GenomicRanges::width(bbTwo)}/2}
                if (fraction_ov > 0.9) {
                    coords <- c(bbOne, bbTwo)
                    mat <- mat |> 
                        dplyr::filter(x >= GenomicRanges::start(coords[2]) & 
                            x <= GenomicRanges::end(coords[2])) |> 
                        dplyr::filter(y >= GenomicRanges::start(coords[1]) & 
                            y <= GenomicRanges::end(coords[1]))
                }
            } 

            ## -- Plot matrix
            p <- ggMatrix(mat, cols = cmap, limits = limits, grid = show_grid) +
                plotFun +
                p_loops + 
                p_borders + 
                ggplot2::labs(
                    x = unique(mat$seqnames1),
                    y = "Genome coordinates"
                ) + 
                ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6), position = 'top') +
                ggplot2::scale_y_reverse(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6))

            ## -- Add tracks 
            if (length(p_tracks$layers)) {
                p <- p + 
                ggplot2::theme(
                    plot.margin = grid::unit(c(5.5, 5.5, 5.5 + 10*ntracks, 5.5), "points")
                ) + patchwork::inset_element(
                    p_tracks, 
                    left = 0, 
                    bottom = -0.1*ntracks, 
                    right = 1, 
                    top = 0, 
                    align_to = 'panel', 
                    clip = FALSE
                ) 
            }
            
        }
        
        else { ##### HORIZONTAL TRIANGULAR MATRIX
            ## -- Convert gis to table and extract x/y
            df <- gis |>
                tibble::as_tibble() |>
                dplyr::mutate(
                    diag = center2 - center1, 
                    x = center1 + (center2 - center1)/2
                ) |> 
                tidyr::drop_na(score) |> 
                dplyr::filter(diag <= maxDistance)

            ## -- Clamp scores to limits
            df <- dplyr::mutate(df, score = scales::oob_squish(score, c(m, M)))

            ## -- Plot matrix
            p <- ggHorizontalMatrix(df, cols = cmap, limits = limits) +
                plotFun +
                p_loops + 
                p_borders + 
                ggplot2::labs(
                    x = "Genomic location",
                    y = "Distance"
                )
        }
    }

    else { ##### MULTI-CHROMOSOMES MAP

        ## -- Abort if more than 1 chr. is plotted AND maxDistance is provided
        if (!is.null(maxDistance)) stop(
            "Horizontal, multi-chromosomes maps are not supported."
        )

        chroms <- gis |>
            tidyr::as_tibble() |>
            dplyr::group_by(seqnames2) |>
            dplyr::slice_max(order_by = end2, n = 1) |>
            dplyr::select(seqnames2, end2) |>
            dplyr::distinct()
        chroms$cumlength <- cumsum(c(0, chroms$end2)[seq_len(nrow(chroms))])
        chroms$end2 <- NULL
        mat <- gis |>
            tibble::as_tibble() |>
            dplyr::left_join(chroms, by = c(seqnames1 = "seqnames2")) |>
            dplyr::rename(cumlength_x = cumlength) |>
            dplyr::left_join(chroms, by = c(seqnames2 = "seqnames2")) |>
            dplyr::rename(cumlength_y = cumlength) |>
            dplyr::mutate(
                x = floor(end1 - (end1 - start1) / 2) + cumlength_x,
                y = floor(end2 - (end2 - start2) / 2) + cumlength_y
            )
        mat <- dplyr::mutate(
            mat, 
            score = ifelse(score > M, M, ifelse(score < m, m, score))
        )
        mat <- rbind(
            mat, 
            mat |> 
                dplyr::mutate(x2 = y, y = x, x = x2) |> 
                dplyr::select(-x2)
        )

        ## -- Plot matrix
        p <- ggMatrix(mat, cols = cmap, limits = limits, grid = show_grid) +
            plotFun +
            ggplot2::labs(
                x = "Genome coordinates",
                y = "Genome coordinates"
            ) + 
            ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6), position = 'top') +
            ggplot2::scale_y_reverse(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6))
        if (chrom_lines) {
            p <- p +
                ggplot2::geom_hline(
                    yintercept = chroms$cumlength[-1], 
                    colour = "black", alpha = 0.75, linewidth = 0.15
                ) +
                ggplot2::geom_vline(
                    xintercept = chroms$cumlength[-1], 
                    colour = "black", alpha = 0.75, linewidth = 0.15
                )
        }

        ## -- Add tracks 
        if (length(p_tracks$layers)) {
            p <- p + 
            ggplot2::theme(
                plot.margin = grid::unit(c(5.5, 5.5, 5.5 + 10*ntracks, 5.5), "points")
            ) + patchwork::inset_element(
                p_tracks, 
                left = 0, 
                bottom = -0.1*ntracks, 
                right = 1, 
                top = 0, 
                align_to = 'panel', 
                clip = FALSE
            ) 
        }
        
    }

    p
})

#' @rdname plotMatrix
#' @export

setMethod("plotMatrix", "matrix", function(
    x, 
    scale = 'log10', 
    limits = NULL, 
    dpi = 500, 
    rasterize = TRUE, 
    cmap = NULL  
) {    

    mat <- x
    has_negative_scores <- any(mat < 0, na.rm = TRUE)

    ## -- Scale score
    if (scale == 'log10') {
        mat <- log10(mat)
    } 
    else if (scale == 'log2') {
        mat <- log2(mat)
    }
    else if (scale == 'exp0.2') {
        mat <- mat^0.2
    }
    else if (scale == 'linear') {
        mat <- mat
    }

    ## -- Set matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    }
    else {
        mat0 <- mat
        sdiag(mat0, 0) <- NA
        .scores <- mat0
        .scores <- .scores[!is.na(.scores)]
        .scores <- .scores[!is.infinite(.scores)]
        M <- max(.scores)
        m <- min(.scores)
        limits <- c(m, M)
    }

    ## -- Choose color map 
    if (is.null(cmap)) {
        if (has_negative_scores) {
            cmap <- bgrColors()
        }
        else {
            cmap <- coolerColors()
        }
    }
    
    # -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(
            raster.dpi = dpi
        )
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    ## -- Convert gis to table and extract x/y
    colnames(mat) <- paste0("bin", seq_len(ncol(mat)))
    mat <- tibble::as_tibble(mat) |>
        dplyr::mutate(x = seq_len(ncol(mat))) |>
        tidyr::pivot_longer(-x, names_to = 'y', values_to = 'score') |> 
        dplyr::mutate(y = gsub('bin', '', y) |> as.numeric()) |>
        tidyr::drop_na(score)

    ## -- Clamp scores to limits
    mat <- dplyr::mutate(mat, score = scales::oob_squish(score, c(m, M)))

    ## -- Plot matrix
    p <- ggMatrix(mat, cols = cmap, limits = limits) +
        plotFun +
        ggplot2::labs(
            x = "Bin number",
            y = "Bin number"
        ) + 
        ggplot2::scale_x_continuous(expand = c(0, 0), position = 'top') +
        ggplot2::scale_y_reverse(expand = c(0, 0))

    p
})

#' @rdname plotMatrix
#' @export

setMethod("plotMatrix", "AggrHiCExperiment", function(
    x, 
    use.scores = 'balanced', 
    scale = 'log10', 
    maxDistance = NULL, 
    loops = NULL, 
    borders = NULL, 
    limits = NULL, 
    dpi = 500, 
    rasterize = TRUE, 
    chrom_lines = TRUE, 
    show_grid = FALSE, 
    cmap = NULL, 
    caption = TRUE  
 ) {
    pairs <- topologicalFeatures(x, 'targets')
    is1D <- all(S4Vectors::first(pairs) == S4Vectors::second(pairs))
    p <- plotMatrix(
        interactions(x), 
        use.scores = use.scores, 
        scale = scale, 
        maxDistance = maxDistance, 
        loops = loops, 
        borders = borders, 
        limits = limits, 
        dpi = dpi, 
        rasterize = rasterize, 
        symmetrical = ifelse(is1D, TRUE, FALSE), 
        chrom_lines = chrom_lines, 
        show_grid = show_grid, 
        cmap = cmap  
    )
    if (caption) {
        p <- p + ggplot2::labs(
            caption = paste(
                sep = '\n',
                paste0('file: ', fileName(x)), 
                paste0('resolution: ', resolution(x)), 
                paste0('focus: ', focus(x)),
                paste0('scores: ', use.scores),
                paste0('scale: ', scale)
            )
        )
    }
    return(p)
})

#' @rdname plotMatrix
#' @export

setMethod("montage", "AggrHiCExperiment", function(
    x, 
    use.scores = 'balanced', 
    scale = 'log10', 
    limits = NULL, 
    dpi = 500, 
    rasterize = TRUE, 
    cmap = NULL
 ) {
    pairs <- topologicalFeatures(x, 'targets')
    is1D <- ifelse(is(pairs, 'GRanges'), TRUE, FALSE)
    sli <- slices(x, use.scores)
    snips <- lapply(seq_len(dim(sli)[3]), function(K){
        mat <- sli[, , K]
        if (is1D) mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
        plotMatrix(
            mat, 
            scale = scale, 
            limits = limits, 
            dpi = dpi, 
            rasterize = rasterize, 
            cmap = cmap
        ) + 
            ggplot2::theme_void() + 
            ggplot2::theme(legend.position = "none") + 
            ggplot2::theme(plot.margin = margin(.002, .002, .002, .002, "npc"))
    })
    p <- patchwork::wrap_plots(snips)
})

ggMatrix <- function(mat, ticks = TRUE, grid = FALSE, cols = coolerColors(), limits) {
    p <- ggplot2::ggplot(mat, ggplot2::aes(x, y, fill = score))
    p <- p + ggplot2::scale_fill_gradientn(
        colors = cols,
        na.value = "#FFFFFF",
        limits = limits
    ) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(.25, "npc"), barwidth = 0.5, frame.colour = "black")) + 
        ggplot2::coord_fixed() +
        ggthemeHiContacts(ticks, grid)
    p
}

ggHorizontalMatrix <- function(df, cols = coolerColors(), limits) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x, diag, fill = score))
    r <- 1/sqrt(2)/sqrt(2)
    p <- p + ggplot2::scale_fill_gradientn(
        colors = cols,
        na.value = "#FFFFFF",
        limits = limits
    ) +
        ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(.25, "npc"), barwidth = 0.5, frame.colour = "black")) + 
        ggplot2::coord_fixed(ratio = r) +
        ggthemeHiContacts(ticks = TRUE, grid = FALSE)
    p
}

################################################################################
#                                                                              #
#                                 P(s)                                         #
#                                                                              #
################################################################################

#' Plotting a P(s) distance law
#' 
#' @name plotPs
#' @aliases plotPs
#' @aliases plotPsSlope
#' 
#' @param x the output data.frame of `distanceLaw` function
#' @param mapping aes to pass on to ggplot2
#' @param xlim xlim
#' @param ylim ylim
#' @return ggplot
#'
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @examples 
#' ## Single P(s)
#' 
#' contacts_yeast <- contacts_yeast()
#' ps <- distanceLaw(contacts_yeast)
#' plotPs(ps, ggplot2::aes(x = binned_distance, y = norm_p))
#' 
#' ## Comparing several P(s)
#' 
#' contacts_yeast <- contacts_yeast()
#' contacts_yeast_eco1 <- contacts_yeast_eco1()
#' ps_wt <- distanceLaw(contacts_yeast)
#' ps_wt$sample <- 'WT'
#' ps_eco1 <- distanceLaw(contacts_yeast_eco1)
#' ps_eco1$sample <- 'eco1'
#' ps <- rbind(ps_wt, ps_eco1)
#' plotPs(ps, ggplot2::aes(x = binned_distance, y = norm_p, group = sample, color = sample))
#' plotPsSlope(ps, ggplot2::aes(x = binned_distance, y = slope, group = sample))
NULL

#' @rdname plotPs
#' @export

plotPs <- function(x, mapping, xlim = c(5000, 4.99e5), ylim = c(1e-8, 1e-4)) {
    p <- ggplot2::ggplot(x, mapping = mapping) + 
        ggplot2::geom_line() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::annotation_logticks() + 
        ggplot2::labs(x = "Genomic distance", y = "Contact frequency")
    if (!is.null(ylim)) {
        p <- p + ggplot2::scale_y_log10(
            limits = ylim, 
            expand = c(0, 0), 
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
    }
    else {
        p <- p + ggplot2::scale_y_log10(
            expand = c(0, 0), 
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
    }
    if (!is.null(xlim)) {
        p <- p + ggplot2::scale_x_log10(
            limits = xlim, 
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        )
    }
    else {
        p <- p + ggplot2::scale_x_log10(
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        )
    }
    p
}

#' @rdname plotPs
#' @export

plotPsSlope <- function(x, mapping, xlim = c(5000, 4.99e5), ylim = c(-3, 0)) {
    gg <- ggplot2::ggplot(x, mapping = mapping) + 
        ggplot2::geom_line() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_y_continuous(
            limits = ylim, 
            expand = c(0, 0)
        ) +
        ggplot2::scale_x_log10(
            limits = xlim, 
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        ) + 
        ggplot2::annotation_logticks(sides = "b") + 
        ggplot2::labs(x = "Genomic distance", y = "Slope of P(s)")
    gg
}

################################################################################
#                                                                              #
#                                 virtual 4C                                   #
#                                                                              #
################################################################################

#' Plotting virtual 4C profiles
#' 
#' @name plot4C
#' @aliases plot4C
#' 
#' @param x GRanges, generally the output of `virtual4C()`
#' @param mapping aes to pass on to ggplot2 (default: 
#' `ggplot2::aes(x = center, y = score, col = seqnames)`)
#' @return ggplot
#' 
#' @import tibble
#' @importFrom scales unit_format
#' @examples 
#' contacts_yeast <- contacts_yeast()
#' v4C <- virtual4C(contacts_yeast, GenomicRanges::GRanges('II:490000-510000'))
#' plot4C(v4C)
NULL

#' @rdname plot4C
#' @export

plot4C <- function(x, mapping = ggplot2::aes(x = center, y = score, col = seqnames)) {
    x <- tibble::as_tibble(x)
    p <- ggplot2::ggplot(x, mapping = mapping) + 
        # ggplot2::geom_line() + 
        ggplot2::geom_area(position = "identity", alpha = 0.2) + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::labs(x = "Position", y = "Contacts with viewpoint") + 
        ggplot2::scale_x_continuous(
            expand = c(0, 0), 
            labels = scales::unit_format(unit = "M", scale = 1e-6)
        )
    p
}

################################################################################
#                                                                              #
#                                 Scalograms                                   #
#                                                                              #
################################################################################

#' Plotting scalograms
#' 
#' @name plotScalogram
#' @aliases plotScalogram
#' 
#' @param x GRanges, the output of `scalogram()`
#' @param ylim Range of distances to use for y-axis in scalograms
#' @return ggplot
#' 
#' @import tibble
#' @examples 
#' contacts_yeast <- HiCExperiment::contacts_yeast()
#' pairsFile(contacts_yeast) <- HiContactsData::HiContactsData(
#'   'yeast_wt', format = 'pairs.gz'
#' )
#' scalo <- scalogram(contacts_yeast['II'])
#' plotScalogram(scalo)
NULL

#' @rdname plotScalogram
#' @export

plotScalogram <- function(x, ylim = c(5e2, 1e5)) {
    x_ <- dplyr::mutate(x, dist_quantile = scales::oob_squish(dist_quantile, c(ylim[[1]], ylim[[2]]))) |>
        dplyr::mutate(y0 = ylim[[1]])
    p <- ggplot2::ggplot() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) + 
        ggplot2::scale_y_log10(
            expand = c(0, 0), 
            limits = c(ylim[[1]], ylim[[2]])
        ) + 
        ggplot2::annotation_logticks(sides = 'l') + 
        ggplot2::labs(x = 'Position along chr.', y = 'Distance of interactions')
    if (length(unique(x_$prob)) == 3) {
        x_$prob <- as.numeric(as.character(x_$prob))
        low <- sort(unique(x_$prob))[1]
        mid <- sort(unique(x_$prob))[2]
        high <- sort(unique(x_$prob))[3]
        x_$prob <- dplyr::case_when(
            x_$prob == low ~ 'low', 
            x_$prob == mid ~ 'mid', 
            x_$prob == high ~ 'high'
        )
        diff <- {high - low} * 100
        p <- p + 
            ggplot2::geom_ribbon(
                data = tidyr::pivot_wider(x_[x_$prob != 'mid',], names_from = prob, values_from = dist_quantile), 
                mapping = ggplot2::aes(x = binned_pos, ymin = low, ymax = high, fill = y0), 
                col = 'black', linewidth = 0.2, linetype = 'dashed'
            ) + 
            ggplot2::geom_line(
                data = x_[x_$prob == 'mid',], 
                ggplot2::aes(x = binned_pos, y = dist_quantile, col = prob)
            ) + 
            ggplot2::scale_fill_gradient(low = '#b5b5b54d', high = '#b5b5b54d') + 
            ggplot2::scale_color_manual(values = '#000000') + 
            ggplot2::labs(fill = paste0(diff, '% (', low, '% - ', high, '%)\nof interactions\nwithin this range')) + 
            ggplot2::labs(col = paste0(mid*100, '% of interaction distance')) + 
            ggplot2::theme(legend.text = ggplot2::element_blank(), legend.position = 'bottom')
    } 
    else {
        p <- p + 
            ggplot2::geom_ribbon(data = x_, ggplot2::aes(
                x = binned_pos, ymax = dist_quantile,
                group = prob, fill = prob
            )) + 
            ggplot2::aes(ymin = y0) + 
            ggplot2::scale_fill_brewer(palette = 'Spectral') 
    }
    if (length(unique(x_$chr)) > 1) {
        p <- p + ggplot2::facet_wrap(~chr)
    }
    p
}

################################################################################
#                                                                              #
#                                 Saddle plots                                 #
#                                                                              #
################################################################################

#' Plotting saddle plots
#' 
#' @name plotSaddle
#' @aliases plotSaddle
#' 
#' @param x a HiCExperiment object with a stored `eigens` metadata
#' @param nbins Number of bins to use to discretize the eigenvectors
#' @param limits limits for color map being used
#' @param plotBins Whether to plot the distribution of bins on top of the plot
#' @param BPPARAM a BiocParallel registered method 
#' @return ggplot
NULL

#' @rdname plotSaddle
#' @export

plotSaddle <- function(x, nbins = 50, limits = c(-1, 1), plotBins = FALSE, BPPARAM = BiocParallel::bpparam()) {
    if (is.null(metadata(x)$eigens)) stop(
        "No eigen vector found in metadata. Run getCompartments(x) first."
    )
    eigens <- metadata(x)$eigens

    ## -- Filter and bin regions by their eigenvector score
    eigens$eigen_bin <- cut(
        eigens$eigen, breaks = stats::quantile(
            eigens$eigen[eigens$eigen != 0], 
            probs = seq(0.025, 0.975, length.out = nbins+1)
        ), 
        include.lowest = TRUE
    ) |> as.numeric()

    ## -- Compute over-expected score in each pair of eigen bins
    BiocParallel::bpprogressbar(BPPARAM) <- TRUE
    df <- BiocParallel::bplapply(
        BPPARAM = BPPARAM, 
        as.character(unique(seqnames(eigens))), 
        function(chr) {.saddle(x[chr], eigens)}
    ) |> dplyr::bind_rows()
    dat <- df |> 
        dplyr::group_by(eigen_bin1, eigen_bin2) |> 
        dplyr::summarize(score = mean(detrended), .groups = 'drop') |> 
        dplyr::mutate(
            x = eigen_bin1 / nbins,
            y = eigen_bin2 / nbins, 
            squished_score = scales::oob_squish(score, limits)
        )

    ## -- Make saddle plot 
    p1 <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) + 
        ggrastr::geom_tile_rast(ggplot2::aes(fill = squished_score)) + 
        ggplot2::scale_y_reverse(labels = scales::percent, expand = c(0, 0)) + 
        ggplot2::scale_x_continuous(
            labels = scales::percent, expand = c(0, 0), 
            position = 'top'
        ) + 
        ggplot2::coord_fixed() + 
        ggplot2::scale_fill_gradientn(
            colors = bgrColors(),
            na.value = "#FFFFFF",
            limits = limits
        ) + 
        ggplot2::labs(
            x = '', 
            y = 'Genomic bins ranked by eigenvector value', 
            fill = "Interaction frequency\n(log2 over expected)"
        ) + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(legend.position = 'bottom') +
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(axis.title.x = ggplot2::element_blank()) 
    if (!plotBins) return(p1)
    
    ## -- Make side barplot 
    dat2 <- tibble::as_tibble(eigens) |> 
        dplyr::group_by(eigen_bin) |> 
        dplyr::summarize(mean_eigen = mean(eigen)) |> 
        dplyr::mutate(x = eigen_bin / nbins)
    p2 <- ggplot2::ggplot(dat2, ggplot2::aes(x = x, y = mean_eigen)) + 
        ggplot2::geom_col() + 
        ggplot2::scale_x_continuous(
            labels = NULL, expand = c(0, 0)
        ) + 
        ggplot2::labs(
            x = '', 
            y = 'Eigen', 
            fill = "Interaction frequency\n(log2 over expected)"
        ) + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(legend.position = 'bottom') +
        ggplot2::theme(panel.border = ggplot2::element_blank()) + 
        ggplot2::theme(panel.grid = ggplot2::element_blank()) + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
        ggplot2::theme(axis.line.x = ggplot2::element_blank()) 
    
    p <- patchwork::wrap_plots(p2, p1, ncol = 1, heights = c(1, 3))

}

.saddle <- function(x_chr, eigens) {
    dx <- detrend(x_chr)
    ints <- interactions(dx)
    an <- anchors(ints)
    ints2 <- ints[IRanges::`%over%`(an[[1]], eigens) & IRanges::`%over%`(an[[2]], eigens)]
    an2 <- anchors(ints2)
    df <- as(ints2, 'data.frame')
    df$eigen_bin1 <- eigens[S4Vectors::subjectHits(findOverlaps(an2[[1]], eigens))]$eigen_bin
    df$eigen_bin2 <- eigens[S4Vectors::subjectHits(findOverlaps(an2[[2]], eigens))]$eigen_bin
    drop_na(df, eigen_bin1, eigen_bin2, detrended) |> 
        dplyr::select(eigen_bin1, eigen_bin2, detrended)
}

################################################################################
#                                                                              #
#                                 ggplot2 extra                                #
#                                                                              #
################################################################################

ggthemeHiContacts <- function(ticks = TRUE, grid = FALSE) {
    t <- ggplot2::theme_bw() +
        ggplot2::theme(text = ggplot2::element_text(size = 8))
    if (ticks) t <- t + ggplot2::theme(axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.2))
    if (grid) {
        t <- t + 
            ggplot2::theme(
                panel.grid.minor = ggplot2::element_line(linewidth = 0.025, colour = "#00000052"), 
                panel.grid.major = ggplot2::element_line(linewidth = 0.025, colour = "#00000052")
            )
            t
    } else {
        t <- t + 
            ggplot2::theme(
                panel.grid.minor = ggplot2::element_blank(), 
                panel.grid.major = ggplot2::element_blank()
            )
            t
    }
}

#' Matrix palettes
#' 
#' @name palettes
#' 
#' @return A vector of colours carefully
#' picked for Hi-C contact heatmaps
#' 
#' @examples
#' bwrColors()
#' bbrColors()
#' bgrColors()
#' afmhotrColors()
#' coolerColors()
#' rainbowColors()
NULL

#' @rdname palettes
#' @export

bwrColors <- function() {
    c("#1659b1", "#4778c2", "#ffffff", "#b13636", "#6C150E")
}

#' @rdname palettes
#' @export

bbrColors <- function() {
    c("#1659b1", "#4778c2", "#a9c3e7", "#ffffff", 
        "#e2adad", "#b13636", "#6C150E"
    )
}

#' @rdname palettes
#' @export

bgrColors <- function() {
    rev(c(
        "#BD202D",
        "#DE614D",
        "#F19374",
        "#F6BA9E",
        "#EAD5CB",
        "#CEDBEB",
        "#AFC5E6",
        "#8FA7D6",
        "#677CBD",
        "#495BA9"
    ))
}

#' @rdname palettes
#' @export

afmhotrColors <- function() {
    c(
        "#ffffff", 
        "#f8f5c3", 
        "#f4ee8d", 
        "#f6be35", 
        "#ee7d32",
        "#c44228", 
        "#821d19", 
        "#381211", 
        "#050606"
    )
}

#' @rdname palettes
#' @export

coolerColors <- function() {
    rev(c(
        "#1A0A10",
        "#7A1128",
        "#B01F29",
        "#D42027",
        "#ED3024",
        "#F15C34",
        "#F78E40",
        "#FAAA4B",
        "#FFCA67",
        "#FDE188",
        "#FFF2A9",
        "#FCF9CE",
        "#FFFEF9"
    ))
}

#' @rdname palettes
#' @export

rainbowColors <- function() {
    rev(c(
        "#FDEE00",
        "#FEC00F",
        "#F37421",
        "#ED1E24",
        "#F26768",
        "#F9B7B8",
        "#FFFCFC",
        "#C6E8ED",
        "#8BD1EC",
        "#56ABDF",
        "#3864AF",
        "#3A57A7",
        "#8D509F"
    ))
}
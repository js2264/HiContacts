#' plotMatrix
#'
#' @param x x
#' @param use.assay use.assay
#' @param scale scale
#' @param limits limits
#' @param dpi dpi
#' @param rasterize rasterize
#' @param symmetrical symmetrical
#' @param chrom_lines chrom_lines
#' @param cmap color map
#'
#' @import ggrastr
#' @import ggplot2
#' @import InteractionSet
#' @import tibble
#' @import dplyr
#' @importFrom GenomicRanges seqnames
#' @export

plotMatrix <- function(x, use.assay = 'balanced', scale = 'log10', loops = NULL, borders = NULL, limits = NULL, dpi = 500, rasterize = TRUE, symmetrical = TRUE, chrom_lines = TRUE, cmap = NULL) {
    `%>%` <- tidyr::`%>%`
    `%over%` <- IRanges::`%over%`
    
    if (!missing(use.assay))
        gis <- assay(x, use.assay)
    else {
        gis <- assay(x)
    }

    ## -- Put metric to plot in `score` column
    has_negative_scores <- any(gis$score < 0, na.rm = TRUE)

    ## -- Scale score
    if (scale == 'log10') {
        gis$score <- log10(gis$score)
    } 
    else if (scale == 'exp0.2') {
        gis$score <- gis$score^0.2
    }

    ## -- Set matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    }
    else {
        an <- anchors(x)
        diff_an <- an[[1]] != an[[2]]
        .scores <- gis$score[diff_an]
        .scores <- .scores[!is.na(.scores)]
        .scores <- .scores[!is.infinite(.scores)]
        M <- max(.scores)
        m <- min(.scores)
        limits <- c(m, M)
    }

    ## -- Choose color map 
    if (is.null(cmap)) {
        if (has_negative_scores) {
            cmap <- bwr_colors
        }
        else {
            cmap <- afmhotr_colors
        }
    }
    
    # -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(raster.dpi = dpi, width = resolution(x), height = resolution(x))
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    ## -- If loops are provided, filter them and add
    if (!is.null(loops)) {
        filtered_loops <- as_tibble(loops[anchors(loops)[[1]] %over% regions(gis) & anchors(loops)[[2]] %over% regions(gis)])
        p_loops <- geom_point(
            data = filtered_loops, 
            aes(x = start1+(end1-start1)/2, y = start2+(end2-start2)/2), 
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
    if (!is.null(borders)) {
        filtered_borders <- as_tibble(borders[borders %over% regions(gis)])
        p_borders <- geom_point(
            data = filtered_borders, 
            aes(x = start+(end-start)/2, y = start+(end-start)/2),
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

    # -- Check number of chromosomes that were extracted
    nseqnames <- length(unique(as.vector(GenomicRanges::seqnames(InteractionSet::anchors(gis)[[1]]))))

    if (nseqnames == 1) { ## Single chromosome coordinates to plot (easy scenario)

        ## -- Convert gis to table and extract x/y
        mat <- gis %>%
            tibble::as_tibble() %>%
            dplyr::mutate(
                x = floor(end1 - (end1 - start1) / 2),
                y = floor(end2 - (end2 - start2) / 2)
            )

        ## -- Clamp scores to limits
        mat <- dplyr::mutate(mat, score = scales::oob_squish(score, c(m, M)))

        ## -- Add lower triangular matrix scores (if symetrical)
        if (symmetrical) {
            mat <- rbind(mat, mat %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% dplyr::select(-x2))
        }

        ## -- Plot matrix
        p <- ggMatrix(mat, cols = cmap, limits = limits) +
            plotFun +
            p_loops + 
            p_borders + 
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
            dplyr::left_join(chroms, by = c(seqnames1 = "seqnames2")) %>%
            dplyr::rename(cumlength_x = cumlength) %>%
            dplyr::left_join(chroms, by = c(seqnames2 = "seqnames2")) %>%
            dplyr::rename(cumlength_y = cumlength) %>%
            dplyr::mutate(
                x = floor(end1 - (end1 - start1) / 2) + cumlength_x,
                y = floor(end2 - (end2 - start2) / 2) + cumlength_y
            )
        mat <- dplyr::mutate(mat, score = ifelse(score > M, M, ifelse(score < m, m, score)))
        mat <- rbind(mat, mat %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% dplyr::select(-x2))

        ## -- Plot matrix
        p <- ggMatrix(mat, cols = cmap, limits = limits) +
            plotFun +
            ggplot2::labs(
                x = "Genome coordinates",
                y = "Genome coordinates"
            )
        if (chrom_lines) {
            p <- p +
                ggplot2::geom_hline(yintercept = chroms$cumlength[-1], colour = "black", alpha = 0.75, size = 0.15) +
                ggplot2::geom_vline(xintercept = chroms$cumlength[-1], colour = "black", alpha = 0.75, size = 0.15)
        }
    }

    p
}

#' ggmatrix
#'
#' @param mat mat
#' @param ticks ticks
#' @param cols cols
#' @param limits limits
#'
#' @import ggplot2
#' @import scales
#' @export

ggMatrix <- function(mat, ticks = TRUE, cols = afmhotr_colors, limits) {
    p <- ggplot2::ggplot(mat, ggplot2::aes(x, y, fill = score))
    p <- p + ggplot2::scale_fill_gradientn(
        colors = cols,
        na.value = "#FFFFFF",
        limits = limits
    ) +
        ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6), position = 'top') +
        ggplot2::scale_y_reverse(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(5, "cm"), barwidth = 0.5, frame.colour = "black")) +
        ggtheme_HiContacts()
    p
}

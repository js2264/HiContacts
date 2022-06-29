#' plotMatrix
#'
#' @param gis gis
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

plotMatrix <- function(gis, colour_by = 'score', scale = 'log10', loops = NULL, borders = NULL, limits = NULL, dpi = 500, rasterize = TRUE, symmetrical = TRUE, chrom_lines = TRUE, cmap = NULL) {
    `%>%` <- tidyr::`%>%`

    ## -- Check that `colour_by` metric is present 
    if (!colour_by %in% colnames(GenomicRanges::mcols(gis))) {
        stop(glue::glue("`{colour_by}` metric not present in the provided GInteractions object. Aborting now."))
    } 

    ## -- Put metric to plot in `score` column
    gis$score <- mcols(gis)[, colour_by]
    has_negative_scores <- any(mcols(gis)[, colour_by] < 0, na.rm = TRUE)

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
        M <- max(gis$score[[1]][diff_an], na.rm = TRUE, is.finite = TRUE)
        m <- min(gis$score[[1]][diff_an], na.rm = TRUE, is.finite = TRUE)
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
        plotFun <- ggrastr::geom_tile_rast(raster.dpi = dpi)
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

#' plotTriangularMatrix
#'
#' @param gis gis
#' @param limits limits
#' @param dist_max dist_max
#' @param cmap cmap
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @export

plotTriangularMatrix <- function(gis, limits = NULL, dist_max = NULL, cmap = afmhotr_colors) {
    `%>%` <- tidyr::`%>%`

    ## -- Convert gis to table and extract x/y
    mat <- gis %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
            x = floor(end1 - (end1 - start1) / 2),
            y = floor(end2 - (end2 - start2) / 2)
        )

    ## -- Truncate to the max distance
    binsize <- unique(mat$start2 - mat$start1)[2]
    if (is.null(dist_max)) {
        dist_max <- diff(range(mat$start2 - mat$start1)) * 0.2
    }
    else {
        if (dist_max > diff(range(mat$x))) dist_max <- diff(range(mat$x))
        if (dist_max < 1) {
            dist_max <- diff(range(mat$start2 - mat$start1)) * dist_max
        }
    }
    max_bin_dist <- floor(dist_max / binsize)
    mat <- dplyr::filter(mat, abs(bin2 - bin1) <= max_bin_dist)

    ## -- Matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    } else {
        M <- max(mat$score, na.rm = TRUE)
        m <- min(mat$score, na.rm = TRUE)
        limits <- c(m, M)
    }

    ## -- Clamp scores to limits
    mat <- dplyr::mutate(mat, score = ifelse(score > M, M, ifelse(score < m, m, score)))

    ## -- Add lower triangular matrix scores
    mat <- rbind(mat, mat %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% dplyr::select(-x2))

    ## -- Get coordinates for rotated squares
    mat_ <- mat %>%
        ## -- Filter upper diagonal
        dplyr::filter(y >= x) %>%
        ## -- Recompute center of each bin and distance to the diagonal
        dplyr::mutate(
            x_ = x + (y - x) / 2,
            dist = (y - x) / {
                sqrt(2) * sqrt(2)
            },
            ID = 1:dplyr::n()
        ) %>%
        ## -- Group by bins
        dplyr::group_by(ID) %>%
        ## -- Compute each corner of each of the 45deg-tilted bins
        dplyr::mutate(
            A_x = x_ - binsize / 2,
            A_y = dist,
            B_x = x_,
            B_y = dist + binsize / 2,
            C_x = x_ + binsize / 2,
            C_y = dist,
            D_x = x_,
            D_y = dist - binsize / 2
        ) %>%
        tidyr::pivot_longer(cols = dplyr::starts_with(c("A_", "B_", "C_", "D_")), names_to = "pt", values_to = "value_coord") %>%
        tidyr::extract(pt, c("pt", "coord_x"), "(.)_(.)") %>%
        tidyr::spread(coord_x, value_coord) %>%
        dplyr::mutate(y = ifelse(y < 0, 0, ifelse(y > max_bin_dist * binsize / 2, max_bin_dist * binsize / 2, y)))

    ## -- Compute triangular truncated border
    M_y <- max(mat_$y)
    uptri <- tibble::tibble(
        x = c(
            min(mat_$x),
            min(mat_$x) + max_bin_dist * binsize / 2,
            max(mat_$x) - max_bin_dist * binsize / 2,
            max(mat_$x)
        ),
        y = c(0, M_y, M_y, 0)
    )

    ## -- Plot matrix
    p <- ggtiltedmatrix(mat_, cols = cmap, limits = limits) +
        ggplot2::geom_polygon() +
        ggplot2::labs(
            x = unique(mat$seqnames1),
            y = ""
        ) +
        ggplot2::geom_polygon(
            data = uptri,
            ggplot2::aes(x, y),
            col = "black",
            fill = NA,
            size = 0.2,
            inherit.aes = FALSE
        )

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
        ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::scale_y_reverse(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(5, "cm"), barwidth = 0.5, frame.colour = "black")) +
        ggtheme_HiContacts()
    p
}

#' APA
#'
#' @param file file
#' @param coords coords
#' @param res res
#' @param limits limits
#' @param dpi dpi
#' @param rasterize rasterize
#' @param symmetrical symmetrical
#' @param BPPARAM BPPARAM
#' @param scale scale
#' @param cmap cmap
#'
#' @import InteractionSet
#' @import ggrastr
#' @import ggplot2
#' @import BiocParallel
#' @import zeallot
#' @import tibble
#' @import dplyr
#' @importFrom GenomicRanges width
#' @export

APA <- function(file, coords, res = NULL, limits = NULL, dpi = 500, rasterize = TRUE, symmetrical = TRUE, BPPARAM = BiocParallel::bpparam(), scale = FALSE, cmap = afmhotr_colors) {
    `%>%` <- tidyr::`%>%`

    ## -- Coerce loops into coords and coords2
    if (class(coords) == "GInteractions") {
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
        breaks <- seq(
            {
                -wi / 2
            } - res0,
            {
                wi / 2
            } + res0,
            length.out = wi / res0 + 3
        )
    }

    ## -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(raster.dpi = dpi)
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    if (!isLoops) {
        
        stop("Loop APA plots are not supported yet.")
    
    }
    else {
        ## -- Load entire matrix in memory
        gis <- cool2gi(file, res = res, coords = NULL)

        ## -- Calculate scores over each coords
        nc <- length(coords)
        br <- seq(1, nc+100, by = 100)
        coords_split <- coords %>%
            plyranges::mutate(
                chunk = as.numeric(cut(seq(1, nc), breaks = br, include.lowest = TRUE))
            )
        mats <- BiocParallel::bplapply(BPPARAM = BPPARAM, unique(coords_split$chunk), function(K) {
            coords_sub <- coords_split[coords_split$chunk == K]
            coords_sub %>% 
                plyranges::mutate(ID = seq(1, length(.))) %>% 
                plyranges::select(-score) %>%
                plyranges::join_overlap_left(gis) %>% 
                plyranges::group_by(ID) %>% 
                tibble::as_tibble() %>% 
                dplyr::left_join(
                    tibble::as_tibble(anchors) %>% 
                        dplyr::mutate(bin1 = seq(1, length(anchors)), center1 = end - (end-start)/2) %>% 
                        dplyr::select(bin1, center1),
                    by = 'bin1'
                ) %>% 
                dplyr::left_join(
                    tibble::as_tibble(anchors) %>% 
                        dplyr::mutate(bin2 = seq(1, length(anchors)), center2 = end - (end-start)/2) %>% 
                        dplyr::select(bin2, center2),
                    by = 'bin2'
                ) %>% 
                dplyr::mutate(
                    x = center1,
                    y = center2,
                    x = breaks[as.numeric(cut(x - midpoint, breaks, include.lowest = TRUE)) + 1],
                    y = breaks[as.numeric(cut(y - midpoint2, breaks, include.lowest = TRUE)) + 1]
                ) %>% 
                dplyr::group_by(x, y) %>%
                dplyr::filter(
                    x >= min(breaks[2:(length(breaks) - 1)]), 
                    x <= max(breaks[2:(length(breaks) - 1)]), 
                    y >= min(breaks[2:(length(breaks) - 1)]), 
                    y <= max(breaks[2:(length(breaks) - 1)])
                ) %>%
                dplyr::summarize(score = mean(score, na.rm = TRUE), n = n()) %>%
                dplyr::mutate(seqnames1 = "Aggr. coordinates")
        }) %>%
            dplyr::bind_rows()
        
        ## -- Calculate scores over each coords
        mats <- mats %>% 
            dplyr::group_by(x, y) %>%
            dplyr::summarize(score = mean(score, na.rm = TRUE)) %>%
            dplyr::mutate(seqnames1 = "Aggr. coordinates")

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
        mats <- dplyr::mutate(mats, score = ifelse(score > M, M, ifelse(score < m, m, score)))

        ## -- Add lower triangular matrix scores
        if (symmetrical) {
            mats <- rbind(mats, mats %>% dplyr::mutate(x2 = y, y = x, x = x2) %>%
                dplyr::select(-x2))
        }

        ## -- Plot matrix
        p <- ggMatrix(mats, cols = cmap, limits = limits) +
            plotFun +
            ggplot2::labs(
                x = unique(mats$seqnames1),
                y = unique(mats$seqnames1)
            )

        p
    }
}

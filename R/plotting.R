################################################################################
#                                                                              #
#                                 Matrix                                       #
#                                                                              #
################################################################################

#' Plotting a contact matrix
#'
#' @rdname contacts-plot
#' 
#' @param x x
#' @param use.scores use.scores
#' @param scale scale
#' @param limits limits
#' @param loops loops
#' @param borders borders
#' @param dpi dpi
#' @param rasterize rasterize
#' @param symmetrical symmetrical
#' @param chrom_lines chrom_lines
#' @param cmap color map
#' @return ggplot
#'
#' @import ggplot2
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
#' @importFrom InteractionSet anchors
#' @importFrom scales oob_squish
#' @export
#' @examples 
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' plotMatrix(
#'     contacts_yeast, 
#'     use.scores = 'balanced', 
#'     scale = 'log10', 
#'     limits = c(-4, -1)
#' )

plotMatrix <- function(
    x, 
    use.scores = NULL, 
    scale = 'log10', 
    loops = NULL, 
    borders = NULL, 
    limits = NULL, 
    dpi = 500, 
    rasterize = TRUE, 
    symmetrical = TRUE, 
    chrom_lines = TRUE, 
    cmap = NULL
) {
    `%over%` <- IRanges::`%over%`
    
    ## -- Extract scores
    if (!is.null(use.scores)) {
        gis <- interactions(x)
        gis$score <- scores(x, use.scores)
    }
    else {
        if ("balanced" %in% names(scores(x))) {
            gis <- interactions(x)
            gis$score <- scores(x, "balanced")
        } 
        else {
            gis <- interactions(x)
            gis$score <- scores(x, 1)
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
        an <- anchors(x)
        diff_an <- an[['first']] != an[['second']]
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
            cmap <- bwrColors()
        }
        else {
            cmap <- afmhotrColors()
        }
    }
    
    # -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(
            raster.dpi = dpi, width = resolution(x), height = resolution(x)
        )
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    ## -- If loops are provided, filter them and add
    if (!is.null(loops)) {
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
    if (!is.null(borders)) {
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

    # -- Check number of chromosomes that were extracted
    nseqnames <- length(unique(as.vector(
        GenomicRanges::seqnames(InteractionSet::anchors(gis)[['first']])
    )))

    if (nseqnames == 1) { ## Single chromosome coordinates to plot

        ## -- Convert gis to table and extract x/y
        mat <- gis |>
            tibble::as_tibble() |>
            dplyr::mutate(
                x = floor(end1 - (end1 - start1) / 2),
                y = floor(end2 - (end2 - start2) / 2)
            ) |> 
            drop_na(score)

        ## -- Clamp scores to limits
        mat <- dplyr::mutate(mat, score = scales::oob_squish(score, c(m, M)))

        ## -- Add lower triangular matrix scores (if symetrical)
        if (symmetrical) {
            mat <- rbind(
                mat, 
                mat |> 
                    dplyr::mutate(x2 = y, y = x, x = x2) |> 
                    dplyr::select(-x2)
            )
            if (!is_symmetrical(x)) {
                coords <- unlist(S4Vectors::zipup(char2pair(focus(x))))
                mat <- mat |> 
                    filter(x >= GenomicRanges::start(coords[1]) & 
                        x <= GenomicRanges::end(coords[1])) |> 
                    filter(y >= GenomicRanges::start(coords[2]) & 
                        y <= GenomicRanges::end(coords[2]))
            }
        } 

        ## -- Plot matrix
        p <- ggMatrix(mat, cols = cmap, limits = limits) +
            plotFun +
            p_loops + 
            p_borders + 
            ggplot2::labs(
                x = unique(mat$seqnames1),
                y = "Genome coordinates", 
                caption = paste(
                    sep = '\n',
                    paste0('file: ', fileName(x)), 
                    paste0('res: ', resolution(x))
                )
            )
    }

    else {
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
        p <- ggMatrix(mat, cols = cmap, limits = limits) +
            plotFun +
            ggplot2::labs(
                x = "Genome coordinates",
                y = "Genome coordinates", 
                caption = paste(
                    sep = '\n',
                    paste0('file: ', fileName(x)), 
                    paste0('res: ', resolution(x))
                )
            )
        if (chrom_lines) {
            p <- p +
                ggplot2::geom_hline(
                    yintercept = chroms$cumlength[-1], 
                    colour = "black", alpha = 0.75, size = 0.15
                ) +
                ggplot2::geom_vline(
                    xintercept = chroms$cumlength[-1], 
                    colour = "black", alpha = 0.75, size = 0.15
                )
        }
    }

    p
}

#' @rdname contacts-plot
#' 
#' @param mat mat
#' @param ticks ticks
#' @param cols cols
#' @param limits limits
#' @return ggplot
#'
#' @import ggplot2
#' @importFrom scales unit_format

ggMatrix <- function(mat, ticks = TRUE, cols = afmhotrColors(), limits) {
    p <- ggplot2::ggplot(mat, ggplot2::aes(x, y, fill = score))
    p <- p + ggplot2::scale_fill_gradientn(
        colors = cols,
        na.value = "#FFFFFF",
        limits = limits
    ) +
        ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6), position = 'top') +
        ggplot2::scale_y_reverse(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(5, "cm"), barwidth = 0.5, frame.colour = "black")) + 
        coord_fixed() +
        ggthemeHiContacts()
    p
}

################################################################################
#                                                                              #
#                                 P(s)                                         #
#                                                                              #
################################################################################

#' Plotting a P(s) distance law
#' 
#' @rdname Ps-plot
#' 
#' @param ... ...
#' @param xlim xlim
#' @param ylim ylim
#' @return ggplot
#'
#' @import ggplot2
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @export
#' @examples 
#' ## Single P(s)
#' 
#' contacts_yeast <- contacts_yeast()
#' ps <- getPs(contacts_yeast)
#' plotPs(ps, ggplot2::aes(x = binned_distance, y = norm_p))
#' 
#' ## Comparing several P(s)
#' 
#' contacts_yeast <- contacts_yeast()
#' contacts_yeast_eco1 <- contacts_yeast_eco1()
#' ps_wt <- getPs(contacts_yeast)
#' ps_wt$sample <- 'WT'
#' ps_eco1 <- getPs(contacts_yeast_eco1)
#' ps_eco1$sample <- 'eco1'
#' ps <- rbind(ps_wt, ps_eco1)
#' plotPs(ps, ggplot2::aes(x = binned_distance, y = norm_p, group = sample, color = sample))

plotPs <- function(..., xlim = c(5000, 4.99e5), ylim = c(1e-8, 1e-4)) {
    gg <- ggplot2::ggplot(...) + 
        ggplot2::geom_line() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_y_log10(
            limits = ylim, 
            expand = c(0, 0), 
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        ggplot2::scale_x_log10(
            limits = xlim, 
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        ) + 
        ggplot2::annotation_logticks() + 
        ggplot2::labs(x = "Genomic distance", y = "Contact frequency")
    gg
}

#' @rdname Ps-plot
#' 
#' @return ggplot
#' 
#' @import ggplot2
#' @export
#' @examples 
#' plotPsSlope(ps, ggplot2::aes(x = binned_distance, y = slope))

plotPsSlope <- function(..., xlim = c(5000, 4.99e5), ylim = c(-3, 0)) {
    gg <- ggplot2::ggplot(...) + 
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
#' @rdname virtual4C-plot
#' 
#' @param x GRanges, generally the output of `virtual4C()`
#' @param mapping aes to pass on to ggplot2
#' @return ggplot
#' 
#' @import ggplot2
#' @import tibble
#' @importFrom scales unit_format
#' @export
#' @examples 
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' v4C <- virtual4C(contacts_yeast, GenomicRanges::GRanges('II:490000-510000'))
#' plot4C(v4C, ggplot2::aes(x = center, y = score))

plot4C <- function(x, mapping) {
    x <- tibble::as_tibble(x)
    p <- ggplot2::ggplot(x, mapping) + 
        ggplot2::geom_line() + 
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
#                                 ggplot2 extra                                #
#                                                                              #
################################################################################

#' ggplot2-related functions
#' 
#' @rdname ggplot2-extra
#'
#' @param ticks ticks
#' @return a custom ggplot2 theme
#' 
#' @import ggplot2

ggthemeHiContacts <- function(ticks = TRUE) {
    t <- ggplot2::theme_bw() +
        ggplot2::theme(
            text = ggplot2::element_text(size = 8),
            panel.grid.minor = ggplot2::element_line(size = 0.025, colour = "#00000052"),
            panel.grid.major = ggplot2::element_line(size = 0.05, colour = "#00000052")
        )
    if (ticks) t <- t + ggplot2::theme(axis.ticks = ggplot2::element_line(colour = "black", size = 0.2))
    t
}

#' Matrix palettes
#'
#' @return ggplot
#'
#' @rdname palettes
#' 
#' @import ggplot2
#' @importFrom scales unit_format
#' @export
#' @examples
#' bwrColors()

bwrColors <- function() {
    c("#1659b1", "#4778c2", "#ffffff", "#b13636", "#6C150E")
}

#' @rdname palettes
#' 
#' @export
#' @examples
#' afmhotrColors()

afmhotrColors <- function() {
    c("#ffffff", "#f8f5c3", "#f4ee8d", "#f6be35", "#ee7d32",
        "#c44228", "#821d19", "#381211", "#050606"
    )
}

#' @rdname palettes
#' 
#' @export
#' @examples
#' bbrColors()

bbrColors <- function() {
    c("#1659b1", "#4778c2", "#a9c3e7", "#ffffff", 
        "#e2adad", "#b13636", "#6C150E"
    )
}

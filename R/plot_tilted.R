gi2mat <- function(gis, limits = NULL, truncate_tip = 0.5) {

    `%>%` <- magrittr::`%>%`

    ## -- Convert gis to table and extract x/y
    mat <- gis %>% 
        tibble::as_tibble() %>%
        dplyr::mutate(
            x = floor(end1 - (end1 - start1)/2), 
            y = floor(end2 - (end2 - start2)/2)
        ) 

    # -- Define plotting approach
    binsize <- GenomicRanges::width(InteractionSet::regions(gis)[1])
    max_bin_dist <- (max(mat$bin2 - mat$bin1) + 1) * truncate_tip

    ## -- Matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    } else {
        M <- max(gis$score[abs(gis$bin1 - gis$bin2) >= 0 & abs(gis$bin1 - gis$bin2) <= max_bin_dist], na.rm = TRUE)
        m <- min(gis$score[abs(gis$bin1 - gis$bin2) >= 0 & abs(gis$bin1 - gis$bin2) <= max_bin_dist], na.rm = TRUE)
        limits <- c(m, M)
    }

    ## -- Clamp scores to limits
    mat <- dplyr::mutate(mat, score = ifelse(score > M, M, ifelse(score < m, m, score)))

    ## -- Add lower triangular matrix scores
    mat <- rbind(mat, mat %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% dplyr::select(-x2))

    ## -- Filter matrix to cut out the tip of the triangle
    mat <- mat[abs(mat$bin2 - mat$bin1) <= max_bin_dist, ]
    
    mat_ <- mat %>% 
        ## -- Filter upper diagonal
        dplyr::filter(y >= x) %>% 
        ## -- Recompute center of each bin and distance to the diagonal
        dplyr::mutate(
            x_ = x + (y - x)/2, 
            dist = (y - x) / {sqrt(2)*sqrt(2)}, 
            ID = 1:dplyr::n()
        ) %>% 
        ## -- Group by bins
        dplyr::group_by(ID) %>% 
        ## -- Compute each corner of each of the 45deg-tilted bins
        dplyr::mutate(
            A_x = x_ - binsize/2,
            A_y = dist, 
            B_x = x_,
            B_y = dist + binsize/2, 
            C_x = x_ + binsize/2,
            C_y = dist, 
            D_x = x_,
            D_y = dist - binsize/2
        ) %>% 
        tidyr::pivot_longer(cols = dplyr::starts_with(c('A_', 'B_', 'C_', 'D_')), names_to = 'pt', values_to = 'value_coord') %>%
        tidyr::extract(pt, c("pt", "coord_x"), "(.)_(.)") %>%
        tidyr::spread(coord_x, value_coord) %>% 
        dplyr::mutate(y = ifelse(y < 0, 0, y)) %>% 
        ## -- Rescale y for proper scaling when plotting
        dplyr::mutate(y = y * { truncate_tip / { sqrt(2) * sqrt(2) } })

}

plotMatrixList <- function(ls, ...) {

    `%>%` <- magrittr::`%>%`

    ## -- Convert gis to table and extract x/y, for each element in `ls`
    mat <- ls[[2]] %>% 
        tibble::as_tibble() %>%
        dplyr::mutate(
            x = floor(end1 - (end1 - start1)/2), 
            y = floor(end2 - (end2 - start2)/2)
        ) 
    max_bin_dist <- (max(mat$bin2 - mat$bin1) + 1) * truncate_tip
    mat_ <- lapply(names(ls), function(n) {
        gis <- ls[[n]]
        m <- gi2mat(gis) %>% dplyr::mutate(item = n)
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

ggtiltedmatrix <- function(mat_, ticks = TRUE, cols = afmhotr_colors, truncate_tip, nmatrices = 1) {
   
    p <- ggplot2::ggplot(mat_, ggplot2::aes(
        x, y,
        group = ID, 
        fill = score
    ))

    p <- p + ggplot2::scale_fill_gradientn(
        colors = cols, 
        na.value = '#FFFFFF'
    ) + 
        ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) + 
        ggplot2::scale_y_continuous(expand = c(0, 0)) + 
        ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(5, 'cm'), barwidth = 0.5, frame.colour = "black")) + 
        ggtheme_coolerr() + 
        ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(), 
            aspect.ratio = 1 / ( sqrt(2)*sqrt(2) / truncate_tip / {nmatrices/sqrt(2)})
        )

    p
}

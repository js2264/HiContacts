plotMatrix <- function(gis, limits = NULL, dpi = 500) {

    ## -- Format matrix
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    } else {
        M <- log10(max(gis$score[abs(gis$bin1 - gis$bin2) >= 0], na.rm = TRUE))
        m <- log10(min(gis$score[abs(gis$bin1 - gis$bin2) >= 0], na.rm = TRUE))
        limits <- c(m, M)
    }

    mat <- gis %>% 
        as_tibble() %>%
        mutate(
            x = floor(end1 - (end1 - start1)/2), 
            y = floor(end2 - (end2 - start2)/2), 
            score = log10(score)
        ) 
    mat <- mutate(mat, score = ifelse(score > M, M, ifelse(score < m, m, score)))
    mat <- rbind(mat, mat %>% mutate(x2 = y, y = x, x = x2) %>% select(-x2))

    ## -- Plot matrix
    p <- ggplot(mat, aes(x, y, fill = score)) + 
        ggrastr::geom_tile_rast(raster.dpi = dpi) + 
        # geom_tile() + 
        scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) + 
        scale_y_reverse(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) + 
        guides(fill = guide_colorbar(barheight = unit(5, 'cm'), barwidth = 0.5, frame.colour = "black")) + 
        ggtheme_matrix() + 
        labs(
            x = unique(mat$seqnames1),
            y = unique(mat$seqnames1)
        ) + scale_fill_gradientn(
            colors = c('#FFFFF2', '#FEFBE0', '#FCF6BD', '#F9F198', '#FFF073', '#FECC50', '#F6A32C', '#F0801A', '#DD5D12', '#B83917', '#6C150E', '#430F11', '#1D0809', '#000000'), 
            na.value = '#FFFFFF', 
            limits = c(m, M)
        )

}

ggtheme_matrix <- function() {
    theme_bw() + 
    theme(
        text = element_text(size=8), 
        panel.grid.minor = element_blank(), 
        aspect.ratio = 1, 
        panel.grid.major = element_line(size = 0.05, colour = '#00000042'), 
        axis.ticks = element_line(colour = "black", size = 0.25), 
        panel.background = element_rect(fill = NA),
        panel.ontop = TRUE
    )
}
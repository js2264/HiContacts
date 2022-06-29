getScalogram <- function(gis, nbins = 250, ylim = c(1e5, 1e8), probs = c(0, 0.15, 0.3, 0.45, 0.6, 0.75, 1)) {
    df <- tibble(
        dist = gis$distance, 
        chr = as.vector(seqnames(anchors(gis)[[1]])), 
        pos1 = start(anchors(gis)[[1]]),
        pos2 = start(anchors(gis)[[2]])
    ) %>% 
        drop_na() %>% 
        pivot_longer(-c(chr, dist), names_to = 'pos_', values_to = 'pos') %>%
        select(-pos_) %>%
        filter(dist > 0, pos > 1) %>% 
        mutate(binned_pos = {seq(1, max(pos), length.out = nbins)}[findInterval(pos, vec = seq(1, max(pos), length.out = nbins))]) %>% 
        group_by(chr, binned_pos) %>% 
        group_modify(~ {
            quantile(.x$dist, probs = probs) %>%
                tibble::enframe(name = "prob", value = "dist_quantile") %>% 
                mutate(prob = str_replace(prob, '%', '') %>% as.numeric() %>% `/`(100) %>% factor(levels = rev(probs)))
        }) %>%
        mutate(dist_quantile = scales::oob_squish(dist_quantile, c(ylim[[1]], ylim[[2]]))) %>% 
        mutate(y0 = ylim[[1]])
    return(df)
}

ggScalogram <- function(..., ylim = c(1e4, 1e7)) {
    gg <- ggplot(...) + 
        geom_ribbon() + 
        theme_minimal() + 
        theme(panel.border = element_rect(fill = NA)) + 
        theme(panel.grid.minor = element_blank()) +
        scale_fill_brewer(palette = 'Spectral') + 
        scale_x_continuous(
            expand = c(0, 0)
        ) + 
        scale_y_log10(
            expand = c(0, 0), 
            limits = c(ylim[[1]], ylim[[2]])
        ) + 
        annotation_logticks(sides = 'l') + 
        labs(x = 'Position along chr.', y = 'Distance of interactions') + 
        facet_wrap(~chr, scales = 'free')
    gg
}

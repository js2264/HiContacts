#' getScalogram
#' 
#' In development.
#' 
#' @param x x
#' @param nbins nbins
#' @param ylim ylim
#' @param probs probs
#' @param aes aes
#' @return a tibble
#' 
#' @importFrom dplyr group_modify
#' @rdname scalograms

getScalogram <- function(x, nbins = 250, ylim = c(1e5, 1e8), probs = c(0, 0.15, 0.3, 0.45, 0.6, 0.75, 1)) {
    pairsFile <- pairsFile(x)
    if (is.null(pairsFile)) {
        stop("Please provide a pairsFile for `x`. Aborting now.")
        # message("pairsFile not specified. The P(s) curve will be an approximation.")
        # pairs <- scores(x, 'raw')
        # scalo <- FALSE
        # return(scalo)
    }
    else {
        message("Importing pairs file ", pairsFile, " in memory. This may take a while...")
        pairs <- pairs2gi(pairsFile)
        an <- anchors(pairs)
        scalo <- tibble(
            dist = pairs$distance, 
            chr = as.vector(seqnames(an[[1]])), 
            pos1 = start(an[[1]]),
            pos2 = start(an[[2]])
        ) %>% 
            drop_na() %>% 
            pivot_longer(-c(chr, dist), names_to = 'pos_', values_to = 'pos') %>%
            select(-pos_) %>%
            filter(dist > 0, pos > 1) %>% 
            mutate(binned_pos = {seq(1, max(pos, na.rm = TRUE), length.out = nbins)}[findInterval(pos, vec = seq(1, max(pos, na.rm = TRUE), length.out = nbins))]) %>% 
            group_by(chr, binned_pos) %>% 
            dplyr::group_modify(~ {
                quantile(.x$dist, probs = probs) %>%
                    tibble::enframe(name = "prob", value = "dist_quantile") %>% 
                    mutate(prob = str_replace(prob, '%', '') %>% 
                        as.numeric() %>% 
                        `/`(100) %>% 
                        factor(levels = rev(probs))
                    )
            }) 
        return(scalo)
    }
}

#' ggScalogram
#' 
#' In development.
#' 
#' @return ggplot
#' 
#' @import ggplot2
#' @importFrom scales oob_squish
#' @rdname scalograms

ggScalogram <- function(x, aes, ylim = c(1e4, 1e7)) {
    x <- x %>%
        filter(chr == 'chr3') %>% 
        mutate(dist_quantile = scales::oob_squish(dist_quantile, c(ylim[[1]], ylim[[2]]))) %>% 
        mutate(y0 = ylim[[1]])
    gg <- ggplot(x, aes(x = binned_pos, ymin = y0, ymax = dist_quantile, group = prob, fill = prob)) + 
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

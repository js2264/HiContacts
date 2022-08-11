#' getPs
#'
#' @param pairs pairs in GenomicInteractions format, e.g. imported from a `.pairs` file with `pair2gi()`
#' @param by_chr by_chr
#' @param filtered_chr filtered_chr
#'
#' @import tibble
#' @import dplyr
#' @export

getPs <- function(pairs, by_chr = FALSE, filtered_chr = c('XII', 'chrXII', 'chr12', '12', 'Mito', 'MT', 'chrM')) {
    df <- tibble::tibble(
        chr = as.vector(GenomeInfoDb::seqnames(InteractionSet::anchors(pairs)[[1]])),
        distance = pairs$distance
    ) %>% 
        tidyr::drop_na() %>% 
        dplyr::filter(!chr %in% filtered_chr) %>% 
        dplyr::mutate(binned_distance = HiContacts::breaks$break_pos[findInterval(distance, vec = HiContacts::breaks$break_pos)])
    if (by_chr) {
        df <- dplyr::group_by(df, chr, binned_distance)
    } 
    else {
        df <- dplyr::group_by(df, binned_distance)
    }
    d <- dplyr::tally(df, name = 'ninter') %>%
        dplyr::mutate(p = ninter/sum(ninter)) %>% 
        dplyr::left_join(HiContacts::breaks, by = c('binned_distance' = 'break_pos')) %>% 
        dplyr::mutate(norm_p = p / binwidth)
    if (by_chr) {
        d <- dplyr::group_by(d, chr)
    } 
    else {
        d <- d
    }
    ps <- dplyr::group_split(d) %>% 
        purrr::map(function(x) {
            dplyr::mutate(x, norm_p_unity = norm_p / {dplyr::slice(x, which.min(abs(x$binned_distance - 100000))) %>% dplyr::pull(norm_p)}) %>% 
            dplyr::mutate(slope = (log10(lead(norm_p)) - log10(norm_p)) / (log10(lead(binned_distance)) - log10(binned_distance))) %>% 
            dplyr::mutate(slope = c(0, predict(loess(slope ~ binned_distance, span = 0.5, data = .))))
        }) %>% 
        dplyr::bind_rows()
    if (by_chr) {
        ps <- dplyr::select(ps, chr, binned_distance, p, norm_p, norm_p_unity, slope)
    } 
    else {
        ps <- dplyr::select(ps, binned_distance, p, norm_p, norm_p_unity, slope)
    }
    return(ps)
}

#' ggPs
#'
#' @param ... ...
#' @param xlim xlim
#' @param ylim ylim
#'
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @export

ggPs <- function(..., xlim = c(5000, 4.99e5), ylim = c(1e-8, 1e-4)) {
    gg <- ggplot(...) + 
        geom_line() + 
        theme_minimal() + 
        theme(panel.border = element_rect(fill = NA)) + 
        theme(panel.grid.minor = element_blank()) +
        scale_y_log10(
            limits = ylim, 
            expand = c(0, 0), 
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_x_log10(
            limits = xlim, 
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        ) + 
        annotation_logticks() + 
        labs(x = "Genomic distance", y = "Contact frequency")
    gg
}

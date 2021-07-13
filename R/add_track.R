#' @import IRanges
#' @import zeallot
#' @import tibble
#' @import dplyr
#' @import GenomicRanges
#' @import ggplot2
#' @import scales
#' @import egg
#' @export

addTracks <- function(p, range, annotations = NULL, profiles = NULL) {
    `%over%` <- IRanges::`%over%`
    `%<-%` <- zeallot::`%<-%`
    c(coords_chr, coords_start, coords_end) %<-% splitCoords(range)

    ## -- Get annotations
    p_annotations <- lapply(names(annotations), function(n) {
        gr <- annotations[[n]]
        gr <- gr[gr %over% GenomicRanges::GRanges(range)]
        df <- tibble::as_tibble(gr) %>%
            dplyr::mutate(start = ifelse(start < coords_start, coords_start, start)) %>%
            dplyr::mutate(end = ifelse(end > coords_end, coords_end, end))
        ggplot2::ggplot(df) +
            ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 1), col = "#444444", size = 0.2) +
            ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
            ggplot2::scale_y_continuous(expand = c(0, 0)) +
            ggtheme_coolerr_tracks()
    })

    ## -- Get annotations
    p_profiles <- lapply(names(profiles), function(n) {
        tr <- profiles[[n]]
        tr <- tr[GenomicRanges::GRanges(range)][[1]]
        tr <- tibble::tibble(x = coords_start + cumsum(tr@lengths) - 1, y = tr@values) %>%
            dplyr::filter(x >= coords_start, x <= coords_end)
        ggplot2::ggplot(tr) +
            ggplot2::geom_path(ggplot2::aes(x, y)) +
            ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
            ggplot2::scale_y_continuous(expand = c(0, 0)) +
            ggtheme_coolerr_tracks()
    })

    ## -- Combine tracks and heatmap
    p_final <- egg::ggarrange(
        plots = c(list(p), p_annotations, p_profiles),
        draw = FALSE,
        ncol = 1,
        heights = c(1, rep(0.1, length(annotations) + length(profiles)))
    )

    return(p_final)
}

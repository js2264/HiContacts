#' ggmatrix
#'
#' @param mat mat
#' @param ticks ticks
#' @param cols cols
#'
#' @import ggplot2
#' @import scales
#' @export

ggmatrix <- function(mat, ticks = TRUE, cols = afmhotr_colors) {
    p <- ggplot2::ggplot(mat, ggplot2::aes(x, y, fill = score))
    p <- p + ggplot2::scale_fill_gradientn(
        colors = cols,
        na.value = "#FFFFFF"
    ) +
        ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::scale_y_reverse(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(5, "cm"), barwidth = 0.5, frame.colour = "black")) +
        ggtheme_coolerr()
    p
}

#' ggtiltedmatrix
#'
#' @param mat_ mat_
#' @param ticks ticks
#' @param cols cols
#' @param truncate_tip truncate_tip
#' @param nmatrices nmatrices
#'
#' @import ggplot2
#' @import scales
#' @export

ggtiltedmatrix <- function(mat_, ticks = TRUE, cols = afmhotr_colors, truncate_tip, nmatrices = 1) {
    p <- ggplot2::ggplot(mat_, ggplot2::aes(
        x, y,
        group = ID,
        fill = score
    ))

    p <- p + ggplot2::scale_fill_gradientn(
        colors = cols,
        na.value = "#FFFFFF"
    ) +
        ggplot2::scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "M", scale = 1e-6)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barheight = ggplot2::unit(5, "cm"), barwidth = 0.5, frame.colour = "black")) +
        ggtheme_coolerr() +
        ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            aspect.ratio = 1 / (sqrt(2) * sqrt(2) / truncate_tip / {
                nmatrices / sqrt(2)
            })
        )

    p
}

#' ggtheme_coolerr
#'
#' @param ticks ticks
#'
#' @import ggplot2
#' @export

ggtheme_coolerr <- function(ticks = TRUE) {
    t <- ggplot2::theme_bw() +
        ggplot2::theme(
            text = ggplot2::element_text(size = 8),
            panel.background = ggplot2::element_rect(fill = NA),
            panel.ontop = FALSE,
            panel.grid.minor = ggplot2::element_line(size = 0.025, colour = "#00000052"),
            panel.grid.major = ggplot2::element_line(size = 0.05, colour = "#00000052"),
            aspect.ratio = 1
        )
    if (ticks) t <- t + ggplot2::theme(axis.ticks = ggplot2::element_line(colour = "black", size = 0.2))
    t
}

#' ggtheme_coolerr_tracks
#'
#' @import ggplot2
#' @export

ggtheme_coolerr_tracks <- function() {
    ggplot2::theme(
        text = ggplot2::element_text(size = 8),
        panel.background = ggplot2::element_rect(fill = NA),
        panel.ontop = FALSE,
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
    )
}

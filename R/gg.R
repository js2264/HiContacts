#' ggtheme_HiContacts
#'
#' @param ticks ticks
#' @rdname ggplot2-HiContacts
#' 
#' @import ggplot2
#' @export

ggtheme_HiContacts <- function(ticks = TRUE) {
    t <- ggplot2::theme_bw() +
        ggplot2::theme(
            text = ggplot2::element_text(size = 8),
            panel.grid.minor = ggplot2::element_line(size = 0.025, colour = "#00000052"),
            panel.grid.major = ggplot2::element_line(size = 0.05, colour = "#00000052")
        )
    if (ticks) t <- t + ggplot2::theme(axis.ticks = ggplot2::element_line(colour = "black", size = 0.2))
    t
}

#' @rdname ggplot2-HiContacts
#'
#' @import ggplot2
#' @export

ggtheme_HiContacts_tracks <- function() {
    ggplot2::theme(
        text = ggplot2::element_text(size = 8),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
    )
}

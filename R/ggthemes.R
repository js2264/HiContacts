#' ggplot2-related functions
#' 
#' @rdname ggplot2-HiContacts
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

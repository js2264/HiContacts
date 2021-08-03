#' addTads
#'
#' @param p p
#' @param tads tads
#' @param coords coords
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRanges
#' @import tibble
#' @import dplyr
#' @import purrr
#' @import ggplot2
#' @export

addTads <- function(p, tads, coords) {

    if (is.character(coords)) {
        coords <- GenomicRanges::GRanges(coords)
    }

    tads_sub <- IRanges::subsetByOverlaps(tads, coords, type = 'within')
    tads_df <- tibble::as_tibble(tads_sub) %>% 
        dplyr::mutate(
            an1 = start, 
            peak_x = start,
            peak_y = end,
            an2 = end, 
            ID = seq(1, n())
        ) %>% 
        dplyr::group_by(ID) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(
            function(df) {
                tibble::tibble(
                    x = c(df$an1, df$peak_x, df$an2),
                    y = c(df$an1, df$peak_y, df$an2), 
                    ID = df$ID
                )
            }
        )
    
    p + ggplot2::geom_path(
        data = tads_df, inherit.aes = FALSE, ggplot2::aes(x, y, group = ID)
    )

}

#' addLoops
#'
#' @param p p
#' @param loops loops
#' @param coords coords
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRanges
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @export

addLoops <- function(p, loops, coords) {

    if (is.character(coords)) {
        coords <- GenomicRanges::GRanges(coords)
    }

    loops_sub <- IRanges::subsetByOverlaps(loops, coords, type = 'within')
    loops_df <- tibble::as_tibble(loops_sub) %>% 
        dplyr::mutate(
            x = first.end - (first.end - first.start)/2,
            y = second.end - (second.end - second.start)/2
        )
    
    p + ggplot2::geom_point(
        data = loops_df, inherit.aes = FALSE, ggplot2::aes(x, y),
        pch = 21, fill = NA, size = 4
    )

}
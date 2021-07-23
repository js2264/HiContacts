#' getPs
#'
#' @param pairs pairs
#'
#' @export

getPs <- function(pairs) {
    ps <- pairs %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
            distance = start2 - start1,
            binned_distance = cut(distance, breaks$break_pos, include.lowest = TRUE),
            binned_distance = breaks$break_pos[-nrow(breaks)][as.numeric(binned_distance)]
        ) %>%
        dplyr::group_by(binned_distance) %>%
        dplyr::mutate(cnt = 1) %>%
        dplyr::summarize(ninter = sum(cnt)) %>%
        dplyr::mutate(p = ninter / sum(ninter)) %>%
        dplyr::left_join(breaks, by = c("binned_distance" = "break_pos")) %>%
        dplyr::mutate(p = p / binwidth, distance = binned_distance) %>%
        dplyr::select(distance, p)
    return(ps)
}

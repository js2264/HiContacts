lsCool <- function(file) {
    `%>%` <- magrittr::`%>%`
    x <- rhdf5::h5ls(file)
    message("\nPossible paths:\n", paste0('\t', x$group, '/', x$name, '\n') %>% unique() %>% stringr::str_replace('//', ''))
    x
}

peekCool <- function(file, path, res = NULL) {
    path <- ifelse(is.null(res), glue::glue('/{path}'), glue::glue('/resolutions/{res}/{path}'))
    res <- as.vector(rhdf5::h5read(file, name = path))
    if (is.list(res)) {
        lapply(
            res, 
            head
        )
    } 
    else {
        head(res)
    }
}

fetchCool <- function(file, path, res = NULL, idx = NULL, ...) {
    path <- ifelse(is.null(res), glue::glue('/{path}'), glue::glue('/resolutions/{res}/{path}'))
    res <- as.vector(rhdf5::h5read(file, name = path, index = list(idx), ...))
    res
}

splitCoords <- function(coords) {
    chr <- stringr::str_replace(coords, ':.*', '')
    chr <- ifelse(length(chr) == 0, NA, chr)
    start <- suppressWarnings(as.numeric(stringr::str_replace_all(coords, '.*:|-.*', '')))
    start <- ifelse(length(start) == 0, NA, start)
    end <- suppressWarnings(as.numeric(stringr::str_replace(coords, '.*-', '')))
    end <- ifelse(length(end) == 0, NA, end)
    list(
        'chr' = chr, 
        'start' = start,
        'end' = end
    )
}

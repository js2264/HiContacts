#' @import magrittr
#' @import rhdf5
#' @import stringr
#' @export

lsCoolFiles <- function(file, full.list = FALSE) {
    `%>%` <- magrittr::`%>%`
    x <- rhdf5::h5ls(file)
    message("\nPossible paths:\n", paste0("\t", x$group, "/", x$name, "\n") %>% unique() %>% stringr::str_replace("//", ""))
    if (full.list) x
}

#' @import tools
#' @import magrittr
#' @import rhdf5
#' @export

lsCoolResolutions <- function(file, full.list = FALSE) {
    if (tools::file_ext(file) != "mcool") stop("Provided file is not .mcool multi-resolution map. Aborting now.")
    `%>%` <- magrittr::`%>%`
    x <- rhdf5::h5ls(file)
    rez <- unique(grep(gsub("/resolutions/", "", x$group), pattern = "/", invert = TRUE, value = TRUE))
    message("\nAvailable resolutions:\n", paste0(rez, collapse = ", "))
    if (full.list) x
}

#' @import glue
#' @import rhdf5
#' @export

peekCool <- function(file, path, res = NULL) {
    path <- ifelse(is.null(res), glue::glue("/{path}"), glue::glue("/resolutions/{res}/{path}"))
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

#' @import glue
#' @import rhdf5
#' @export

fetchCool <- function(file, path, res = NULL, idx = NULL, ...) {
    path <- ifelse(is.null(res), glue::glue("/{path}"), glue::glue("/resolutions/{res}/{path}"))
    as.vector(rhdf5::h5read(file, name = path, index = list(idx), ...))
}

#' @import GenomicRanges
#' @import stringr
#' @export

splitCoords <- function(coords) {
    if (class(coords) == "GRanges") {
        chr <- as.vector(GenomicRanges::seqnames(coords))
        start <- GenomicRanges::start(coords)
        end <- GenomicRanges::end(coords)
        list(
            "chr" = chr,
            "start" = start,
            "end" = end
        )
    }
    else {
        chr <- stringr::str_replace(coords, ":.*", "")
        chr <- ifelse(length(chr) == 0, NA, chr)
        start <- suppressWarnings(as.numeric(stringr::str_replace_all(coords, ".*:|-.*", "")))
        start <- ifelse(length(start) == 0, NA, start)
        end <- suppressWarnings(as.numeric(stringr::str_replace(coords, ".*-", "")))
        end <- ifelse(length(end) == 0, NA, end)
        list(
            "chr" = chr,
            "start" = start,
            "end" = end
        )
    }
}

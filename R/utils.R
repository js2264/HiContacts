#' lsCoolFiles
#'
#' @param file file
#' @param full.list full.list
#'
#' @import rhdf5
#' @import stringr
#' @import tidyr
#' @export

lsCoolFiles <- function(file, full.list = FALSE) {
    `%>%` <- tidyr::`%>%`
    x <- rhdf5::h5ls(file)
    message("\nPossible paths:\n", paste0("\t", x$group, "/", x$name, "\n") %>% unique() %>% stringr::str_replace("//", ""))
    if (full.list) x
}

#' lsCoolResolutions
#'
#' @param file file
#' @param full.list full.list
#'
#' @import tools
#' @import rhdf5
#' @import tidyr
#' @export

lsCoolResolutions <- function(file, full.list = FALSE) {
    if (tools::file_ext(file) != "mcool") stop("Provided file is not .mcool multi-resolution map. Aborting now.")
    `%>%` <- tidyr::`%>%`
    x <- rhdf5::h5ls(file)
    rez <- gsub("/resolutions/", "", x$group) %>%
        grep(., , pattern = "/", invert = TRUE, value = TRUE) %>%
        unique() %>%
        as.numeric() %>%
        sort() %>%
        as.character()
    message("\nAvailable resolutions:\n", paste0(rez, collapse = ", "))
    if (full.list) x
}

#' peekCool
#'
#' @param file file
#' @param path path
#' @param res res
#'
#' @importFrom glue glue
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

#' fetchCool
#'
#' @param file file
#' @param path path
#' @param res res
#' @param idx idx
#' @param ... ...
#'
#' @importFrom glue glue
#' @import rhdf5
#' @export

fetchCool <- function(file, path, res = NULL, idx = NULL, ...) {
    path <- ifelse(is.null(res), glue::glue("/{path}"), glue::glue("/resolutions/{res}/{path}"))
    as.vector(rhdf5::h5read(file, name = path, index = list(idx), ...))
}

#' splitCoords
#'
#' @param coords coords
#'
#' @import stringr
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
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

#' getHicStats
#'
#' @param hicstuff_log hicstuff_log
#'
#' @import stringr
#' @export

getHicStats <- function(hicstuff_log) {
    stats <- list()

    ## -- Parse log file
    nFiltered <- readLines(hicstuff_log) %>%
        grep("INFO :: .* pairs kept", ., value = TRUE) %>%
        stringr::str_replace_all("INFO ::| pairs.*", "") %>%
        as.numeric()
    nDangling <- readLines(hicstuff_log) %>%
        grep("INFO :: .* pairs discarded", ., value = TRUE) %>%
        stringr::str_replace_all(".* Uncuts: |, Weirds.*", "") %>%
        as.numeric() # "Uncut"
    nSelf <- readLines(hicstuff_log) %>%
        grep("INFO :: .* pairs discarded", ., value = TRUE) %>%
        stringr::str_replace_all(".* Loops: |, Uncuts.*", "") %>%
        as.numeric() # "loop"
    nDumped <- readLines(hicstuff_log) %>%
        grep("INFO :: .* pairs discarded", ., value = TRUE) %>%
        stringr::str_replace_all(".* Weirds:", "") %>%
        as.numeric() # "Weird"
    nDups <- readLines(hicstuff_log) %>%
        grep("INFO :: .* PCR duplicates have", ., value = TRUE) %>%
        stringr::str_replace_all(".*out \\(| \\/ .*", "") %>%
        as.numeric()
    nCis <- readLines(hicstuff_log) %>%
        grep("INFO :: Proportion of inter contacts", ., value = TRUE) %>%
        stringr::str_replace_all(".*intra:|, inter:.*", "") %>%
        as.numeric()
    nTrans <- readLines(hicstuff_log) %>%
        grep("INFO :: Proportion of inter contacts", ., value = TRUE) %>%
        stringr::str_replace_all(".* inter: |\\)", "") %>%
        as.numeric()

    stats[["Threshold for dangling pairs (uncut)"]] <- readLines(hicstuff_log) %>%
        grep("INFO :: Filtering with thresholds", ., value = TRUE) %>%
        stringr::str_replace_all(".*thresholds: | loops=.*", "") %>%
        stringr::str_replace(".*=", "") %>%
        as.numeric()
    stats[["Threshold for self pairs (loop)"]] <- readLines(hicstuff_log) %>%
        grep("INFO :: Filtering with thresholds", ., value = TRUE) %>%
        stringr::str_replace_all(".*=", "") %>%
        as.numeric()
    stats[["nFragments"]] <- readLines(hicstuff_log) %>%
        grep("INFO :: .* mapped with Q", ., value = TRUE) %>%
        stringr::str_replace_all(".*/|\\)", "") %>%
        as.numeric() %>%
        `/`(2)
    stats[["nFragments"]] <- readLines(hicstuff_log) %>%
        grep("INFO :: .* mapped with Q", ., value = TRUE) %>%
        stringr::str_replace_all(".*/|\\)", "") %>%
        as.numeric() %>%
        `/`(2)
    stats[["nPairs"]] <- nFiltered + nDangling + nSelf + nDumped
    stats[["pctPairs"]] <- round(stats[["nPairs"]] / stats[["nFragments"]], 4) * 100
    stats[["nFiltered"]] <- nFiltered
    stats[["pctFiltered"]] <- round(stats[["nFiltered"]] / stats[["nPairs"]], 4) * 100
    stats[["nDangling"]] <- nDangling # "Uncut"
    stats[["pctDangling"]] <- round(stats[["nDangling"]] / stats[["nPairs"]], 4) * 100
    stats[["nSelf"]] <- nSelf # "loop"
    stats[["pctSelf"]] <- round(stats[["nSelf"]] / stats[["nPairs"]], 4) * 100
    stats[["nDumped"]] <- nDumped # "Weird"
    stats[["pctDumped"]] <- round(stats[["nDumped"]] / stats[["nPairs"]], 4) * 100
    stats[["nUnique"]] <- nFiltered - nDups
    stats[["pctUnique"]] <- round(stats[["nUnique"]] / stats[["nPairs"]], 4) * 100
    stats[["CisTransRatio"]] <- round(nCis / nTrans, 4)
    stats[["CisTransPct"]] <- round(nCis / (nCis + nTrans), 4) * 100

    return(stats)
}

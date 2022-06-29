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
    if (class(coords)[1] == "GRanges") {
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

formatCoords <- function(coords) {
    if (class(coords)[1] == "GRanges") {
        chr <- as.vector(GenomicRanges::seqnames(coords))
        start <- GenomicRanges::start(coords)
        end <- GenomicRanges::end(coords)
        paste0(chr, ':', format(start, big.mark=","), '-', format(end, big.mark=","))
    }
    else {
        chr <- stringr::str_replace(coords, ":.*", "")
        start <- suppressWarnings(as.numeric(stringr::str_replace_all(coords, ".*:|-.*", "")))
        end <- suppressWarnings(as.numeric(stringr::str_replace(coords, ".*-", "")))
        if (is.na(start)) {
            return(chr)
        }
        paste0(chr, ':', format(start, big.mark=","), '-', format(end, big.mark=","))
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

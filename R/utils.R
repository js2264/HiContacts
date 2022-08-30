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

#' formatCoords
#'
#' @param coords coords
#'
#' @import stringr
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @export

formatCoords <- function(coords) {
    if (class(coords)[1] == "GRanges") {
        chr <- as.vector(GenomicRanges::seqnames(coords))
        start <- GenomicRanges::start(coords)
        end <- GenomicRanges::end(coords)
        paste0(chr, ':', format(start, big.mark=","), '-', format(end, big.mark=","))
    }
    else {
        if (grepl('x', coords)) {
            coords
        }
        else {
            chr <- stringr::str_replace(coords, ":.*", "")
            start <- suppressWarnings(as.numeric(stringr::str_replace_all(coords, ".*:|-.*", "")))
            end <- suppressWarnings(as.numeric(stringr::str_replace(coords, ".*-", "")))
            if (is.na(start)) {
                return(chr)
            }
            paste0(chr, ':', format(start, big.mark=",", scientific = FALSE), '-', format(end, big.mark=",", scientific = FALSE))
        }
    }
}

char2pair <- function(char) {
    if (methods::is(char, 'Pairs')) {
        return(char)
    }
    splitst <- stringr::str_split(char, ' x ')[[1]]
    S4Vectors::Pairs(
        GenomicRanges::GRanges(splitst[[1]]), 
        GenomicRanges::GRanges(splitst[[2]])
    )
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

fullContactInteractions <- function(chr, start, end, binning) {
    full_anchors <- GRanges(
        seqnames = chr, 
        IRanges::IRanges(
            start = seq(start, end-1, by = binning),
            width = binning
        )
    )
    GInteractions(
        full_anchors[rep(seq_along(full_anchors), length(full_anchors))], 
        full_anchors[rep(seq_along(full_anchors), each = length(full_anchors))], 
        full_anchors
    )
}

`sdiag<-` <- function(A, k = 0, value) {
    p <- ncol(A)
    n <- nrow(A)
    if (k>p-1||-k > n-1) return()
    if (k >= 0) {
        i <- 1:n
        j <- (k+1):p
    } 
    else {
        i <- (-k+1):n
        j <- 1:p
    }
    if (length(i)>length(j)) i <- i[1:length(j)] else j <- j[1:length(i)]
    ii <- i + (j-1) * n 
    A[ii] <- value
    A
} 

sort_pairs <- function(pairs) {
    p <- S4Vectors::zipup(pairs)
    p_sorted <- GRangesList(lapply(p, function(gr) {
        sort(gr)
    }))
    S4Vectors::zipdown(p_sorted)
}

as_GInteractions <- function(df) {
    gi <- GInteractions(
        anchor1 = GRanges(
            df$seqnames1, IRanges::IRanges(df$start1, df$end1)
        ),
        anchor2 = GRanges(
            df$seqnames2, IRanges::IRanges(df$start2, df$end2)
        )
    )
    if ('score' %in% colnames(df)) gi$score <- df$score
    gi
}
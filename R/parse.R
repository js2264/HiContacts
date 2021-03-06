#' getAnchors
#'
#' @param file file
#' @param resolution resolution
#' @param balanced balanced
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom IRanges IRanges
#' @export

getAnchors <- function(file, resolution = NULL, balanced = "cooler") {
    bins <- fetchCool(file, "bins", resolution)
    anchors <- GenomicRanges::GRanges(
        bins$chr,
        IRanges::IRanges(bins$start + 1, bins$end),
        seqinfo = cool2seqinfo(file, resolution)
    )
    names(anchors) <- paste(GenomicRanges::seqnames(anchors), GenomicRanges::start(anchors), GenomicRanges::end(anchors), sep = "_")
    if ("weight" %in% names(bins) & {
        balanced == "cooler" | balanced == TRUE
    }) {
        anchors$weight <- bins$weight
    } else {
        weight <- 1
    }
    return(anchors)
}

#' getCounts
#'
#' This was adapted from `dovetail-genomics/coolR`
#'
#' @param file file
#' @param coords coords
#' @param anchors anchors
#' @param coords2 coords2
#' @param resolution resolution
#'
#' @import zeallot
#' @import tidyr
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom glue glue
#' @importFrom S4Vectors subjectHits
#' @export

getCounts <- function(file,
                      coords,
                      anchors,
                      coords2 = NULL,
                      resolution = NULL) {
    
    `%<-%` <- zeallot::`%<-%`

    ## Get chunks to parse
    if (length(coords) <= 1) { ## coords = NULL or 1 GRanges

        ## Process coordinates
        c(coords_chr, coords_start, coords_end) %<-% splitCoords(coords)
        if (any(is.na(coords_start) & !is.na(coords_chr))) {
            coords_start <- rep(1, length(coords_start))
            coords_end <- GenomeInfoDb::seqlengths(anchors)[coords_chr]
        }

        ## Check that queried chr. exists
        if (any(!coords_chr %in% as.vector(GenomicRanges::seqnames(anchors)) & !is.na(coords_chr))) {
            sn <- paste0(unique(as.vector(GenomicRanges::seqnames(anchors))), collapse = ", ")
            stop(glue::glue("{coords_chr} not in file. Available seqnames: {sn}"))
        }

        if (is.na(coords_chr)) {
            chunks <- NULL
        } 
        else {
            start <- GenomicRanges::GRanges(coords_chr, IRanges::IRanges(coords_start + 1, width = 1))
            end <- GenomicRanges::GRanges(coords_chr, IRanges::IRanges(coords_end, width = 1))
            start_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(start, anchors))
            end_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(end, anchors))
            bin1_idx <- fetchCool(file, path = "indexes/bin1_offset", resolution, idx = seq(start_idx, end_idx))
            slice <- sum(bin1_idx[-1] - bin1_idx[-length(bin1_idx)])
            chunks <- seq(
                bin1_idx[1] + 1,
                bin1_idx[1] + 1 + slice
            )
        }
    }

    else { ## If there are several GRanges provided (for APA)
        
        coords <- subsetByOverlaps(coords, as(cool2seqinfo(file, resolution), 'GRanges'), type = 'within') ## Filter coords to exclude those out of mcool 
        coords <- GenomicRanges::reduce(coords)
        c(coords_chr, coords_start, coords_end) %<-% splitCoords(coords)

        ## Check that queried chr. exists
        if (any(!coords_chr %in% as.vector(GenomicRanges::seqnames(anchors)) & !is.na(coords_chr))) {
            sn <- paste0(unique(as.vector(GenomicRanges::seqnames(anchors))), collapse = ", ")
            stop(glue::glue("{coords_chr} not in file. Available seqnames: {sn}"))
        }

        start <- GenomicRanges::GRanges(coords_chr, IRanges::IRanges(coords_start + 1, width = 1))
        end <- GenomicRanges::GRanges(coords_chr, IRanges::IRanges(coords_end, width = 1))
        start_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(start, anchors))
        end_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(end, anchors))
        chunks <- lapply(seq_along(start_idx), function(K) {
            bin1_idx <- fetchCool(file, path = "indexes/bin1_offset", resolution, idx = seq(start_idx[K], end_idx[K]))
            slice <- sum(bin1_idx[-1] - bin1_idx[-length(bin1_idx)])
            seq(
                bin1_idx[1] + 1,
                bin1_idx[1] + 1 + slice
            )
        }) %>% unlist()

    }

    ## Reading the chunks from the cool file
    df <- tidyr::tibble(
        bin1_id = fetchCool(file, path = "pixels/bin1_id", resolution, idx = chunks) + 1,
        bin2_id = fetchCool(file, path = "pixels/bin2_id", resolution, idx = chunks) + 1,
        count = fetchCool(file, path = "pixels/count", resolution, idx = chunks)
    )

    ## ------- Case when coords2 are provided...
    if (is.null(coords2)) {
        df <- df[df$bin2_id %in% df$bin1_id, ]
        return(df)
    }
    else {
        c(coords_chr2, coords_start2, coords_end2) %<-% splitCoords(coords2)
        if (is.na(coords_start2) & !is.na(coords_chr2)) {
            coords_start2 <- 1
            coords_end2 <- GenomeInfoDb::seqlengths(anchors)[coords_chr2]
        }
        else if (is.na(coords_chr2)) {
            chunks <- NULL
        } 
        else {
            start <- GenomicRanges::GRanges(coords_chr2, IRanges::IRanges(coords_start2 + 1, width = 1))
            end <- GenomicRanges::GRanges(coords_chr2, IRanges::IRanges(coords_end2, width = 1))
            start_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(start, anchors))
            end_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(end, anchors))
            bin2_idx <- fetchCool(file, path = "indexes/bin1_offset", res, idx = seq(start_idx, end_idx))
            slice <- sum(bin2_idx[-1] - bin2_idx[-length(bin2_idx)])
            chunks <- seq(
                bin2_idx[1] + 1,
                bin2_idx[1] + 1 + slice
            )
            df2 <- tidyr::tibble(
                bin1_id = fetchCool(file, path = "pixels/bin1_id", res, idx = chunks),
                bin2_id = fetchCool(file, path = "pixels/bin2_id", res, idx = chunks),
                count = fetchCool(file, path = "pixels/count", res, idx = chunks)
            )
            df <- df[df$bin2_id %in% df2$bin1_id, ]
            return(df)
        }
    }

}

#' fetchCool
#'
#' @param file file
#' @param path path
#' @param resolution resolution
#' @param idx idx
#' @param ... ...
#'
#' @importFrom glue glue
#' @import rhdf5
#' @export

fetchCool <- function(file, path, resolution = NULL, idx = NULL, ...) {
    path <- ifelse(is.null(resolution), glue::glue("/{path}"), glue::glue("/resolutions/{resolution}/{path}"))
    as.vector(rhdf5::h5read(file, name = path, index = list(idx), ...))
}

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

lsCoolResolutions <- function(file, full.list = FALSE, silent = FALSE) {
    if (tools::file_ext(file) != "mcool") stop("Provided file is not .mcool multi-resolution map. Aborting now.")
    `%>%` <- tidyr::`%>%`
    x <- rhdf5::h5ls(file)
    rez <- gsub("/resolutions/", "", x$group) %>%
        grep(., , pattern = "/", invert = TRUE, value = TRUE) %>%
        unique() %>%
        as.numeric() %>%
        sort() %>%
        as.character()
    if (!silent) message("\nAvailable resolutions:\n", paste0(rez, collapse = ", "))
    if (full.list) {
        x
    } else {
        invisible(as.integer(rez))
    }
}

#' peekCool
#'
#' @param file file
#' @param path path
#' @param resolution resolution
#'
#' @importFrom glue glue
#' @import rhdf5
#' @export

peekCool <- function(file, path, resolution = NULL) {
    path <- ifelse(is.null(resolution), glue::glue("/{path}"), glue::glue("/resolutions/{resolution}/{path}"))
    resolution <- as.vector(rhdf5::h5read(file, name = path))
    if (is.list(resolution)) {
        lapply(
            resolution,
            head
        )
    }
    else {
        head(resolution)
    }
}

#' cool2seqinfo
#'
#' @param file file
#' @param resolution resolution
#'
#' @importFrom GenomeInfoDb Seqinfo
#' @export

cool2seqinfo <- function(file, resolution = NULL) {
    chroms <- fetchCool(file, "chroms", resolution)
    seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames = as.vector(chroms$name),
        seqlengths = as.vector(chroms$length)
    )
    return(seqinfo)
}

#' cool2gi
#'
#' @param file file
#' @param coords coords
#' @param coords2 coords2
#' @param resolution resolution
#'
#' @import InteractionSet
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges resize
#' @export

cool2gi <- function(file, coords = NULL, coords2 = NULL, resolution = NULL) {
    anchors <- getAnchors(file, resolution)
    cnts <- getCounts(file, coords = coords, anchors = anchors, coords2 = coords2, resolution = resolution)
    gi <- InteractionSet::GInteractions(
        anchors[cnts$bin1_id],
        anchors[cnts$bin2_id],
        count = cnts$count
    )
    gi$bin1 <- cnts$bin1_id
    gi$bin2 <- cnts$bin2_id
    if (!is.null(gi$anchor1.weight) & !is.null(gi$anchor2.weight)) {
        gi$score <- gi$count * gi$anchor1.weight * gi$anchor2.weight
    } else {
        gi$score <- gi$count
    }
    InteractionSet::regions(gi)$chr <- GenomicRanges::seqnames(InteractionSet::regions(gi))
    InteractionSet::regions(gi)$center <- GenomicRanges::start(GenomicRanges::resize(InteractionSet::regions(gi), fix = "center", width = 1))
    return(gi)
}

#' gi2cm
#'
#' @param gi gi
#'
#' @import InteractionSet
#' @importFrom GenomicRanges mcols
#' @export

gi2cm <- function(gi) {
    InteractionSet::inflate(
        gi,
        rows = 1:length(InteractionSet::regions(gi)),
        columns = 1:length(InteractionSet::regions(gi)),
        fill = GenomicRanges::mcols(gi)[['score']]
    )
}

cm2matrix <- function(cm, replace_NA = NA) {
    m <- Matrix::as.matrix(cm)
    m[is.na(m)] <- replace_NA
    m
}

#' readPairs
#'
#' @param file pairs file: `<readname>\t<chr1>\t<start1>\t<chr2>\t<start2>`
#'
#' @import tibble
#' @import dplyr
#' @export

pairs2gi <- function(
    file,
    chr1.field = 2, 
    start1.field = 3, 
    chr2.field = 4, 
    start2.field = 5, 
    nThread = 16, 
    nrows = Inf
) {

    ## Use zgrep if pairs file is zipped (.gz)
    if (grepl('.gz$', file)) {
        grep_cmd <- "zgrep"
    }
    else {
        grep_cmd <- "grep"
    }

    headers <- data.table::fread(glue::glue("{grep_cmd} -v '^#' {file}"), nrows = 1)
    
    anchors1 <- data.table::fread(
        glue::glue("{grep_cmd} -v '^#' {file}"), 
        header = FALSE, 
        drop = {1:ncol(headers)}[!{1:ncol(headers) %in% c(chr1.field, start1.field)}], 
        nThread = nThread, 
        col.names = c('chr', 'start'), 
        nrows = nrows
    ) 
    anchors2 <- data.table::fread(
        glue::glue("{grep_cmd} -v '^#' {file}"), 
        header = FALSE, 
        drop = {1:ncol(headers)}[!{1:ncol(headers) %in% c(chr2.field, start2.field)}], 
        nThread = nThread, 
        col.names = c('chr', 'start'), 
        nrows = nrows
    ) 
    anchor_one <- GRanges(
        anchors1[['chr']],
        IRanges(anchors1[['start']], width = 1)
    )
    anchor_two <- GRanges(
        anchors2[['chr']],
        IRanges(anchors2[['start']], width = 1)
    )

    gi <- GenomicInteractions::GenomicInteractions(anchor_one, anchor_two)
    gi$distance <- GenomicInteractions::calculateDistances(gi) 

    return(gi)

}

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
        bin_id = seq_along(bins$chr), 
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

#' getCounts2
#'
#' Function to extract counts for a non-square matrix (@ a pair of coordinates).
#' 
#' This was adapted from `dovetail-genomics/coolR`
#'
#' @param file file
#' @param pair pair (e.g. S4Vectors::Pairs(GRanges("II:200000-300000"), GRanges("II:70000-100000")))
#' @param anchors anchors
#' @param resolution resolution
#' @param BPPARAM BPPARAM
#'
#' @import methods
#' @import zeallot
#' @import tidyr
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom glue glue
#' @importFrom S4Vectors subjectHits
#' @export

getCounts2 <- function(file,
                      pair,
                      anchors,
                      resolution = NULL, 
                      BPPARAM = BiocParallel::bpparam()) {
    
    `%<-%` <- zeallot::`%<-%`

    # Make sure pair is sorted (first is first, second is after)
    pair <- sort_pairs(pair)

    # For each pair, get the coords for first and second.
    if (length(pair) > 1) {
        coords_list <- lapply(pair, function(x) {
            coords <- unlist(S4Vectors::zipup(x))
            splitCoords(coords)
        })
    }
    else {
        coords <- unlist(S4Vectors::zipup(pair))
        coords_list <- list(splitCoords(coords))
    }

    ## Check that queried chr. exists
    chrs <- unlist(lapply(coords_list, '[[', 'chr'))
    if (any(!chrs %in% as.vector(GenomicRanges::seqnames(anchors)) & !is.na(chrs))) {
        sn <- paste0(unique(as.vector(GenomicRanges::seqnames(anchors))), collapse = ", ")
        stop(glue::glue("Some chr. are not available. Available seqnames: {sn}"))
    }

    df <- BiocParallel::bplapply(seq_along(coords_list), function(i) {
        # For a single pair of GRanges: Create start and end GRanges...
        .x <- coords_list[[i]]
        start <- GenomicRanges::GRanges(.x$chr, IRanges::IRanges(.x$start + 1, width = 1))
        end <- GenomicRanges::GRanges(.x$chr, IRanges::IRanges(.x$end, width = 1))
        start_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(start, anchors))
        end_idx <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(end, anchors))
        bin1_idx_1 <- fetchCool(file, path = "indexes/bin1_offset", resolution, idx = seq(start_idx[1], end_idx[1]))
        bin1_idx_2 <- fetchCool(file, path = "indexes/bin1_offset", resolution, idx = seq(start_idx[2], end_idx[2]))
        chunks_1 <- seq(
            bin1_idx_1[1] + 1,
            bin1_idx_1[1] + 1 + {sum(bin1_idx_1[-1] - bin1_idx_1[-length(bin1_idx_1)])}
        )
        chunks_2 <- seq(
            bin1_idx_2[1] + 1,
            bin1_idx_2[1] + 1 + {sum(bin1_idx_2[-1] - bin1_idx_2[-length(bin1_idx_2)])}
        )
        valid_bin2 <- unique(fetchCool(file, path = "pixels/bin1_id", resolution, idx = chunks_2))
        tidyr::tibble(
            bin1_id = fetchCool(file, path = "pixels/bin1_id", resolution, idx = chunks_1),
            bin2_id = fetchCool(file, path = "pixels/bin2_id", resolution, idx = chunks_1),
            count = fetchCool(file, path = "pixels/count", resolution, idx = chunks_1)
        ) %>% 
            dplyr::filter(bin2_id %in% valid_bin2) %>% 
            dplyr::mutate(pair = paste0(
                as.character(S4Vectors::first(pair[i])), 
                '_', 
                as.character(S4Vectors::second(pair[i]))
            ))
    }, BPPARAM = BPPARAM) %>% dplyr::bind_rows()

    return(df)
}

#' Function to extract counts for a square matrix from an mcool file

getCounts <- function(file,
                      coords,
                      anchors,
                      resolution = NULL) {
    
    `%<-%` <- zeallot::`%<-%`
    `%within%` <- IRanges::`%within%`
    check_cool_format(file, resolution)

    ## Process coordinates
    c(coords_chr, coords_start, coords_end) %<-% splitCoords(coords)
    
    # If only chr. names are provided, find their start and stop
    if (any(is.na(coords_start) & !is.na(coords_chr))) {
        coords_start <- rep(1, length(coords_start))
        coords_end <- GenomeInfoDb::seqlengths(anchors)[coords_chr]
    }

    ## Check that queried chr. exist
    if (any(!coords_chr %in% as.vector(GenomicRanges::seqnames(anchors)) & !is.na(coords_chr))) {
        sn <- paste0(unique(as.vector(GenomicRanges::seqnames(anchors))), collapse = ", ")
        stop(glue::glue("{coords_chr} not in file. Available seqnames: {sn}"))
    }

    ## Find out which chunks of the mcool to recover
    if (is.na(coords_chr)) {
        chunks <- NULL
    } 
    else {
        gr <- GRanges(coords_chr, IRanges::IRanges(coords_start, coords_end))
        sub <- anchors %within% gr
        s <- anchors$bin_id[sub]
        bin_idx <- fetchCool(file, path = "indexes/bin1_offset", resolution, idx = s)
        chunks <- seq(min(bin_idx)+1, max(bin_idx)+1, by = 1)
    }

    ## Reading the chunks from the cool file
    df <- tidyr::tibble(
        bin1_id = fetchCool(file, path = "pixels/bin1_id", resolution, idx = chunks) + 1,
        bin2_id = fetchCool(file, path = "pixels/bin2_id", resolution, idx = chunks) + 1,
        count = fetchCool(file, path = "pixels/count", resolution, idx = chunks)
    )

    ## Filter to only get interesting bins
    df <- df[df$bin2_id %in% df$bin1_id, ]

    return(df)

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
    check_cool_format(file, resolution)
    path <- ifelse(is.null(resolution), glue::glue("/{path}"), glue::glue("/resolutions/{resolution}/{path}"))
    as.vector(rhdf5::h5read(file, name = path, index = list(idx), ...))
}

#' lsCoolFiles
#'
#' @param file file
#'
#' @import rhdf5
#' @import stringr
#' @import tidyr
#' @import GenomicInteractions
#' @export

lsCoolFiles <- function(file) {
    `%>%` <- tidyr::`%>%`
    x <- rhdf5::h5ls(file) %>% 
        mutate(path = paste0(group, "/", name)) %>% 
        pull(path) %>% 
        unique() %>% 
        stringr::str_replace("//", "/")
    len <- length(x)
    if (len > 10) {
        mess <- c(
            paste0(x[seq_len(5)], "\n"), 
            paste0("... (", len-10, " more paths)\n"), 
            paste0(x[(len-5+1):len], "\n")
        )
        message(mess)
    }
    else {
        message(x)
    }
    invisible(x)
}

#' lsCoolResolutions
#'
#' @param file file
#' @param silent silent
#'
#' @import tools
#' @import rhdf5
#' @import tidyr
#' @export

lsCoolResolutions <- function(file, silent = FALSE) {
    `%>%` <- tidyr::`%>%`
    if (is_cool(file)) {
        x <- rhdf5::h5ls(file)
        bin_ends <- peekCool(file, '/bins/end', n = 2)
        rez <- bin_ends[2] - bin_ends[1]
    }
    if (is_mcool(file)) {
        x <- rhdf5::h5ls(file)
        rez <- gsub("/resolutions/", "", x$group) %>%
            grep(., , pattern = "/", invert = TRUE, value = TRUE) %>%
            unique() %>%
            as.numeric() %>%
            sort() %>%
            as.character()
    }
    if (!silent) message(S4Vectors::coolcat("resolutions(%d): %s", rez))
    invisible(as.integer(rez))
}

#' peekCool
#'
#' @param file file
#' @param path path
#' @param resolution resolution
#' @param n n
#'
#' @importFrom Matrix head
#' @importFrom glue glue
#' @import rhdf5
#' @export

peekCool <- function(file, path, resolution = NULL, n = 10) {
    check_cool_format(file, resolution)
    path <- ifelse(is.null(resolution), glue::glue("/{path}"), glue::glue("/resolutions/{resolution}/{path}"))
    resolution <- as.vector(rhdf5::h5read(file, name = path))
    if (is.list(resolution)) {
        lapply(resolution, Matrix::head, n = n)
    }
    else {
        Matrix::head(resolution, n = n)
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
    check_cool_format(file, resolution)
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
#' @param coords NULL, character, or GRanges. Can also be a Pairs object of paired GRanges.
#' @param resolution resolution
#'
#' @import InteractionSet
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges resize
#' @export

cool2gi <- function(file, coords = NULL, resolution = NULL) {
    `%<-%` <- zeallot::`%<-%`
    check_cool_format(file, resolution)
    
    # Check if the provided coords are GRanges or Pairs
    is_pair <- is(coords, 'Pairs')

    # Get anchors from mcool
    anchors <- getAnchors(file, resolution)

    # Get raw counts for bins from mcool
    if (!is_pair) {
        c(coords_chr, coords_start, coords_end) %<-% splitCoords(coords)
        cnts <- getCounts(file, coords = coords, anchors = anchors, resolution = resolution)
    }
    else if (is_pair) {
        coords_list <- lapply(coords, function(x) {
            coords <- unlist(S4Vectors::zipup(x))
            splitCoords(coords)
        })
        coords_chr <- unlist(lapply(coords_list, '[[', 'chr'))
        cnts <- getCounts2(file, pair = coords, anchors = anchors, resolution = resolution)
    }

    # Associate raw counts for bins to corresponding anchors
    gi <- InteractionSet::GInteractions(
        anchors[cnts$bin1_id],
        anchors[cnts$bin2_id],
        count = cnts$count
    )
    gi$bin1 <- cnts$bin1_id
    gi$bin2 <- cnts$bin2_id
    
    if (!is.null(coords_chr) & all(!is.na(coords_chr)) & !is_pair) {
        # Make sure no extra GInteractions is pulled from cool (happends e.g. when fetching whole chrs.)
        # DEFINITELY HACKY HERE
        sub <- seqnames(InteractionSet::anchors(gi)[[1]]) == coords_chr & 
            seqnames(InteractionSet::anchors(gi)[[2]]) == coords_chr
        gi <- gi[sub]
        regs <- unique(c(
            InteractionSet::anchors(gi)[[1]], 
            InteractionSet::anchors(gi)[[2]]
        ))
        names(regs) <- paste(
            GenomicRanges::seqnames(regs), 
            GenomicRanges::start(regs), 
            GenomicRanges::end(regs), 
            sep = "_"
        )
        InteractionSet::replaceRegions(gi) <- regs
    }

    # Get balanced counts if they exist
    if (!is.null(gi$anchor1.weight) & !is.null(gi$anchor2.weight)) {
        gi$score <- gi$count * gi$anchor1.weight * gi$anchor2.weight
    } 
    else {
        gi$score <- gi$count
    }

    # Add extra info
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
        rows = seq_along(InteractionSet::regions(gi)),
        columns = seq_along(InteractionSet::regions(gi)),
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
#' @importFrom data.table fread
#' @importFrom glue glue
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicInteractions GenomicInteractions
#' @importFrom GenomicInteractions calculateDistances
#' @importFrom IRanges IRanges
#' @import tibble
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
        drop = {seq_len(ncol(headers))}[!{seq_len(ncol(headers)) %in% c(chr1.field, start1.field)}], 
        nThread = nThread, 
        col.names = c('chr', 'start'), 
        nrows = nrows
    ) 
    anchors2 <- data.table::fread(
        glue::glue("{grep_cmd} -v '^#' {file}"), 
        header = FALSE, 
        drop = {seq_len(ncol(headers))}[!{seq_len(ncol(headers)) %in% c(chr2.field, start2.field)}], 
        nThread = nThread, 
        col.names = c('chr', 'start'), 
        nrows = nrows
    ) 
    anchor_one <- GenomicRanges::GRanges(
        anchors1[['chr']],
        IRanges::IRanges(anchors1[['start']], width = 1)
    )
    anchor_two <- GenomicRanges::GRanges(
        anchors2[['chr']],
        IRanges::IRanges(anchors2[['start']], width = 1)
    )
    gi <- GenomicInteractions::GenomicInteractions(anchor_one, anchor_two)
    gi$distance <- GenomicInteractions::calculateDistances(gi) 
    return(gi)
}

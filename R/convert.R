#' cool2seqinfo
#'
#' @param file file
#' @param res res
#'
#' @importFrom GenomeInfoDb Seqinfo
#' @export

cool2seqinfo <- function(file, res = NULL) {
    chroms <- fetchCool(file, "chroms", res)
    seqinfo <- GenomeInfoDb::Seqinfo(
        seqnames = as.vector(chroms$name),
        seqlengths = as.vector(chroms$length)
    )
    return(seqinfo)
}

#' cool2gi
#'
#' @param file file
#' @param balanced balanced
#' @param coords coords
#' @param coords2 coords2
#' @param res res
#'
#' @import InteractionSet
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges resize
#' @export

cool2gi <- function(file, balanced = "cooler", coords = NULL, coords2 = NULL, res = NULL) {
    anchors <- getAnchors(file, res, balanced = balanced)
    cnts <- getCounts(file, coords = coords, anchors = anchors, coords2 = coords2, res = res)
    gi <- InteractionSet::GInteractions(
        anchors[cnts$bin1_id],
        anchors[cnts$bin2_id],
        count = cnts$count
    )
    gi$bin1 <- cnts$bin1_id
    gi$bin2 <- cnts$bin2_id

    if ({
        !is.null(gi$anchor1.weight) & !is.null(gi$anchor2.weight)
    } & balanced == "cooler") {
        gi$score <- log10(gi$count * gi$anchor1.weight * gi$anchor2.weight)
    }
    else if (balanced == "ICE") {
        gi <- iceGis(gi)
        gi$score <- log10(gi$score + 1)
    }
    else {
        gi$score <- log10(gi$count + 1)
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
#' @export

gi2cm <- function(gi) {
    InteractionSet::inflate(
        gi,
        rows = 1:length(InteractionSet::regions(gi)),
        columns = 1:length(InteractionSet::regions(gi)),
        fill = gi$count
    )
}

#' gi2mat
#'
#' @param gis gis
#' @param limits limits
#' @param truncate_tip truncate_tip
#'
#' @import tibble
#' @import dplyr
#' @import InteractionSet
#' @import tidyr
#' @importFrom GenomicRanges width
#' @export

gi2mat <- function(gis, limits = NULL, truncate_tip = 0.5) {
    `%>%` <- tidyr::`%>%`

    ## -- Convert gis to table and extract x/y
    mat <- gis %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
            x = floor(end1 - (end1 - start1) / 2),
            y = floor(end2 - (end2 - start2) / 2)
        )

    # -- Define plotting approach
    binsize <- GenomicRanges::width(InteractionSet::regions(gis)[1])
    max_bin_dist <- (max(mat$bin2 - mat$bin1) + 1) * truncate_tip

    ## -- Matrix limits
    if (!is.null(limits)) {
        m <- limits[1]
        M <- limits[2]
    } else {
        M <- max(gis$score[abs(gis$bin1 - gis$bin2) >= 0 & abs(gis$bin1 - gis$bin2) <= max_bin_dist], na.rm = TRUE)
        m <- min(gis$score[abs(gis$bin1 - gis$bin2) >= 0 & abs(gis$bin1 - gis$bin2) <= max_bin_dist], na.rm = TRUE)
        limits <- c(m, M)
    }

    ## -- Clamp scores to limits
    mat <- dplyr::mutate(mat, score = ifelse(score > M, M, ifelse(score < m, m, score)))

    ## -- Add lower triangular matrix scores
    mat <- rbind(mat, mat %>% dplyr::mutate(x2 = y, y = x, x = x2) %>% dplyr::select(-x2))

    ## -- Filter matrix to cut out the tip of the triangle
    mat <- mat[abs(mat$bin2 - mat$bin1) <= max_bin_dist, ]

    mat_ <- mat %>%
        ## -- Filter upper diagonal
        dplyr::filter(y >= x) %>%
        ## -- Recompute center of each bin and distance to the diagonal
        dplyr::mutate(
            x_ = x + (y - x) / 2,
            dist = (y - x) / {
                sqrt(2) * sqrt(2)
            },
            ID = 1:dplyr::n()
        ) %>%
        ## -- Group by bins
        dplyr::group_by(ID) %>%
        ## -- Compute each corner of each of the 45deg-tilted bins
        dplyr::mutate(
            A_x = x_ - binsize / 2,
            A_y = dist,
            B_x = x_,
            B_y = dist + binsize / 2,
            C_x = x_ + binsize / 2,
            C_y = dist,
            D_x = x_,
            D_y = dist - binsize / 2
        ) %>%
        tidyr::pivot_longer(cols = dplyr::starts_with(c("A_", "B_", "C_", "D_")), names_to = "pt", values_to = "value_coord") %>%
        tidyr::extract(pt, c("pt", "coord_x"), "(.)_(.)") %>%
        tidyr::spread(coord_x, value_coord) %>%
        dplyr::mutate(y = ifelse(y < 0, 0, y)) %>%
        ## -- Rescale y for proper scaling when plotting
        dplyr::mutate(y = y * {
            truncate_tip / {
                sqrt(2) * sqrt(2)
            }
        })
}

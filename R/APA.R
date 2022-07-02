
#' APA
#'
#' @param file file
#' @param coords coords
#' @param res res
#' @param limits limits
#' @param dpi dpi
#' @param rasterize rasterize
#' @param symmetrical symmetrical
#' @param BPPARAM BPPARAM
#' @param scale scale
#' @param cmap cmap
#'
#' @import InteractionSet
#' @import ggrastr
#' @import ggplot2
#' @import BiocParallel
#' @import zeallot
#' @import tibble
#' @import dplyr
#' @importFrom GenomicRanges width
#' @export

APA <- function(file, coords, res = NULL, limits = NULL, dpi = 500, rasterize = TRUE, symmetrical = TRUE, BPPARAM = BiocParallel::bpparam(), scale = FALSE, cmap = afmhotr_colors) {
    `%>%` <- tidyr::`%>%`

    ## -- Coerce loops into coords and coords2
    if (class(coords) == "GInteractions") {
        isLoops <- TRUE
        loops <- coords
        coords <- InteractionSet::anchors(loops)[[1]]
        coords2 <- InteractionSet::anchors(loops)[[2]]
    }
    else {
        isLoops <- FALSE
        coords2 <- NULL
    }

    ## -- Check that coords are all the same width
    if (length(unique(GenomicRanges::width(coords))) != 1) {
        stop("Please provide GRanges that are all the same width. Aborting now.")
    }
    else {
        wi <- unique(GenomicRanges::width(coords))
        if (is.null(res)) {
            res0 <- unique(GenomicRanges::width(cool2gi(file, res = res, coords = coords[1], coords2 = coords[1]))[[1]])
        }
        else {
            res0 <- res
        }
        breaks <- seq(
            {
                -wi / 2
            } - res0,
            {
                wi / 2
            } + res0,
            length.out = wi / res0 + 3
        )
    }

    ## -- Define plotting approach
    if (rasterize) {
        plotFun <- ggrastr::geom_tile_rast(raster.dpi = dpi)
    }
    else {
        plotFun <- ggplot2::geom_tile()
    }

    if (!isLoops) {
        
        stop("Loop APA plots are not supported yet.")
    
    }
    else {
        ## -- Load entire matrix in memory
        gis <- cool2gi(file, res = res, coords = NULL)

        ## -- Calculate scores over each coords
        nc <- length(coords)
        br <- seq(1, nc+100, by = 100)
        coords_split <- coords %>%
            plyranges::mutate(
                chunk = as.numeric(cut(seq(1, nc), breaks = br, include.lowest = TRUE))
            )
        mats <- BiocParallel::bplapply(BPPARAM = BPPARAM, unique(coords_split$chunk), function(K) {
            coords_sub <- coords_split[coords_split$chunk == K]
            coords_sub %>% 
                plyranges::mutate(ID = seq(1, length(.))) %>% 
                plyranges::select(-score) %>%
                plyranges::join_overlap_left(gis) %>% 
                plyranges::group_by(ID) %>% 
                tibble::as_tibble() %>% 
                dplyr::left_join(
                    tibble::as_tibble(anchors) %>% 
                        dplyr::mutate(bin1 = seq(1, length(anchors)), center1 = end - (end-start)/2) %>% 
                        dplyr::select(bin1, center1),
                    by = 'bin1'
                ) %>% 
                dplyr::left_join(
                    tibble::as_tibble(anchors) %>% 
                        dplyr::mutate(bin2 = seq(1, length(anchors)), center2 = end - (end-start)/2) %>% 
                        dplyr::select(bin2, center2),
                    by = 'bin2'
                ) %>% 
                dplyr::mutate(
                    x = center1,
                    y = center2,
                    x = breaks[as.numeric(cut(x - midpoint, breaks, include.lowest = TRUE)) + 1],
                    y = breaks[as.numeric(cut(y - midpoint2, breaks, include.lowest = TRUE)) + 1]
                ) %>% 
                dplyr::group_by(x, y) %>%
                dplyr::filter(
                    x >= min(breaks[2:(length(breaks) - 1)]), 
                    x <= max(breaks[2:(length(breaks) - 1)]), 
                    y >= min(breaks[2:(length(breaks) - 1)]), 
                    y <= max(breaks[2:(length(breaks) - 1)])
                ) %>%
                dplyr::summarize(score = mean(score, na.rm = TRUE), n = n()) %>%
                dplyr::mutate(seqnames1 = "Aggr. coordinates")
        }) %>%
            dplyr::bind_rows()
        
        ## -- Calculate scores over each coords
        mats <- mats %>% 
            dplyr::group_by(x, y) %>%
            dplyr::summarize(score = mean(score, na.rm = TRUE)) %>%
            dplyr::mutate(seqnames1 = "Aggr. coordinates")

        ## -- Matrix limits
        if (!is.null(limits)) {
            m <- limits[1]
            M <- limits[2]
        }
        else {
            M <- max(mats$score, na.rm = TRUE)
            m <- min(mats$score, na.rm = TRUE)
            limits <- c(m, M)
        }

        ## -- Clamp scores to limits
        mats <- dplyr::mutate(mats, score = ifelse(score > M, M, ifelse(score < m, m, score)))

        ## -- Add lower triangular matrix scores
        if (symmetrical) {
            mats <- rbind(mats, mats %>% dplyr::mutate(x2 = y, y = x, x = x2) %>%
                dplyr::select(-x2))
        }

        ## -- Plot matrix
        p <- ggMatrix(mats, cols = cmap, limits = limits) +
            plotFun +
            ggplot2::labs(
                x = unique(mats$seqnames1),
                y = unique(mats$seqnames1)
            )

        p
    }
}

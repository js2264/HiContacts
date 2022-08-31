#' getPs
#'
#' @param x A `contacts` object
#' @param by_chr by_chr
#' @param filtered_chr filtered_chr
#' @return a tibble
#'
#' @import tibble
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr left_join
#' @importFrom dplyr group_split
#' @importFrom dplyr slice
#' @importFrom dplyr pull
#' @importFrom dplyr bind_rows
#' @importFrom dplyr select
#' @importFrom dplyr tally
#' @importFrom dplyr arrange
#' @importFrom dplyr lead
#' @rdname Ps
#' @export
#' @examples 
#' library(HiContacts)
#' data(contacts_yeast)
#' getPs(contacts_yeast)

getPs <- function(
    x, 
    by_chr = FALSE, 
    filtered_chr = c('XII', 'chrXII', 'chr12', '12', 'Mito', 'MT', 'chrM')
) {
    pairsFile <- pairsFile(x)
    if (is.null(pairsFile)) {
        # stop("Please provide a pairsFile for `x`. Aborting now.")
        message("pairsFile not specified. The P(s) curve will be an approximation.")
        pairs <- assay(x, 'raw')
        df <- tibble::tibble(
            chr = as.vector(GenomeInfoDb::seqnames(InteractionSet::anchors(pairs)[[1]])),
            distance = InteractionSet::pairdist(pairs, type = 'gap'),
            n = pairs$score
        ) %>% 
            tidyr::drop_na() %>% 
            dplyr::filter(!chr %in% filtered_chr) %>% 
            dplyr::mutate(binned_distance = breaks$break_pos[findInterval(distance, vec = breaks$break_pos, all.inside = TRUE)])
        if (by_chr) {
            df <- dplyr::group_by(df, chr, binned_distance)
        } 
        else {
            df <- dplyr::group_by(df, binned_distance)
        }
        d <- dplyr::summarize(df, ninter = sum(n)) %>%
            dplyr::mutate(p = ninter/sum(ninter)) %>% 
            dplyr::left_join(breaks, by = c('binned_distance' = 'break_pos')) %>% 
            dplyr::mutate(norm_p = p / binwidth)
        if (by_chr) {
            d <- dplyr::group_by(d, chr)
        } 
        else {
            d <- d
        }
        ps <- dplyr::group_split(d) %>% 
            purrr::map(function(x) {
                x %>% 
                    dplyr::mutate(
                        norm_p_unity = norm_p / 
                            {dplyr::slice(x, which.min(abs(x$binned_distance - 100000))) %>% dplyr::pull(norm_p)}
                    ) %>% 
                    dplyr::mutate(
                        slope = (log10(dplyr::lead(norm_p)) - log10(norm_p)) / 
                            (log10(dplyr::lead(binned_distance)) - log10(binned_distance))
                    ) %>% 
                    dplyr::mutate(
                        slope = c(0, predict(
                            loess(slope ~ binned_distance, span = 0.5, data = .)
                        ))
                    )
            }) %>% 
            dplyr::bind_rows()
        if (by_chr) {
            ps <- dplyr::select(ps, chr, binned_distance, p, norm_p, norm_p_unity, slope) %>% 
                dplyr::arrange(chr, binned_distance)
        } 
        else {
            ps <- dplyr::select(ps, binned_distance, p, norm_p, norm_p_unity, slope) %>% 
                dplyr::arrange(binned_distance)
        }
        return(ps)
    }
    else {
        message("Importing pairs file ", pairsFile, " in memory. This may take a while...")
        pairs <- pairs2gi(pairsFile)
        df <- tibble::tibble(
            chr = as.vector(GenomeInfoDb::seqnames(InteractionSet::anchors(pairs)[[1]])),
            distance = pairs$distance
        ) %>% 
            tidyr::drop_na() %>% 
            dplyr::filter(!chr %in% filtered_chr) %>% 
            dplyr::mutate(binned_distance = breaks$break_pos[findInterval(distance, vec = breaks$break_pos, all.inside = TRUE)])
        if (by_chr) {
            df <- dplyr::group_by(df, chr, binned_distance)
        } 
        else {
            df <- dplyr::group_by(df, binned_distance)
        }
        d <- dplyr::tally(df, name = 'ninter') %>%
            dplyr::mutate(p = ninter/sum(ninter)) %>% 
            dplyr::left_join(breaks, by = c('binned_distance' = 'break_pos')) %>% 
            dplyr::mutate(norm_p = p / binwidth)
        if (by_chr) {
            d <- dplyr::group_by(d, chr)
        } 
        else {
            d <- d
        }
        ps <- dplyr::group_split(d) %>% 
            purrr::map(function(x) {
                dplyr::mutate(x, norm_p_unity = norm_p / {dplyr::slice(x, which.min(abs(x$binned_distance - 100000))) %>% dplyr::pull(norm_p)}) %>% 
                dplyr::mutate(slope = (log10(dplyr::lead(norm_p)) - log10(norm_p)) / (log10(dplyr::lead(binned_distance)) - log10(binned_distance))) %>% 
                dplyr::mutate(slope = c(0, predict(loess(slope ~ binned_distance, span = 0.5, data = .))))
            }) %>% 
            dplyr::bind_rows()
        if (by_chr) {
            ps <- dplyr::select(ps, chr, binned_distance, p, norm_p, norm_p_unity, slope) %>% 
                dplyr::arrange(binned_distance)
        } 
        else {
            ps <- dplyr::select(ps, binned_distance, p, norm_p, norm_p_unity, slope) %>% 
                dplyr::arrange(binned_distance)
        }
        return(ps)
    }
}

#' plotPs
#'
#' @param ... ...
#' @param xlim xlim
#' @param ylim ylim
#' @return ggplot
#'
#' @import ggplot2
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @rdname Ps
#' @export
#' @examples 
#' 
#' ## Single P(s)
#' 
#' library(HiContacts)
#' data(contacts_yeast)
#' ps <- getPs(contacts_yeast)
#' plotPs(ps, aes(x = binned_distance, y = norm_p))
#' 
#' ## Comparing several P(s)
#' 
#' library(HiContacts)
#' data(contacts_yeast)
#' data(contacts_yeast_eco1)
#' ps_wt <- getPs(contacts_yeast)
#' ps_wt$sample <- 'WT'
#' ps_eco1 <- getPs(contacts_yeast_eco1)
#' ps_eco1$sample <- 'eco1'
#' ps <- bind_rows(ps_wt, ps_eco1)
#' plotPs(ps, aes(x = binned_distance, y = norm_p, group = sample, color = sample))

plotPs <- function(..., xlim = c(5000, 4.99e5), ylim = c(1e-8, 1e-4)) {
    gg <- ggplot2::ggplot(...) + 
        ggplot2::geom_line() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_y_log10(
            limits = ylim, 
            expand = c(0, 0), 
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        ggplot2::scale_x_log10(
            limits = xlim, 
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        ) + 
        ggplot2::annotation_logticks() + 
        ggplot2::labs(x = "Genomic distance", y = "Contact frequency")
    gg
}

#' plotPsSlope
#'
#' @return ggplot
#' 
#' @import ggplot2
#' @rdname Ps
#' @export
#' @examples 
#' library(HiContacts)
#' data(contacts_yeast)
#' ps <- getPs(contacts_yeast)
#' plotPsSlope(ps, aes(x = binned_distance, y = slope))

plotPsSlope <- function(..., xlim = c(5000, 4.99e5), ylim = c(-3, 0)) {
    gg <- ggplot2::ggplot(...) + 
        ggplot2::geom_line() + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::scale_y_continuous(
            limits = ylim, 
            expand = c(0, 0)
        ) +
        ggplot2::scale_x_log10(
            limits = xlim, 
            expand = c(0, 0), 
            breaks = c(1, 10, 100, 1000, 10000, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09, 1e+10),
            labels = c('1', '10', '100', '1kb', '10kb', '100kb', '1Mb', '10Mb', '100Mb', '1Gb', '10Gb')
        ) + 
        ggplot2::annotation_logticks(sides = "b") + 
        ggplot2::labs(x = "Genomic distance", y = "Slope of P(s)")
    gg
}

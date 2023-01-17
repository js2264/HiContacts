#' Compute a scalogram of contacts
#' 
#' @name scalogram
#' 
#' @param x A `HiCExperiment` object
#' @param dist_min Minimum distance for interactions to be considered.
#' @param nbins Number of bins to divide each chromosome
#' @param probs Quantiles of interactions 
#' @return a tibble
#' 
#' @return a tibble
#'
#' @importFrom tibble tibble
#' @importFrom tidyr drop_na
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr group_modify
#' 
#' @examples 
#' contacts_yeast <- HiCExperiment::contacts_yeast()
#' pairsFile(contacts_yeast) <- HiContactsData::HiContactsData(
#'   'yeast_wt', format = 'pairs.gz'
#' )
#' scalo <- scalogram(contacts_yeast['II'])
#' scalo
NULL

#' @rdname scalogram
#' @export

scalogram <- function(x, dist_min = 0, nbins = 250, probs = c(0.25, 0.5, 0.75)) {
    pairsFile <- HiCExperiment::pairsFile(x)
    if (is.null(pairsFile)) {
        stop("No PairsFile is associated with the provided HiCExperiment object.")
    }
    message("Importing pairs file ", pairsFile, " in memory. This may take a while...")
    pairs <- BiocIO::import(pairsFile, format = 'pairs')
    an_ <- HiCExperiment::anchors(pairs)
    scalo <- tibble::tibble(
        dist = pairs$distance, 
        chr = as.vector(GenomicRanges::seqnames(an_[[1]])), 
        pos1 = GenomicRanges::start(an_[[1]]),
        pos2 = GenomicRanges::start(an_[[2]])
    ) %>% 
        tidyr::drop_na() %>% 
        tidyr::pivot_longer(-c(chr, dist), names_to = 'pos_', values_to = 'pos') %>%
        dplyr::select(-pos_) %>%
        dplyr::filter(dist > dist_min, pos > 1) 
    seq_ <- seq(1, max(scalo$pos, na.rm = TRUE), length.out = nbins)
    df <- scalo %>% 
        dplyr::mutate(binned_pos = seq_[findInterval(pos, vec = seq_)]) %>% 
        dplyr::group_by(chr, binned_pos) %>% 
        dplyr::group_modify(~ {
            tibble::tibble(
                prob = factor(probs, rev(probs)), 
                dist_quantile = quantile(.x$dist, probs = probs)
            )
        }) 
    return(df)
}


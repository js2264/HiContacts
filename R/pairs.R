#' readPairs
#'
#' @param file pairs file: <readname>\t<chr1>
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
    
    headers <- data.table::fread(glue::glue("grep -v '^#' {file}"), nrows = 1)
    
    anchors1 <- data.table::fread(
        glue::glue("grep -v '^#' {file}"), 
        header = FALSE, 
        drop = {1:ncol(headers)}[!{1:ncol(headers) %in% c(chr1.field, start1.field)}], 
        nThread = nThread, 
        col.names = c('chr', 'start'), 
        nrows = nrows
    ) 
    anchors2 <- data.table::fread(
        glue::glue("grep -v '^#' {file}"), 
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

#' getPs
#'
#' @param pairs pairs in GenomicInteractions format, e.g. imported from a `.pairs` file with `pair2gi()`
#' @param by_chr by_chr
#' @param filtered_chr filtered_chr
#'
#' @import tibble
#' @import dplyr
#' @export

getPs <- function(pairs, by_chr = FALSE, filtered_chr = c('XII', 'chrXII', 'Mito', 'MT')) {
    df <- tibble(
        chr = as.vector(seqnames(anchors(pairs)[[1]])),
        distance = pairs$distance
    ) %>% 
        drop_na() %>% 
        filter(!chr %in% filtered_chr) %>% 
        mutate(binned_distance = coolerr::breaks$break_pos[findInterval(distance, vec = coolerr::breaks$break_pos)])
    if (by_chr) {
        df <- group_by(df, chr, binned_distance)
    } 
    else {
        df <- group_by(df, binned_distance)
    }
    d <- tally(df, name = 'ninter') %>%
        mutate(p = ninter/sum(ninter)) %>% 
        left_join(coolerr::breaks, by = c('binned_distance' = 'break_pos')) %>% 
        mutate(norm_p = p / binwidth)
    if (by_chr) {
        d <- group_by(d, chr)
    } 
    else {
        d <- d
    }
    ps <- group_split(d) %>% 
        map(function(x) {mutate(x, norm_p_unity = norm_p / {dplyr::slice(x, which.min(abs(x$binned_distance - 100000))) %>% pull(norm_p)})}) %>% 
        bind_rows()
    if (by_chr) {
        ps <- select(ps, chr, binned_distance, p, norm_p, norm_p_unity)
    } 
    else {
        ps <- select(ps, binned_distance, p, norm_p, norm_p_unity)
    }
    return(ps)
}

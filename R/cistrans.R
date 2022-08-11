cis_trans <- function(x) {
    if (!is.null(focus(x))) {
        stop('Please provide a contact matrix over the entire genome. Aborting now.')
    }
    cnts <- assay(x, 'raw') %>% 
        as_tibble() %>% 
        relocate(c(seqnames1, seqnames2))
    cnts_dup <- cnts %>% 
        dplyr::rename(seqnames1 = seqnames2, seqnames2 = seqnames1) %>% 
        relocate(c(seqnames1, seqnames2))
    cnts <- rbind(cnts, cnts_dup)
    res <- cnts %>% 
        group_by(seqnames1, seqnames2) %>% 
        summarize(n = sum(score)) %>% 
        mutate(type = ifelse(seqnames1 == seqnames2, 'cis', 'trans')) %>% 
        group_by(seqnames1, type) %>% 
        dplyr::rename(chr = seqnames1) %>%
        summarize(n = sum(n)) %>% 
        pivot_wider(names_from = type, values_from = n) %>% 
        mutate(
            n_total = sum(cis + trans), 
            cis_pct = cis/n_total, 
            trans_pct = trans/n_total
        )
    return(res)
}

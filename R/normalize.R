iceGis <- function(gi) {
    icedn <- reticulate::import('iced.normalization')
    cm <- gi2cm(gi)
    cm
    m <- cm@matrix
    m_toremove <- order(rowSums(is.na(m)), decreasing = TRUE)[1:{nrow(m)*0.00}]
    m_toremove <- NA
    m_norm <- log10(icedn$ICE_normalization(m[ !{1:nrow(m) %in% m_toremove} , !{1:nrow(m) %in% m_toremove} ]))
    new_cm <- cm[!{1:nrow(m) %in% m_toremove} , !{1:nrow(m) %in% m_toremove}]
    InteractionSet::as.matrix(new_cm) <- m_norm
    is <- InteractionSet::deflate(new_cm)
    new_gis <- InteractionSet::GInteractions(
        InteractionSet::anchors(is)[[1]],
        InteractionSet::anchors(is)[[2]],
        score = assay(is)[,1]
    ) 
    new_gis$bin1 <- tidyr::as_tibble(new_gis) %>% 
        dplyr::left_join(
            tidyr::as_tibble(gi), 
            by = c("seqnames1", "start1", "end1", "width1", "strand1", "seqnames2", "start2", "end2", "width2", "strand2")
        ) %>% 
        dplyr::pull(bin1)
    new_gis$bin2 <- tidyr::as_tibble(new_gis) %>% 
        dplyr::left_join(
            tidyr::as_tibble(gi), 
            by = c("seqnames1", "start1", "end1", "width1", "strand1", "seqnames2", "start2", "end2", "width2", "strand2")
        ) %>% 
        dplyr::pull(bin2)
    return(new_gis)
}
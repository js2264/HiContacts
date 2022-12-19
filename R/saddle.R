################################################################################
#                                                                              #
#                                 Saddle plots                                 #
#                                                                              #
################################################################################

plotSaddle <- function(x, nbins = 51, limits = c(-1, 1), BPPARAM = BiocParallel::bpparam()) {
    if (is.null(metadata(x)$eigens)) stop(
        "No eigen vector found in metadata. Run getCompartments(x) first."
    )
    eigens <- metadata(x)$eigens
    BPPARAM <- BiocParallel::SerialParam()

    ## -- Filter and bin regions by their eigenvector score
    filtered_eigens <- eigens[eigens$eigen != 0]
    filtered_eigens <- filtered_eigens[filtered_eigens$eigen >= quantile(filtered_eigens$eigen, 0.025) & filtered_eigens$eigen <= quantile(filtered_eigens$eigen, 0.975)]
    filtered_eigens$eigen_bin <- cut(
        filtered_eigens$eigen, breaks = quantile(
            filtered_eigens$eigen, probs = seq(0, 1, length.out = nbins+1)
        )
    ) |> as.numeric()

    ## -- Compute over-expected score in each pair of eigen bins
    df <- BiocParallel::bplapply(
        BPPARAM = BPPARAM, 
        seqnames(seqinfo(filtered_eigens)), 
        function(chr) {.saddle(x[chr], filtered_eigens)}
    ) |> dplyr::bind_rows()
    dat <- df |> 
        dplyr::group_by(eigen_bin1, eigen_bin2) |> 
        dplyr::summarize(score = mean(detrended), .groups = 'drop') |> 
        dplyr::mutate(
            x = eigen_bin1 / nbins,
            y = eigen_bin2 / nbins, 
            squished_score = scales::oob_squish(score, limits)
        )

    ## -- Make saddle plot 
    p1 <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) + 
        ggrastr::geom_tile_rast(ggplot2::aes(fill = squished_score)) + 
        ggplot2::scale_y_reverse(labels = scales::percent, expand = c(0, 0)) + 
        ggplot2::scale_x_continuous(
            labels = scales::percent, expand = c(0, 0), 
            position = 'top'
        ) + 
        ggplot2::coord_fixed() + 
        ggplot2::scale_fill_gradientn(
            colors = bgrColors(),
            na.value = "#FFFFFF",
            limits = limits
        ) + 
        ggplot2::labs(
            x = '', 
            y = 'Genomic bins ranked by eigenvector value', 
            fill = "Interaction frequency\n(log2 over expected)"
        ) + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(legend.position = 'bottom') +
        ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)) + 
        ggplot2::theme(axis.title.x = ggplot2::element_blank()) 
    
    ## -- Make side barplot 
    dat2 <- tibble::as_tibble(filtered_eigens) |> 
        dplyr::group_by(eigen_bin) |> 
        dplyr::summarize(mean_eigen = mean(eigen)) |> 
        dplyr::mutate(x = eigen_bin / nbins)
    p2 <- ggplot2::ggplot(dat2, ggplot2::aes(x = x, y = mean_eigen)) + 
        ggplot2::geom_col() + 
        ggplot2::scale_x_continuous(
            labels = NULL
        ) + 
        ggplot2::labs(
            x = '', 
            y = 'Eigen', 
            fill = "Interaction frequency\n(log2 over expected)"
        ) + 
        ggplot2::theme_minimal() + 
        ggplot2::theme(legend.position = 'bottom') +
        ggplot2::theme(panel.border = ggplot2::element_blank()) + 
        ggplot2::theme(panel.grid = ggplot2::element_blank()) + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
        ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
        ggplot2::theme(axis.line.x = ggplot2::element_blank()) 
    
    p <- cowplot::plot_grid(p2, p1, ncol = 1, rel_heights = c(1, 4), axis = 'lr', align = 'h')

}

.saddle <- function(x_chr, filtered_eigens) {
    dx <- detrend(x_chr)
    df <- as(dx, 'data.frame')
    df$eigen_bin1 <- df |>
        dplyr::left_join(as.data.frame(filtered_eigens), by = c('bin_id1' = 'bin_id')) |> 
        dplyr::pull(eigen_bin)
    df$eigen_bin2 <- df |>
        dplyr::left_join(as.data.frame(filtered_eigens), by = c('bin_id2' = 'bin_id')) |> 
        dplyr::pull(eigen_bin)
    df <- drop_na(df, eigen_bin1, eigen_bin2, detrended) |> 
        dplyr::select(eigen_bin1, eigen_bin2, detrended)
    return(df)
}

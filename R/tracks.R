#' @title Aligning tracks with contact maps
#' @name tracks
#' @rdname tracks
#' @description
#' 
#' Aligning tracks with contact maps
#'
#' @param x A `HiCExperiment` object over a full genome
#' @param bw a bigwig track imported as 'Rle'
#' @param track.name Name to give to the track in regions metadata
#' @return A `HiCExperiment` object with 2 added columns in `regions(x)`

alignTrack <- function(x, bw, track.name = 'track') {
    re <- regions(interactions(x))
    mc <- S4Vectors::mcols(re)
    mc[[track.name]] <- bw[re]
    mc[[paste0(track.name, '.mean')]] <- BiocGenerics::mean(
        mc[[track.name]], na.rm = TRUE
    )
    S4Vectors::mcols(regions(interactions(x))) <- mc
    return(x)
}


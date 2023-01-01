#' @title Generic functions
#' 
#' @name AllGenerics
#' @aliases plotMatrix
#' 
#' @description
#' 
#' Generics functions created in HiContacts package. 
#' 
#' @param x Can be a `HiCExperiment`, a GInteractions or an array/matrix object. 
#' @param use.scores use.scores
#' @param scale scale
#' @param limits limits
#' @param max.distance maximum distance. If provided, the heatmap is plotted 
#'   horizontally. 
#' @param loops loops
#' @param borders borders
#' @param dpi dpi
#' @param rasterize rasterize
#' @param symmetrical symmetrical
#' @param chrom_lines chrom_lines
#' @param show_grid show_grid
#' @param cmap color map
NULL

setGeneric("plotMatrix", function(
    x, 
    use.scores = NULL, 
    scale = 'log10', 
    max.distance = NULL, 
    loops = NULL, 
    borders = NULL, 
    limits = NULL, 
    dpi = 500, 
    rasterize = TRUE, 
    symmetrical = TRUE, 
    chrom_lines = TRUE, 
    show_grid = FALSE, 
    cmap = NULL  
) {standardGeneric("plotMatrix")})

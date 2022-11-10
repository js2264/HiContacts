setClassUnion("GRangesOrGInteractions", members = c("GRanges", "GInteractions"))
setClassUnion("GRangesOrPairsOrcharacterOrNULL", members = c("GRanges", "Pairs", "character", "NULL"))
setClassUnion("numericOrcharacter", members = c("numeric", "character"))
setClassUnion("characterOrNULL", members = c("character", "NULL"))

#' @title `Contacts` S4 class and methods
#'
#' @slot fileName Path of (m)cool file
#' @slot focus Chr. coordinates for which interaction counts are extracted 
#'   from the .(m)cool file.
#' @slot resolutions Resolutions available in the .(m)cool file.
#' @slot resolution Current resolution
#' @slot interactions Genomic Interactions extracted from the .(m)cool object
#' @slot scores Available interaction scores. 
#' @slot topologicalFeatures Topological features associated with the dataset 
#'   (e.g. loops (\<Pairs\>), borders (\<GRanges\>), 
#'   viewpoints (\<GRanges\>), etc...)
#' @slot pairsFile Path to the .pairs file associated with the .(m)cool file
#' @slot metadata metadata associated with the .(m)cool file.
#' 
#' @import methods
#' @importClassesFrom S4Vectors Pairs
#' @importClassesFrom S4Vectors Annotated
#' @importFrom methods setClass
#' @export
#' @rdname Contacts 

methods::setClass("Contacts", 
    contains = c("Annotated"), 
    slots = c(
        fileName = "character",
        focus = "characterOrNULL", 
        resolutions = "numeric", 
        resolution = "numeric", 
        interactions = "GInteractions",
        scores = "SimpleList", 
        topologicalFeatures = "SimpleList",
        pairsFile = "characterOrNULL",
        metadata = "list"
    )
)

#' @title Finding loops in contact map
#' @name getLoops
#' @rdname getLoops
#' @description
#' 
#' Find loops using chromosight.  
#' 
#' This function is actually provided by the `HiCool` package rather than 
#' the `HiContacts` package. `HiCool` provides a self-managed
#' `conda` environment, and this limits 
#' 
#' @param ... Parameters passed to `HiCool::getLoops()`.
#' @export

getLoops <- function(...){
    stop("`getLoops` is provided by the `HiCool` package. Please install HiCool and use `HiCool::getLoops()`.")
}

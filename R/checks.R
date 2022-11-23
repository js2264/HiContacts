#' @title Checks functions
#' 
#' @name checks
#' @rdname checks
#' 
#' @description 
#' 
#' Useful functions to validate the nature/structure of (m)cool files or 
#' `HiCExperiment` objects.
#'  All these check functions should return a logical.
#' 
#' @param x A `HiCExperiment` object
#' @param ... `HiCExperiment` objects
#' @return Logical
NULL

#' @rdname checks

is_symmetrical <- function(x) {
    if (is.null(focus(x))) {
        return(TRUE)
    }
    if (grepl(' x ', focus(x))) {
        return(FALSE)
    }
    else {
        return(TRUE)
    }
}

#' @rdname checks

is_comparable <- function(...) {
    are_HiCExperiment(...)
    err <- c()
    if (!is_same_seqinfo(...)) {
        err <- c(err, "seqinfos")
    }
    if (!is_same_resolution(...)) {
        err <- c(err, "resolutions")
    }
    if (!is_same_bins(...)) {
        err <- c(err, "bins")
    }
    if (!is_same_regions(...)) {
        err <- c(err, "regions")
    }
    if (length(err) > 0) {
        mess <- paste0("Provided `HiCExperiment` have different ", paste(err, collapse = ' & '), '.')
        stop(mess)
    }
    TRUE
}

#' @rdname checks

are_HiCExperiment <- function(...) {
    args <- list(...)
    if (!all(unlist(lapply(args, is, 'HiCExperiment')))) {
        stop("Provided arguments are not all `HiCExperiment` objects. 
        Please only use `HiCExperiment` objects with this function.")
    }
    TRUE
}

#' @rdname checks

is_same_seqinfo <- function(...) {
    contacts_list <- list(...)
    all(unlist(lapply(contacts_list, function(x) {
        identical(seqinfo(contacts_list[[1]]), seqinfo(x))
    })))
}

#' @rdname checks

is_same_resolution <- function(...) {
    contacts_list <- list(...)
    all(unlist(lapply(contacts_list, function(x) {
        identical(resolution(contacts_list[[1]]), resolution(x))
    })))
}

#' @rdname checks

is_same_bins <- function(...) {
    contacts_list <- list(...)
    all(unlist(lapply(contacts_list, function(x) {
        identical(bins(contacts_list[[1]]), bins(x))
    })))
}

#' @rdname checks

is_same_regions <- function(...) {
    contacts_list <- list(...)
    all(unlist(lapply(contacts_list, function(x) {
        identical(regions(contacts_list[[1]]), regions(x))
    })))
}


#' Checks functions
#' 
#' Useful functions to validate the nature/structure of (m)cool files or 
#' `contacts` objects.
#' 
#' @param path Path of a (m)cool file
#' @param contacts A `contacts` object
#' @param resolution Resolution
#' @param pair Pairs object with length of 1
#' @param ... `contacts` object
#' @return Logical
#' 
#' @rdname checks

check_cool_file <- function(path) {
    if (!file.exists(path) | {
        isTRUE(nzchar(Sys.readlink(path), keepNA = TRUE)) & 
        !file.exists(Sys.readlink(path))
    }) {
        stop('File not found. Aborting now')
    }
    if (!{is_cool(path) | is_mcool(path)}) {
        stop('Provided file is not a .cool/.mcool file.\n  Aborting now.')
    }
    TRUE
}

#' @rdname checks

check_resolution <- function(contacts, resolution) {
    available_res <- resolutions(contacts)
    if (!resolution %in% available_res) 
        stop("Resolution not stored in the cool file. Aborting now.")
    TRUE
}

#' @rdname checks

check_cool_format <- function(path, resolution) {
    if (is_mcool(path)) {
        if (is.null(resolution)) {
            stop("File is in .mcool format, a resolution must be provided. Aborting now.")
        }
        if (!resolution %in% lsCoolResolutions(path)) {
            stop("Resolution not stored in cool file. Aborting now.")
        }
    }
    if (is_cool(path)) {
        if (!is.null(resolution)) {
            stop("File is in .cool format, please do not specify any resolution. Aborting now.")
        }
    }
    TRUE
}

#' @rdname checks

is_mcool <- function(path) {
    x <- lsCoolFiles(path)
    all(grepl('^/resolutions', x))
}

#' @rdname checks

is_cool <- function(path) {
    x <- lsCoolFiles(path)
    all(!grepl('^/resolutions', x))
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

#' @rdname checks

is_comparable <- function(...) {
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
        mess <- paste0("Provided contacts have different ", paste(err, collapse = ' & '), '.')
        stop(mess)
    }
    TRUE
}

#' @rdname checks

is_square <- function(pair) {
    w1 <- GenomicRanges::width(S4Vectors::first(pair))
    w2 <- GenomicRanges::width(S4Vectors::second(pair))
    if (w1 != w2) {
        stop("Provided pair is not square.")
    }
    TRUE
}

#' @rdname checks
#' @export
#' @examples
#' library(HiContacts)
#' contacts_yeast <- contacts_yeast()
#' is_symmetrical(contacts_yeast)

is_symmetrical <- function(contacts) {
    if (is.null(focus(contacts))) {
        return(TRUE)
    }
    if (grepl(' x ', focus(contacts))) {
        return(FALSE)
    }
    else {
        return(TRUE)
    }
}

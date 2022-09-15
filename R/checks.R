#' @name Checks functions
#' @title Checks functions
#' 
#' @param cool_path Path of a (m)cool file
#' @param contacts A `contacts` object
#' @param resolution Resolution
#' @param pair Pairs object with length of 1
#' @param ... `contacts` object
#' @return Logical
#' 
#' @rdname checks

check_cool_file <- function(cool_path) {
    if (!file.exists(cool_path) | {
        isTRUE(nzchar(Sys.readlink(cool_path), keepNA = TRUE)) & 
        !file.exists(Sys.readlink(cool_path))
    }) {
        stop('File not found. Aborting now')
    }
    if (!{is_cool(cool_path) | is_mcool(cool_path)}) {
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

check_cool_format <- function(cool_path, resolution) {
    if (is_mcool(cool_path)) {
        if (is.null(resolution)) {
            stop("File is in .mcool format, a resolution must be provided. Aborting now.")
        }
        if (!resolution %in% lsCoolResolutions(cool_path, verbose = FALSE)) {
            stop("Resolution not stored in cool file. Aborting now.")
        }
    }
    if (is_cool(cool_path)) {
        if (!is.null(resolution)) {
            stop("File is in .cool format, please do not specify any resolution. Aborting now.")
        }
    }
    TRUE
}

#' @rdname checks

is_mcool <- function(cool_path) {
    x <- lsCoolFiles(cool_path)
    all(grepl('^/resolutions', x))
}

#' @rdname checks

is_cool <- function(cool_path) {
    x <- lsCoolFiles(cool_path)
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

is_centered <- function(contacts) {
    if (is.character(focus(contacts))) {
        if (grepl(' x ', focus(contacts))) {
            return(TRUE)
        }
        else {
            return(FALSE)
        }
    }
    else if (methods::is(focus(contacts), 'Pairs')) {
        return(TRUE)
    }
}

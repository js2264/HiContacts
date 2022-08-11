check_cool_file <- function(cool_path) {
    if (!file.exists(cool_path)) {
        stop('File not found. Aborting now')
    }
    if (!{is_cool(cool_path) | is_mcool(cool_path)}) {
        stop('Provided file is not a .cool/.mcool file.\n  Aborting now.')
    }
    TRUE
}

check_resolution <- function(contacts, resolution) {
    available_res <- resolutions(contacts)
    if (!resolution %in% available_res) 
        stop("New resolution not stored in the cool file. Aborting now.")
    TRUE
}

check_cool_format <- function(cool_path, resolution) {
    if (is_mcool(cool_path)) {
        if (is.null(resolution)) {
            stop("File is in .mcool format, a resolution must be provided. Aborting now.")
        }
    }
    if (is_cool(cool_path)) {
        if (!is.null(resolution)) {
            stop("File is in .cool format, please do not specify any resolution. Aborting now.")
        }
    }
    TRUE
}

is_mcool <- function(cool_path) {
    tools::file_ext(cool_path) == 'mcool'
}

is_cool <- function(cool_path) {
    tools::file_ext(cool_path) == 'cool'
}

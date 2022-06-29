check_cool_file <- function(cool_path) {
    ext <- gsub('.*\\.', '', cool_path)
    if (!file.exists(cool_path)) {
        stop('File not found. Aborting now')
    }
    if (!ext %in% c('mcool', 'cool')) {
        stop('Provided file is not a .cool/.mcool file.\n  Aborting now')
    }
    TRUE
}

check_resolution <- function(contacts, resolution) {
    available_res <- resolutions(contacts)
    if (!resolution %in% available_res) 
        stop("New resolution not stored in the cool file")
    TRUE
}

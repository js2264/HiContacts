#' Example datasets provided in `HiContacts` & `HiContactsData`
#' 
#' @format An object of class \code{"GRanges"}.
#' @docType data
#' @usage data(centros_yeast)
#' @source HiContacts
#' @rdname datasets
#' @examples
#' data(centros_yeast)
#' centros_yeast
"centros_yeast"

#' @rdname datasets
#' 
#' @format An object of class \code{"contacts"}.
#' @importFrom HiContactsData HiContactsData
#' @export
#' @examples
#' contacts_yeast()
#' contacts_yeast_eco1()
#' full_contacts_yeast()

contacts_yeast <- function() {
    fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    contacts(fpath, 'II', resolution = 1000)
}

#' @rdname datasets
#' @export

contacts_yeast_eco1 <- function() {
    fpath <- HiContactsData::HiContactsData('yeast_eco1', 'mcool')
    contacts(fpath, 'II', resolution = 1000)
}

#' @rdname datasets
#' @export

full_contacts_yeast <- function() {
    data(centros_yeast)
    fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    x <- contacts(fpath, resolution = 16000)
    topologicalFeatures(x, 'centromeres') <- centros_yeast
    return(x)
}

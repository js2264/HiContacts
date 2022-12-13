#' Example datasets provided in `HiContacts` & `HiContactsData`
#' 
#' @format An object of class \code{"GRanges"}.
#' @docType data
#' @usage data(centros_yeast)
#' @source HiContacts
#' @rdname datasets
#' @importClassesFrom GenomicRanges GRanges
#' @examples
#' data(centros_yeast)
#' centros_yeast
"centros_yeast"

#' @rdname datasets
#' 
#' @format An object of class \code{"Contacts"}.
#' @export
#' @examples
#' contacts_yeast()
#' contacts_yeast_eco1()
#' full_contacts_yeast()

contacts_yeast <- function() {
    .Deprecated(
        "HiCExperiment::contacts_yeast", 
        msg = "`contacts_yeast` is deprecated; see '?HiCExperiment::contacts_yeast' instead."
    )
    HiCExperiment::contacts_yeast()
}

#' @rdname datasets
#' @export

contacts_yeast_eco1 <- function() {
    .Deprecated(
        "HiCExperiment::contacts_yeast_eco1", 
        msg = "`contacts_yeast_eco1` is deprecated; see '?HiCExperiment::contacts_yeast_eco1' instead."
    )
    HiCExperiment::contacts_yeast_eco1()
}

#' @rdname datasets
#' @export

full_contacts_yeast <- function() {
    .Deprecated(
        "HiCExperiment::full_contacts_yeast", 
        msg = "`full_contacts_yeast` is deprecated; see '?HiCExperiment::full_contacts_yeast' instead."
    )
    HiCExperiment::full_contacts_yeast()
}

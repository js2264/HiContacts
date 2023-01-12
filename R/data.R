#' Example datasets provided in `HiContacts` & `HiContactsData`
#' 
#' @docType data
#' @usage contacts_yeast()
#' @usage contacts_yeast_eco1()
#' @usage full_contacts_yeast()
#' @source HiContacts
#' @rdname datasets
#' @importClassesFrom GenomicRanges GRanges
#' 
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
    HiCExperiment::contacts_yeast(full = TRUE)
}

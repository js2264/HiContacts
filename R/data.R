#' Data provided in HiContacts
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

#' @format An object of class \code{"contacts"}.
#' @docType data
#' @usage data(contacts_yeast)
#' @source HiContactsData
#'   fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#'   contacts_yeast <- contacts(fpath, 'II', resolution = 1000)
#'   usethis::use_data(contacts_yeast, overwrite = TRUE)
#' @rdname datasets
#' @examples
#' data(contacts_yeast)
#' contacts_yeast
"contacts_yeast"

#' @format An object of class \code{"contacts"}.
#' @docType data
#' @usage data(contacts_yeast_eco1)
#' @source HiContactsData
#'   fpath <- HiContactsData::HiContactsData('yeast_eco1', 'mcool')
#'   contacts_yeast_eco1 <- contacts(fpath, 'II', resolution = 1000)
#'   usethis::use_data(contacts_yeast_eco1, overwrite = TRUE)
#' @rdname datasets
#' @examples
#' data(contacts_yeast_eco1)
#' contacts_yeast_eco1
"contacts_yeast_eco1"

#' @format An object of class \code{"contacts"}.
#' @docType data
#' @usage data(full_contacts_yeast)
#' @source HiContactsData
#'   data(centros_yeast)
#'   fpath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#'   full_contacts_yeast <- contacts(fpath, resolution = 16000)
#'   features(full_contacts_yeast, 'centromeres') <- centros_yeast
#'   usethis::use_data(full_contacts_yeast, overwrite = TRUE)
#' @rdname datasets
#' @examples
#' data(full_contacts_yeast)
#' full_contacts_yeast
"full_contacts_yeast"

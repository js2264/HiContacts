#' @rdname Contacts
#' 
#' @name setAs
#' @docType methods
#' @aliases setAs,Contacts-method
#'
#' @export
#' @examples 
#' as(contacts_yeast, 'GInteractions')
#' as(contacts_yeast, 'ContactMatrix')
#' as(contacts_yeast, 'matrix')[seq_len(10), seq_len(10)]
#' as(contacts_yeast, 'data.frame')

setAs("Contacts", "GInteractions", function(from) interactions(from))
setAs("Contacts", "ContactMatrix", function(from) {
    if ('balanced' %in% names(scores(from))) {
        x <- interactions(from)
        x$score <- scores(from, 'balanced')
        gi2cm(x)
    } 
    else {
        x <- interactions(from)
        x$score <- scores(from, 1)
        gi2cm(x)
    }
})
setAs("Contacts", "matrix", function(from) {
    as(from, "ContactMatrix") |> cm2matrix()
})
setAs("Contacts", "data.frame", function(from) {
    x <- interactions(from)
    x <- as.data.frame(x)
    x <- x[, !colnames(x) %in% c("chr1", "chr2", "bin_id1.1", "bin_id2.1")]
    for (n in names(scores(from))) {
        x[[n]] <- scores(from, n)
    }
    return(x)
})

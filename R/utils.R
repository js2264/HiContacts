`sdiag<-` <- function(A, k = 0, value) {
    p <- ncol(A)
    n <- nrow(A)
    if (k>p-1||-k > n-1) return()
    if (k >= 0) {
        i <- seq_len(n)
        j <- (k+1):p
    } 
    else {
        i <- (-k+1):n
        j <- seq_len(p)
    }
    if (length(i)>length(j)) i <- i[seq_along(j)] else j <- j[seq_along(i)]
    ii <- i + (j-1) * n 
    A[ii] <- value
    A
} 

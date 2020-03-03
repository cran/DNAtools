#' Collapse m/p output to vector
#' 
#' Collapse a m/p-matrix from dbCompare/dbExpect to a vector.
#' 
#' Collapse a m/p-matrix from dbCompare/dbExpect to a vector with entry i being
#' the sum of all entries from m/p-matrix satisfying 2*m+p=i.
#' 
#' @param x Either a object of class 'dbcompare' (result from dbCompare) or
#' 'matrix'.
#' @return A vector of length 2*max(m)+1 with entries begin the sum of entries
#' i in m/p-matrix satisfying i=2*m+p.
#' @author Torben Tvedebrink
#' @examples
#' 
#'   \dontrun{
#'   data(dbExample)
#'   res <- dbCompare(dbExample, hit=5, trace=TRUE)
#'   dbCollapse(res) ## same as dbCompare(dbExample, hit=5, trace=TRUE, collapse=TRUE)
#'   }
#' 
#' @export dbCollapse
dbCollapse <- function(x) {
  if (inherits(x, "dbcompare")) 
    mpmatrix <- x$m else if (!inherits(x, "matrix")) 
    stop("Input must be of class 'matrix' or 'dbcompare'") else mpmatrix <- x
  res <- sapply(1:((nL <- (ncol(mpmatrix) - 1)) * 2 + 1), function(i) sum(diag(mpmatrix[(MP <- mpcollapse(i - 
    1, nL))$m2 + 1, , drop = FALSE][, MP$m1 + 1, drop = FALSE]), na.rm = TRUE))
  names(res) <- 0:(length(res) - 1)
  res
}

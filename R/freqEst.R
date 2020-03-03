#' Simple allele frequency estimation
#' 
#' Estimates allele frequencies from a database with DNA profiles
#' 
#' Computes the allele frequencies for a given database.
#' 
#' @param x A database of the form ['id','locus1 allele1','locus1
#' allele2',...,'locusN allele 1','locusN allele2'].
#' @return Returns a list of probability vectors - one vector for each locus.
#' @author James Curran and Torben Tvedebrink
#' @examples
#' 
#'   data(dbExample)
#'   freqEst(dbExample)
#' 
#' @export freqEst
freqEst <- function(x) {
  idcol <- seq(from = 2, by = 2, len = (ncol(x) - 1)/2)
  x1 <- do.call("cbind", lapply(x[, idcol], paste))
  x2 <- do.call("cbind", lapply(x[, idcol + 1], paste))
  names(x2) <- names(x1)
  xx <- as.data.frame(rbind(x1, x2))
  lapply(xx, function(allele) table(allele)/length(allele))
}

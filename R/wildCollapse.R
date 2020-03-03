wildCollapse <- function(x) {
  maxAlleles <- nrow(x) - 1
  wildCol <- rep(NA, maxAlleles + 1)
  for (s in 0:maxAlleles) wildCol[s + 1] <- sum(diag(as.matrix(x[0:s + 1, rev(0:s + 1)])))
  names(wildCol) <- 0:maxAlleles
  wildCol
}

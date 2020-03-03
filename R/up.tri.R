## Extracts the upper left triangle of a quadratic matrix
up.tri <- function(x, diag = TRUE, droplast = FALSE) {
  x <- as.matrix(x)
  res <- (row(x) + col(x)) - 1 <= ncol(x)
  if (!diag) 
    res[(row(x) + col(x) - 1) == ncol(x)] <- FALSE
  if (droplast) 
    res[nrow(res), 1] <- FALSE
  res
}

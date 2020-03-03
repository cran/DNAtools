## Creates a matrix/list of cell names
dbCats <- function(nloci, vector = FALSE) {
  res <- outer(0:nloci, 0:nloci, function(i, j) paste(i, j, sep = "/"))
  res[!up.tri(res)] <- NA
  if (vector) 
    res <- t(res)[up.tri(res)]
  res
}

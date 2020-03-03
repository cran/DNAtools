## Input: Variance matrix (or list of them). Output: Collapsed versions
varCollapse <- function(x, nl) {
  if (is.list(x)) {
    colx <- replicate(length(x), matrix(0, 2 * nl + 1, 2 * nl + 1), simplify = FALSE)
    names(colx) <- names(x)
  } else colx <- matrix(0, 2 * nl + 1, 2 * nl + 1)
  for (i in 0:nl) {
    for (j in 0:(nl - i)) {
      for (l in 0:nl) {
        for (k in 0:(nl - l)) {
          if (is.list(x)) {
          for (n in 1:length(x)) {
            colx[[n]][2 * i + j + 1, 2 * l + k + 1] <- colx[[n]][2 * i + j + 1, 2 * 
            l + k + 1] + x[[n]][m2v(i = i + 1, j = j + 1, L = nl), m2v(i = l + 1, 
            j = k + 1, L = nl)]
          }
          } else colx[2 * i + j + 1, 2 * l + k + 1] <- colx[2 * i + j + 1, 2 * l + k + 
          1] + x[m2v(i = i + 1, j = j + 1, L = nl), m2v(i = l + 1, j = k + 1, L = nl)]
        }
      }
    }
  }
  colx
}

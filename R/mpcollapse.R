## Collapses a (m,p)-matrix to a (2*m+p)-vector
mpcollapse <- function(mallele, nloci) {
  tmp <- data.frame(m2 = (m <- floor(mallele/2):0), m1 = mallele - 2 * m)
  tmp[tmp$m1 <= nloci & tmp$m2 <= nloci, ]
}


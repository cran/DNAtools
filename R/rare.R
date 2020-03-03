
## The recursive step of the expected value function. See eq. (2) in Tvedebrink et al.
rare <- function(q) {
  S <- ncol(q)
  M <- replicate(S, matrix(0, S + 1, S + 1), simplify = FALSE)
  M[[1]][1, 1] <- q[1, 1]  # P_{0,1}
  M[[1]][1, 2] <- q[2, 1]  # P_{1,1}
  M[[1]][2, 1] <- q[3, 1]  # P_{2,1}
  for (s in 2:S) {
    for (m in 1:(s + 1)) {
      for (p in 1:(s - m + 2)) {
        if (m == 1 & p > 1) 
          M[[s]][m, p] <- q[1, s] * M[[s - 1]][m, p] + q[2, s] * M[[s - 1]][m, p - 1] else if (m > 1 & p == 1) 
          M[[s]][m, p] <- q[1, s] * M[[s - 1]][m, p] + q[3, s] * M[[s - 1]][m - 1, p] else if (m == 1 & p == 1) 
          M[[s]][m, p] <- q[1, s] * M[[s - 1]][m, p] else M[[s]][m, p] <- q[1, s] * M[[s - 1]][m, p] + q[2, s] * M[[s - 1]][m, p - 
          1] + q[3, s] * M[[s - 1]][m - 1, p]
      }
    }
  }
  x <- M[[S]]
  dimnames(x) <- list(0:S, 0:S)
  x
}

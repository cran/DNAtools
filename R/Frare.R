## The recursive step of the expected value function. See eq. (2) in Tvedebrink et al.
Frare <- function(q) {
  S <- ncol(q)
  allCats <- as.vector(t(outer(0:S, 0:S, FUN = paste, sep = "/")))
  M <- replicate(S, matrix(0, (S + 1)^2, (S + 1)^2, dimnames = list(Genuine = paste(1:((S + 
    1)^2), allCats, sep = ":"), Wildcard = paste(1:((S + 1)^2), allCats, sep = ":"))), simplify = FALSE)
  M[[1]][1 * (S + 1) + (1 + 0), 0 * (S + 1) + (1 + 0)] <- q[1, 1]  # P_{1/0/0/0} {m/p/Fm/Fp}
  M[[1]][0 * (S + 1) + (1 + 1), 0 * (S + 1) + (1 + 0)] <- q[2, 1]  # P_{0/1/0/0}
  M[[1]][0 * (S + 1) + (1 + 0), 1 * (S + 1) + (1 + 0)] <- q[3, 1]  # P_{0/0/1/0}
  M[[1]][0 * (S + 1) + (1 + 0), 0 * (S + 1) + (1 + 1)] <- q[4, 1]  # P_{0/0/0/1}
  M[[1]][0 * (S + 1) + (1 + 0), 0 * (S + 1) + (1 + 0)] <- q[5, 1]  # P_{0/0/0/0}
  for (s in 2:S) {
    for (m in 1:(s + 1)) {
      for (p in 1:(s - m + 2)) {
        for (fm in 1:(s + 1)) {
          for (fp in 1:(s - fm + 2)) {
          if (m > 1) 
            M[[s]][(m - 1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] <- M[[s]][(m - 
            1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] + q[1, s] * M[[s - 1]][(m - 
            2) * (S + 1) + p, (fm - 1) * (S + 1) + fp]
          if (p > 1) 
            M[[s]][(m - 1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] <- M[[s]][(m - 
            1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] + q[2, s] * M[[s - 1]][(m - 
            1) * (S + 1) + (p - 1), (fm - 1) * (S + 1) + fp]
          if (fm > 1) 
            M[[s]][(m - 1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] <- M[[s]][(m - 
            1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] + q[3, s] * M[[s - 1]][(m - 
            1) * (S + 1) + p, (fm - 2) * (S + 1) + fp]
          if (fp > 1) 
            M[[s]][(m - 1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] <- M[[s]][(m - 
            1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] + q[4, s] * M[[s - 1]][(m - 
            1) * (S + 1) + p, (fm - 1) * (S + 1) + (fp - 1)]
          M[[s]][(m - 1) * (S + 1) + p, (fm - 1) * (S + 1) + fp] <- M[[s]][(m - 1) * 
            (S + 1) + p, (fm - 1) * (S + 1) + fp] + q[5, s] * M[[s - 1]][(m - 1) * 
            (S + 1) + p, (fm - 1) * (S + 1) + fp]
          }
        }
      }
    }
  }
  x <- M[[S]]
  remCats <- as.vector(t(outer(0:S, 0:S, FUN = function(x, y, n = S) (x + y) <= n)))
  dimnames(x) <- list(Genuine = allCats, Wildcard = allCats)
  x[remCats, remCats]
}

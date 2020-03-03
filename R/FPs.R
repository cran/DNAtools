## Computes the P_{m/p/Fm/Fp} for a given locus R and r not used
FPs <- function(p, t, k = rep(0, 3), r = 0, R = 0) {
  s2 <- sum(p^2)
  s3 <- sum(p^3)
  s4 <- sum(p^4)
  d <- (1 + 2 * t) * (1 + t)
  pmatch <- 2 * (t^2 * (1 - t) * (1 - s2) + 2 * t * (1 - t)^2 * (s2 - s3) + (1 - t)^3 * (s2^2 - 
    s4))
  pfmatch <- 6 * t^3 + t^2 * (1 - t) * (8 + 3 * s2) + t * (1 - t)^2 * (12 * s2 - 6 * s3) + 
    (1 - t)^3 * (4 * s3 - 3 * s4)
  p2 <- pmatch + pfmatch
  ppartial <- 4 * (t * (1 - t)^2 * (1 - 3 * s2 + 2 * s3) + (1 - t)^3 * (s2 - 2 * s3 + 2 * 
    s4 - s2^2))
  pfpartial <- t^2 * (1 - t) * (1 - s2) + t * (1 - t)^2 * (2 - 4 * s2 + 2 * s3) + (1 - t)^3 * 
    (2 * s2 - 4 * s3 + 3 * s4 - 1 * s2^2)
  p1 <- ppartial + pfpartial
  pmismatch <- (1 - t)^3 * (1 - 6 * s2 + 8 * s3 - 6 * s4 + 3 * s2^2)
  p0 <- pmismatch
  if (all(k == 0)) 
    res <- c(p0, p1, p2)/d else res <- c(k[3] * p0/d, k[2] * (1 - t)^2 * (1 - 3 * s2 + 2 * s3)/(1 + t) + k[3] * p1/d, 
    k[1] + k[2] * (2 * t^2 + 3 * t * (1 - t) + (1 - t)^2 * (3 * s2 - 2 * s3))/(1 + t) + 
      k[3] * p2/d)
  res
}

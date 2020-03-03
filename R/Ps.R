## Computes the P_{0/0}, P_{0/1}, P_{1/0} for a given locus R and r not used
Ps <- function(p, t, k = rep(0, 3), r = 0, R = 0) {
  s1 <- sum(p)
  s2 <- sum(p^2)
  s3 <- sum(p^3)
  s4 <- sum(p^4)
  d <- (1 + t) * (1 + 2 * t)
  p0 <- t^2 * (1 - t) * (1 - s2) + 2 * t * (1 - t)^2 * (1 - 2 * s2 + s3) + (1 - t)^3 * (1 - 
    4 * s2 + 4 * s3 + 2 * s2^2 - 3 * s4)
  p1 <- 8 * t^2 * (1 - t) * (1 - s2) + 4 * t * (1 - t)^2 * (1 - s3) + 4 * (1 - t)^3 * (s2 - 
    s3 - s2^2 + s4)
  p2 <- 6 * t^3 + t^2 * (1 - t) * (2 + 9 * s2) + 2 * t * (1 - t)^2 * (2 * s2 + s3) + (1 - 
    t)^3 * (2 * s2^2 - s4)
  if (all(k == 0)) 
    res <- c(p0, p1, p2)/d else res <- c(k[3] * p0/d, k[2] * (1 - t) * (1 - s2) + k[3] * p1/d, k[1] + k[2] * (t + (1 - 
    t) * s2) + k[3] * p2/d)
  res
}

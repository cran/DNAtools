context("NoA: Regression tests comparing to old version of package")

# res_3cntr_locuswise
# res_3cntr_locuswise_NoA
#source("include-test_noa.R")
# Loads/source()s helper-noa.R

res_3cntr_locuswise <- DNAtools::Pnm_all(m = 3, theta = 0, probs = freqs, locuswise = TRUE)

res_3cntr_conv <- DNAtools::Pnm_all(m = 3, theta = 0, probs = freqs, locuswise = FALSE)
res_3cntr_NoA <- DNAtools::convolve(res_3cntr_locuswise)

test_that("PnmAll: numContrib = 3, theta = 0", {
  expect_equal(OLD_res_3cntr_locuswise, res_3cntr_locuswise)
})

test_that("convolve: numContrib = 3, theta = 0", {
  expect_equal(OLD_res_3cntr_NoA, res_3cntr_NoA)
  expect_equal(res_3cntr_conv, res_3cntr_NoA)
})

# A battery...
for (i in seq_along(test_big_cache)) {
  #print(i)
  x <- test_big_cache[[i]]
  
  if (x$m >= 2) {
    next
  }

  res_locuswise <- DNAtools::Pnm_all(m = x$m, theta = x$theta, probs = x$freqs, locuswise = TRUE)
  res_NoA <- DNAtools::convolve(res_locuswise)
  res_conv <- DNAtools::Pnm_all(m = x$m, theta = x$theta, probs = x$freqs, locuswise = FALSE)
  
  test_that(paste0("Battery #", i, ": PnmAll: numContrib = ", x$m, ", theta = ", x$theta), {
    expect_equal(res_locuswise, x$OLD_res_locuswise)
  })
  
  test_that(paste0("Battery #", i, ": convolve: numContrib = ", x$m, ", theta = ", x$theta), {
    expect_equal(res_NoA, x$OLD_res_NoA)
    expect_equal(res_NoA, res_conv)
  })
}

test_that("Pnm_locus", {
  f <- c(0.0447047384895345, 0.169317566829043, 0.0614798043637582, 
         0.163182339345515, 0.191425934897317, 0.0164065573454084, 0.122803731911374, 
         0.051369281234042, 0.1338414148308, 0.0454686307532083)
  expect_equal(sum(f), 1)
  
  for (m in 1:10) {
    p <- Pnm_locus(m = m, theta = 0, alleleProbs = f)
    expect_equal(sum(p), 1, info = paste0("Pnm_locus with f for m = ", m))
  }
  
  
  f3 <- c(f, f, f)
  f3 <- f3 / sum(f3)
  for (m in 1:10) {
    p <- Pnm_locus(m = m, theta = 0, alleleProbs = f3)
    expect_equal(sum(p), 1, info = paste0("Pnm_locus with f3 for m = ", m))
  }
})

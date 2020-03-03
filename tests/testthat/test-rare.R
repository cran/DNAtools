context("Rare")

################################################################################

test_that("Rare", {
  expect_equal(RPs(p = freqs[[1]], t = 0, r = 0, R = 0, k = rep(0, 3)), 
               c(0.595219266518844, 0.37643791086126, 0.0283428226198924))
})

data(dbExample)


test_that("dbCompare rare", {
  dbExample_rare <- dbExample
  #dbExample_rare[c(1, 2), 4] <- 99
  dbExample_rare[c(1, 2, 999), c(1, 2, 4)] <- 99
  dbExample_rare[c(1, 999), c(2, 8)] <- 99

  res <- dbCompare(dbExample, hit = 5, trace = FALSE, threads = 3,
                   Rallele = TRUE)
  
  expect_equal(res$m, structure(c(102L, 206L, 165L, 72L, 22L, 6L, 0L, 0L, 0L, 0L, 0L, 
                                  1368L, 2114L, 1477L, 556L, 149L, 19L, 2L, 0L, 0L, 0L, 0L, 7122L, 
                                  10013L, 5710L, 1821L, 360L, 44L, 3L, 0L, 0L, 0L, 0L, 21878L, 
                                  26084L, 12566L, 3250L, 493L, 41L, 0L, 0L, 0L, 0L, 0L, 44189L, 
                                  43656L, 17049L, 3361L, 379L, 26L, 0L, 0L, 0L, 0L, 0L, 59463L, 
                                  47418L, 14642L, 2135L, 156L, 5L, 0L, 0L, 0L, 0L, 0L, 54601L, 
                                  34320L, 7570L, 719L, 34L, 0L, 0L, 0L, 0L, 0L, 0L, 34203L, 15463L, 
                                  2220L, 116L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 13571L, 4145L, 310L, 
                                  0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 3281L, 472L, 0L, 0L, 0L, 0L, 
                                  0L, 0L, 0L, 0L, 0L, 353L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
                                  0L), .Dim = c(11L, 11L), .Dimnames = list(match = c("0", "1", 
                                                                                      "2", "3", "4", "5", "6", "7", "8", "9", "10"), partial = c("0", 
                                                                                                                                                 "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))))
})


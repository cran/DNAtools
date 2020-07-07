context("dbCompare: Regression tests comparing to old version of package")

# OLD_db_comp_dbExample_5
#source("include-test_compare.R")
# Loads/source()s helper-compare.R

data(dbExample)

################################################################################

test_that("dbCompare: dbExample, hit = 5, threads = 1", {
  db_comp_dbExample_5_t1 <- dbCompare(dbExample, hit = 5, trace = FALSE, threads = 1)
  
  expect_equal(OLD_db_comp_dbExample_5, db_comp_dbExample_5_t1)
})

test_that("dbCompare: dbExample, hit = 5, threads = 3", {
  db_comp_dbExample_5_t3 <- dbCompare(dbExample, hit = 5, trace = FALSE, threads = 3)
  
  expect_equal(OLD_db_comp_dbExample_5, db_comp_dbExample_5_t3)
})


test_that("dbVariance", {
  tmp_res <- dbVariance(probs = freqs[1:3], theta = 0, n = 100)
  #dput(apply(apply(tmp_res, 1, abs), 1, sum))
  expect_equal(apply(apply(tmp_res, 1, abs), 1, sum), 
               c(`0/0` = 12025.4906980002, `0/1` = 5811.84997479591, `0/2` = 9364.23536992672, 
                 `0/3` = 4135.34056563665, `1/0` = 594.248274779467, `1/1` = 2305.17380020229, 
                 `1/2` = 1205.22554346765, `2/0` = 119.443531637947, `2/1` = 110.427603014486, 
                 `3/0` = 3.2459841091192), 
               tolerance = 1e-7)
})

test_that("dbExpect", {
  tmp_res <- dbExpect(probs = freqs[1:3], theta = 0)
  #dput(tmp_res)
  expect_equal(tmp_res,
               structure(c(0.213951897063799, 0.0299140991222018, 0.00139207464114651, 
                           2.15593082248657e-05, 0.401463915258025, 0.0373935880117773, 
                           0.000869392353447335, NA, 0.251015506206786, 0.0116813630967216, 
                           NA, NA, 0.0522966049378706, NA, NA, NA), .Dim = c(4L, 4L), .Dimnames = list(
                             match = c("0", "1", "2", "3"), partial = c("0", "1", "2", 
                                                                        "3"))),
               tolerance = 1e-7)
})

test_that("dbExpect", {
  loci <- 3
  tmp_res <- dbCompare(dbExample[, 1:(2*loci + 1)], hit = 5, trace = FALSE, threads = 1)
  tmp_opt <- optim.relatedness(tmp_res$m, probs = freqs[1:loci],
                               theta.step = 1e-2,
                               solnp.ctrl = list(tol = 10^(-6), 
                                                 rho = 1, 
                                                 delta = 1e-7, trace = FALSE))
  
  #dput(tmp_opt)
  expect_equal(tmp_opt$value,
               structure(list(theta = c(0, 0.01, 0.02, 0.03), value = c(24334.542747719, 
                                                                        5688.37755223684, 1561.82240248446, 426.733772127926)), row.names = c(NA, 
                                                                                                                                              4L), class = "data.frame"), 
               tolerance = 1e-7)
})


test_that("dbCompare: threaded vs non-threaded: small data", {
  for (h in c(5, 7)) {
    a1 <- dbCompare(dbExample, hit = h, trace = FALSE, threads = 1)
    a2 <- dbCompare(dbExample, hit = h, trace = FALSE, threads = 2)
    a3 <- dbCompare(dbExample, hit = h, trace = FALSE, threads = 3)
    a4 <- dbCompare(dbExample, hit = h, trace = FALSE, threads = 4)
    
    expect_equal(a1$m, a2$m, info = paste0("hit = ", h))
    expect_equal(a1$m, a3$m, info = paste0("hit = ", h))
    expect_equal(a1$m, a4$m, info = paste0("hit = ", h))
  }
})

test_that("dbCompare: threaded vs non-threaded: medium data", {
  db_medium <- rbind(dbExample, dbExample)
  
  for (h in c(5, 7)) {
    a1 <- dbCompare(db_medium, hit = h, trace = FALSE, threads = 1)
    a3 <- dbCompare(db_medium, hit = h, trace = FALSE, threads = 3)
    
    # Ordering
    a1$hits <- a1$hits[order(a1$hits$id1, a1$hits$id2), ]
    a3$hits <- a3$hits[order(a3$hits$id1, a3$hits$id2), ]
    rownames(a1$hits) <- NULL
    rownames(a3$hits) <- NULL
    
    a1$m
    a3$m
    a1$m - a3$m

    expect_equal(a1, a3, info = paste0("hit = ", h))
  }
})

if (FALSE) {
  db_big <- rbind(dbExample, dbExample, dbExample, dbExample)
  db_big <- rbind(db_big, db_big)
  nrow(db_big)
  
  # rbenchmark::benchmark(
  #   threads_4 = dbCompare(db_big, hit = 5, trace = FALSE, threads = 4),
  #   single = dbCompare(db_big, hit = 5, trace = FALSE, threads = 1),
  #   replications = 2
  # )
  # 
  # microbenchmark::microbenchmark(
  #   threads_4 = dbCompare(db_big, hit = 5, trace = FALSE, threads = 4),
  #   single = dbCompare(db_big, hit = 5, trace = FALSE, threads = 1),
  #   times = 2
  # )
}

################################################################################

if (FALSE) {
  test_that("dbCompare: dbExample, hit = 5, threads = 2", {
    db_comp_dbExample_5_t2 <- dbCompare(dbExample, hit = 5, trace = FALSE, threads = 2)
    OLD_db_comp_dbExample_5_sorted <- OLD_db_comp_dbExample_5
  
    # Sort:
    db_comp_dbExample_5_t2$hits <- db_comp_dbExample_5_t2$hits[with(db_comp_dbExample_5_t2$hits, order(id1, id2, match, partial)), ]
    OLD_db_comp_dbExample_5_sorted$hits <- OLD_db_comp_dbExample_5_sorted$hits[with(OLD_db_comp_dbExample_5_sorted$hits, order(id1, id2, match, partial)), ]
    # Remove row names
    rownames(db_comp_dbExample_5_t2$hits) <- NULL
    rownames(OLD_db_comp_dbExample_5_sorted$hits) <- NULL
  
    expect_equal(OLD_db_comp_dbExample_5_sorted, db_comp_dbExample_5_t2)
  })
    
  
  test_that("dbCompare: !!!RcppParallel!!! dbExample, hit = 5, threads = 2", {
    db_comp_dbExample_5_t2 <- dbCompare(dbExample, hit = 5, trace = FALSE, threads = 666)
    OLD_db_comp_dbExample_5_sorted <- OLD_db_comp_dbExample_5
    
    # Sort:
    db_comp_dbExample_5_t2$hits <- db_comp_dbExample_5_t2$hits[with(db_comp_dbExample_5_t2$hits, order(id1, id2, match, partial)), ]
    OLD_db_comp_dbExample_5_sorted$hits <- OLD_db_comp_dbExample_5_sorted$hits[with(OLD_db_comp_dbExample_5_sorted$hits, order(id1, id2, match, partial)), ]
    # Remove row names
    rownames(db_comp_dbExample_5_t2$hits) <- NULL
    rownames(OLD_db_comp_dbExample_5_sorted$hits) <- NULL
    
    expect_equal(OLD_db_comp_dbExample_5_sorted, db_comp_dbExample_5_t2)
    #expect_equal(1, 2)
  })
  
  
}



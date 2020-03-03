context("Simulate")

test_that("genTypeRec()", {
  tmp_res <- genTypeRec(x = c(1), t = 0, n = 10)
  expect_equal(1, unique(unname(unlist(tmp_res))))
  
  tmp_res <- genTypeRec(x = c(1), t = 0.1, n = 10)
  expect_equal(1, unique(unname(unlist(tmp_res))))
})

test_that("genRypeRec()", {
  tmp_res <- genRypeRec(x = c(1), t = 0, k = c(0.25, 0.25, 0.5), n = 10)
  expect_equal(1, unique(unname(unlist(tmp_res))))
  
  tmp_res <- genRypeRec(x = c(1), t = 0.1, k = c(0.25, 0.25, 0.5), n = 10)
  expect_equal(1, unique(unname(unlist(tmp_res))))
})

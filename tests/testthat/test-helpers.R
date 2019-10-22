test_that("enforce_symm works", {
  mat_asym <- matrix(1:4, 2, 2)
  
  expect_equal(enforce_symm(mat_asym, method = "lower"),
               matrix(c(1, 2, 2, 4), 2, 2))
  expect_equal(enforce_symm(mat_asym, method = "upper"),
               matrix(c(1, 3, 3, 4), 2, 2))
  expect_equal(enforce_symm(mat_asym, method = "average"),
               matrix(c(1, 2.5, 2.5, 4), 2, 2))
})

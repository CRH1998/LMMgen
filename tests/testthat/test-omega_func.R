test_that("omega_func() calculates the variance/covariance matrix", {
  expect_equal(dim(omega_func(semi_def_matrix = list(matrix(1:9,3,3)), sigma2_vec = c(1))), dim(matrix(1:9,3,3)))
})

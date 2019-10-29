test_that("gcoda_sub works", {
  set.seed(1)
  n <- 100000
  pi <- 0.5
  mus <- c(0, 1)
  sigmas <- c(0.1, 0.1)
  inds <- rbinom(n = n, size = 1, prob = 0.5)
  samples <- 
    rnorm(n = n, mean = mus[1], sigmas[1]) * inds +
    rnorm(n = n, mean = mus[2], sigmas[2]) * (1 - inds)
  samples <- exp(matrix(samples, ncol = 1))
  
  df_marginal_mix <- get_marginal_mix(data = samples,
                                      mean0 = mus[1],
                                      var0 = sigmas[1]^2)
  
  testthat::expect_true(abs(df_marginal_mix$pi - pi) < 0.01)
  testthat::expect_true(abs(df_marginal_mix$mu0 - mus[1]) < 0.01)
  testthat::expect_true(abs(df_marginal_mix$sigma0 - sigmas[1]) < 0.01)
  testthat::expect_true(abs(df_marginal_mix$sigma_overall - 
                              sd(log(samples))) < 0.01)
})

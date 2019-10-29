# test_that("gcoda_sub works", {
#   mat_test <- enforce_symm(matrix(1:100, nrow = 100, ncol = 100))
#   
#   for(lambda in seq(10, 90, by = 10)) {
#     fit_gcoda <- gcoda:::gcoda_sub(A = mat_test, lambda = 50)
#     fit_CMTmicrobiome <- gcoda_sub(S = mat_test, lambda = 50)
#     
#     expect_equal(fit_gcoda$iSig, fit_CMTmicrobiome$Omega)
#     expect_equal(fit_gcoda$nloglik, fit_CMTmicrobiome$negLogLik)
#   }
# })

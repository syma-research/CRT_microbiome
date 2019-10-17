pick_rho_pencopula <- function(data, rholist, K) {
  folds <- sample.int(n = K, size = nrow(data), replace = TRUE)
  
  doParallel::registerDoParallel(cores = K)
  logLik <- foreach::`%dopar%`(
    foreach::foreach(k = seq_len(K),
                     .combine = "+"),
    {
      data_train <- data[folds != k, ]
      s_train <- get_s(data = data_train)
      
      l_fit_glasso <- lapply(rholist, 
                             function(rho) {
                               glasso::glasso(s = s_train,
                                              rho = rho)
                             })
      
      data_test <- data[folds == k, ]
      s_test <- get_s(data = data_test)
      vapply(seq_along(rholist),
             function(i) {
               wi <- l_fit_glasso[[i]]$wi
               log(det(wi)) - sum(diag(wi %*% s_test))
             },
             0.0) * nrow(s_test)
    })
  doParallel::stopImplicitCluster()
  
  data.frame(rho = rholist, logLik = logLik)
}

iRho <- function(rho) sinpi(rho/6) * 2

get_s <- function(data, method = "spearman", random = TRUE, R = 50) {
  s_s <- cor2(x = data, method = method, random = random, R = R)
  iRho(s_s)
}

lh <- function(s, w) {}

gcoda <- function(data, lambda_list, 
                  gamma_EBIC = 0.5,
                  K_CV = 5,
                  ncores = 6) {
  if(any(data <= 0))
    stop("Data cannot have zeros!")
  
  df_eval <- data.frame(lambda = lambda_list,
                        df = NA_real_,
                        AIC = NA_real_,
                        BIC = NA_real_,
                        EBIC = NA_real_)
  # CLT for compositional data
  data <- log(data) - rowMeans(log(data))
  s_data <- var(data)
  
  doParallel::registerDoParallel(cores = ncores)
  l_fits_evals <- foreach::`%dopar%`(
    foreach::foreach(lambda = lambda_list),
    {
      # With parallelization cannot utilize bigger lambda fits on the path
      # to speed up convergence
      # FIXME?
      fit_gcoda <- gcoda_sub(S = s_data, 
                             lambda = lambda,
                             threshold = 1e-6)
      
      df <- (sum(fit_gcoda != 0) - ncol(data)) / 2
      negLogLik <- negLogLik_gcoda(S = s_data, Omega = fit_gcoda)
      AIC <- negLogLik * nrow(s_data) + 2 * df 
      BIC <- negLogLik * nrow(s_data) + log(nrow(data)) * df 
      EBIC <- BIC + 4 * gamma_EBIC * log(ncol(data)) * df
      return(list(fit = fit_gcoda,
                  evals = c(df, AIC, BIC, EBIC)))
    })
  doParallel::stopImplicitCluster()
  
  l_fits <- lapply(l_fits_evals, function(x) x$fit)
  df_eval[, c("df", "AIC", "BIC", "EBIC")] <- 
    t(vapply(l_fits_evals,
             function(x) x$evals,
             rep(0.0, 4)))
  
  if(!is.null(K_CV)) {
    folds <- sample.int(n = K_CV, size = nrow(data), replace = TRUE)
    doParallel::registerDoParallel(cores = ncores)
    ## FIXME
    # parallelization not optimized
    # ideally should parallel over combinations of k and lambda
    # instead of just over k
    negLogLik_CV <- foreach::`%dopar%`(
      foreach::foreach(k = seq_len(K_CV),
                       .combine = "+"),
      {
        data_train <- data[folds != k, ]
        s_train <- var(data_train)
        
        l_fit_gcoda <- lapply(lambda_list, 
                              function(lambda)
                                gcoda_sub(S = s_train, 
                                          lambda = lambda,
                                          threshold = 1e-6)
        )
        
        data_test <- data[folds == k, ]
        s_test <- var(data_test)
        return(vapply(seq_along(lambda_list),
                      function(i)
                        negLogLik_gcoda(S = s_test, Omega = l_fit_gcoda[[i]]) * 
                        nrow(s_test),
                      0.0))
      })
    doParallel::stopImplicitCluster()
    
    df_eval$negLogLik_CV <- negLogLik_CV / nrow(data)
  }
  
  return(list(fits = l_fits, df_eval = df_eval))
}

gcoda_sub <- function(S, Omega = NULL, 
                      lambda = 0.1, 
                      tol_err = 1e-4, 
                      k_max = 100,
                      symm = TRUE,
                      threshold = NULL) {
  p <- ncol(S)
  if(is.null(Omega)) {
    Omega <- diag(p)
  }
  
  err <- 1
  k <- 0
  fval_cur <- Inf
  
  while(err > tol_err && k < k_max) {
    
    Stilde <- calc_Stilde(S = S, Omega = Omega)
    Omega_new <- glasso_wrapper(Stilde, lambda, 
                                symm = symm, threshold = threshold)
    
    fval_new <- obj_gcoda(Omega = Omega_new, S = S, lambda = lambda)
    xerr <- max(abs(Omega_new - Omega) / (abs(Omega_new) + 1))
    err <- min(xerr, abs(fval_cur - fval_new)/(abs(fval_new) + 1))
    
    k <- k + 1
    Omega <- Omega_new
    fval_cur <- fval_new
  }
  
  if(k >= k_max)
    warning("Maximum Iteration:", k_max, 
            "&& Relative error:", err, "!\n")
  
  return(Omega)
}

calc_Stilde <- function(S, Omega) {
  p <- ncol(S)
  # Omega %*% 1
  rowsum_Omega <- rowSums(Omega)
  sum_Omega <- sum(rowsum_Omega)
  rowsum_Omega_norm <- rowsum_Omega / sum_Omega
  # One column of S %*% Omega %*% 11^T / (1^T %*% Omega %*% 1)
  oneCol_S_Omega_norm <- rowSums(S * rep(rowsum_Omega_norm, each = p))  
  return(S - oneCol_S_Omega_norm - rep(oneCol_S_Omega_norm, each = p) + 
    sum(rowsum_Omega_norm * oneCol_S_Omega_norm) + 1 / sum_Omega)
}

obj_gcoda <- function(S, Omega, lambda) {
  pen <- lambda * sum(abs(Omega))
  return(negLogLik_gcoda(S, Omega) + pen)
}

negLogLik_gcoda <- function(S, Omega) {
  p <- ncol(S)
  rowsum_Omega <- rowSums(Omega)
  sum_Omega <- sum(rowsum_Omega)
  negLogLik <- -log(det(Omega)) + sum(Omega * S) + log(sum_Omega) - 
    sum(rowsum_Omega * rowSums(S * rep(rowsum_Omega, each = p))) / sum_Omega
  return(negLogLik)
}

get_mean_logisticMVN <- function(data) {
  if(any(data <= 0))
    stop("Data cannot have zeros!")
  
  mu <- apply(log(mat_X_p_pseudocount), 1, mean)
  return(mu)
}
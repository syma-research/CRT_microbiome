choose_lambda <- function(data, lambda_list, 
                          methods = c("CV", "EBIC"),
                          K_CV = 5,
                          gamma_EBIC = 0.5,
                          ncores = 6) {
  
  df_return <- data.frame(lambda = lambda_list,
                          negLogLik_CV = NA_real_,
                          df = NA_real_,
                          AIC = NA_real_,
                          BIC = NA_real_,
                          EBIC = NA_real_)
  
  ind_methods <- c("CV" = FALSE, "EBIC" = FALSE)
  ind_methods[] <- names(ind_methods) %in% methods
  if(!any(ind_methods))
    stop("At least one evaluation metric needs to be specified!")
  
  if(ind_methods["CV"]) {
    folds <- sample.int(n = K_CV, size = nrow(data), replace = TRUE)
    doParallel::registerDoParallel(cores = ncores)
    ## FIXME
    # parallelization not optimized
    # ideally should parallel over combinations of k and lambda
    # instead of just over k
    logLik <- foreach::`%dopar%`(
      foreach::foreach(k = seq_len(K_CV),
                       .combine = "+"),
      {
        data_train <- data[folds != k, ]
        s_train <- get_s(data = data_train)
        
        l_fit_glasso <- lapply(lambda_list, 
                               function(lambda) {
                                 glasso_wrapper(S = s_train, 
                                                lambda = lambda)
                               })
        
        data_test <- data[folds == k, ]
        s_test <- get_s(data = data_test)
        vapply(seq_along(lambda_list),
               function(i) {
                 wi <- l_fit_glasso[[i]]
                 logLik_S(S = s_test, n = nrow(data_test), Theta = wi)
               },
               0.0)
      })
    doParallel::stopImplicitCluster()
    
    df_return$negLogLik_CV <- -logLik
  }
  
  if(ind_methods["EBIC"]) {
    s_data <- get_s(data = data)
    doParallel::registerDoParallel(cores = ncores)
    ICs <- foreach::`%dopar%`(
      foreach::foreach(lambda = lambda_list,
                       .combine = rbind),
      {
        fit_glasso <- glasso_wrapper(S = s_data, 
                                     lambda = lambda)
        df <- (sum(0 + (abs(fit_glasso) > 1e-6)) - ncol(data)) / 2
        logLik <- logLik_S(S = s_data, n = nrow(data), Theta = fit_glasso)
        AIC <- -logLik * 2 + 2 * df 
        BIC <- -logLik * 2 + log(nrow(data)) * df 
        EBIC <- BIC + 4 * gamma_EBIC * log(ncol(data)) * df
        c(df, AIC, BIC, EBIC)
      })
    doParallel::stopImplicitCluster()
    
    df_return[, c("df", "AIC", "BIC", "EBIC")] <- ICs
  }
  
  return(df_return)
}

iRho <- function(rho) sinpi(rho/6) * 2

get_s <- function(data, method = "spearman", 
                  random = TRUE, sim = FALSE, 
                  R = 1000) {
  s_s <- cor2(x = data, method = method, random = random,
              sim = sim, R = R)
  iRho(s_s)
}

glasso_wrapper <- function(S, lambda, source = "glasso",
                           simplify = FALSE,
                           symm = TRUE) {
  icov <- diag(1/(diag(S) + lambda))
  
  if(simplify) 
    z <- which(rowSums(abs(S) > lambda) > 1)
  else 
    z <- seq_len(nrow(S))
  q <- length(z)
  if (q > 0) {
    if(source == "huge") {
      out.glasso <- huge::huge.glasso(x = S[z, z, drop = FALSE],
                                      lambda = lambda,
                                      verbose = FALSE)
      icov[z, z] <- out.glasso$icov[[1]]
    }
    if(source == "hugec") {
      out.glasso <- .Call("_huge_hugeglasso",
                          S[z, z, drop = FALSE],
                          lambda,
                          FALSE,
                          FALSE,
                          FALSE,
                          PACKAGE = "huge")
      icov[z, z] <- out.glasso$icov[[1]]
    }
    if(source == "glasso") {
      out.glasso <- glasso::glasso(s = S[z, z, drop = FALSE],
                                   rho = lambda)
      icov[z, z] <- out.glasso$wi
    }
  }
  
  if(symm)
    icov <- enforce_symm(icov)
  return(icov)
}

sequential_gcoda <- function(sample0, ind_feature, 
                             params,
                             sd_proposal = 1,
                             m = 100,
                             debug_dir = NULL) {
  if(!is.null(debug_dir)) {
    dir.create(debug_dir, recursive = TRUE)
    debug_dir <- normalizePath(debug_dir)
  }
  
  samples <- 
    list(list(sample = sample0,
              logLik = 
                log_dgcoda(x = sample0,
                           mu = params$mu,
                           Omega = params$Omega),
              success = TRUE
    ))
  
  if(!is.null(debug_dir)) {
    save(samples,
         file = paste0(debug_dir, 
                       "/samples.RData"))
  }
  
  for(i_step in seq_len(m)) {
    samples[[i_step + 1]] <- 
      one_step_gcoda(current = samples[[i_step]],
                     ind_feature = ind_feature,
                     params = params,
                     sd_proposal = sd_proposal,
                     debug_dir = debug_dir)
    
    if(!is.null(debug_dir))
      save(samples,
           file = paste0(debug_dir, 
                         "/samples.RData"))
  }
  
  return(
    Reduce("rbind", 
           lapply(
             samples[-1],
             function(x) x$sample))
  )
}

log_dgcoda <- function(x, mu, Omega) {
  A <- (log(x) - mu) %*% t((log(x) - mu))
  p <- ncol(A)
  Omega_O <- rowSums(Omega)
  S_Omega <- sum(Omega_O)
  nloglik <- sum(Omega_O * A) - 
    sum(Omega_O * rowSums(A * rep(Omega_O, each = p))) / S_Omega +
    sum(log(x))
  return(-nloglik)
}

one_step_gcoda <- 
  function(current,
           ind_feature,
           params,
           sd_proposal,
           debug_dir) {
    x_feature_star <- 
      SparseDOSSA2:::expit(
        rnorm(1,
              mean = SparseDOSSA2:::logit(current$sample[ind_feature]),
              sd = sd_proposal))

    if(!is.null(debug_dir))
      save(x_feature_star, file = paste0(debug_dir, "/proposal.RData"))
    
    x_star <- current$sample
    x_star[ind_feature] <- x_feature_star
    x_star[!ind_feature] <- 
      (1 - x_feature_star) * 
      x_star[!ind_feature] / 
      sum(x_star[!ind_feature])
    
    logLik_star <- 
      log_dgcoda(x = x_star, 
                 mu = params$mu,
                 Omega = params$Omega)
    
    log_paccept <- logLik_star - current$logLik -
      log(current$sample[ind_feature]) - log(1 - current$sample[ind_feature]) +
      log(x_feature_star) + log(1 - x_feature_star)
    
    if(log(runif(1)) < log_paccept)
      return(list(sample = x_star,
                  logLik = logLik_star,
                  success = TRUE))
    else
      return(list(sample = current$sample,
                  logLik = current$logLik,
                  success = FALSE))
  }

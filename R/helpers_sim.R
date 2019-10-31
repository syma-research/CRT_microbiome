enforce_ZI <- function(data, pi) {
  if(ncol(data) != length(pi))
    stop("Number of features in data and pi must agree!")
  
  data_mod <- 
    vapply(seq_len(ncol(data)),
           function(i) 
             enforce_ZI_one(data[, i, drop = TRUE],
                            pi[i]),
           data[, 1, drop = TRUE])
  
  data_mod <- t(apply(data_mod, 1, function(x) x / sum(x)))
  dimnames(data_mod) <- dimnames(data)
  
  return(data_mod)
}

enforce_ZI_one <- function(x, pi) {
  if(any(x <= 0))
    stop("x shouldn't have any zeros!")
  x[rank(x) / length(x) <= 1- pi] <- 0
  return(x)
}

rcopulasso <- function(n, mean, sd, pi, sigma, norm = TRUE) {
  if(length(mean) != length(sd) |
     length(mean) != length(pi) |
     length(mean) != nrow(sigma))
    stop("Parameter dimensions must agree!")
  
  # sample marginals
  mat_marginals <- 
    vapply(seq_len(length(mean)),
           function(i)
             rZILogitN_one(n = n,
                           mean = mean[i],
                           sd = sd[i],
                           pi = pi[i]),
           rep(0.0, n))
  # arrange from smallest to largest for shuffling
  mat_marginals <- 
    apply(mat_marginals, 2, function(x) x[order(x)])
  
  # sample ranks
  mat_rank <- 
    mvtnorm::rmvnorm(n = n, sigma = sigma)
  mat_rank <- apply(mat_rank, 2, rank)
  
  mat_samples <- 
    vapply(seq_len(length(mean)),
           function(i)
             mat_marginals[, i, drop = TRUE][mat_rank[, i, drop = TRUE]],
           rep(0.0, n))
  
  if(norm)
    mat_samples <- apply(mat_samples, 1, function(x) x / sum(x))
  
  return(mat_samples)
}

rZILogitN_one <- function(n, mean, sd, pi) {
  return(expit(rnorm(n = n, mean = mean, sd = sd)) *
           rbinom(n = n, size = 1, prob = pi))
}

rgcoda <- function(n, mean0, mean1, sd0, sd1, pi, sigma) {
  if(length(mean0) != length(mean1) |
     length(mean0) != length(sd0) |
     length(mean0) != length(sd1) |
     length(mean0) != length(pi) |
     length(mean0) != nrow(sigma)
  )
    stop("Parameter dimensions must agree!")
  
  # sample marginals
  l_marginals <- 
    lapply(seq_len(length(mean0)),
           function(i)
             rMixLogN_one(n = n,
                          mean0 = mean0[i], mean1 = mean1[i],
                          sd0 = sd0[i], sd1 = sd1[i],
                          pi = pi[i]))
  mat_marginals <- vapply(l_marginals, 
                          function(x) x$abd,
                          rep(0.0, n))
  mat_ind <- vapply(l_marginals, 
                    function(x) x$ind,
                    rep(1, n))
  mat_order <- apply(mat_marginals, 2, order)
  
  # sample ranks
  mat_rank <- 
    mvtnorm::rmvnorm(n = n, sigma = sigma)
  mat_rank <- apply(mat_rank, 2, rank)
  
  mat_samples <- 
    vapply(seq_len(length(mean0)),
           function(i)
             mat_marginals[, i, drop = TRUE][mat_order[, i, drop = TRUE]][
               mat_rank[, i, drop = TRUE]],
           rep(0.0, n))
  mat_ind <- 
    vapply(seq_len(length(mean0)),
           function(i)
             mat_ind[, i, drop = TRUE][mat_order[, i, drop = TRUE]][
               mat_rank[, i, drop = TRUE]],
           rep(0.0, n))
  
  
  mat_samples <- apply(mat_samples, 1, function(x) x / sum(x))
  
  return(list(abd = mat_samples,
              ind = t(mat_ind)))
}

rMixLogN_one <- function(n, mean0, mean1, sd0, sd1, pi) {
  ind <- rbinom(n = n, size = 1, prob = pi)
  return(list(abd = exp(rnorm(n = n, mean = mean0, sd = sd0) * (1 - ind) +
                          rnorm(n = n, mean = mean1, sd = sd1) * ind),
              ind = ind))
}
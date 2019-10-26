mcmc_gcoda <- function(p, index,
                       params, R = 10000) {
  if(length(index) != length(p) |
     !is.logical(index) |
     sum(index) != 1)
    stop("Index does not follow format!")
  if(any(p <= 0))
    stop("p must be positive!")
  
  Omega <- params$Omega
  sigma <- sqrt(solve(Omega)[index, index, drop = TRUE])
  mu <- params$mu
  rowsum_Omega <- rowSums(Omega)
  sum_Omega <- sum(rowsum_Omega)
  A <- Omega - rowsum_Omega %*% t(rowsum_Omega) / sum_Omega
  
  logp <- log(p)
  l_samples <- list()
  l_samples[[1]] <- list(logp = logp,
                         index = index,
                         sigma = sigma,
                         mu = mu,
                         A = A,
                         logLik = log_dmvnorm(logp - mu, A) - 
                           log(1 - exp(logp[index])))
  
  for(r in (seq_len(R) + 1)) {
    l_samples[[r]] <- 
      one_step_gcoda(logp = l_samples[[r - 1]]$logp,
                     index = index,
                     sigma = sigma, 
                     mu = mu, 
                     A = A,
                     logLik = l_samples[[r - 1]]$logLik)
  }
  
  return(l_samples[-1])
}

one_step_gcoda <- function(logp, index, sigma, mu, A, logLik) {
  logp_star <- logp
  while(TRUE) {
    logp_star[index] <- logp[index] + rnorm(n = 1, sd = sigma)
    if(logp_star[index] < 0) break
  }
  logp_star[!index] <- 
    logp[!index] + 
    log(1 - exp(logp_star[index])) - 
    log(1 - exp(logp[index]))
  logLik_star <- log_dmvnorm(logp_star - mu, A) - log(1 - exp(logp_star[index]))
  prop <- exp(logLik_star - logLik) * 
    pnorm(-logp[index], sd = sigma) / 
    pnorm(-logp_star[index], sd = sigma)
    
  if(runif(1) < prop)
    return(list(accept = TRUE,
                prop = prop,
                logp = logp_star,
                logLik = logLik_star))
  else
    return(list(accept = FALSE,
                prop = prop,
                logp = logp,
                logLik = logLik))
}

mcmc_copulasso <- function(p, index,
                           params,
                           R = 10000) {
  if(length(index) != length(p) |
     !is.logical(index) |
     sum(index) != 1)
    stop("Index does not follow format!")
  A <- params$Omega - diag(length(p))
  u <- p_to_u(p, mu = params$mu, sigma = params$sigma, pi = params$pi)
  g <- qnorm(u)
  
  l_samples <- list()
  l_samples[[1]] <- list(g = g,
                         u = u,
                         p = p,
                         loglik_g = log_dmvnorm(g, A),
                         loglik_u = sum(log(u_to_d(u, params$pi))))
  for(r in (seq_len(R) + 1)) {
    l_samples[[r]] <- 
      one_step_copulasso(g = l_samples[[r - 1]]$g,
                         u = l_samples[[r - 1]]$u,
                         p = l_samples[[r - 1]]$p,
                         index = index,
                         mu = params$mu, sigma = params$sigma, pi = params$pi, 
                         A = A,
                         loglik_g = l_samples[[r - 1]]$loglik_g, 
                         loglik_u = l_samples[[r - 1]]$loglik_u)
  }
  
  return(l_samples[-1])
}

mcmc2_copulasso <- function(p, index,
                            params,
                            steps = 10,
                            R = 1000,
                            ncores = 6) {
  if(length(index) != length(p) |
     !is.logical(index) |
     sum(index) != 1)
    stop("Index does not follow format!")
  A <- params$Omega - diag(length(p))
  u <- p_to_u(p, mu = params$mu, sigma = params$sigma, pi = params$pi)
  g <- qnorm(u)
  
  ll_paths <- list()
  
  # first, walk steps to generate "center" of MCMC sampling
  ll_paths[[1]] <- vector("list", steps + 1)
  ll_paths[[1]][[steps + 1]] <- list(g = g,
                                     u = u,
                                     p = p,
                                     loglik_g = log_dmvnorm(g, A),
                                     loglik_u = sum(log(u_to_d(u, params$pi))))
  for(i_step in seq_len(steps)) {
    ll_paths[[1]][[steps + 1 - i_step]] <- 
      one_step_copulasso(
        g = ll_paths[[1]][[steps + 2 - i_step]]$g,
        u = ll_paths[[1]][[steps + 2 - i_step]]$u,
        p = ll_paths[[1]][[steps + 2 - i_step]]$p,
        index = index,
        mu = params$mu, sigma = params$sigma, pi = params$pi, A = A,
        loglik_g = ll_paths[[1]][[steps + 2 - i_step]]$loglik_g, 
        loglik_u = ll_paths[[1]][[steps + 2 - i_step]]$loglik_u)
  }
  sample_center <- ll_paths[[1]][[1]]
  ll_paths[[1]] <- ll_paths[[1]][-1]
  
  # Now start from the center, simultaneously sample new samples
  doParallel::registerDoParallel(cores = ncores)
  ll_paths_new <- foreach::`%dopar`(
    foreach::foreach(r = seq_len(R)),
    {
      l_samples <- list()
      l_samples[[1]] <- sample_center
      for(i_step in seq_len(steps)) {
        l_samples[[i_step + 1]] <- 
          one_step_copulasso(
            g = l_samples[[i_step]][[steps + 2 - i_step]]$g,
            u = l_samples[[i_step]][[steps + 2 - i_step]]$u,
            p = l_samples[[i_step]][[steps + 2 - i_step]]$p,
            index = index,
            mu = params$mu, sigma = params$sigma, pi = params$pi, A = A,
            loglik_g = l_samples[[i_step]][[steps + 2 - i_step]]$loglik_g, 
            loglik_u = l_samples[[i_step]][[steps + 2 - i_step]]$loglik_u)
      }
      l_samples <- l_samples[-1]
      return(l_samples)
    })
  doParallel::stopImplicitCluster()
  
  ll_paths <- c(ll_paths, ll_paths_new)
  return(ll_paths[-1])
}

one_step_copulasso <- function(g, u, p, index,
                               mu, sigma, pi, A,
                               loglik_g, loglik_u) {
  
  g_star <- g
  u_star <- u
  p_star <- p
  g_star[index] <- g[index] + rnorm(n = 1)
  u_star[index] <- pnorm(g_star[index])
  p_star[index] <- u_to_p(u_star[index], 
                          mu = mu[index], 
                          sigma = sigma[index], 
                          pi = pi[index])
  p_star[!index] <- p[!index] * (1 - p_star[index]) / (1 - p[index])
  u_star[!index] <- p_to_u(p_star[!index], 
                           mu = mu[!index], 
                           sigma = sigma[!index], 
                           pi = pi[!index])
  g_star[!index] <- qnorm(u_star[!index])
  loglik_g_star <- log_dmvnorm(g_star, A)
  loglik_u_star <- sum(log(u_to_d(u_star, pi)))
  prop <- exp(loglik_g_star - loglik_g + loglik_u_star - loglik_u)
  if(runif(1) < prop)
    return(list(accept = TRUE,
                prop = prop,
                g = g_star,
                u = u_star,
                p = p_star,
                loglik_g = loglik_g_star,
                loglik_u = loglik_u_star))
  else
    return(list(accept = FALSE,
                prop = prop,
                g = g,
                u = u,
                p = p,
                loglik_g = loglik_g,
                loglik_u = loglik_u))
}

p_to_u <- function(p, mu, sigma, pi) {
  to_return <- pi / 2
  to_return[p > 0] <- 
    pi[p > 0] + 
    pnorm(logit(p[p > 0]), 
          mean = mu[p > 0], 
          sd = sigma[p > 0]) * 
    (1 - pi[p > 0])
  
  return(to_return)
}

u_to_p <- function(u, mu, sigma, pi) {
  to_return <- rep(0, length = length(u))
  to_return[u > pi] <- 
    expit(qnorm((u[u > pi] - pi[u > pi]) / 
                  (1 - pi[u > pi]), 
                mean = mu[u >= pi], 
                sd = sigma[u > pi]))
  
  return(to_return)
}

u_to_d <- function(u, pi) {
  to_return <- pi
  to_return[u > pi] <- 
    dnorm(qnorm(u[u > pi] - pi[u > pi])) * 
    (1 - pi[u > pi])
  
  return(to_return)
}

log_dmvnorm <- function(g, A) {
  return((-t(g) %*% A %*% g / 2)[1, 1, drop = TRUE])
}

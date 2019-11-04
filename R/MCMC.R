# mcmc_gcoda <- function(x, index,
#                        params, R = 10000) {
#   if(length(index) != length(x) |
#      !is.logical(index) |
#      sum(index) != 1)
#     stop("Index does not follow format!")
#   if(any(x <= 0))
#     stop("x must be positive!")
#   
#   Omega <- params$Omega
#   sigma <- sqrt(solve(Omega)[index, index, drop = TRUE])
#   mu <- params$mu
#   rowsum_Omega <- rowSums(Omega)
#   sum_Omega <- sum(rowsum_Omega)
#   A <- Omega - rowsum_Omega %*% t(rowsum_Omega) / sum_Omega
#   
#   logx <- log(x)
#   l_samples <- list()
#   l_samples[[1]] <- list(logx = logx,
#                          index = index,
#                          sigma = sigma,
#                          mu = mu,
#                          A = A,
#                          logLik = log_dmvnorm(logx - mu, A) - 
#                            log(1 - exp(logx[index])))
#   
#   for(r in (seq_len(R) + 1)) {
#     l_samples[[r]] <- 
#       one_step_gcoda(logx = l_samples[[r - 1]]$logx,
#                      index = index,
#                      sigma = sigma, 
#                      mu = mu, 
#                      A = A,
#                      logLik = l_samples[[r - 1]]$logLik)
#   }
#   
#   return(l_samples[-1])
# }
# 
# one_step_gcoda <- function(logx, index, sigma, mu, A, logLik) {
#   logp_star <- logx
#   while(TRUE) {
#     logp_star[index] <- logx[index] + rnorm(n = 1, sd = sigma)
#     if(logp_star[index] < 0) break
#   }
#   logp_star[!index] <- 
#     logx[!index] + 
#     log(1 - exp(logp_star[index])) - 
#     log(1 - exp(logx[index]))
#   logLik_star <- log_dmvnorm(logp_star - mu, A) - log(1 - exp(logp_star[index]))
#   prop <- exp(logLik_star - logLik) * 
#     pnorm(-logx[index], sd = sigma) / 
#     pnorm(-logp_star[index], sd = sigma)
#   
#   if(runif(1) < prop)
#     return(list(accept = TRUE,
#                 prop = prop,
#                 logx = logp_star,
#                 logLik = logLik_star))
#   else
#     return(list(accept = FALSE,
#                 prop = prop,
#                 logx = logx,
#                 logLik = logLik))
# }

mcmc_copulasso <- function(a, index,
                           params,
                           R = 10000) {
  if(length(index) != length(a) |
     !is.logical(index) |
     sum(index) != 1)
    stop("Index does not follow format!")
  
  index_nonzero <- a > 0 | index
  Sigma <- solve(params$Omega)
  Omega_nonzero <- solve(Sigma[index_nonzero, index_nonzero])
  mean_cond <- cond_mean_MVN(Sigma = Sigma, index = index_nonzero)
  var_cond <- cond_var_MVN(Sigma = Sigma, index = index_nonzero)
  
  u <- p_to_u(a[index_nonzero], 
              mu = params$mu[index_nonzero], 
              sigma = params$sigma[index_nonzero], 
              pi = params$pi[index_nonzero])
  g <- qnorm(u)
  
  l_samples <- list()
  l_samples[[1]] <- 
    list(a = a,
         u = u,
         g = g,
         loglik_g = 
           log_dmvnorm(g = g, Omega = Omega_nonzero) +
           log(mvtnorm::pmvnorm(upper = qnorm(p = params$pi[!index_nonzero]),
                                mean = (mean_cond %*% g)[, 1, drop = TRUE],
                                sigma = var_cond)))
  for(r in seq_len(R)) {
    l_samples[[r + 1]] <- 
      one_step_copulasso(
        a = l_samples[[r]]$a,
        u = l_samples[[r]]$u,
        g = l_samples[[r]]$g,
        index = index, index_nonzero = index_nonzero,
        mu = params$mu, sigma = params$sigma, pi = params$pi, 
        Omega_nonzero = Omega_nonzero,
        mean_cond = mean_cond, var_cond = var_cond,
        loglik_g = l_samples[[r]]$loglik_g)
  }
  
  return(l_samples)
}

cond_mean_MVN <- function(Sigma, index) {
  t(solve(Sigma[index, index], Sigma[index, !index]))
}

cond_var_MVN <- function(Sigma, index) {
  Sigma[!index, !index] - 
    Sigma[!index, index] %*% solve(Sigma[index, index], Sigma[index, !index])
}
# mcmc2_copulasso <- function(a, index,
#                             params,
#                             steps = 10,
#                             R = 1000,
#                             ncores = 6) {
#   if(length(index) != length(a) |
#      !is.logical(index) |
#      sum(index) != 1)
#     stop("Index does not follow format!")
#   
#   u <- p_to_u(p = a, mu = params$mu, sigma = params$sigma, pi = params$pi)
#   g <- qnorm(p = u)
#   
#   ll_paths <- list()
#   
#   # first, walk steps to generate "center" of MCMC sampling
#   ll_paths[[1]] <- list()
#   ll_paths[[1]][[1]] <- 
#     list(g = g,
#          u = u,
#          a = a,
#          loglik_g = log_dmvnorm(g = g, 
#                                 Omega = params$Omega))
#   for(i_step in seq_len(steps)) {
#     ll_paths[[1]][[i_step + 1]] <- 
#       one_step_copulasso(
#         a = ll_paths[[1]][[i_step]]$a,
#         u = ll_paths[[1]][[i_step]]$u,
#         g = ll_paths[[1]][[i_step]]$g,
#         index = index,
#         mu = params$mu, sigma = params$sigma, pi = params$pi, 
#         Omega = params$Omega,
#         loglik_g = ll_paths[[1]][[i_step]]$loglik_g)
#   }
#   sample_center <- ll_paths[[1]][[steps + 1]]
#   
#   # Now start from the center, simultaneously sample new samples
#   doParallel::registerDoParallel(cores = ncores)
#   ll_paths_new <- foreach::`%dopar`(
#     foreach::foreach(r = seq_len(R)),
#     {
#       l_samples <- list()
#       l_samples[[1]] <- sample_center
#       for(i_step in seq_len(steps)) {
#         l_samples[[i_step + 1]] <- 
#           one_step_copulasso(
#             a = l_samples[[i_step]]$a,
#             u = l_samples[[i_step]]$u,
#             g = l_samples[[i_step]]$g,
#             index = index,
#             mu = params$mu, sigma = params$sigma, pi = params$pi, 
#             Omega = params$Omega,
#             loglik_g = l_samples[[i_step]]$loglik_g)
#       }
#       l_samples <- l_samples[-1]
#       return(l_samples)
#     })
#   doParallel::stopImplicitCluster()
#   
#   ll_paths <- c(ll_paths, ll_paths_new)
#   return(ll_paths)
# }

one_step_copulasso <- function(a, u, g, 
                               index, index_nonzero,
                               mu, sigma, pi, 
                               Omega_nonzero, mean_cond, var_cond,
                               loglik_g) {
  a_star <- a
  u_star <- u
  g_star <- g
  g_star[index[index_nonzero]] <- g[index[index_nonzero]] + rnorm(n = 1)
  u_star[index[index_nonzero]] <- pnorm(g_star[index[index_nonzero]])
  a_star[index] <- u_to_p(u = u_star[index[index_nonzero]], 
                          mu = mu[index], 
                          sigma = sigma[index], 
                          pi = pi[index])
  a_star[!index] <- a[!index] * (1 - a_star[index]) / (1 - a[index])
  u_star[!index[index_nonzero]] <- p_to_u(p = a_star[!index & index_nonzero], 
                                          mu = mu[!index & index_nonzero], 
                                          sigma = sigma[!index & index_nonzero], 
                                          pi = pi[!index & index_nonzero])
  g_star[!index[index_nonzero]] <- qnorm(u_star[!index[index_nonzero]])
  loglik_g_star <- log_dmvnorm(g = g_star, Omega = Omega_nonzero) +
    log(mvtnorm::pmvnorm(upper = qnorm(p = pi[!index_nonzero]),
                         mean = (mean_cond %*% g_star)[, 1, drop = TRUE],
                         sigma = var_cond))
  prop <- 
    exp(loglik_g_star - loglik_g + 
          (length(g) - 2) * (log(1 - a_star[index]) - log(1 - a[index])))
  if(runif(1) < prop)
    return(list(accept = TRUE,
                prop = prop,
                a = a_star,
                u = u_star,
                g = g_star,
                loglik_g = loglik_g_star))
  else
    return(list(accept = FALSE,
                prop = prop,
                a = a,
                u = u,
                g = g,
                loglik_g = loglik_g))
}

p_to_u <- function(p, mu, sigma, pi) {
  if(sum(p == 0) > 1)
    stop("Something's wrong!")
  to_return <- runif(n = length(p)) * pi
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

log_dmvnorm <- function(g, Omega) {
  return((-t(g) %*% Omega %*% g / 2)[1, 1, drop = TRUE])
}

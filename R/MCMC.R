logit <- function(x) log(x) - log(1 - x)

expit <- function(x) exp(x) / (1 + exp(x))

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
  return((-t(g) %*% A %*% g / 2)[1, 1])
}

one_step <- function(g, u, p,
                     A, 
                     mu, sigma, pi,
                     loglik_g, loglik_u) {
  g_star <- g
  u_star <- u
  p_star <- p
  g_star[1] <- g[1] + rnorm(n = 1)
  u_star[1] <- pnorm(g_star[1])
  p_star[1] <- u_to_p(u_star[1], mu = mu[1], sigma = sigma[1], pi = pi[1])
  p_star[-1] <- p[-1] * (1 - p_star[-1]) / (1 - p[-1])
  u_star[-1] <- p_to_u(p_star[-1], mu = mu[-1], sigma = sigma[-1], pi = pi[-1])
  g_star[-1] <- qnorm(u_star[-1])
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

metro <- function(p, mu, sigma, pi, Sigma, 
                  R = 1000, K = 3) {
  A <- solve(Sigma) - diag(ncol(Sigma))
  u <- p_to_u(p, mu = mu, sigma = sigma, pi = pi)
  g <- qnorm(u)
  
  doParallel::registerDoParallel(cores = K)
  l_chains <- foreach::`%dopar%`(
    foreach::foreach(k = 1:K),
    {
      # Initialize p
      g[1] <- g[1] + rnorm(n = 1)
      u[1] <- pnorm(g[1])
      p[1] <- u_to_p(u[1], mu = mu[1], sigma = sigma[1], pi = pi[1])
      samples <- list()
      samples[[1]] <- list(g = g,
                           u = u,
                           p = p,
                           loglik_g = log_dmvnorm(g, A),
                           loglik_u = sum(log(u_to_d(u, pi))))
      for(r in 2:(R + 1)) {
        samples[[r]] <- 
          one_step(g = samples[[r - 1]]$g,
                   u = samples[[r - 1]]$u,
                   p = samples[[r - 1]]$p,
                   A = A, mu = mu, sigma = sigma, pi = pi,
                   loglik_g = samples[[r - 1]]$loglik_g, 
                   loglik_u = samples[[r - 1]]$loglik_u)
      }
      return(samples)
    })
  doParallel::stopImplicitCluster()
  
  return(l_chains)
}

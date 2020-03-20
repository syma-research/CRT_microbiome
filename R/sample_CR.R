sample_CR <- function(data, ind_feature, 
                      params_a, params_x,
                      m = 100, space_size = 10,
                      debug_dir = NULL) {
  samples_CR <- future.apply::future_lapply(
    seq_len(nrow(data)),
    function(i_sample) {
      x <- data[i_sample, ]
      d <- sample.int(n = m + 1, size = 1)
      samples_forward <- 
        samples_backward <- 
        list(list(sample = x,
                  logLik = SparseDOSSA2:::dx(x = x,
                                             pi0 = params_a$pi0,
                                             mu = params_a$mu,
                                             sigma = params_a$sigma,
                                             Omega = params_a$Omega,
                                             log.p = TRUE) -
                    dZILogit(x = x[ind_feature],
                             pi0 = params_x$pi0[ind_feature],
                             mu = params_x$mu[ind_feature],
                             sigma = params_x$sigma[ind_feature],
                             log.p = TRUE),
                  success = TRUE
                    ))
      if(d < m + 1)
        for(i_step in seq_len((m + 1 - d) * space_size)) {
          samples_forward[[i_step + 1]] <- 
            one_step_CR(current = samples_forward[[i_step]],
                        ind_feature = ind_feature,
                        params_a = params_a,
                        params_x = params_x)
        }
      if(d > 1)
          for(i_step in seq_len((d - 1) * space_size)) {
            samples_backward[[i_step + 1]] <- 
              one_step_CR(current = samples_backward[[i_step]],
                          ind_feature = ind_feature,
                          params_a = params_a,
                          params_x = params_x)
          }
      return(list(samples = c(samples_forward[-1],
                              rev(samples_backward[-1])),
                  d = d))
    }
  )
  
  if(!is.null(debug_dir))
    save(samples_CR, file = debug_dir)
  
  samples_for_testing <- seq_len(m) %>% 
    purrr::map(function(i_dataset) {
      seq_len(nrow(data)) %>% 
        sapply(function(i_sample) {
          samples_CR[[i_sample]]$samples[[i_dataset * space_size]]$sample
        }) %>% 
        t()
    })
  
  return(samples_for_testing)
}

one_step_CR <- function(current, ind_feature, 
                          params_a, params_x) {
  x_feature_star <- rbinom(n = 1, size = 1,
                           prob = 1 - params_x$pi0[ind_feature]) *
    SparseDOSSA2:::expit(rnorm(n = 1, 
                               mean = params_x$mu[ind_feature],
                               sd = params_x$sigma[ind_feature]))
  x_star <- current$sample
  x_star[ind_feature] <- x_feature_star
  x_star[!ind_feature] <- 
    (1 - x_feature_star) * 
    x_star[!ind_feature] / 
    sum(x_star[!ind_feature])
  
  logLik_star <- SparseDOSSA2:::dx(x = x_star, 
                                   pi0 = params_a$pi0,
                                   mu = params_a$mu,
                                   sigma = params_a$sigma,
                                   Omega = params_a$Omega,
                                   log.p = TRUE) -
    dZILogit(x = x_feature_star,
             pi0 = params_x$pi0[ind_feature],
             mu = params_x$mu[ind_feature],
             sigma = params_x$sigma[ind_feature],
             log.p = TRUE)
  
  if(log(runif(1)) < logLik_star - current$logLik)
    return(list(sample = x_star,
                logLik = logLik_star,
                success = TRUE))
  else
    return(list(sample = current$sample,
                logLik = current$logLik,
                success = FALSE))
}

dZILogit <- function(x, pi0, mu, sigma,
                     log.p = FALSE) {
  if(x == 0)
    logLik <- log(pi0)
  else
    logLik <- log(1 - pi0) - (SparseDOSSA2:::logit(x) - mu)^2 / sigma^2 / 2
  
  if(log.p)
    return(logLik)
  else
    return(exp(logLik))
}

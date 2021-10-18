sample_CR <- function(data, ind_feature, 
                      params,
                      sampling_method = "parallel",
                      proposal = "unif", 
                      pi0_proposal = 0.1,
                      sd_proposal = 1,
                      m = 100, space_size = 10,
                      debug_dir = NULL) {
  if(!is.null(debug_dir)) {
    dir.create(debug_dir, recursive = TRUE)
    debug_dir <- normalizePath(debug_dir)
  }
  
  samples_CR <- future.apply::future_lapply(
    seq_len(nrow(data)),
    function(i_sample) {
      x <- data[i_sample, ]
      if(sampling_method == "sequential") {
        d <- sample.int(n = m + 1, size = 1)
        samples_forward <- 
          samples_backward <- 
          list(list(sample = x,
                    logLik = 
                      log_drenorm(x = x,
                                  ind_feature = ind_feature,
                                  pi0 = params$pi0,
                                  mu = params$mu,
                                  sigma = params$sigma,
                                  Omega = params$Omega,
                                  Sigma = params$Sigma),
                    success = TRUE
          ))
        if(!is.null(debug_dir)) {
          save(samples_forward,
               file = paste0(debug_dir, 
                             "/sample_", i_sample,
                             "_samples_forward.RData"))
          save(samples_backward,
               file = paste0(debug_dir, 
                             "/sample_", i_sample,
                             "_samples_backward.RData"))
        }
        if(d < m + 1)
          for(i_step in seq_len((m + 1 - d) * space_size)) {
            samples_forward[[i_step + 1]] <- 
              one_step_CR(current = samples_forward[[i_step]],
                          ind_feature = ind_feature,
                          params = params,
                          proposal = proposal, 
                          pi0_proposal = pi0_proposal,
                          sd_proposal = sd_proposal)
            if(!is.null(debug_dir)) {
              save(samples_forward,
                   file = paste0(debug_dir, 
                                 "/sample_", i_sample,
                                 "_samples_forward.RData"))
            }
          }
        if(d > 1)
          for(i_step in seq_len((d - 1) * space_size)) {
            samples_backward[[i_step + 1]] <- 
              one_step_CR(current = samples_backward[[i_step]],
                          ind_feature = ind_feature,
                          params = params,
                          proposal = proposal, 
                          pi0_proposal = pi0_proposal,
                          sd_proposal = sd_proposal)
            if(!is.null(debug_dir)) {
              save(samples_backward,
                   file = paste0(debug_dir, 
                                 "/sample_", i_sample,
                                 "_samples_backward.RData"))
            }
          }
        return(list(samples = c(samples_forward[-1],
                                rev(samples_backward[-1])),
                    d = d))
      }
      if(sampling_method == "parallel") {
        sample_center <- list(
          sample = x,
          logLik = log_drenorm(x = x,
                               ind_feature = ind_feature,
                               pi0 = params$pi0,
                               mu = params$mu,
                               sigma = params$sigma,
                               Omega = params$Omega,
                               Sigma = params$Sigma),
          success = TRUE)
        if(!is.null(debug_dir)) {
          sample_to_center <- 
            c(sample_center,
              list(i_step = 0))
          save(sample_to_center,
               file = paste0(debug_dir, 
                             "/sample_", i_sample,
                             "_sample_to_center.RData"))
        }
        # walk back to find center of the "star"
        for(i_step in seq_len(space_size)) {
          sample_center <- one_step_CR(current = sample_center,
                                       ind_feature = ind_feature,
                                       params = params,
                                       proposal = proposal, 
                                       pi0_proposal = pi0_proposal,
                                       sd_proposal = sd_proposal)
          if(!is.null(debug_dir)) {
            sample_to_center <- 
              c(sample_center,
                list(i_step = i_step))
            save(sample_to_center,
                 file = paste0(debug_dir, 
                               "/sample_", i_sample,
                               "_sample_to_center.RData"))
          }
        }
        
        ll_samples <- list()
        for(i_dataset in seq_len(m)) {
          l_samples <- list(sample_center)
          for(i_step in seq_len(space_size)) {
            l_samples[[i_step + 1]] <- 
              one_step_CR(current = l_samples[[i_step]],
                          ind_feature = ind_feature,
                          params = params,
                          proposal = proposal, 
                          pi0_proposal = pi0_proposal,
                          sd_proposal = sd_proposal)
            if(!is.null(debug_dir)) {
              save(l_samples,
                   file = paste0(debug_dir, 
                                 "/sample_", i_sample,
                                 "_m_", i_dataset,
                                 "_samples.RData"))
            }
          }
          ll_samples[[i_dataset]] <- l_samples[-1]
        }
        return(list(samples = Reduce("c", ll_samples),
                    sample_center = sample_center))
      }
    })
  
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

sequential_CR <- function(sample0, ind_feature, 
                          params,
                          proposal = "unif", 
                          pi0_proposal = 0.1,
                          sd_proposal = 1,
                          m = 100,
                          debug_dir = NULL) {
  if(!is.null(debug_dir)) {
    dir.create(debug_dir, recursive = TRUE)
    debug_dir <- normalizePath(debug_dir)
  }
  
  samples_forward <- 
    samples_backward <- 
    list(list(sample = sample0,
              logLik = 
                log_drenorm(x = sample0,
                            ind_feature = ind_feature,
                            pi0 = params$pi0,
                            mu = params$mu,
                            sigma = params$sigma,
                            Omega = params$Omega,
                            Sigma = params$Sigma),
              success = TRUE
    ))
  
  if(!is.null(debug_dir)) {
    save(samples_forward,
         file = paste0(debug_dir, 
                       "/samples_forward.RData"))
    save(samples_backward,
         file = paste0(debug_dir, 
                       "/samples_backward.RData"))
  }
  
  # break point
  d <- sample.int(n = m, size = 1)
  if(d < m) {
    for(i_step in seq_len(m - d)) {
      samples_forward[[i_step + 1]] <- 
        one_step_CR(current = samples_forward[[i_step]],
                    ind_feature = ind_feature,
                    params = params,
                    proposal = proposal, 
                    pi0_proposal = pi0_proposal,
                    sd_proposal = sd_proposal,
                    debug_dir = debug_dir)
      
      if(!is.null(debug_dir))
        save(samples_forward,
             file = paste0(debug_dir, 
                           "/samples_forward.RData"))
    }
  }
  if(d > 1) {
    for(i_step in seq_len(d - 1)) {
      samples_backward[[i_step + 1]] <- 
        one_step_CR(current = samples_backward[[i_step]],
                    ind_feature = ind_feature,
                    params = params,
                    proposal = proposal, 
                    pi0_proposal = pi0_proposal,
                    sd_proposal = sd_proposal,
                    debug_dir = debug_dir)
      
      if(!is.null(debug_dir))
        save(samples_backward,
             file = paste0(debug_dir, 
                           "/samples_backward.RData"))
    }
  }
  
  return(
    Reduce("rbind", 
           lapply(c(rev(samples_backward[-1]),
                    samples_forward[-1]),
                  function(x) x$sample))
  )
}

sequential_CR2 <- function(sample0, ind_feature, 
                           params,
                           proposal = "unif", 
                           pi0_proposal = 0.1,
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
                log_drenorm(x = sample0,
                            ind_feature = ind_feature,
                            pi0 = params$pi0,
                            mu = params$mu,
                            sigma = params$sigma,
                            Omega = params$Omega,
                            Sigma = params$Sigma),
              success = TRUE
    ))
  
  if(!is.null(debug_dir)) {
    save(samples,
         file = paste0(debug_dir, 
                       "/samples.RData"))
  }
  
  for(i_step in seq_len(m)) {
    samples[[i_step + 1]] <- 
      one_step_CR(current = samples[[i_step]],
                  ind_feature = ind_feature,
                  params = params,
                  proposal = proposal, 
                  pi0_proposal = pi0_proposal,
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

one_step_CR <- function(current, ind_feature, 
                        params,
                        proposal,
                        pi0_proposal = 0.1,
                        sd_proposal = 1,
                        debug_dir = NULL) {
  if(proposal == "unif")
    x_feature_star <- 
      rbinom(n = 1, size = 1, prob = 1 - pi0_proposal) *
      runif(1, 0, 1)
  if(proposal == "ZILogit")
    x_feature_star <- 
      rbinom(n = 1, size = 1, prob = 1 - pi0_proposal) *
      ifelse(current$sample[ind_feature] == 0,
             runif(1, 0, 1),
             SparseDOSSA2:::expit(
               rnorm(1,
                     mean = SparseDOSSA2:::logit(current$sample[ind_feature]),
                     sd = sd_proposal))
      )
  
  if(!is.null(debug_dir))
    save(x_feature_star, file = paste0(debug_dir, "/proposal.RData"))
  
  x_star <- current$sample
  x_star[ind_feature] <- x_feature_star
  x_star[!ind_feature] <- 
    (1 - x_feature_star) * 
    x_star[!ind_feature] / 
    sum(x_star[!ind_feature])
  
  logLik_star <- 
    log_drenorm(x = x_star, 
                ind_feature = ind_feature,
                pi0 = params$pi0,
                mu = params$mu,
                sigma = params$sigma,
                Omega = params$Omega,
                Sigma = params$Sigma)
  
  log_paccept <- logLik_star - current$logLik
  if(x_feature_star == 0 & current$sample[ind_feature] != 0)
    log_paccept <- 
    log_paccept +
    log(1 - pi0_proposal) -
    log(pi0_proposal)
  if(x_feature_star != 0 & current$sample[ind_feature] == 0)
    log_paccept <- 
    log_paccept -
    log(1 - pi0_proposal) +
    log(pi0_proposal)
  if(x_feature_star != 0 & current$sample[ind_feature] != 0 &
     proposal == "ZILogit")
    log_paccept <- 
    log_paccept - 
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

log_drenorm <- function(x, ind_feature,
                        pi0, mu, sigma, Omega, Sigma) {
  
  log_dx <- SparseDOSSA2:::dx(x = x,
                              pi0 = pi0,
                              mu = mu,
                              sigma = sigma,
                              Omega = Omega,
                              Sigma = Sigma,
                              log.p = TRUE)
  
  if(x[ind_feature] != 0 & x[ind_feature] != 1)
    log_dx <- log_dx + (sum(x > 0) - 2) * log(1 - x[ind_feature])
  
  return(log_dx)
}

dZILogit <- function(x, pi0, mu, sigma,
                     log.p = FALSE) {
  if(x == 0)
    logLik <- log(pi0)
  else
    logLik <- 
      log(1 - pi0) + 
      dnorm(x = SparseDOSSA2:::logit(x),
            mean = mu,
            sd = sigma,
            log = TRUE) -
      log(x) - log(1 - x) # jacobian
  
  if(log.p)
    return(logLik)
  else
    return(exp(logLik))
}

get_marginal <- function(data) {
  mat_marginal <-
    vapply(seq_len(ncol(data)),
           function(i) 
             c(pi = pi_nonzero(data[, i, drop = TRUE]),
               mu = mean_nonzero(data[, i, drop = TRUE]),
               sigma = sd_nonzero(data[, i, drop = TRUE])),
           c(0.0, 0.0, 0.0))
  df_marginal <- data.frame(t(mat_marginal))
  rownames(df_marginal) <- colnames(data)
  return(df_marginal)
}

mean_nonzero <- function(x) {
  x <- x[x != 0]
  return(mean(log(x)))
}

sd_nonzero <- function(x) {
  x <- x[x != 0]
  if(length(x > 1))
    return(sd(log(x)))
  else
    return(0)
}

pi_nonzero <- function(x) {
  mean(x > 0)
}

get_marginal_mix <- function(data, df_marginal = get_marginal(data), 
                             mean0, var0) {
  mat_marginal_mix <- 
    vapply(seq_len(ncol(data)),
           function(i) {
             msg <- capture.output(
               fit_EM_mix <- mixtools::normalmixEM2comp(
                 x = log(data[, i, drop = TRUE]),
                 lambda = 0.5,
                 mu = c(mean0, df_marginals$mu[i]),
                 sigsqrd = c(var0, (df_marginals$sigma[i])^2),
                 verb = FALSE))
             if(fit_EM_mix$mu[1] > fit_EM_mix$mu[2])
               stop("Something went wrong!")
             return(c(pi = fit_EM_mix$lambda[2], 
                      mu0 = fit_EM_mix$mu[1],
                      mu1 = fit_EM_mix$mu[2],
                      sigma0 = fit_EM_mix$sigma[1],
                      sigma1 = fit_EM_mix$sigma[2]))},
           c(0.0, 0.0, 0.0, 0.0, 0.0))
  df_marginal_mix <- data.frame(t(mat_marginal_mix))
  
  # add overall mean and sd
  df_marginal_mix$mu_overall <- 
    df_marginal_mix$mu0 * (1 - df_marginal_mix$pi) +
    df_marginal_mix$mu1 * df_marginal_mix$pi
  df_marginal_mix$sigma_overall <- 
    sqrt(
      df_marginal_mix$sigma0^2 * (1 - df_marginal_mix$pi) +
      df_marginal_mix$sigma1^2 * df_marginal_mix$pi +
      df_marginal_mix$pi * (1 - df_marginal_mix$pi) * 
        (df_marginal_mix$mu1 - df_marginal_mix$mu0)^2
    )
  
  rownames(df_marginal_mix) <- colnames(data)
  return(df_marginal_mix)
}
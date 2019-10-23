get_marginals <- function(data) {
  mat_marginal <-
    vapply(seq_len(ncol(data)),
           function(i) 
             c(pi_nonzero(data[, i, drop = TRUE]),
               mean_nonzero(data[, i, drop = TRUE]),
               sd_nonzero(data[, i, drop = TRUE])),
           c(0.0, 0.0, 0.0))
  df_marginal <- data.frame(t(mat_marginal))
  dimnames(df_marginal) <- list(colnames(data),
                                c("pi", "mu", "sigma"))
  return(df_marginal)
}

mean_nonzero <- function(x) {
  x <- x[x != 0]
  return(mean(logit(x)))
}

sd_nonzero <- function(x) {
  x <- x[x != 0]
  if(length(x > 1))
    return(sd(logit(x)))
  else
    return(0)
}

pi_nonzero <- function(x) {
  mean(x > 0)
}

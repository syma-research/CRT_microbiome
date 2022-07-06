fit_delta <- function(delta, z) {
  if(!is.logical(delta))
    stop("delta should be logical vector.")
  if(sum(!delta) <= 2)
    return("all_one")
  if(sum(delta) <= 2)
    return("all_zero")
  
  delta <- factor(delta, levels = c(FALSE, TRUE))
  fit_control <- 
    caret::trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 5)
  tune_grid <- expand.grid(
    alpha = c(1),
    lambda = 10^seq(-4, 1, by = 0.5)
  )
  vb <- capture.output(
    fitted_delta <- 
      caret::train(x = z, 
                   y = delta,
                   method = "glmnet", 
                   family = "binomial",
                   metric = "Accuracy",
                   trControl = fit_control,
                   tuneGrid = tune_grid))
  return(fitted_delta)
}

predict_delta <- function(fitted_delta, z) {
  if("character" %in% class(fitted_delta)) {
    if(fitted_delta == "all_one")
      return(rep(1, nrow(z)))
    if(fitted_delta == "all_zero")
      return(rep(0, nrow(z)))
  }
  
  return(caret::predict.train(fitted_delta, newdata = z, type = "prob")[, 2])
}

fit_mean <- function(x, z, method = "lasso") {
  if(length(x) == 0)
    return(NULL)
  if(length(x) <= 2)
    return(mean(x))
  
  fit_control <- 
    caret::trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 5)
  if(method == "lasso") {
    tune_grid <- expand.grid(
      alpha = c(1),
      lambda = 10^seq(-4, 1, by = 0.5)
    )
    vb <- capture.output(
      fitted_mean <- 
        caret::train(x = z, 
                     y = x,
                     method = "glmnet",
                     family = "gaussian",
                     metric = "RMSE",
                     trControl = fit_control,
                     tuneGrid = tune_grid))
  } else if(method == "gbt") {
    tune_grid <- expand.grid(
      eta = c(0.1, 0.2, 0.4),
      max_depth = c(2, 4, 6), 
      colsample_bytree = c(0.6, 0.8, 1),
      min_child_weight = 1,
      subsample = 1,
      nrounds = c(5, 10, 25),
      gamma = 0
    )
    vb <- capture.output(
      fitted_mean <- 
        caret::train(x = z, 
                     y = x,
                     method = "xgbTree", 
                     metric = "RMSE",
                     trControl = fit_control,
                     tuneGrid = tune_grid))
  }
  fitted_mean$R2 <- 1 - min(fitted_mean$results$RMSE) / var(x)
  return(fitted_mean)
}

predict_mean <- function(fitted_mean, z) {
  if(is.null(fitted_mean))
    return(numeric(0))
  if("numeric" %in% class(fitted_mean))
    return(rep(fitted_mean, nrow(z)))
  
  return(caret::predict.train(fitted_mean, newdata = z))  
}

fit_var <- function(res2, z, method = "lasso") {
  if(length(res2) == 0)
    return(NULL)
  if(length(res2) <= 2)
    return(mean(res2))
  
  fit_control <- 
    caret::trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 5)
  if(method == "lasso") {
    tune_grid <- expand.grid(
      alpha = c(1),
      lambda = 10^seq(-4, 1, by = 0.5)
    )
    vb <- capture.output(
      fitted_var <- 
        caret::train(x = z, 
                     y = res2,
                     method = "glmnet",
                     family = "gaussian",
                     metric = "RMSE",
                     trControl = fit_control,
                     tuneGrid = tune_grid))
  } else if(method == "gbt") {
    tune_grid <- expand.grid(
      eta = c(0.1, 0.2, 0.4),
      max_depth = c(2, 4, 6), 
      colsample_bytree = c(0.6, 0.8, 1),
      min_child_weight = 1,
      subsample = 1,
      nrounds = c(5, 10, 25),
      gamma = 0
    )
    vb <- capture.output(
      fitted_var <- 
        caret::train(x = z, 
                     y = res2,
                     method = "xgbTree", 
                     metric = "RMSE",
                     trControl = fit_control,
                     tuneGrid = tune_grid))
  }
  fitted_var$R2 <- 1 - min(fitted_var$results$RMSE) / var(res2)
  return(fitted_var)
}

predict_var <- function(fitted_var, z, cutoff = 0.01) {
  if(is.null(fitted_var))
    return(numeric(0))
  if("numeric" %in% class(fitted_var))
    return(rep(fitted_var, nrow(z)))
  
  var_pred <- caret::predict.train(fitted_var, newdata = z)
  var_pred[var_pred <= cutoff] <- cutoff
  return(var_pred)  
}

test_oneX <- function(x, y, z, epsilon, m = 1e4,
                      method = "component",
                      method_mean = "lasso",
                      method_var = "lasso",
                      seed = 1) {
  set.seed(seed)
  
  n <- length(x)
  
  # fit y vs. z model
  fit_y <-  
    glmnet::cv.glmnet(x = z,
                      y = y,
                      family = "gaussian",
                      alpha = 1)
  yhat <- glmnet:::predict.cv.glmnet(
    fit_y, 
    newx = z,
    s = "lambda.min")[, 1]
  # else {
  #   fit_control <- 
  #     caret::trainControl(
  #       method = "repeatedcv",
  #       number = 5,
  #       repeats = 5)
  #   tune_grid <- expand.grid(
  #     eta = c(0.1, 0.2, 0.4),
  #     max_depth = c(2, 4, 6), 
  #     colsample_bytree = c(0.6, 0.8, 1),
  #     min_child_weight = 1,
  #     subsample = 1,
  #     nrounds = c(5, 10, 25),
  #     gamma = 0)
  #   # fit y vs. z model
  #   fit_y <- 
  #     caret::train(x = z, 
  #                  y = y,
  #                  method = "xgbTree", 
  #                  metric = "RMSE",
  #                  trControl = fit_control,
  #                  tuneGrid = tune_grid)
  #   yhat <- caret::predict.train(fit_y)
  # }
  y_res <- y - yhat
  
  if(method == "component") {
    # fit x vs z model
    delta <- x > 0
    fitted_delta <- fit_delta(delta = delta, z = z)
    fitted_mean <- fit_mean(x = log(x[delta] + epsilon), z = z[delta, , drop = FALSE], 
                            method = method_mean)
    x_res_nonzero <- 
      log(x[delta] + epsilon) - 
      predict_mean(fitted_mean = fitted_mean, 
                   z = z[delta, , drop = FALSE])
    fitted_var <- fit_var(res2 = 
                            x_res_nonzero^2, 
                          z = z[delta, , drop = FALSE],
                          method = method_var)
    x_res_nonzero_stand <- x_res_nonzero / 
      sqrt(predict_var(fitted_var = fitted_var, z = z[delta, , drop = FALSE]))
    
    if(sum(delta) <= 2)
      return(list(p = 1,
                  ps_full = c(1, 1, 1, 1),
                  sds = c(NA_real_, NA_real_, NA_real_),
                  fitted_delta = fitted_delta,
                  fitted_mean = fitted_mean,
                  fitted_var = fitted_var,
                  code_delta = 1))
    
    predicted_delta <- predict_delta(fitted_delta = fitted_delta, z = z)
    predicted_mean <- predict_mean(fitted_mean = fitted_mean, z = z)
    predicted_var_noModel <- rep(mean(x_res_nonzero^2), n)
    predicted_var <- predict_var(fitted_var = fitted_var, z = z)
    
    mean_x <- (1 - predicted_delta) * log(epsilon) + 
      predicted_delta * predicted_mean
    var_x_noModel <- 
      (1 - predicted_delta) * (log(epsilon) - mean_x)^2 + 
      predicted_delta * ((predicted_mean - mean_x)^2 + predicted_var_noModel)
    var_x <- 
      (1 - predicted_delta) * (log(epsilon) - mean_x)^2 + 
      predicted_delta * ((predicted_mean - mean_x)^2 + predicted_var)
    
    x_res <- log(x + epsilon) - mean_x
    deltatilde <- matrix(rbinom(n = n * m, size = 1, prob = predicted_delta),
                         nrow = n)
    xtilde_nonzero_noModel <- seq(1, m) %>% 
      vapply(
        function(i_m)
          sample(x_res_nonzero, replace = TRUE, size = n),
        rep(0.0, n)) +
      predicted_mean
    xtilde_nonzero <- seq(1, m) %>% 
      vapply(
        function(i_m)
          sample(x_res_nonzero_stand, replace = TRUE, size = n),
        rep(0.0, n)) * sqrt(predicted_var) +
      predicted_mean
    xtilde_noModel <- 
      deltatilde * xtilde_nonzero_noModel  +
      (1 - deltatilde) * log(epsilon)
    xtilde <- 
      deltatilde * xtilde_nonzero  +
      (1 - deltatilde) * log(epsilon)
    xtilde_res_noModel <- xtilde_noModel - mean_x
    xtilde_res <- xtilde - mean_x
    
    sd_x_res <- sqrt(mean(x_res^2))
    sd_xtilde_res_noModel <- sqrt(mean(xtilde_res_noModel^2))
    sd_xtilde_res <- sqrt(mean(xtilde_res^2))
    
    # not assuming heterogeneity
    stat_original <- sum(y_res * x_res)
    stat_perm <- apply(y_res * xtilde_res_noModel, 2, sum)
    p1 <- (sum(abs(stat_original) <= abs(stat_perm)) + 1) /
      (m + 1)
    # standardize
    stat_original <- stat_original / sd_x_res
    stat_perm <- stat_perm / sd_xtilde_res_noModel
    p2 <- (sum(abs(stat_original) <= abs(stat_perm)) + 1) /
      (m + 1)
    
    # user heterogeneous sampling
    stat_original <- sum(y_res * x_res)
    stat_perm <- apply(y_res * xtilde_res, 2, sum)
    p3 <- (sum(abs(stat_original) <= abs(stat_perm)) + 1) /
      (m + 1)
    # standardize
    stat_original <- stat_original / sd_x_res
    stat_perm <- stat_perm / sd_xtilde_res
    p4 <- (sum(abs(stat_original) <= abs(stat_perm)) + 1) /
      (m + 1)
    
    # # heterogeneity only due to zero component
    # stat_original <- sum(y_res * x_res / var_x_noModel)
    # stat_perm <- apply(y_res * xtilde_res_noModel / var_x_noModel, 2, sum)
    # p2 <- (sum(abs(stat_original) <= abs(stat_perm)) + 1) /
    #   (m + 1)
    # 
    # # also model heterogeneity in zero component
    # stat_original <- sum(y_res * x_res / var_x)
    # stat_perm <- apply(y_res * xtilde_res / var_x, 2, sum)
    # p3 <- (sum(abs(stat_original) <= abs(stat_perm)) + 1) /
    #   (m + 1)
    
    return(list(p = p4,
                ps_full = c(p1, p2, p3, p4),
                sds = c(sd_x_res, sd_xtilde_res_noModel, sd_xtilde_res),
                fitted_delta = fitted_delta,
                fitted_mean = fitted_mean,
                fitted_var = fitted_var,
                code_delta = 0))
  } else if(method == "naive"){
    if(sum(x > 0) <= 2)
      return(list(p = 1, fitted_naive = NULL))
    
    # first transform x
    x_trans <- log(x + epsilon)
    fitted_naive <- fit_mean(x = x_trans, z = z, method = method_mean)
    
    mean_x <- predict_mean(fitted_mean = fitted_naive, z = z)
    x_res <- x_trans - mean_x
    
    xtilde_res <- seq(1, m) %>% 
      vapply(function(i_m) {
        sample(x_res, replace = FALSE, size = n)
      },
      rep(0.0, n))
    
    stat_original <- sum(y_res * x_res)
    stat_perm <- apply(y_res * xtilde_res, 2, sum)
    p <- (sum(abs(stat_original) <= abs(stat_perm)) + 1) /
      (m + 1)
    
    return(list(p = p,
                fitted_naive = fitted_naive))
  }
}


test_oneX_clean <- function(
  x, y, z, epsilon, m = 1e4,
  family = "gaussian",
  method = "component",
  method_mean = "lasso",
  method_var = "lasso",
  seed = 1) {
  set.seed(seed)
  
  n <- length(x)
  
  # fit y vs. z model
  fit_y <-  
    glmnet::cv.glmnet(x = z,
                      y = y,
                      family = family,
                      alpha = 1)
  yhat <- glmnet:::predict.cv.glmnet(
    fit_y, 
    newx = z,
    s = "lambda.min")[, 1]
 
  # fit x vs z model
  delta <- x > 0
  fitted_delta <- fit_delta(delta = delta, z = z)
  fitted_mean <- fit_mean(x = log(x[delta] + epsilon), z = z[delta, , drop = FALSE], 
                          method = method_mean)
  x_res_nonzero <- 
    log(x[delta] + epsilon) - 
    predict_mean(fitted_mean = fitted_mean, 
                 z = z[delta, , drop = FALSE])
  fitted_var <- fit_var(res2 = 
                          x_res_nonzero^2, 
                        z = z[delta, , drop = FALSE],
                        method = method_var)
  x_res_nonzero_stand <- x_res_nonzero / 
    sqrt(predict_var(fitted_var = fitted_var, z = z[delta, , drop = FALSE]))
  
  if(sum(delta) <= 2)
    return(1)

  predicted_delta <- predict_delta(fitted_delta = fitted_delta, z = z)
  predicted_mean <- predict_mean(fitted_mean = fitted_mean, z = z)
  predicted_var <- predict_var(fitted_var = fitted_var, z = z)
  
  mean_x <- (1 - predicted_delta) * log(epsilon) + 
    predicted_delta * predicted_mean
  var_x <- 
    (1 - predicted_delta) * (log(epsilon) - mean_x)^2 + 
    predicted_delta * ((predicted_mean - mean_x)^2 + predicted_var)
  
  x_res <- log(x + epsilon) - mean_x
  deltatilde <- matrix(rbinom(n = n * m, size = 1, prob = predicted_delta),
                       nrow = n)
  while(TRUE) {
    ind_allzero <- apply(deltatilde == 0, 2, all)
    if(all(!ind_allzero))
      break
    deltatilde[, ind_allzero] <- rbinom(n = n * sum(ind_allzero), 
                                        size = 1, 
                                        prob = predicted_delta)
    
  }
  
  xtilde_nonzero <- seq(1, m) %>% 
    vapply(
      function(i_m)
        sample(x_res_nonzero_stand, replace = TRUE, size = n),
      rep(0.0, n)) * sqrt(predicted_var) +
    predicted_mean
  xtilde <- 
    deltatilde * xtilde_nonzero  +
    (1 - deltatilde) * log(epsilon)
  xtilde_res <- xtilde - mean_x
  
  beta_original <- coef(glm(y ~ {x_res / sd(x_res)}, family = family, offset = yhat))[2]
  beta_perm <- apply(xtilde_res, 
                     2,
                     function(i_xtilde_res) {
                       coef(glm(y ~ {i_xtilde_res / sd(i_xtilde_res)}, family = family, offset = yhat))[2]
                     })
  p <- (sum(abs(beta_original) <= abs(beta_perm)) + 1) /
    (m + 1)
  
  return(p)
}


microTensor <- 
  function(x, y, covariates = NULL, 
           family = "gaussian",
           method = "component", 
           method_mean = "lasso", method_var = "lasso",
           epsilon, m = 1e4, seed = 1,
           debug_path = NULL) {
    
    ps <- c()
    for(i_feature in seq(1, ncol(x)))  {
      i_x <- x[, i_feature]
      i_z <- x[, -i_feature, drop = FALSE]
      i_z <- t(apply(i_z, 1, tss_withzero))
      i_z <- log(i_z + epsilon)
      if(!is.null(covariates))
        i_z <- cbind(i_z, covariates)
      
      i_fit <- 
        test_oneX_clean(
          y = y, x = i_x, z = i_z,
          epsilon = epsilon,
          m = m,
          family = family,
          method = method,
          method_mean = method_mean,
          method_var = method_var,
          seed = seed * i_feature)
      
      ps <- c(ps, i_fit)
      if(!is.null(debug_path) & i_feature %% 20 == 0) {
        save(ps, file = debug_path)
      }
    }
    
    return(ps)
  }


# stat_oneX_nodewise <- function(x, y, z, epsilon, 
#                                m = 99,
#                                method_mean = "lasso",
#                                method_var = "lasso") {
#   n <- length(x)
#   log_z <- log(z + epsilon)
#   
#   # fit x vs z model
#   delta <- x > 0
#   if(sum(delta) <= 2) {
#     return(c(pval = 1, stat = 0))
#   }
#   
#   fitted_delta <- fit_delta(delta = delta, z = log_z)
#   fitted_mean <- fit_mean(x = log(x[delta]) - log(1 - x[delta]), 
#                           z = log_z[delta, , drop = FALSE], 
#                           method = method_mean)
#   x_res_nonzero <- 
#     log(x[delta]) - log(1 - x[delta]) -
#     predict_mean(fitted_mean = fitted_mean, 
#                  z = log_z[delta, , drop = FALSE])
#   fitted_var <- fit_var(res2 = 
#                           x_res_nonzero^2, 
#                         z = log_z[delta, , drop = FALSE],
#                         method = method_var)
#   x_res_nonzero_stand <- x_res_nonzero / 
#     sqrt(predict_var(fitted_var = fitted_var, z = log_z[delta, , drop = FALSE]))
#   
#   predicted_delta <- predict_delta(fitted_delta = fitted_delta, z = log_z)
#   predicted_mean <- predict_mean(fitted_mean = fitted_mean, z = log_z)
#   # predicted_var_noModel <- rep(mean(x_res_nonzero^2), n)
#   predicted_var <- predict_var(fitted_var = fitted_var, z = log_z)
#   
#   deltatilde <- matrix(rbinom(n = n * m, size = 1, prob = predicted_delta),
#                        nrow = n)
#   # xtilde_nonzero_noModel <- seq(1, m) %>% 
#   #   vapply(
#   #     function(i_m)
#   #       sample(x_res_nonzero, replace = TRUE, size = n),
#   #     rep(0.0, n)) %>% 
#   #   {exp(. + predicted_mean)}
#   xtilde_nonzero <- seq(1, m) %>% 
#     vapply(
#       function(i_m)
#         sample(x_res_nonzero_stand, replace = TRUE, size = n),
#       rep(0.0, n))  %>% 
#     {. * sqrt(predicted_var) + predicted_mean} %>% 
#     {exp(.) / (1 + exp(.))}
#   
#   xtilde <- 
#     deltatilde * xtilde_nonzero
#   
#   # fit y model
#   fit_y <-  
#     glmnet::cv.glmnet(x = log(cbind(x, z) + epsilon),
#                       y = y,
#                       family = "gaussian",
#                       alpha = 1)
#   stat_original <- 
#     glmnet:::coef.cv.glmnet(
#       fit_y, 
#       s = "lambda.min")[2, 1]
#   stat_tilde <- seq(1, m) %>% 
#     vapply(function(i_m) {
#       fit_y <-  
#         glmnet::cv.glmnet(x = log(cbind(xtilde[, i_m], 
#                                         z * (1 - xtilde[, i_m]) / (1 - x)) + 
#                                     epsilon),
#                           y = y,
#                           family = "gaussian",
#                           alpha = 1)
#       glmnet:::coef.cv.glmnet(
#         fit_y, 
#         s = "lambda.min")[2, 1]
#     },
#     0.0)
#   
#   pval <- (sum(abs(stat_original) <= abs(stat_tilde)) + 1) /
#     (m + 1)
#   stat <- max(c(abs(stat_original), abs(stat_tilde))) -
#     median(c(abs(stat_original), abs(stat_tilde)))
#   
#   return(c(pval = pval, stat = stat))
# }
# 
# seqstep_select <- function(pvals, stats, 
#                            q = 0.2,
#                            c_const = 0.2501) {
#   features <- seq(1, length(pvals))[order(-stats)]
#   pvals <- pvals[order(-stats)]
#   ind_select <- seq(1, length(features)) %>% 
#     purrr::map_lgl(function(k) 
#       ((1 + sum(pvals[seq(1, k)] > c_const)) / 
#          max(c(sum(pvals[seq(1, k)] <= c_const), 1))) <= 
#         (1 - c_const) / c_const * q
#     )
#   if(all(!ind_select))
#     return(rep(FALSE, length(features)))
#   else {
#     k_max <- max(seq(1, length(features))[ind_select])
#     return(seq(1, length(features)) %in% 
#              features[seq(1, k_max)[pvals[seq(1, k_max)] <= c_const]])
#   }
# }
# 
# stat_oneX_nodewise_dCRT <- function(
#   x, y, z, epsilon, m = 1e4,
#   method_mean = "lasso",
#   method_var = "lasso") {
#   
#   n <- length(x)
#   
#   # fit y vs. z model
#   fit_y <-  
#     glmnet::cv.glmnet(x = z,
#                       y = y,
#                       family = "gaussian",
#                       alpha = 1)
#   yhat <- glmnet:::predict.cv.glmnet(
#     fit_y, 
#     newx = z,
#     s = "lambda.min")[, 1]
#   y_res <- y - yhat
#   
#   # fit x vs z model
#   delta <- x > 0
#   if(sum(delta) <= 2)
#     return(c(pval = 1, stat = 0))
#   fitted_delta <- fit_delta(delta = delta, z = z)
#   
#   fitted_mean <- fit_mean(x = log(x[delta] + epsilon), z = z[delta, , drop = FALSE], 
#                           method = method_mean)
#   x_res_nonzero <- 
#     log(x[delta] + epsilon) - 
#     predict_mean(fitted_mean = fitted_mean, 
#                  z = z[delta, , drop = FALSE])
#   fitted_var <- fit_var(res2 = 
#                           x_res_nonzero^2, 
#                         z = z[delta, , drop = FALSE],
#                         method = method_var)
#   x_res_nonzero_stand <- x_res_nonzero / 
#     sqrt(predict_var(fitted_var = fitted_var, z = z[delta, , drop = FALSE]))
#   
#   predicted_delta <- predict_delta(fitted_delta = fitted_delta, z = z)
#   predicted_mean <- predict_mean(fitted_mean = fitted_mean, z = z)
#   predicted_var <- predict_var(fitted_var = fitted_var, z = z)
#   
#   mean_x <- (1 - predicted_delta) * log(epsilon) + 
#     predicted_delta * predicted_mean
#   var_x <- 
#     (1 - predicted_delta) * (log(epsilon) - mean_x)^2 + 
#     predicted_delta * ((predicted_mean - mean_x)^2 + predicted_var)
#   
#   x_res <- log(x + epsilon) - mean_x
#   deltatilde <- matrix(rbinom(n = n * m, size = 1, prob = predicted_delta),
#                        nrow = n)
#   xtilde_nonzero <- seq(1, m) %>% 
#     vapply(
#       function(i_m)
#         sample(x_res_nonzero_stand, replace = TRUE, size = n),
#       rep(0.0, n)) * sqrt(predicted_var) +
#     predicted_mean
#   xtilde <- 
#     deltatilde * xtilde_nonzero  +
#     (1 - deltatilde) * log(epsilon)
#   xtilde_res <- xtilde - mean_x
#   
#   sd_x_res <- sqrt(mean(x_res^2))
#   sd_xtilde_res <- sqrt(mean(xtilde_res^2))
#   
#   # user heterogeneous sampling
#   stat_original <- sum(y_res * x_res) / sd_x_res
#   stat_tilde <- apply(y_res * xtilde_res, 2, sum) / sd_xtilde_res
#   
#   pval <- (sum(abs(stat_original) <= abs(stat_tilde)) + 1) /
#     (m + 1)
#   stat <- max(c(abs(stat_original), abs(stat_tilde))) -
#     median(c(abs(stat_original), abs(stat_tilde)))
#   
#   return(c(pval = pval, stat = stat))
# }

simulate_y <- function(x, effect_size, n_signal, epsilon, 
                       family = "gaussian",
                       seed = 0) {
  if(n_signal != round(n_signal))
    stop("n_signal must be an integer!")
  # if(n_signal <= 1)
  #   stop("Needs at least two signals!")
  
  set.seed(seed)
  ind_TP <- seq(1, ncol(x)) %in%
    sample(seq(1, ncol(x)), size = n_signal, replace = FALSE)
  effects <- runif(min = -2, max = 2, n = n_signal) * effect_size
  
  xTP_renorm <- x[, ind_TP, drop = FALSE] %>%
    apply(1, tss_withzero) %>%
    t()
  # log_xTP_renorm_cont <-
  #   log(xTP_renorm[, ceiling(seq(1, n_signal/2)), drop = FALSE] + epsilon)
  # xTP_bin <- 
  #   (xTP_renorm[, seq(ceiling(seq(1, n_signal/2)) + 1, n_signal), drop = FALSE] > 0) * 1
  # x_TP <- cbind(log_xTP_renorm_cont, xTP_bin)
  x_TP <- log(xTP_renorm + epsilon)
  x_TP <-
    apply(x_TP, 2, function(x) x - mean(x))
  
  link_ey <- (x_TP %*% effects)[, 1]
  if(family == "gaussian") {
    y <- link_ey + rnorm(n = nrow(x), mean = 0, sd = 1)
  }
  else if (family == "binomial") {
    y <- rbinom(n = nrow(x), size = 1, prob = exp(link_ey) / (1 + exp(link_ey)))
  }
  snr <- var(link_ey) / 1
  
  return(list(ind_TP = ind_TP,
              effects = effects,
              x_TP = x_TP,
              link_ey = link_ey,
              y = y,
              snr = snr,
              seed = seed))
}

test_oneX_binaryFix <- function(
  x, y, z, epsilon, m = 1e4,
  family = "binomial",
  method = "component",
  method_mean = "lasso",
  method_var = "lasso",
  seed = 1) {
  set.seed(seed)
  
  n <- length(x)
  
  # fit y vs. z model
  fit_y <-  
    glmnet::cv.glmnet(x = z,
                      y = y,
                      family = family,
                      alpha = 1)
  yhat <- glmnet:::predict.cv.glmnet(
    fit_y, 
    newx = z,
    s = "lambda.min")[, 1]
  
  # fit x vs z model
  delta <- x > 0
  fitted_delta <- fit_delta(delta = delta, z = z)
  fitted_mean <- fit_mean(x = log(x[delta] + epsilon), z = z[delta, , drop = FALSE], 
                          method = method_mean)
  x_res_nonzero <- 
    log(x[delta] + epsilon) - 
    predict_mean(fitted_mean = fitted_mean, 
                 z = z[delta, , drop = FALSE])
  fitted_var <- fit_var(res2 = 
                          x_res_nonzero^2, 
                        z = z[delta, , drop = FALSE],
                        method = method_var)
  x_res_nonzero_stand <- x_res_nonzero / 
    sqrt(predict_var(fitted_var = fitted_var, z = z[delta, , drop = FALSE]))
  
  if(sum(delta) <= 2)
    return(1)
  
  predicted_delta <- predict_delta(fitted_delta = fitted_delta, z = z)
  predicted_mean <- predict_mean(fitted_mean = fitted_mean, z = z)
  predicted_var <- predict_var(fitted_var = fitted_var, z = z)
  
  mean_x <- (1 - predicted_delta) * log(epsilon) + 
    predicted_delta * predicted_mean
  var_x <- 
    (1 - predicted_delta) * (log(epsilon) - mean_x)^2 + 
    predicted_delta * ((predicted_mean - mean_x)^2 + predicted_var)
  
  x_res <- log(x + epsilon) - mean_x
  deltatilde <- matrix(rbinom(n = n * m, size = 1, prob = predicted_delta),
                       nrow = n)
  while(TRUE) {
    ind_allzero <- apply(deltatilde == 0, 2, all)
    if(all(!ind_allzero))
      break
    deltatilde[, ind_allzero] <- rbinom(n = n * sum(ind_allzero), 
                                        size = 1, 
                                        prob = predicted_delta)
    
  }
  
  xtilde_nonzero <- seq(1, m) %>% 
    vapply(
      function(i_m)
        sample(x_res_nonzero_stand, replace = TRUE, size = n),
      rep(0.0, n)) * sqrt(predicted_var) +
    predicted_mean
  xtilde <- 
    deltatilde * xtilde_nonzero  +
    (1 - deltatilde) * log(epsilon)
  xtilde_res <- xtilde - mean_x
  
  # use logistf fit
  # beta_original <- coef(logistf::logistf(y ~ {x_res / sd(x_res)}, offset = yhat))[2]
  # beta_perm <- apply(xtilde_res, 
  #                    2,
  #                    function(i_xtilde_res) {
  #                      coef(logistf::logistf(y ~ {i_xtilde_res / sd(i_xtilde_res)}, offset = yhat))[2]
  #                    })
  # p1 <- (sum(abs(beta_original) <= abs(beta_perm)) + 1) /
  #   (m + 1)
  
  # use linear fit
  yhat_linear <- exp(yhat) / (1 + exp(yhat))
  beta_original <- coef(lm(y ~ {x_res / sd(x_res)}, offset = yhat_linear))[2]
  beta_perm <- apply(xtilde_res, 
                     2,
                     function(i_xtilde_res) {
                       coef(lm(y ~ {i_xtilde_res / sd(i_xtilde_res)}, offset = yhat_linear))[2]
                     })
  p2 <- (sum(abs(beta_original) <= abs(beta_perm)) + 1) /
    (m + 1)
  
  # return(c(p1, p2))
  return(p2)
}

microTensor_binaryFix <- 
  function(x, y, covariates = NULL, 
           family = "binomial",
           method = "component", 
           method_mean = "lasso", method_var = "lasso",
           epsilon, m = 1e4, seed = 1,
           debug_path = NULL) {
    
    # ps <- matrix(nrow = 0, ncol = 2)
    ps <- c()
    for(i_feature in seq(1, ncol(x)))  {
      i_x <- x[, i_feature]
      i_z <- x[, -i_feature, drop = FALSE]
      i_z <- t(apply(i_z, 1, tss_withzero))
      i_z <- log(i_z + epsilon)
      if(!is.null(covariates))
        i_z <- cbind(i_z, covariates)
      
      i_fit <- 
        test_oneX_binaryFix(
          y = y, x = i_x, z = i_z,
          epsilon = epsilon,
          m = m,
          family = family,
          method = method,
          method_mean = method_mean,
          method_var = method_var,
          seed = seed * i_feature)
      
      # ps <- rbind(ps, i_fit)
      ps <- c(ps, i_fit)
      if(!is.null(debug_path) & i_feature %% 20 == 0) {
        save(ps, file = debug_path)
      }
    }
    
    return(ps)
  }

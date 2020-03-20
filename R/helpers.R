logit <- function(x) log(x) - log(1 - x)

expit <- function(x) 
  exp(x) / (1 + exp(x))

cor2 <- function(x, y = x, 
                 method = "spearman", 
                 random = TRUE, 
                 sim = FALSE,
                 R = 30,
                 ncores = 6) {
  if(any(is.na(x)) | any(is.na(y)))
    stop("Missing values currently not supported!")
  if (nrow(x) != nrow(y)) 
    stop("Number of rows between x and y must agree!")
  if (!(all(x >= 0) & all(y >= 0))) 
    stop("Only suitable for feature abundance (i.e., non-negative) table!")
  method <- match.arg(method, c("spearman"))
  if (!random | (all(x > 0) & all(y > 0))) {
    return(cor(x, y, method = method))
  }
  if (!sim) {
    x_fill <- vapply(seq_len(ncol(x)),
                     function(p) {
                       rank(x[, p, drop = TRUE])
                     },
                     rep(0.0, nrow(x)))
    y_fill <- vapply(seq_len(ncol(y)),
                     function(p) {
                       rank(y[, p, drop = TRUE])
                     },
                     rep(0.0, nrow(y)))
    return(vapply(seq_len(ncol(x)), 
                  function(p_x)
                    vapply(seq_len(ncol(y)),
                           function(p_y)
                             rank_spearman(x_fill[, p_x, drop = TRUE],
                                           y_fill[, p_y, drop = TRUE]),
                           0.0),
                  rep(0.0, ncol(y))))
  }
  else {
    ind_x <- x != 0
    ind_y <- y != 0
    minMeasure_x <- min(setdiff(x, 0))/2
    minMeasure_y <- min(setdiff(y, 0))/2
    x_fill <- x
    y_fill <- y
    doParallel::registerDoParallel(cores = ncores)
    sum_cor <- foreach::`%dopar%`(
      foreach::foreach(r = seq_len(R),
                       .combine = "+"),
      {
        x_fill[!ind_x] <- runif(n = sum(!ind_x), min = -minMeasure_x, 
                                max = minMeasure_x)
        y_fill[!ind_y] <- runif(n = sum(!ind_y), min = -minMeasure_y, 
                                max = minMeasure_y)
        return(cor(x_fill, y_fill, method = method))
      }
    )
    doParallel::stopImplicitCluster()
    average_cor <- sum_cor / R
    # to ensure symmetry
    return(enforce_symm(average_cor, method = "average"))
  }
}

rank_spearman <- function(x, y) {
  if(length(x) != length(y))
    stop("Number of observations for x and y must agree!")
  n <- length(x)
  sum(x * y) / (n - 1) / n / (n + 1) * 12 - (n + 1) / (n - 1) * 3
}

lower_tri <- function(x, warning = TRUE) {
  if(!isSymmetric(x) & warning) 
    warning("x is not symmetric!")
  
  x[lower.tri(x)]
}

upper_tri <- function(x, warning = TRUE) {
  if(!isSymmetric(x) & warning) 
    warning("x is not symmetric!")
  
  x[upper.tri(x)]
}

enforce_symm <- function(x, method = "upper") {
  if(nrow(x) != ncol(x)) 
    stop("x does not appear to be a covariance matrix!")
  x_out <- x
  if(!isSymmetric(x_out)) {
    if(method == "average") {
      lower_averaged <- (lower_tri(x_out, warning = FALSE) + 
                           lower_tri(t(x_out), warning = FALSE)) / 2
      x_out[lower.tri(x_out)] <- lower_averaged
      x_out[upper.tri(x_out)] <- 
        t(x_out)[upper.tri(x_out)]
    }
    if(method == "lower")
      x_out[upper.tri(x_out)] <- upper_tri(t(x_out), warning = FALSE)
    if(method == "upper")
      x_out[lower.tri(x_out)] <- lower_tri(t(x_out), warning = FALSE)
  }
  
  return(x_out)
}

simplify_feature <- function(x) {
  mat_x <- Reduce("rbind", strsplit(x, "|", fixed = TRUE))
  mat_x_ind <- gsub("NA", "", gsub(".*\\_\\_", "", mat_x)) != ""
  x_simple <- vapply(
    seq_along(x),
    function(i) {
      rev(mat_x[i, , drop = TRUE][mat_x_ind[i, , drop = TRUE]])[1]
    },
    "")
  if(anyDuplicated(x_simple))
    stop("Duplicates after simplifying!")
  
  return(x_simple)
}

longify_abd <- function(x, 
                        feature = "feature", 
                        sample = "sample",
                        abundance = "abd") {
  tb_long <- x %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(feature) %>% 
    tidyr::pivot_longer(-tidyselect::contains(feature), 
                        names_to = sample, values_to = abundance)
}
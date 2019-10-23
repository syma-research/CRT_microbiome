plot_cors <- function(cor1, cor2, labels = c("cor1", "cor2")) {
  df_plot <- data.frame(cor1 = lower_tri(cor1),
                        cor2 = lower_tri(cor2))
  
  p <- ggplot(aes(x = cor1, y = cor2), data = df_plot) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    xlab(labels[1]) + ylab(labels[2]) 
}

log_with_zero <- function(x, small_val = 1e-16) log10(x + small_val) 

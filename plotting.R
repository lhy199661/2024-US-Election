plot_2024_results <- function(Y,
                              monitor_result,
                              filtered,
                              time,
                              train_size,
                              r) {
  
  l <- monitor_result$results
  
  res <- fit_mu_online(Y, k = 15, verbose = TRUE)
  
  r <- unlist(
    lapply(res$mu_history, function(x) tail(x,1))
  )
  t <- as.numeric(predict(fit_mu_all(Y[1:10000])))
  r[1:10000] <- t
  r <- c(r,NA)
  
  df <- data.frame(
    time = time,
    Y = Y,
    Upper = l$Upper_cal,
    r = r
  )
  
  df_filtered <- data.frame(
    time = time[filtered],
    Y = Y[filtered]
  )
  
  p <- ggplot(df, aes(x = time)) +
    
    geom_line(
      aes(y = Y),
      color = "gray60"
    ) +
    
    geom_line(
      aes(y = Upper),
      color = "red",
      linewidth = 0.7
    ) +
    
    geom_line(
      aes(y = r),
      color = "green3",
      linewidth = 0.7
    ) +
    
    geom_point(
      data = df_filtered,
      aes(y = Y),
      shape = 4,
      color = "blue",
      size = 0.8 * 2.5
    ) +
    
    geom_vline(
      xintercept = time[train_size],
      color = "magenta",
      linetype = "dashed"
    ) +
    
    labs(
      x = "Date",
      y = "Sum Of Probability"
    ) +
    
    theme_bw()
  
  print(p)
}


plot_2020_results <- function(Y,
                              monitor_result,
                              filtered,
                              time,
                              train_size,
                              r) {
  
  l <- monitor_result$results
  
  res <- fit_mu_online(Y, k = 15, verbose = TRUE)
  r <- unlist(
    lapply(res$mu_history, function(x) tail(x,1))
  )
  t <- as.numeric(predict(fit_mu_all(Y[1:10000])))
  r[1:10000] <- t
  r <- c(r,NA)
  r <- append(r, rep(NA,1079), after=11593)
  
  Y <- append(Y, rep(NA,1079), after=11593)
  
  
  Ul <- l$Upper_cal
  Ul <- append(Ul, rep(NA,1079), after=11593)
  
  Ul_break <- Ul
  
  filtered <- ifelse(filtered < 11594,
                     filtered,
                     filtered + 1079)
  
  df_Yp <- data.frame(
    time = time,
    Y = Y
  )
  
  df_covid <- data.frame(
    time = time,
    covid = r
  )
  
  df_Ul <- data.frame(
    time = time,
    Ul = Ul_break
  )
  
  df_points <- data.frame(
    time = time[filtered],
    Y = Y[filtered]
  )
  
  q <- ggplot() +
    geom_line(
      data = df_Yp,
      aes(x = time, y = Y),
      color = "gray60"
    ) +
    
    geom_line(
      data = df_covid,
      aes(x = time, y = covid),
      color = "green3"
    ) +
    
    geom_line(
      data = df_Ul,
      aes(x = time, y = Ul),
      color = "red"
    ) +
    
    geom_point(
      data = df_points,
      aes(x = time, y = Y),
      shape = 4,
      color = "blue",
      size = 1
    ) +
    
    geom_vline(
      xintercept = time[train_size],
      color = "magenta",
      linetype = "dashed"
    ) +
    
    geom_vline(
      xintercept = time[11593],
      color = "orange",
      linetype = "dashed"
    ) +
    
    geom_vline(
      xintercept = time[12673],
      color = "orange",
      linetype = "dashed"
    ) +
    
    theme_bw() +
    
    labs(
      x = "Date",
      y = "Sum Of Probability"
    )
  
  print(q)
}

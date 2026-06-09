# =========================================================
# plotting.R
# Plotting functions for 2020 and 2024 results
# =========================================================

plot_2024_results <- function(Y,
                              monitor_result,
                              filtered,
                              time,
                              train_size) {
  
  l <- monitor_result$results
  
  idx <- l$original_index
  
  df <- data.frame(
    time = time[idx],
    Y = l$Z,
    Upper = l$Upper,
    Mu_hat = l$Mu_hat
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
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    
    geom_line(
      aes(y = Mu_hat),
      color = "green3",
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    
    geom_point(
      data = df_filtered,
      aes(x = time, y = Y),
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
  invisible(p)
}

plot_2020_results <- function(Y,
                              monitor_result,
                              filtered,
                              time,
                              train_size,
                              covid_insert_after = 11593,
                              covid_gap = 1079,
                              covid_start_full = 11594,
                              covid_end_full = 12672) {
  
  l <- monitor_result$results
  
  # --------------------------------------------------
  # Map effective-series results back to no-COVID index
  # --------------------------------------------------
  
  Mu_nocovid <- rep(NA_real_, length(Y))
  Upper_nocovid <- rep(NA_real_, length(Y))
  PredMean_nocovid <- rep(NA_real_, length(Y))
  
  Mu_nocovid[l$original_index] <- l$Mu_hat
  Upper_nocovid[l$original_index] <- l$Upper
  PredMean_nocovid[l$original_index] <- l$PredMean
  
  # --------------------------------------------------
  # Insert COVID gap to match full 2020 time vector
  # --------------------------------------------------
  
  Y_plot <- append(
    Y,
    rep(NA, covid_gap),
    after = covid_insert_after
  )
  
  Mu_plot <- append(
    Mu_nocovid,
    rep(NA, covid_gap),
    after = covid_insert_after
  )
  
  Upper_plot <- append(
    Upper_nocovid,
    rep(NA, covid_gap),
    after = covid_insert_after
  )
  
  PredMean_plot <- append(
    PredMean_nocovid,
    rep(NA, covid_gap),
    after = covid_insert_after
  )
  
  # --------------------------------------------------
  # Signal points are original no-COVID indices
  # Convert them to full-data indices for plotting
  # --------------------------------------------------
  
  filtered_plot <- ifelse(
    filtered <= covid_insert_after,
    filtered,
    filtered + covid_gap
  )
  
  filtered_plot <- filtered_plot[
    !is.na(filtered_plot) &
      filtered_plot >= 1 &
      filtered_plot <= length(time)
  ]
  
  df_plot <- data.frame(
    time = time,
    Y = Y_plot,
    Mu_hat = Mu_plot,
    Upper = Upper_plot,
    PredMean = PredMean_plot
  )
  
  df_points <- data.frame(
    time = time[filtered_plot],
    Y = Y_plot[filtered_plot]
  )
  
  q <- ggplot(df_plot, aes(x = time)) +
    
    geom_line(
      aes(y = Y),
      color = "gray60"
    ) +
    
    geom_line(
      aes(y = Upper),
      color = "red",
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    
    geom_line(
      aes(y = Mu_hat),
      color = "green3",
      linewidth = 0.7,
      na.rm = TRUE
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
      xintercept = time[covid_start_full],
      color = "orange",
      linetype = "dashed"
    ) +
    
    geom_vline(
      xintercept = time[covid_end_full],
      color = "orange",
      linetype = "dashed"
    ) +
    
    theme_bw() +
    
    labs(
      x = "Date",
      y = "Sum Of Probability"
    )
  
  print(q)
  invisible(q)
}
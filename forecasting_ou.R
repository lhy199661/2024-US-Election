# =========================
# Run OU Forecast
# =========================

run_ou_forecast <- function(Y, train_size, B) {
  
  bootstrap_forecast_calibrated(
    Y,
    train_size = train_size,
    B = 200
  )
}

# =========================
# Find Signal Points
# =========================

get_signal_points <- function(monitor_result) {
  
  l <- monitor_result$results
  
  l_true <- l[l$Viol_cal == TRUE, ]
  
  change_points <- find_first_changes(
    l_true$time,
    l_true$Z
  )
  
  filtered <- change_points$time
  
  na.omit(filtered)
}
# =========================================================
# forecasting_ou.R
# OU forecasting wrapper and signal extraction
# =========================================================

run_ou_forecast <- function(Y,
                            train_size,
                            B = 200,
                            alpha = 0.05,
                            k_gam = 15,
                            cores_limit = NULL,
                            verbose = TRUE,
                            skip_flat = TRUE,
                            time_index = NULL) {
  
  bootstrap_forecast(
    Z = Y,
    B = B,
    alpha = alpha,
    train_size = train_size,
    k_gam = k_gam,
    cores_limit = cores_limit,
    verbose = verbose,
    skip_flat = skip_flat,
    time_index = time_index
  )
}

get_signal_points <- function(monitor_result) {
  get_signal_points_from_results(
    monitor_result = monitor_result
  )
}

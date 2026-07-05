# ============================================================
# Model_checking.R
# Simulation, Model Selection, Reliability, and Diagnostics
# ============================================================

source("functions.R")
source("data_2020.R")
source("data_2024.R")
source("forecasting_ou.R")
source("plotting.R")


library(ggplot2)
library(mgcv)
library(gridExtra)

output_dir <- file.path("outputs", "model_checking_figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================
# 1. Monte Carlo Simulation
# ============================================================

simulate_trending_ou_exact <- function(n,
                                       theta,
                                       sigma2,
                                       omega2,
                                       seed) {
  
  set.seed(seed)
  
  time_index <- seq_len(n)
  h <- 1
  
  M_t <- 5 +
    0.001 * time_index +
    0.5 * sin(2 * pi * time_index / 200)
  
  Mprime_t <- numeric(n)
  Mprime_t[1] <- M_t[2] - M_t[1]
  Mprime_t[n] <- M_t[n] - M_t[n - 1]
  Mprime_t[2:(n - 1)] <- (M_t[3:n] - M_t[1:(n - 2)]) / 2
  
  mu_t <- M_t + Mprime_t / theta
  
  theta <- max(theta, 1e-6)
  sigma2 <- max(sigma2, 1e-6)
  omega2 <- max(omega2, 1e-6)
  
  a <- exp(-theta * h)
  q_var <- sigma2 / (2 * theta) * (1 - exp(-2 * theta * h))
  
  Q <- numeric(n)
  Z <- numeric(n)
  
  Q[1] <- M_t[1]
  Z[1] <- Q[1] + rnorm(1, 0, sqrt(omega2))
  
  for (t in 2:n) {
    
    Q[t] <- a * Q[t - 1] +
      (1 - a) * mu_t[t] +
      rnorm(1, 0, sqrt(q_var))
    
    Z[t] <- Q[t] +
      rnorm(1, 0, sqrt(omega2))
  }
  
  list(
    Z = Z,
    Q = Q,
    M_t = M_t,
    mu_t = mu_t,
    true_params = list(
      theta = theta,
      sigma2 = sigma2,
      omega2 = omega2
    )
  )
}

estimate_ou_params_sim <- function(Z,
                                   M_hat,
                                   h_vec = NULL,
                                   time_index = NULL,
                                   init_theta = NULL,
                                   max_theta = 2,
                                   max_var = 5) {
  
  n <- length(Z)
  
  if (is.null(h_vec)) {
    h_vec <- rep(1, n)
  }
  
  if (is.null(time_index)) {
    time_index <- seq_len(n)
  }
  
  Mprime <- num_derivative(
    M_hat,
    time_index = time_index
  )
  
  if (is.null(init_theta)) {
    h_trans <- pmax(h_vec[-1], 1e-8)
    dZ <- diff(Z)
    Xreg <- (M_hat[-n] - Z[-n]) * h_trans
    
    fit0 <- tryCatch(
      lm(dZ ~ 0 + Xreg),
      error = function(e) NULL
    )
    
    init_theta <- if (!is.null(fit0)) {
      as.numeric(coef(fit0)[1])
    } else {
      NA_real_
    }
    
    if (!is.finite(init_theta) || init_theta <= 0) {
      init_theta <- 0.1
    }
  }
  
  init_var <- max(var(diff(Z), na.rm = TRUE), 1e-4)
  
  init_params <- c(
    init_theta,
    log(init_var / 2),
    log(init_var / 2)
  )
  
  opt <- optim(
    par = init_params,
    fn = log_likelihood_ou_recursive,
    method = "L-BFGS-B",
    lower = c(
      1e-6,
      log(1e-6),
      log(1e-6)
    ),
    upper = c(
      max_theta,
      log(max_var),
      log(max_var)
    ),
    control = list(maxit = 1000),
    Z = Z,
    M_hat = M_hat,
    Mprime = Mprime,
    h_vec = h_vec,
    time_index = time_index
  )
  
  list(
    theta = opt$par[1],
    sigma2 = exp(opt$par[2]),
    omega2 = exp(opt$par[3]),
    convergence = opt$convergence
  )
}

estimate_from_simulation <- function(seed,
                                     n = 1000,
                                     train_size = 500,
                                     theta_true = 0.15,
                                     sigma2_true = 0.55,
                                     omega2_true = 0.20,
                                     use_true_M = FALSE) {
  
  sim <- simulate_trending_ou_exact(
    n = n,
    theta = theta_true,
    sigma2 = sigma2_true,
    omega2 = omega2_true,
    seed = seed
  )
  
  Z_train <- sim$Z[1:train_size]
  time_train <- seq_len(train_size)
  h_train <- rep(1, train_size)
  
  if (use_true_M) {
    M_train <- sim$M_t[1:train_size]
  } else {
    fit <- fit_mu(
      Z_train,
      k = 15,
      time_index = time_train
    )
    M_train <- fit$M_hat
  }
  
  est <- estimate_ou_params_sim(
    Z = Z_train,
    M_hat = M_train,
    h_vec = h_train,
    time_index = time_train,
    init_theta = NULL,
    max_theta = 1,
    max_var = 5
  )
  
  data.frame(
    theta = est$theta,
    sigma2 = est$sigma2,
    omega2 = est$omega2,
    convergence = est$convergence
  )
}

set.seed(34351)

# ------------------------------------------------------------
# Simulate data
# ------------------------------------------------------------

sim_data <- simulate_trending_ou_exact(
  n = 1000,
  theta = 0.15,
  sigma2 = 0.55,
  omega2 = 0.20,
  seed = 123
)

Y_sim <- sim_data$Z

time_sim <- seq_along(Y_sim)
time_index <- time_sim

# ------------------------------------------------------------
# Run forecasting model
# ------------------------------------------------------------

monitor_result_sim <- run_ou_forecast(
  Y = Y_sim,
  train_size = 500,
  B = 200,
  skip_flat = TRUE
)

# ------------------------------------------------------------
# Extract signal points
# ------------------------------------------------------------

filtered_sim <- get_signal_points(
  monitor_result_sim
)

cat(
  "Number of simulation signals:",
  length(filtered_sim),
  "\n"
)

# ------------------------------------------------------------
# Plot simulation result
# Use plot_sim_results because no Covid-gap correction is needed
# ------------------------------------------------------------

p_sim <- plot_sim_results(
  Y = Y_sim,
  monitor_result = monitor_result_sim,
  filtered = filtered_sim,
  time = time_sim,
  train_size = 500
)

ggsave(
  filename = file.path(output_dir, "Figure_simulation_forecast.eps"),
  plot = p_sim,
  device = "eps",
  width = 8,
  height = 5
)

# ------------------------------------------------------------
# Run for 1000 times
# ------------------------------------------------------------

n_rep <- 1000

mc_results <- do.call(
  rbind,
  lapply(
    1:n_rep,
    estimate_from_simulation,
    use_true_M = FALSE
  )
)

mc_results <- mc_results[
  mc_results$convergence == 0,
]

cat("Monte Carlo means:\n")
print(
  colMeans(
    mc_results[, c("theta", "sigma2", "omega2")],
    na.rm = TRUE
  )
)

cat("Monte Carlo standard deviations:\n")
print(
  apply(
    mc_results[, c("theta", "sigma2", "omega2")],
    2,
    sd,
    na.rm = TRUE
  )
)

# ============================================================
# 2. Parameter Distribution Plot
# ============================================================

dist_df <- data.frame(
  value = c(
    mc_results$theta,
    mc_results$sigma2,
    mc_results$omega2
  ),
  parameter = factor(
    rep(
      c("theta", "sigma2", "omega2"),
      each = nrow(mc_results)
    ),
    levels = c("theta", "sigma2", "omega2"),
    labels = c(
      expression(theta),
      expression(sigma^2),
      expression(omega^2)
    )
  ),
  true_value = c(
    rep(0.15, nrow(mc_results)),
    rep(0.55, nrow(mc_results)),
    rep(0.20, nrow(mc_results))
  )
)

p_dist <- ggplot(
  dist_df,
  aes(x = value)
) +
  geom_density(
    color = "black",
    linewidth = 0.7
  ) +
  geom_vline(
    aes(xintercept = true_value),
    color = "red",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = label_parsed
  ) +
  theme_bw() +
  labs(
    x = NULL,
    y = "Density"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure_parameter_distribution.eps"
  ),
  plot = p_dist,
  device = "eps",
  width = 7,
  height = 4.5
)

# ============================================================
# 3. Unbiasedness Check
# ============================================================

run_unbiased_check <- function(param_name,
                               true_values,
                               base_theta = 0.15,
                               base_sigma2 = 0.55,
                               base_omega2 = 0.20,
                               n_each = 100,
                               n = 1000,
                               train_size = 500) {
  
  out <- list()
  
  for (v in true_values) {
    
    estimates <- numeric(n_each)
    
    for (r in 1:n_each) {
      
      theta_use <- base_theta
      sigma2_use <- base_sigma2
      omega2_use <- base_omega2
      
      if (param_name == "theta") theta_use <- v
      if (param_name == "sigma2") sigma2_use <- v
      if (param_name == "omega2") omega2_use <- v
      
      sim <- simulate_trending_ou_exact(
        n = n,
        theta = theta_use,
        sigma2 = sigma2_use,
        omega2 = omega2_use,
        seed = 10000 + r + round(v * 1000)
      )
      
      Z_train <- sim$Z[1:train_size]
      M_train <- sim$M_t[1:train_size]   # use true M(t)
      
      time_train <- seq_len(train_size)
      h_train <- rep(1, train_size)
      
      est <- estimate_ou_params_sim(
        Z = Z_train,
        M_hat = M_train,
        h_vec = h_train,
        time_index = time_train,
        init_theta = theta_use,
        max_theta = 2,
        max_var = 5
      )
      
      estimates[r] <- est[[param_name]]
    }
    
    out[[as.character(v)]] <- data.frame(
      parameter = param_name,
      true_value = v,
      estimated_value = mean(estimates, na.rm = TRUE)
    )
  }
  
  do.call(rbind, out)
}

theta_grid <- seq(
  0.10,
  0.50,
  length.out = 10
)

sigma2_grid <- seq(
  0.15,
  0.75,
  length.out = 10
)

omega2_grid <- seq(
  0.15,
  0.75,
  length.out = 10
)

unbiased_df <- rbind(
  run_unbiased_check("theta", theta_grid),
  run_unbiased_check("sigma2", sigma2_grid),
  run_unbiased_check("omega2", omega2_grid)
)

unbiased_df$parameter <- factor(
  unbiased_df$parameter,
  levels = c("theta", "sigma2", "omega2"),
  labels = c(
    expression(theta),
    expression(sigma^2),
    expression(omega^2)
  )
)

p_unbiased <- ggplot(
  unbiased_df,
  aes(
    x = true_value,
    y = estimated_value
  )
) +
  geom_point(
    color = "steelblue",
    size = 2
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "red",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = label_parsed
  ) +
  theme_bw() +
  labs(
    x = "True Value",
    y = "Estimated Value"
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure_parameter_unbiasedness.eps"
  ),
  plot = p_unbiased,
  device = "eps",
  width = 7,
  height = 4.5
)

# ============================================================
# 4. Reliability Diagram
# ============================================================

make_reliability_data <- function(monitor_result,
                                  probs = c(0.10, 0.25, 0.50, 0.75, 0.90, 0.95)) {
  
  l <- monitor_result$results
  draws_list <- monitor_result$bootstrap_draws
  
  out <- lapply(probs, function(p) {
    
    covered <- sapply(seq_along(draws_list), function(i) {
      
      draws <- draws_list[[i]]
      
      if (is.null(draws) || length(draws) == 0) {
        return(NA)
      }
      
      if (i > nrow(l) || is.na(l$Z[i])) {
        return(NA)
      }
      
      q_p <- as.numeric(
        quantile(
          draws,
          probs = p,
          na.rm = TRUE,
          type = 8
        )
      )
      
      l$Z[i] <= q_p
    })
    
    data.frame(
      predicted = p,
      observed = mean(covered, na.rm = TRUE),
      n_eval = sum(!is.na(covered))
    )
  })
  
  do.call(rbind, out)
}

reliability_df <- make_reliability_data(
  monitor_result = monitor_result_2024,
  probs = c(0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
)

p_reliability <- ggplot(
  reliability_df,
  aes(x = predicted, y = observed)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_line(
    color = "dodgerblue",
    linewidth = 1
  ) +
  geom_point(
    color = "orange",
    size = 2.5
  ) +
  theme_bw() +
  labs(
    x = "Predicted Probability",
    y = "Observed Frequency"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )

ggsave(
  filename = file.path(output_dir, "Figure_reliability_diagram.eps"),
  plot = p_reliability,
  device = "eps",
  width = 6,
  height = 5
)

# ============================================================
# 5. Model Diagnostics: Standardized Kalman Innovations
# ============================================================

get_ou_standardized_residuals <- function(monitor_result,
                                          use_test_only = TRUE) {
  
  l <- monitor_result$results
  
  idx <- !is.na(l$PredMean) &
    !is.na(l$PredVar) &
    is.finite(l$PredMean) &
    is.finite(l$PredVar) &
    l$PredVar > 0
  
  if (use_test_only && !is.null(monitor_result$train_size_effective)) {
    idx <- idx & l$effective_time > monitor_result$train_size_effective
  }
  
  residuals <- rep(NA_real_, nrow(l))
  
  residuals[idx] <- (l$Z[idx] - l$PredMean[idx]) / sqrt(l$PredVar[idx])
  
  residuals[is.finite(residuals)]
}

resid_2020 <- get_ou_standardized_residuals(
  monitor_result = monitor_result_2020,
  use_test_only = TRUE
)

resid_2024 <- get_ou_standardized_residuals(
  monitor_result = monitor_result_2024,
  use_test_only = TRUE
)

save_diagnostics_plot <- function(residuals,
                                  title_prefix,
                                  filename) {
  
  postscript(
    file = file.path(output_dir, filename),
    width = 8,
    height = 4.5,
    horizontal = FALSE,
    paper = "special"
  )
  
  par(mfrow = c(1, 2))
  
  qqnorm(
    residuals,
    main = paste0(title_prefix, " Q-Q Plot"),
    ylab = "Standardized residuals"
  )
  qqline(
    residuals,
    col = "red"
  )
  
  acf(
    residuals,
    main = paste0(title_prefix, " ACF Plot"),
    lag.max = 25
  )
  
  dev.off()
}

save_diagnostics_plot(
  residuals = resid_2024,
  title_prefix = "2024 U.S. Presidential Election",
  filename = "Figure_2024_model_diagnostics.eps"
)

save_diagnostics_plot(
  residuals = resid_2020,
  title_prefix = "2020 U.S. Presidential Election",
  filename = "Figure_2020_model_diagnostics.eps"
)

cat("Model checking figures saved to:", output_dir, "\n")
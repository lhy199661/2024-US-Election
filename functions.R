# ============================================================
# functions.R
# US Election P2P Betting Strategy
# Trending OU with Additive Noise + Bradley-Terry Model
# Irregular-gap Kalman version
# ============================================================

library(mgcv)
library(nlme)
library(doParallel)
library(parallel)
library(foreach)
library(car)
library(pROC)
library(Metrics)
library(MuMIn)
library(rngtools)
library(doRNG)
library(ggplot2)
library(ResourceSelection)

# ============================================================
# 1. GAM smoother for M(t)
# ============================================================

fit_mu <- function(Z, k = 15, time_index = NULL) { 
  
  if (is.null(time_index)) {
    time_index <- seq_along(Z)
  }
  
  df <- data.frame(
    Z = Z,
    t = time_index
  )
  
  k_valid <- min(k, max(3, nrow(df) - 1))
  
  fit <- gam(
    Z ~ s(t, bs = "cs", k = k_valid),
    data = df
  )
  
  M_hat <- as.numeric(predict(fit))
  M_hat <- ifelse(is.na(M_hat), mean(Z, na.rm = TRUE), M_hat)
  
  return(list(
    fit = fit,
    M_hat = M_hat,
    time_index = time_index
  ))
}

num_derivative <- function(x, time_index = NULL) {
  
  n <- length(x)
  
  if (is.null(time_index)) {
    time_index <- seq_len(n)
  }
  
  d <- numeric(n)
  
  if (n <= 2) {
    
    h <- diff(time_index)
    h <- ifelse(h <= 0, 1, h)
    
    d <- diff(x) / h
    d <- c(d, tail(d, 1))
    
  } else {
    
    h1 <- time_index[2] - time_index[1]
    hn <- time_index[n] - time_index[n - 1]
    
    h1 <- ifelse(h1 <= 0, 1, h1)
    hn <- ifelse(hn <= 0, 1, hn)
    
    d[1] <- (x[2] - x[1]) / h1
    d[n] <- (x[n] - x[n - 1]) / hn
    
    for (i in 2:(n - 1)) {
      h_mid <- time_index[i + 1] - time_index[i - 1]
      h_mid <- ifelse(h_mid <= 0, 1, h_mid)
      d[i] <- (x[i + 1] - x[i - 1]) / h_mid
    }
  }
  
  d <- ifelse(is.na(d) | is.infinite(d), 0, d)
  
  return(d)
}

construct_mu <- function(M_hat, theta, time_index = NULL) {
  
  Mprime <- num_derivative(M_hat, time_index = time_index)
  mu_hat <- M_hat + Mprime / theta
  
  return(mu_hat)
}

construct_mu_next <- function(M_hat,
                              M_next,
                              theta,
                              time_index,
                              time_next) {
  
  h_next <- time_next - tail(time_index, 1)
  h_next <- ifelse(h_next <= 0, 1, h_next)
  
  Mprime_next <- (M_next - tail(M_hat, 1)) / h_next
  mu_next <- M_next + Mprime_next / theta
  
  return(as.numeric(mu_next))
}

# ============================================================
# 2. Kalman likelihood and OU estimation
# ============================================================

log_likelihood_ou_recursive <- function(params,
                                        Z,
                                        M_hat,
                                        Mprime = NULL,
                                        h_vec = NULL,
                                        time_index = NULL) {
  
  theta <- as.numeric(params[1])
  log_sigma2 <- as.numeric(params[2])
  log_omega2 <- as.numeric(params[3])
  
  if (any(is.na(c(theta, log_sigma2, log_omega2))) ||
      any(is.infinite(c(theta, log_sigma2, log_omega2)))) {
    return(1e10)
  }
  
  if (theta < 1e-6 ||
      log_sigma2 < -20 || log_sigma2 > 20 ||
      log_omega2 < -20 || log_omega2 > 20) {
    return(1e10)
  }
  
  theta <- max(theta, 1e-6)
  sigma2 <- max(exp(log_sigma2), 1e-6)
  omega2 <- max(exp(log_omega2), 1e-6)
  
  n <- length(Z)
  
  if (is.null(h_vec)) {
    h_vec <- rep(1, n)
  }
  
  if (is.null(Mprime)) {
    Mprime <- num_derivative(M_hat, time_index = time_index)
  }
  
  mu_hat <- M_hat + Mprime / theta
  
  Q_prev <- mu_hat[1]
  P_prev <- sigma2 / (2 * theta)
  
  ll <- 0
  
  for (t in 1:n) {
    
    if (t == 1) {
      
      Q_pred <- mu_hat[t]
      P_pred <- sigma2 / (2 * theta)
      
    } else {
      
      h_t <- max(h_vec[t], 1e-8)
      a_t <- exp(-theta * h_t)
      q_t <- sigma2 / (2 * theta) * (1 - exp(-2 * theta * h_t))
      
      Q_pred <- a_t * Q_prev + (1 - a_t) * mu_hat[t]
      P_pred <- a_t^2 * P_prev + q_t
    }
    
    P_pred <- max(P_pred, 1e-10)
    S_t <- max(P_pred + omega2, 1e-10)
    
    err <- Z[t] - Q_pred
    
    ll <- ll - 0.5 * (
      log(2 * pi * S_t) +
        err^2 / S_t
    )
    
    K_t <- P_pred / S_t
    
    Q_prev <- Q_pred + K_t * err
    P_prev <- max((1 - K_t) * P_pred, 1e-10)
  }
  
  ll <- ifelse(is.na(ll) | is.infinite(ll), -1e10, ll)
  
  return(-ll)
}

estimate_ou_params <- function(Z,
                               M_hat,
                               h_vec = NULL,
                               time_index = NULL) {
  
  n <- length(Z)
  
  if (is.null(h_vec)) {
    h_vec <- rep(1, n)
  }
  
  if (is.null(time_index)) {
    time_index <- seq_len(n)
  }
  
  Mprime <- num_derivative(M_hat, time_index = time_index)
  
  h_trans <- h_vec[-1]
  h_trans <- pmax(h_trans, 1e-8)
  
  dZ <- diff(Z)
  Xreg <- (M_hat[-n] - Z[-n]) * h_trans
  
  fit0 <- lm(dZ ~ 0 + Xreg)
  
  theta_hat <- as.numeric(coef(fit0)[1])
  if (!is.finite(theta_hat) || theta_hat <= 0) {
    theta_hat <- 0.1
  }
  
  se_theta <- summary(fit0)$coefficients[1, 2]
  if (!is.finite(se_theta) || se_theta <= 0) {
    se_theta <- max(theta_hat / 3, 0.01)
  }
  
  res <- resid(fit0)
  sigma2_hat <- max(var(res), 1e-6)
  if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
    sigma2_hat <- max(var(dZ, na.rm = TRUE), 1e-6)
  }
  
  sigma2_init <- sigma2_hat / 2
  omega2_init <- sigma2_hat / 2
  
  theta_lower <- max(theta_hat - 3 * se_theta, 1e-6)
  theta_upper <- max(theta_hat + 3 * se_theta, theta_lower * 10)
  
  var_se_factor <- 3 / sqrt(2 * (n - 1))
  
  sigma2_lower <- max(sigma2_init - var_se_factor * sigma2_hat, 1e-6)
  sigma2_upper <- sigma2_init + var_se_factor * sigma2_hat
  
  omega2_lower <- max(omega2_init - var_se_factor * sigma2_hat, 1e-6)
  omega2_upper <- omega2_init + var_se_factor * sigma2_hat
  
  init_params <- c(
    theta_hat,
    log(sigma2_init),
    log(omega2_init)
  )
  
  lower_bounds <- c(
    theta_lower,
    log(sigma2_lower),
    log(omega2_lower)
  )
  
  upper_bounds <- c(
    theta_upper,
    log(sigma2_upper),
    log(omega2_upper)
  )
  
  opt <- tryCatch({
    optim(
      par = init_params,
      fn = log_likelihood_ou_recursive,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = 1000),
      Z = Z,
      M_hat = M_hat,
      Mprime = Mprime,
      h_vec = h_vec,
      time_index = time_index
    )
  }, error = function(e) NULL)
  
  if (!is.null(opt)) {
    theta_opt <- opt$par[1]
    sigma2_opt <- exp(opt$par[2])
    omega2_opt <- exp(opt$par[3])
  } else {
    theta_opt <- theta_hat
    sigma2_opt <- sigma2_init
    omega2_opt <- omega2_init
  }
  
  mu_hat <- M_hat + Mprime / theta_opt
  
  return(list(
    theta = theta_opt,
    sigma2 = sigma2_opt,
    omega2 = omega2_opt,
    M_hat = M_hat,
    mu_hat = mu_hat
  ))
}

# ============================================================
# 3. Kalman conditional moments with irregular h
# ============================================================

kalman_filter_state <- function(Z,
                                mu_hat,
                                theta,
                                sigma2,
                                omega2,
                                h_vec = NULL) {
  
  n <- length(Z)
  
  if (is.null(h_vec)) {
    h_vec <- rep(1, n)
  }
  
  theta <- max(theta, 1e-6)
  sigma2 <- max(sigma2, 1e-6)
  omega2 <- max(omega2, 1e-6)
  
  Q_filt <- rep(NA_real_, n)
  P_filt <- rep(NA_real_, n)
  
  Q_prev <- mu_hat[1]
  P_prev <- sigma2 / (2 * theta)
  
  for (t in 1:n) {
    
    if (t == 1) {
      
      Q_pred <- mu_hat[t]
      P_pred <- sigma2 / (2 * theta)
      
    } else {
      
      h_t <- max(h_vec[t], 1e-8)
      a_t <- exp(-theta * h_t)
      q_t <- sigma2 / (2 * theta) * (1 - exp(-2 * theta * h_t))
      
      Q_pred <- a_t * Q_prev + (1 - a_t) * mu_hat[t]
      P_pred <- a_t^2 * P_prev + q_t
    }
    
    P_pred <- max(P_pred, 1e-10)
    S_t <- max(P_pred + omega2, 1e-10)
    
    err <- Z[t] - Q_pred
    K_t <- P_pred / S_t
    
    Q_prev <- Q_pred + K_t * err
    P_prev <- max((1 - K_t) * P_pred, 1e-10)
    
    Q_filt[t] <- Q_prev
    P_filt[t] <- P_prev
  }
  
  return(list(
    Q_filt = Q_filt,
    P_filt = P_filt,
    Q_last = Q_prev,
    P_last = P_prev
  ))
}

kalman_conditional_next <- function(Z_past,
                                    z_current,
                                    mu_past,
                                    mu_current,
                                    mu_next,
                                    theta,
                                    sigma2,
                                    omega2,
                                    h_past = NULL,
                                    h_current = 1,
                                    h_next = 1) {
  
  theta <- max(theta, 1e-6)
  sigma2 <- max(sigma2, 1e-6)
  omega2 <- max(omega2, 1e-6)
  
  if (length(Z_past) > 0) {
    
    if (is.null(h_past)) {
      h_past <- rep(1, length(Z_past))
    }
    
    kf <- kalman_filter_state(
      Z = Z_past,
      mu_hat = mu_past,
      theta = theta,
      sigma2 = sigma2,
      omega2 = omega2,
      h_vec = h_past
    )
    
    Q_prev <- kf$Q_last
    P_prev <- kf$P_last
    
    h_current <- max(h_current, 1e-8)
    a_current <- exp(-theta * h_current)
    q_current <- sigma2 / (2 * theta) * (1 - exp(-2 * theta * h_current))
    
    Q_pred_t <- a_current * Q_prev + (1 - a_current) * mu_current
    P_pred_t <- a_current^2 * P_prev + q_current
    
  } else {
    
    Q_pred_t <- mu_current
    P_pred_t <- sigma2 / (2 * theta)
  }
  
  P_pred_t <- max(P_pred_t, 1e-10)
  S_t <- max(P_pred_t + omega2, 1e-10)
  
  K_t <- P_pred_t / S_t
  innovation <- z_current - Q_pred_t
  
  Q_filt_t <- Q_pred_t + K_t * innovation
  P_filt_t <- max((1 - K_t) * P_pred_t, 1e-10)
  
  h_next <- max(h_next, 1e-8)
  a_next <- exp(-theta * h_next)
  q_next <- sigma2 / (2 * theta) * (1 - exp(-2 * theta * h_next))
  
  Q_pred_next <- a_next * Q_filt_t + (1 - a_next) * mu_next
  P_pred_next <- a_next^2 * P_filt_t + q_next
  
  S_next <- max(P_pred_next + omega2, 1e-10)
  
  return(list(
    mean = as.numeric(Q_pred_next),
    variance = as.numeric(S_next),
    sd = sqrt(as.numeric(S_next)),
    Q_current_filtered = as.numeric(Q_filt_t),
    P_current_filtered = as.numeric(P_filt_t),
    innovation_variance = as.numeric(S_t),
    kalman_gain = as.numeric(K_t)
  ))
}

simulate_bootstrap_by_kalman <- function(Z_init,
                                         mu_hat,
                                         theta,
                                         sigma2,
                                         omega2,
                                         h_vec = NULL) {
  
  n <- length(mu_hat)
  
  if (is.null(h_vec)) {
    h_vec <- rep(1, n)
  }
  
  Zb <- numeric(n)
  Zb[1] <- Z_init
  
  if (n == 1) {
    return(Zb)
  }
  
  for (s in 1:(n - 1)) {
    
    if (s == 1) {
      Z_past <- numeric(0)
      mu_past <- numeric(0)
      h_past <- numeric(0)
    } else {
      Z_past <- Zb[1:(s - 1)]
      mu_past <- mu_hat[1:(s - 1)]
      h_past <- h_vec[1:(s - 1)]
    }
    
    cond <- kalman_conditional_next(
      Z_past = Z_past,
      z_current = Zb[s],
      mu_past = mu_past,
      mu_current = mu_hat[s],
      mu_next = mu_hat[s + 1],
      theta = theta,
      sigma2 = sigma2,
      omega2 = omega2,
      h_past = h_past,
      h_current = h_vec[s],
      h_next = h_vec[s + 1]
    )
    
    Zb[s + 1] <- rnorm(
      1,
      mean = cond$mean,
      sd = cond$sd
    )
  }
  
  return(Zb)
}

# ============================================================
# 4. Effective series and live bootstrap forecast
# train_size is defined on the original full series
# ============================================================

make_effective_series <- function(Z,
                                  skip_flat = TRUE,
                                  time_index = NULL) {
  
  if (is.null(time_index)) {
    time_index <- seq_along(Z)
  }
  
  if (length(time_index) != length(Z)) {
    stop("time_index must have the same length as Z.")
  }
  
  if (!skip_flat) {
    idx <- seq_along(Z)
  } else {
    idx <- c(1, which(diff(Z) != 0) + 1)
  }
  
  effective_time_index <- time_index[idx]
  h_vec <- c(1, diff(effective_time_index))
  
  list(
    Z = Z[idx],
    original_index = idx,
    time_index = effective_time_index,
    h_vec = h_vec
  )
}

convert_raw_train_size_to_effective <- function(original_index,
                                                train_size_raw) {
  
  train_size_eff <- sum(original_index <= train_size_raw)
  
  if (train_size_eff < 5) {
    stop("Effective training sample is too small after removing flat periods.")
  }
  
  return(train_size_eff)
}

bootstrap_forecast_live <- function(
    Z,
    B = 200,
    alpha = 0.05,
    train_size = 500,
    k_gam = 15,
    cores_limit = NULL,
    verbose = TRUE,
    skip_flat = TRUE,
    time_index = NULL
) {
  
  eff <- make_effective_series(
    Z,
    skip_flat = skip_flat,
    time_index = time_index
  )
  
  Z_eff <- eff$Z
  original_index <- eff$original_index
  ou_time_index <- eff$time_index
  h_eff <- eff$h_vec
  
  n <- length(Z_eff)
  train_size_raw <- train_size
  
  if (skip_flat) {
    train_size_eff <- convert_raw_train_size_to_effective(
      original_index = original_index,
      train_size_raw = train_size_raw
    )
  } else {
    train_size_eff <- train_size_raw
  }
  
  if (train_size_eff >= n) {
    stop("Effective train_size must be smaller than the number of effective observations.")
  }
  
  if (verbose) {
    message("Raw training index: ", train_size_raw)
    message("Effective training size: ", train_size_eff)
    message("Effective observations: ", n)
  }
  
  PredMean <- rep(NA_real_, n)
  PredVar <- rep(NA_real_, n)
  Upper <- rep(NA_real_, n)
  ProbDecrease <- rep(NA_real_, n)
  Signal <- rep(FALSE, n)
  Viol <- rep(NA, n)
  
  M_hat_store <- rep(NA_real_, n)
  Mu_hat_store <- rep(NA_real_, n)
  
  theta_store <- rep(NA_real_, n)
  sigma2_store <- rep(NA_real_, n)
  omega2_store <- rep(NA_real_, n)
  
  bootstrap_draws_list <- vector("list", n)
  bootstrap_params_list <- vector("list", n)
  
  # Initial training trend for plotting
  Z_train0 <- Z_eff[1:train_size_eff]
  time_train0 <- ou_time_index[1:train_size_eff]
  h_train0 <- h_eff[1:train_size_eff]
  
  gam_train0 <- fit_mu(
    Z_train0,
    k = k_gam,
    time_index = time_train0
  )
  
  M_train0 <- gam_train0$M_hat
  
  params_train0 <- estimate_ou_params(
    Z = Z_train0,
    M_hat = M_train0,
    h_vec = h_train0,
    time_index = time_train0
  )
  
  Mu_train0 <- params_train0$mu_hat
  
  M_hat_store[1:train_size_eff] <- M_train0
  Mu_hat_store[1:train_size_eff] <- Mu_train0
  
  theta_store[1:train_size_eff] <- params_train0$theta
  sigma2_store[1:train_size_eff] <- params_train0$sigma2
  omega2_store[1:train_size_eff] <- params_train0$omega2
  
  # Parallel setup
  if (is.null(cores_limit)) {
    cores <- max(1, parallel::detectCores() - 1)
  } else {
    cores <- max(1, cores_limit)
  }
  
  cl <- parallel::makeCluster(cores)
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }, add = TRUE)
  
  doParallel::registerDoParallel(cl)
  
  parallel::clusterExport(
    cl,
    c(
      "fit_mu",
      "num_derivative",
      "construct_mu",
      "construct_mu_next",
      "estimate_ou_params",
      "log_likelihood_ou_recursive",
      "kalman_filter_state",
      "kalman_conditional_next",
      "simulate_bootstrap_by_kalman"
    ),
    envir = environment()
  )
  
  parallel::clusterEvalQ(cl, {
    library(mgcv)
  })
  
  if (verbose) {
    message("Using ", cores, " cores.")
  }
  
  for (t in train_size_eff:(n - 1)) {
    
    if (verbose && t %% 50 == 0) {
      message("Forecasting effective time t = ", t, " / ", n - 1)
    }
    
    Ztr <- Z_eff[1:t]
    time_tr <- ou_time_index[1:t]
    h_tr <- h_eff[1:t]
    
    time_next <- ou_time_index[t + 1]
    h_next <- h_eff[t + 1]
    
    # Step 1: fit M(t) and estimate parameters
    gam_result <- fit_mu(
      Ztr,
      k = k_gam,
      time_index = time_tr
    )
    
    M0 <- gam_result$M_hat
    
    params <- estimate_ou_params(
      Z = Ztr,
      M_hat = M0,
      h_vec = h_tr,
      time_index = time_tr
    )
    
    mu0 <- params$mu_hat
    
    M_next <- as.numeric(
      predict(
        gam_result$fit,
        newdata = data.frame(t = time_next)
      )
    )
    
    mu_next <- construct_mu_next(
      M_hat = M0,
      M_next = M_next,
      theta = params$theta,
      time_index = time_tr,
      time_next = time_next
    )
    
    # Step 2: Kalman conditional m_t^(0)(Z_t), v_t^(0)
    if (t == 1) {
      Z_past <- numeric(0)
      mu_past <- numeric(0)
      h_past <- numeric(0)
    } else {
      Z_past <- Ztr[1:(t - 1)]
      mu_past <- mu0[1:(t - 1)]
      h_past <- h_tr[1:(t - 1)]
    }
    
    cond0 <- kalman_conditional_next(
      Z_past = Z_past,
      z_current = Ztr[t],
      mu_past = mu_past,
      mu_current = mu0[t],
      mu_next = mu_next,
      theta = params$theta,
      sigma2 = params$sigma2,
      omega2 = params$omega2,
      h_past = h_past,
      h_current = h_tr[t],
      h_next = h_next
    )
    
    PredMean[t + 1] <- cond0$mean
    PredVar[t + 1] <- cond0$variance
    
    M_hat_store[1:t] <- M0
    Mu_hat_store[1:t] <- mu0
    M_hat_store[t + 1] <- M_next
    Mu_hat_store[t + 1] <- mu_next
    
    theta_store[t + 1] <- params$theta
    sigma2_store[t + 1] <- params$sigma2
    omega2_store[t + 1] <- params$omega2
    
    # Step 3: bootstrap
    boot_mat <- foreach::foreach(
      b = 1:B,
      .combine = rbind,
      .errorhandling = "remove",
      .packages = c("mgcv"),
      .options.RNG = 100000 + t
    ) %dorng% {
      
      Zb <- simulate_bootstrap_by_kalman(
        Z_init = Ztr[1],
        mu_hat = mu0,
        theta = params$theta,
        sigma2 = params$sigma2,
        omega2 = params$omega2,
        h_vec = h_tr
      )
      
      gam_b <- fit_mu(
        Zb,
        k = k_gam,
        time_index = time_tr
      )
      
      Mb <- gam_b$M_hat
      
      pb <- estimate_ou_params(
        Z = Zb,
        M_hat = Mb,
        h_vec = h_tr,
        time_index = time_tr
      )
      
      mub <- pb$mu_hat
      
      Mb_next <- as.numeric(
        predict(
          gam_b$fit,
          newdata = data.frame(t = time_next)
        )
      )
      
      mub_next <- construct_mu_next(
        M_hat = Mb,
        M_next = Mb_next,
        theta = pb$theta,
        time_index = time_tr,
        time_next = time_next
      )
      
      if (t == 1) {
        Z_past_actual <- numeric(0)
        mu_past_b <- numeric(0)
        h_past_actual <- numeric(0)
      } else {
        Z_past_actual <- Ztr[1:(t - 1)]
        mu_past_b <- mub[1:(t - 1)]
        h_past_actual <- h_tr[1:(t - 1)]
      }
      
      cond_b <- kalman_conditional_next(
        Z_past = Z_past_actual,
        z_current = Ztr[t],
        mu_past = mu_past_b,
        mu_current = mub[t],
        mu_next = mub_next,
        theta = pb$theta,
        sigma2 = pb$sigma2,
        omega2 = pb$omega2,
        h_past = h_past_actual,
        h_current = h_tr[t],
        h_next = h_next
      )
      
      c(
        draw = rnorm(1, mean = cond_b$mean, sd = cond_b$sd),
        theta = pb$theta,
        sigma2 = pb$sigma2,
        omega2 = pb$omega2
      )
    }
    
    if (is.null(boot_mat) || length(boot_mat) == 0) {
      stop(
        paste0(
          "No valid bootstrap results at effective time t = ",
          t,
          ". Check GAM fitting, MLE, or parallel workers."
        )
      )
    }
    
    boot_mat <- as.data.frame(boot_mat)
    
    boot_mat$draw <- as.numeric(boot_mat$draw)
    boot_mat$theta <- as.numeric(boot_mat$theta)
    boot_mat$sigma2 <- as.numeric(boot_mat$sigma2)
    boot_mat$omega2 <- as.numeric(boot_mat$omega2)
    
    draw_ok <- is.finite(boot_mat$draw)
    draws <- boot_mat$draw[draw_ok]
    
    bootstrap_draws_list[[t + 1]] <- draws
    
    param_ok <- is.finite(boot_mat$theta) &
      is.finite(boot_mat$sigma2) &
      is.finite(boot_mat$omega2)
    
    bootstrap_params_list[[t + 1]] <- boot_mat[
      param_ok,
      c("theta", "sigma2", "omega2"),
      drop = FALSE
    ]
    
    if (length(draws) == 0) {
      stop(
        paste0(
          "No valid bootstrap draws at effective time t = ",
          t,
          ". Check GAM fitting, MLE, or parallel workers."
        )
      )
    }
    
    if (length(draws) < B && verbose) {
      warning(
        paste0(
          "Only ",
          length(draws),
          " valid bootstrap draws out of B = ",
          B,
          " at effective time t = ",
          t,
          "."
        )
      )
    }
    
    Upper[t + 1] <- as.numeric(
      quantile(
        draws,
        probs = 1 - alpha,
        na.rm = TRUE,
        type = 8
      )
    )
    
    # Paper signal:
    # Z_t > U_{t+1}^{alpha}, equivalently P(Z_{t+1} <= Z_t) >= 1 - alpha.
    prob_decrease <- mean(draws <= Z_eff[t], na.rm = TRUE)
    ProbDecrease[t] <- prob_decrease
    
    Signal[t] <- isTRUE(prob_decrease >= 1 - alpha)
    
    # Coverage diagnostic:
    # realised Z_{t+1} outside U_{t+1}^{alpha}
    Viol[t + 1] <- Z_eff[t + 1] > Upper[t + 1]
  }
  
  eval_idx <- which(!is.na(Viol))
  
  coverage <- if (length(eval_idx) > 0) {
    mean(!Viol[eval_idx], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  results_df <- data.frame(
    effective_time = seq_len(n),
    original_index = original_index,
    time_index = ou_time_index,
    h = h_eff,
    Z = Z_eff,
    M_hat = M_hat_store,
    Mu_hat = Mu_hat_store,
    PredMean = PredMean,
    PredVar = PredVar,
    Upper = Upper,
    ProbDecrease = ProbDecrease,
    Signal = Signal,
    Viol = Viol,
    theta = theta_store,
    sigma2 = sigma2_store,
    omega2 = omega2_store
  )
  
  return(list(
    results = results_df,
    bootstrap_draws = bootstrap_draws_list,
    bootstrap_params = bootstrap_params_list,
    coverage = coverage,
    skip_flat = skip_flat,
    train_size_raw = train_size_raw,
    train_size_effective = train_size_eff
  ))
}

bootstrap_forecast_calibrated <- function(
    Z,
    B = 200,
    alpha = 0.05,
    train_size = 500,
    k_gam = 15,
    cores_limit = NULL,
    verbose = TRUE,
    do_calibrate = FALSE,
    skip_flat = TRUE,
    time_index = NULL
) {
  
  if (isTRUE(do_calibrate)) {
    warning("do_calibrate is ignored to avoid in-sample calibration leakage.")
  }
  
  bootstrap_forecast_live(
    Z = Z,
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

bootstrap_forecast <- bootstrap_forecast_live

# Parameter Output
get_bootstrap_param_ci <- function(monitor_result,
                                   row = NULL,
                                   level = 0.95,
                                   method = c("percentile", "normal")) {
  
  method <- match.arg(method)
  
  if (is.null(monitor_result$bootstrap_params)) {
    stop("monitor_result does not contain bootstrap_params.")
  }
  
  if (is.null(row)) {
    nonempty <- which(
      vapply(
        monitor_result$bootstrap_params,
        function(x) {
          !is.null(x) && nrow(x) > 0
        },
        logical(1)
      )
    )
    
    if (length(nonempty) == 0) {
      stop("No non-empty bootstrap parameter estimates found.")
    }
    
    row <- tail(nonempty, 1)
  }
  
  boot_params <- monitor_result$bootstrap_params[[row]]
  
  if (is.null(boot_params) || nrow(boot_params) == 0) {
    stop("No bootstrap parameter estimates at the selected row.")
  }
  
  point_est <- monitor_result$results[
    row,
    c("theta", "sigma2", "omega2"),
    drop = FALSE
  ]
  
  alpha <- 1 - level
  
  boot_se <- apply(
    boot_params,
    2,
    sd,
    na.rm = TRUE
  )
  
  if (method == "percentile") {
    
    ci <- apply(
      boot_params,
      2,
      quantile,
      probs = c(alpha / 2, 1 - alpha / 2),
      na.rm = TRUE,
      type = 8
    )
    
    lower <- ci[1, ]
    upper <- ci[2, ]
    
  } else {
    
    z <- qnorm(1 - alpha / 2)
    
    lower <- as.numeric(point_est) - z * boot_se
    upper <- as.numeric(point_est) + z * boot_se
  }
  
  data.frame(
    parameter = c("theta", "sigma2", "omega2"),
    estimate = as.numeric(point_est),
    bootstrap_se = as.numeric(boot_se),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    level = level,
    method = method,
    row = row,
    B_used = nrow(boot_params)
  )
}

# ============================================================
# 5. Signal extraction
# ============================================================

find_first_changes <- function(time_vec, Z_vec) {
  
  if (length(Z_vec) == 0) {
    return(data.frame(time = numeric(0), Z = numeric(0)))
  }
  
  change_idx <- c(1, which(diff(Z_vec) != 0) + 1)
  
  data.frame(
    time = time_vec[change_idx],
    Z = Z_vec[change_idx]
  )
}

get_signal_points_from_results <- function(monitor_result) {
  
  l <- monitor_result$results
  
  l_signal <- l[
    !is.na(l$Signal) & l$Signal == TRUE,
  ]
  
  if (nrow(l_signal) == 0) {
    return(numeric(0))
  }
  
  change_points <- find_first_changes(
    time_vec = l_signal$original_index,
    Z_vec = l_signal$Z
  )
  
  na.omit(change_points$time)
}

get_signal_points <- get_signal_points_from_results

# ============================================================
# 6. Change-point utilities
# ============================================================

find_next_change <- function(x, pos) {
  
  if (is.na(pos) || pos >= length(x)) {
    return(NA)
  }
  
  start_val <- x[pos]
  idx <- which(x[(pos + 1):length(x)] != start_val)[1]
  
  if (is.na(idx)) {
    return(NA)
  }
  
  pos + idx
}

get_last_change_index <- function(x, targets) {
  
  if (length(targets) == 0) {
    return(integer(0))
  }
  
  vapply(targets, function(i) {
    
    if (is.na(i) || i <= 1) {
      return(NA_integer_)
    }
    
    if (i > length(x)) {
      return(NA_integer_)
    }
    
    prev_idx <- which(x[1:(i - 1)] != x[i])
    
    if (length(prev_idx) == 0) {
      return(NA_integer_)
    }
    
    as.integer(max(prev_idx))
    
  }, integer(1))
}

# ============================================================
# 7. 2020 index mapping
# ============================================================

convert_2020_index_to_full <- function(idx,
                                       covid_start = 11594,
                                       covid_gap = 1079) {
  
  ifelse(
    idx < covid_start,
    idx,
    idx + covid_gap
  )
}

# ============================================================
# 8. Bradley-Terry data and model
# ============================================================

build_bt_dataset <- function(candidate1,
                             candidate2,
                             filtered,
                             Y,
                             index_map = NULL) {
  
  empty_bt_data <- data.frame(
    Candidate1 = integer(0),
    Candidate2 = integer(0),
    diff_Delta = numeric(0),
    diff_signal = numeric(0),
    C1_trade = numeric(0),
    C1_signal = numeric(0),
    C2_trade = numeric(0),
    C2_signal = numeric(0),
    signal_index = integer(0),
    trade_index = integer(0)
  )
  
  if (length(filtered) == 0) {
    return(empty_bt_data)
  }
  
  filtered1 <- vapply(filtered, function(p) {
    find_next_change(Y, p)
  }, numeric(1))
  
  valid <- !is.na(filtered) & !is.na(filtered1)
  
  filtered <- filtered[valid]
  filtered1 <- filtered1[valid]
  
  if (length(filtered) == 0) {
    return(empty_bt_data)
  }
  
  if (is.null(index_map)) {
    candidate_signal_idx <- filtered
    candidate_trade_idx <- filtered1
  } else {
    candidate_signal_idx <- index_map(filtered)
    candidate_trade_idx <- index_map(filtered1)
  }
  
  in_range <- !is.na(candidate_signal_idx) &
    !is.na(candidate_trade_idx) &
    candidate_signal_idx >= 1 &
    candidate_trade_idx >= 1 &
    candidate_signal_idx <= length(candidate1) &
    candidate_signal_idx <= length(candidate2) &
    candidate_trade_idx <= length(candidate1) &
    candidate_trade_idx <= length(candidate2)
  
  filtered <- filtered[in_range]
  filtered1 <- filtered1[in_range]
  candidate_signal_idx <- candidate_signal_idx[in_range]
  candidate_trade_idx <- candidate_trade_idx[in_range]
  
  if (length(filtered) == 0) {
    return(empty_bt_data)
  }
  
  last_change_indices1 <- get_last_change_index(
    candidate1,
    candidate_signal_idx
  )
  
  last_change_indices2 <- get_last_change_index(
    candidate2,
    candidate_signal_idx
  )
  
  has_previous_change <- !is.na(last_change_indices1) &
    !is.na(last_change_indices2)
  
  filtered <- filtered[has_previous_change]
  filtered1 <- filtered1[has_previous_change]
  candidate_signal_idx <- candidate_signal_idx[has_previous_change]
  candidate_trade_idx <- candidate_trade_idx[has_previous_change]
  last_change_indices1 <- last_change_indices1[has_previous_change]
  last_change_indices2 <- last_change_indices2[has_previous_change]
  
  if (length(filtered) == 0) {
    return(empty_bt_data)
  }
  
  C1_trade <- candidate1[candidate_trade_idx]
  C1_signal <- candidate1[candidate_signal_idx]
  C1_prev <- candidate1[last_change_indices1]
  
  C2_trade <- candidate2[candidate_trade_idx]
  C2_signal <- candidate2[candidate_signal_idx]
  C2_prev <- candidate2[last_change_indices2]
  
  Delta1 <- C1_signal - C1_prev
  Delta2 <- C2_signal - C2_prev
  
  utility_1 <- 1 / C1_trade - 1 / C1_signal
  utility_2 <- 1 / C2_trade - 1 / C2_signal
  
  label <- as.integer(utility_1 > utility_2)
  
  bt_data <- data.frame(
    Candidate1 = label,
    Candidate2 = 1 - label,
    diff_Delta = Delta1 - Delta2,
    diff_signal = C1_signal - C2_signal,
    C1_trade = C1_trade,
    C1_signal = C1_signal,
    C2_trade = C2_trade,
    C2_signal = C2_signal,
    signal_index = filtered,
    trade_index = filtered1
  )
  
  bt_data <- bt_data[complete.cases(bt_data), ]
  
  return(bt_data)
}

fit_bt_model <- function(bt_data) {
  
  model <- glm(
    Candidate1 ~ 0 + diff_Delta + diff_signal,
    data = bt_data,
    family = binomial()
  )
  
  probs <- predict(model, type = "response")
  
  return(list(
    model = model,
    probabilities = probs
  ))
}

# ============================================================
# 9. Return and strategy evaluation
# ============================================================

calc_returns_metrics <- function(o_buy,
                                 o_sell,
                                 method = c("absolute", "probability", "odds")) {
  
  method <- match.arg(method)
  
  if (length(o_buy) != length(o_sell)) {
    stop("Lengths must match")
  }
  
  n <- length(o_buy)
  
  if (method == "absolute") {
    R <- o_sell - o_buy
  } else if (method == "probability") {
    R <- 1 / o_sell - 1 / o_buy
  } else {
    R <- (o_sell - o_buy) / o_buy
  }
  
  ER <- mean(R, na.rm = TRUE)
  sdR <- sd(R, na.rm = TRUE)
  SR <- ifelse(sdR > 0, ER / sdR, NA)
  
  list(
    method = method,
    returns = R,
    mean_return = ER,
    sd_return = sdR,
    sharpe = SR,
    n = n
  )
}

evaluate_strategy <- function(bt_pred,
                              bt_test,
                              return_method = "absolute") {
  
  selected_value <- ifelse(
    bt_pred == 1,
    1 / bt_test$C1_trade,
    1 / bt_test$C2_trade
  )
  
  selected_baseline <- ifelse(
    bt_pred == 1,
    1 / bt_test$C1_signal,
    1 / bt_test$C2_signal
  )
  
  profit <- selected_value - selected_baseline
  
  metrics <- calc_returns_metrics(
    o_buy = selected_baseline,
    o_sell = selected_value,
    method = return_method
  )
  
  list(
    profit = profit,
    profit_sum = sum(profit, na.rm = TRUE),
    selected_value = selected_value,
    selected_baseline = selected_baseline,
    metrics = metrics
  )
}

# ============================================================
# 10. Simulation
# ============================================================

simulate_trending_ou <- function(n,
                                 theta,
                                 sigma2,
                                 omega2,
                                 seed) {
  
  set.seed(seed)
  
  h <- 1
  t_vec <- 1:n
  
  mu_t <- 5 +
    0.001 * t_vec +
    0.5 * sin(2 * pi * t_vec / 200)
  
  theta <- max(theta, 1e-6)
  sigma2 <- max(sigma2, 1e-6)
  omega2 <- max(omega2, 1e-6)
  
  Q <- numeric(n)
  Z <- numeric(n)
  
  Q[1] <- mu_t[1]
  
  sigma <- sqrt(sigma2)
  omega <- sqrt(omega2)
  
  Z[1] <- Q[1] + rnorm(1, 0, omega)
  
  for (t in 2:n) {
    dW <- rnorm(1, 0, sqrt(h))
    Q[t] <- Q[t - 1] +
      theta * (mu_t[t - 1] - Q[t - 1]) * h +
      sigma * dW
    
    Z[t] <- Q[t] + rnorm(1, 0, omega)
  }
  
  list(
    Z = Z,
    Q = Q,
    mu_t = mu_t,
    true_params = list(
      theta = theta,
      sigma2 = sigma2,
      omega2 = omega2
    )
  )
}

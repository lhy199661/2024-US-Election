#Uploading the packages
library(mgcv)
library(nlme)
library(doParallel)
library(car)
library(pROC)
library(Metrics)
library(MuMIn)
library(rngtools)
library(doRNG)
library(ggplot2)
library(ResourceSelection)

# When to place bet.
# ---------- Fit smoother ----------
# We use gam to get the smoothing M(t)
fit_mu <- function(Z, k = 15) { 
  df <- data.frame(Z=Z, t=seq_along(Z))
  k_valid <- min(k, max(3, nrow(df)-1))
  
  # Using linear method,
  fit <- lm(Z ~ t, data=df)
  mu_hat <- as.numeric(predict(fit))
  mu_hat <- ifelse(is.na(mu_hat), mean(Z, na.rm=T), mu_hat)
  
  # Return the outcome
  return(list(fit=fit, mu_hat=mu_hat))
}

# ---------- Kalman filtering ---------- 
log_likelihood_ou_recursive <- function(params, Z, mu_hat) {
  theta <- as.numeric(params[1])
  log_sigma2 <- as.numeric(params[2])
  log_omega2 <- as.numeric(params[3])
  
  if (any(is.na(c(theta, log_sigma2, log_omega2))) || any(is.infinite(c(theta, log_sigma2, log_omega2)))) {
    return(1e10)
  }
  
  if (theta < 1e-6 || log_sigma2 < -20 || log_sigma2 > 20 || 
      log_omega2 < -20 || log_omega2 > 20) {
    return(1e10)
  }
  
  theta <- max(theta, 1e-6)
  sigma2 <- max(exp(log_sigma2), 1e-6)
  omega2 <- max(exp(log_omega2), 1e-6)
  
  n <- length(Z)
  h <- 1
  a <- exp(-theta * h)
  
  Q_prev <- mu_hat[1]
  P_prev <- sigma2 / (2*theta)
  
  ll <- 0
  
  # Kalman filtering,
  for (t in 1:n) {
    
    if (t == 1) {
      
      Q_pred <- mu_hat[t]
      P_pred <- sigma2/(2*theta)
      
    } else {
      
      Q_pred <- a * Q_prev + mu_hat[t]*(1 - a)
      
      P_pred <- a^2 * P_prev + 
        (sigma2/(2*theta))*(1 - exp(-2*theta*h))
    }
    
    P_pred <- max(P_pred, 1e-10)
    
    S_t <- P_pred + omega2
    S_t <- max(S_t, 1e-10)
    
    err <- Z[t] - Q_pred
    
    ll <- ll - 0.5 * (
      log(2*pi*S_t) + err^2 / S_t
    )
    
    K_t <- P_pred / S_t
    
    Q_prev <- Q_pred + K_t * err
    
    P_prev <- (1 - K_t) * P_pred
    P_prev <- max(P_prev, 1e-10)
  }
  
  # Calculate negative loglikelihood,
  ll <- ifelse(is.na(ll) | is.infinite(ll), -1e10, ll)
  
  return(-ll)
}

# ---------- Buid the derivative of M(t) ----------
num_derivative <- function(x, h=1) {
  
  n <- length(x)
  d <- numeric(n)
  
  if (n <= 2) {
    
    d <- diff(x)/h
    d <- c(d, tail(d,1))
    
  } else {
    
    d[1] <- (x[2]-x[1])/h
    d[n] <- (x[n]-x[n-1])/h
    
    d[2:(n-1)] <- (x[3:n]-x[1:(n-2)])/(2*h)
  }
  
  d <- ifelse(is.na(d) | is.infinite(d), 0, d)
  
  return(d)
}

# ---------- Estimate parameters ----------
estimate_ou_params <- function(Z, m_hat) {
  
  n <- length(Z)
  h <- 1
  
  mprime <- num_derivative(m_hat, h)
  
  # Step 1: ΔZ_t = θ (μ_t - Z_t) Δt + ε_t
  dZ <- diff(Z)
  Xreg <- (m_hat[-n] - Z[-n]) * h
  
  fit0 <- lm(dZ ~ 0 + Xreg)
  
  theta_hat <- as.numeric(coef(fit0)[1])
  theta_hat <- max(theta_hat, 1e-6)
  
  se_theta <- summary(fit0)$coefficients[1,2]
  se_theta <- ifelse(is.na(se_theta), 0.01, se_theta)
  
  # Step 2
  res <- resid(fit0)
  
  sigma2_hat <- var(res)
  sigma2_hat <- max(sigma2_hat, 1e-6)
  
  sigma2_init <- sigma2_hat / 2
  omega2_init <- sigma2_hat / 2
  
  # Step 3 boundary
  theta_lower <- max(theta_hat - 3*se_theta, 1e-6)
  theta_upper <- theta_hat + 3*se_theta
  
  var_se_factor <- 3 / sqrt(2*(n-1))
  
  sigma2_lower <- max(
    sigma2_init - var_se_factor*sigma2_hat,
    1e-6
  )
  
  sigma2_upper <- sigma2_init + var_se_factor*sigma2_hat
  
  omega2_lower <- max(
    omega2_init - var_se_factor*sigma2_hat,
    1e-6
  )
  
  omega2_upper <- omega2_init + var_se_factor*sigma2_hat
  
  log_sigma2_init <- log(sigma2_init)
  log_omega2_init <- log(omega2_init)
  
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
  
  init_params <- c(
    theta_hat,
    log_sigma2_init,
    log_omega2_init
  )
  
  # Step 4 optimization
  opt <- tryCatch({
    
    optim(
      par = init_params,
      fn = log_likelihood_ou_recursive,
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = 1000),
      Z = Z,
      mu_hat = m_hat
    )
    
  }, error=function(e) NULL)
  
  if (!is.null(opt)) {
    
    theta_opt <- opt$par[1]
    sigma2_opt <- exp(opt$par[2])
    omega2_opt <- exp(opt$par[3])
    
  } else {
    
    theta_opt <- theta_hat
    sigma2_opt <- sigma2_init
    omega2_opt <- omega2_init
  }
  
  # Step 5
  mu_final <- m_hat + mprime/theta_opt
  
  return(list(
    theta = theta_opt,
    sigma2 = sigma2_opt,
    omega2 = omega2_opt,
    mu_hat = mu_final
  ))
}

# ---------- Parameters limitation ----------
clamp <- function(x, min_val, max_val) {
  max(min(x, max_val), min_val)
}

# ---------- Parameters coefficients ----------
cond_coeffs <- function(theta, sigma2, omega2) {
  
  theta <- max(as.numeric(theta), 1e-6)
  sigma2 <- max(as.numeric(sigma2), 1e-6)
  omega2 <- max(as.numeric(omega2), 1e-6)
  
  h <- 1
  
  a <- exp(-theta*h)
  
  P_inf <- max(sigma2/(2*theta), 1e-10)
  
  S_inf <- P_inf + omega2
  S_inf <- max(S_inf, 1e-10)
  
  c <- (P_inf * a) / S_inf
  
  v <- S_inf - (P_inf^2 * a^2) / S_inf
  v <- max(v, 1e-10)
  
  return(list(c=c, v=v, a=a))
}

# ---------- Compute coverage ----------
compute_coverage_at_p <- function(bootstrap_draws_list, Z, p) {
  
  covered <- sapply(seq_along(bootstrap_draws_list), function(t1) {
    
    dr <- bootstrap_draws_list[[t1]]
    
    if (all(is.na(dr)) || is.na(Z[t1]) || length(dr) < 5) { 
      return(NA)
    }
    
    q <- quantile(dr, p, na.rm=T, type=8)
    
    Z[t1] <= q
  })
  
  return(mean(covered, na.rm=T))
}

# ---------- 95% one-side prediction interval ----------
find_calibrated_p <- function(bootstrap_draws_list,
                              Z,
                              target=0.95) {
  
  lower <- 0.5
  upper <- 0.999
  tol <- 1e-3
  max_iter <- 20 
  
  for (i in 1:max_iter) {
    
    mid <- (lower + upper) / 2
    
    cov_mid <- compute_coverage_at_p(
      bootstrap_draws_list,
      Z,
      mid
    )
    
    if (is.na(cov_mid)) break
    
    if (abs(cov_mid - target) < tol) break
    
    if (cov_mid < target) {
      lower <- mid
    } else {
      upper <- mid
    }
  }
  
  pstar <- (lower + upper) / 2
  pstar <- max(min(pstar, 0.999), 0.5)
  
  return(pstar)
}

# ---------- Main function ----------
bootstrap_forecast_calibrated <- function(
    Z,
    B=200,
    alpha=0.05,
    train_size=500,
    k_gam=15,
    cores_limit = NULL,
    verbose = TRUE,
    do_calibrate = TRUE
) {
  
  n <- length(Z)
  
  Upper_raw <- Upper_cal <- PredMean <- rep(NA_real_, n)
  
  Viol_raw <- Viol_cal <- rep(NA, n)
  
  bootstrap_draws_list <- vector("list", n)
  
  # Starting Parallel System,
  if (is.null(cores_limit)) {
    cores <- detectCores() - 1
  } else {
    cores <- cores_limit
  }
  
  cores <- max(1, cores)
  
  cl <- makeCluster(cores)
  
  registerDoParallel(cl)
  
  clusterExport(
    cl,
    c(
      "fit_mu",
      "estimate_ou_params",
      "cond_coeffs",
      "num_derivative",
      "clamp",
      "log_likelihood_ou_recursive"
    )
  )
  
  clusterEvalQ(cl, {
    library(mgcv)
    library(optimx)
  })
  
  if (verbose) {
    message("Using ", cores, " cores for bootstrap replicates.")
  }
  
  start_t <- train_size
  
  for (t in start_t:(n-1)) {
    
    if (verbose && t %% 50 == 0) {
      message("Processing t = ", t, " / ", n-1)
    }
    
    Ztr <- Z[1:t]
    
    gam_result <- fit_mu(Ztr, k=k_gam)
    
    mu0 <- gam_result$mu_hat
    
    mu1 <- as.numeric(
      predict(
        gam_result$fit,
        newdata=data.frame(t=t+1)
      )
    )
    
    mu1 <- ifelse(
      is.na(mu1) | is.infinite(mu1),
      mean(mu0, na.rm=T),
      mu1
    )
    
    para <- estimate_ou_params(Ztr, mu0)
    
    cc <- cond_coeffs(
      para$theta,
      para$sigma2,
      para$omega2
    )
    
    pred_val <- mu1 + cc$c * (Z[t] - mu0[t])
    
    PredMean[t+1] <- ifelse(
      is.na(pred_val) | is.infinite(pred_val),
      mean(Ztr, na.rm=T),
      pred_val
    )
    
    mu0_local <- mu0
    t_local <- t
    Z_t_local <- Z[t]
    cc_c <- cc$c
    cc_v <- cc$v
    
    dr <- foreach(
      b = 1:B,
      .combine = c,
      .errorhandling = "remove",
      .packages = c("mgcv", "optimx")
    ) %dopar% {
      
      set.seed(b)
      
      Zb <- numeric(t_local)
      Zb[1] <- Ztr[1]
      
      for(s in 1:(t_local-1)) {
        
        m <- mu0_local[s+1] +
          cc_c * (Zb[s] - mu0_local[s])
        
        Zb[s+1] <- rnorm(1, m, sqrt(cc_v))
      }
      
      gam_b <- fit_mu(Zb, k=k_gam)
      
      mb <- gam_b$mu_hat
      
      mb1 <- as.numeric(
        predict(
          gam_b$fit,
          newdata=data.frame(t=t_local+1)
        )
      )
      
      mb1 <- ifelse(
        is.na(mb1),
        mean(mb, na.rm=T),
        mb1
      )
      
      pb <- estimate_ou_params(Zb, mb)
      
      cb <- cond_coeffs(
        pb$theta,
        pb$sigma2,
        pb$omega2
      )
      
      rnorm(
        1,
        mb1 + cb$c * (Z_t_local - mb[t_local]),
        sqrt(cb$v)
      )
    }
    
    dr <- dr[!is.na(dr) & is.finite(dr)]
    
    bootstrap_draws_list[[t+1]] <- dr
    
    if (length(dr) >= 5) {
      
      q_raw <- quantile(
        dr,
        1-alpha,
        na.rm=T,
        type=8
      )
      
      Upper_raw[t+1] <- ifelse(
        is.na(q_raw),
        mean(dr) + 2*sd(dr),
        q_raw
      )
      
      Upper_raw[t+1] <- max(
        Upper_raw[t+1],
        min(Z, na.rm=T) - 0.1
      )
      
      Viol_raw[t+1] <- Z[t+1] > Upper_raw[t+1]
      
    } else {
      
      Upper_raw[t+1] <- mean(Ztr, na.rm=T) +
        2*sd(Ztr, na.rm=T)
      
      Viol_raw[t+1] <- Z[t+1] > Upper_raw[t+1]
      
      if (verbose) {
        warning(
          paste0(
            "time point",
            t+1,
            " efficient sampling(",
            length(dr),
            "/",
            B,
            "), prediction interval"
          )
        )
      }
    }
  }
  
  stopCluster(cl)
  
  pstar <- 1 - alpha
  
  if (do_calibrate) {
    
    if (verbose) {
      message("starting calibrated")
    }
    
    pstar <- find_calibrated_p(
      bootstrap_draws_list,
      Z,
      1-alpha
    )
    
    for (t1 in 1:n) {
      
      dr <- bootstrap_draws_list[[t1]]
      
      if (length(dr) >= 5) {
        
        q_cal <- quantile(
          dr,
          pstar,
          na.rm=T,
          type=8
        )
        
        Upper_cal[t1] <- ifelse(
          is.na(q_cal),
          mean(dr) + 2*sd(dr),
          q_cal
        )
        
        Upper_cal[t1] <- max(
          Upper_cal[t1],
          min(Z, na.rm=T) - 0.1
        )
        
        Viol_cal[t1] <- Z[t1] > Upper_cal[t1]
        
      } else if (!is.na(Upper_raw[t1])) {
        
        Upper_cal[t1] <- Upper_raw[t1]
        Viol_cal[t1] <- Viol_raw[t1]
      }
    }
    
  } else {
    
    Upper_cal <- Upper_raw
    Viol_cal <- Viol_raw
  }
  
  cov_raw <- mean(
    !is.na(Viol_raw) & !Viol_raw,
    na.rm=T
  )
  
  cov_cal <- mean(
    !is.na(Viol_cal) & !Viol_cal,
    na.rm=T
  )
  
  results_df <- data.frame(
    time = 1:n,
    Z = Z,
    PredMean = PredMean,
    Upper_raw = Upper_raw,
    Upper_cal = Upper_cal,
    Viol_raw = Viol_raw,
    Viol_cal = Viol_cal
  )
  
  return(list(
    results = results_df,
    p_star = pstar,
    bootstrap_draws = bootstrap_draws_list,
    coverage = list(
      raw=cov_raw,
      calibrated=cov_cal
    )
  ))
}

# ---------- Gam for all training data ----------
fit_mu_all <- function(Z, k = 15) {
  df <- data.frame(Z=Z, t=seq_along(Z))
  gam(Z ~ s(t, bs="cs", k=min(k, max(3, nrow(df)-1))), data=df)
}

# ---------- Find first changing points ----------
find_first_changes <- function(time_vec, Z_vec) {
  
  if (length(Z_vec) == 0) {
    return(data.frame(
      time = numeric(0),
      Z = numeric(0)
    ))
  }
  
  change_idx <- c(
    1,
    which(diff(Z_vec) != 0) + 1
  )
  
  data.frame(
    time = time_vec[change_idx],
    Z = Z_vec[change_idx]
  )
}


# ---------- Find nearest previous change point ----------
get_last_change_index <- function(x, targets) {
  
  sapply(targets, function(i) {
    
    if (i <= 1) {
      return(NA)
    }
    
    prev_idx <- which(
      x[1:(i - 1)] != x[i]
    )
    
    if (length(prev_idx) == 0) {
      return(NA)
    }
    
    max(prev_idx)
  })
}


# ---------- Find next change point ----------
find_next_change <- function(x, pos) {
  
  start_val <- x[pos]
  
  idx <- which(
    x[(pos + 1):length(x)] != start_val
  )[1]
  
  if (is.na(idx)) {
    return(NA)
  }
  
  pos + idx
}

fit_mu_online <- function(Z, k = 15, verbose = FALSE) {
  n <- length(Z)
  mu_history <- vector("list", n)
  mu_next <- rep(NA_real_, n)
  
  for (t in seq(2, n)) {
    tt <- 1:t
    k_use <- min(k, max(3, t - 1))
    
    if (length(unique(tt)) < 3) {
      mu_history[[t]] <- rep(mean(Z[1:t]), t)
      mu_next[t] <- mean(Z[1:t])
      next
    }
    
    suppressWarnings({
      fit_t <- tryCatch(
        gam(Z[1:t] ~ s(tt, bs = "cs", k = k_use)),
        error = function(e) NULL
      )
    })
    
    if (is.null(fit_t)) {
      mu_hat_t <- rep(mean(Z[1:t]), t)
      mu_next[t] <- mean(Z[1:t])
    } else {
      mu_hat_t <- as.numeric(predict(fit_t))
      mu_next[t] <- as.numeric(
        predict(
          fit_t,
          newdata = data.frame(tt = t + 1)
        )
      )
    }
    
    mu_history[[t]] <- mu_hat_t
    
    if (verbose && t %% 50 == 0) {
      message("Online GAM update at t = ", t, "/", n)
    }
  }
  
  list(
    mu_history = mu_history,
    mu_next = mu_next
  )
}

# ---------- Calculate Returns ----------
calc_returns_metrics <- function(o_buy, o_sell, method = c("probability","odds")) {
  
  method <- match.arg(method)
  
  if(length(o_buy) != length(o_sell)) {
    stop("Lengths must match")
  }
  
  n <- length(o_buy)
  
  if(method == "probability") {
    
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

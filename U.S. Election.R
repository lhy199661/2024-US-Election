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
  if (theta < 1e-6 || log_sigma2 < -20 || log_sigma2 > 20 || log_omega2 < -20 || log_omega2 > 20) {
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
      P_pred <- a^2 * P_prev + (sigma2/(2*theta))*(1 - exp(-2*theta*h))
    }
    
    P_pred <- max(P_pred, 1e-10)
    S_t <- P_pred + omega2
    S_t <- max(S_t, 1e-10)
    
    err <- Z[t] - Q_pred
    ll <- ll - 0.5*(log(2*pi*S_t) + err^2 / S_t)
    
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
  
  # Step 1: ΔZ_t = θ (μ_t - Z_t) Δt + ε_t find the initial parameters
  dZ <- diff(Z)
  Xreg <- (m_hat[-n] - Z[-n]) * h
  
  fit0 <- lm(dZ ~ 0 + Xreg)
  
  theta_hat <- as.numeric(coef(fit0)[1])
  theta_hat <- max(theta_hat, 1e-6)
  
  se_theta <- summary(fit0)$coefficients[1,2]
  se_theta <- ifelse(is.na(se_theta), 0.01, se_theta)
  
  # Step 2: find the SE of parameters,
  res <- resid(fit0)
  sigma2_hat <- var(res)
  sigma2_hat <- max(sigma2_hat, 1e-6)
  
  # Initial parameters of sigma2 and omega2,
  sigma2_init <- sigma2_hat / 2
  omega2_init <- sigma2_hat / 2
  
  # Step 3: build the boundary,
  # theta
  theta_lower <- max(theta_hat - 3*se_theta, 1e-6)
  theta_upper <- theta_hat + 3*se_theta
  
  # Find the SE using X^2 distribution,
  var_se_factor <- 3 / sqrt(2*(n-1))
  
  sigma2_lower <- max(sigma2_init - var_se_factor*sigma2_hat, 1e-6)
  sigma2_upper <- sigma2_init + var_se_factor*sigma2_hat
  
  omega2_lower <- max(omega2_init - var_se_factor*sigma2_hat, 1e-6)
  omega2_upper <- omega2_init + var_se_factor*sigma2_hat
  
  # Log optimization
  log_sigma2_init <- log(sigma2_init)
  log_omega2_init <- log(omega2_init)
  
  lower_bounds <- c(theta_lower,
                    log(sigma2_lower),
                    log(omega2_lower))
  
  upper_bounds <- c(theta_upper,
                    log(sigma2_upper),
                    log(omega2_upper))
  
  init_params <- c(theta_hat,
                   log_sigma2_init,
                   log_omega2_init)
  
  # Step 4: L-BFGS-B optimization,
  opt <- tryCatch({
    optim(par = init_params,
          fn = log_likelihood_ou_recursive,
          method = "L-BFGS-B",
          lower = lower_bounds,
          upper = upper_bounds,
          control = list(maxit = 1000),
          Z = Z,
          mu_hat = m_hat)
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
  
  # Step 5: rebuild μ̂(t) = M(t) + M'(t)/θ,
  mu_final <- m_hat + mprime/theta_opt
  
  return(list(theta = theta_opt,
              sigma2 = sigma2_opt,
              omega2 = omega2_opt,
              mu_hat = mu_final))
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
find_calibrated_p <- function(bootstrap_draws_list, Z, target=0.95) {
  lower <- 0.5
  upper <- 0.999
  tol <- 1e-3
  max_iter <- 20 
  
  for (i in 1:max_iter) {
    mid <- (lower + upper) / 2
    cov_mid <- compute_coverage_at_p(bootstrap_draws_list, Z, mid)
    
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
bootstrap_forecast_calibrated <- function(Z, 
                                          B=200, 
                                          alpha=0.05, 
                                          train_size=500, 
                                          k_gam=15,
                                          cores_limit = NULL,
                                          verbose = TRUE,
                                          do_calibrate = TRUE) {
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
  
  # Input useful functions,
  clusterExport(cl, c("fit_mu", "estimate_ou_params", "cond_coeffs", "num_derivative", "clamp", "log_likelihood_ou_recursive"))
  clusterEvalQ(cl, {
    library(mgcv)
    library(optimx)
  })
  
  if (verbose) message("Using ", cores, " cores for bootstrap replicates.")
  
  # Deal with training data,
  start_t <- train_size 
  for (t in start_t:(n-1)) {
    if (verbose && t %% 50 == 0) {
      message("Processing t = ", t, " / ", n-1)
    }
    
    Ztr <- Z[1:t]
    gam_result <- fit_mu(Ztr, k=k_gam)
    mu0 <- gam_result$mu_hat
    mu1 <- as.numeric(predict(gam_result$fit, newdata=data.frame(t=t+1)))
    mu1 <- ifelse(is.na(mu1) | is.infinite(mu1), mean(mu0, na.rm=T), mu1)
    
    # Using LBFGS,
    para <- estimate_ou_params(Ztr, mu0)
    cc <- cond_coeffs(para$theta, para$sigma2, para$omega2)
    
    # Point estimate,
    pred_val <- mu1 + cc$c * (Z[t] - mu0[t])
    PredMean[t+1] <- ifelse(is.na(pred_val) | is.infinite(pred_val), mean(Ztr, na.rm=T), pred_val)
    
    # Bootstrap Sampling,
    mu0_local <- mu0
    t_local <- t
    Z_t_local <- Z[t]
    cc_c <- cc$c
    cc_v <- cc$v
    
    dr <- foreach(b = 1:B, .combine = c, .errorhandling = "remove", .packages = c("mgcv", "optimx")) %dopar% {
      set.seed(b) 
      Zb <- numeric(t_local)
      Zb[1] <- Ztr[1]
      for(s in 1:(t_local-1)) {
        m <- mu0_local[s+1] + cc_c * (Zb[s] - mu0_local[s])
        Zb[s+1] <- rnorm(1, m, sqrt(cc_v))
      }
      gam_b <- fit_mu(Zb, k=k_gam)
      mb <- gam_b$mu_hat
      mb1 <- as.numeric(predict(gam_b$fit, newdata=data.frame(t=t_local+1)))
      mb1 <- ifelse(is.na(mb1), mean(mb, na.rm=T), mb1)
      
      # LBFGS parameters optimization,
      pb <- estimate_ou_params(Zb, mb)
      cb <- cond_coeffs(pb$theta, pb$sigma2, pb$omega2)
      
      rnorm(1, mb1 + cb$c * (Z_t_local - mb[t_local]), sqrt(cb$v))
    }
    
    # Prediction interval,
    dr <- dr[!is.na(dr) & is.finite(dr)]
    bootstrap_draws_list[[t+1]] <- dr
    
    if (length(dr) >= 5) { 
      q_raw <- quantile(dr, 1-alpha, na.rm=T, type=8)
      Upper_raw[t+1] <- ifelse(is.na(q_raw), mean(dr) + 2*sd(dr), q_raw)
      Upper_raw[t+1] <- max(Upper_raw[t+1], min(Z, na.rm=T) - 0.1)
      Viol_raw[t+1] <- Z[t+1] > Upper_raw[t+1]
    } else {
      Upper_raw[t+1] <- mean(Ztr, na.rm=T) + 2*sd(Ztr, na.rm=T)
      Viol_raw[t+1] <- Z[t+1] > Upper_raw[t+1]
      if (verbose) {
        warning(paste0("time point", t+1, "efficient sampling(", length(dr), "/", B, "), prediction interval"))
      }
    }
  }
  
  # Close parallel,
  stopCluster(cl)
  
  # Calculate percentage,
  pstar <- 1 - alpha
  if (do_calibrate) {
    if (verbose) message("starting calibrated")
    pstar <- find_calibrated_p(bootstrap_draws_list, Z, 1-alpha)
    
    for (t1 in 1:n) {
      dr <- bootstrap_draws_list[[t1]]
      if (length(dr) >= 5) {
        q_cal <- quantile(dr, pstar, na.rm=T, type=8)
        Upper_cal[t1] <- ifelse(is.na(q_cal), mean(dr) + 2*sd(dr), q_cal)
        Upper_cal[t1] <- max(Upper_cal[t1], min(Z, na.rm=T) - 0.1)
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
  
  # Calculate covergance,
  cov_raw <- mean(!is.na(Viol_raw) & !Viol_raw, na.rm=T)
  cov_cal <- mean(!is.na(Viol_cal) & !Viol_cal, na.rm=T)
  
  
  # Sort out the output
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
    coverage = list(raw=cov_raw, calibrated=cov_cal)
  ))
}

# ---------- Output the parameter estimation ----------
generate_stat_report <- function(params_history, mu_history, train_size) {
  
  # Test sample 
  params_test <- params_history[(train_size+1):length(params_history)]
  mu_test     <- mu_history[(train_size+1):length(mu_history)]
  
  # Extract scalar params, 
  theta_vec <- sapply(params_test, function(x) if(is.null(x)) NA else x$theta)
  sigma_vec <- sapply(params_test, function(x) if(is.null(x)) NA else x$sigma)
  omega_vec <- sapply(params_test, function(x) if(is.null(x)) NA else x$omega)
  
  mat <- cbind(theta_vec, sigma_vec, omega_vec)
  colnames(mat) <- c("θ̂","σ̂","ω̂")
  
  stats_param <- data.frame(
    Parameter = colnames(mat),
    Point_Estimate = apply(mat, 2, function(x) tail(na.omit(x),1)),
    Mean = apply(mat, 2, mean, na.rm=TRUE),
    SD   = apply(mat, 2, sd, na.rm=TRUE),
    CI_low = apply(mat, 2, function(x){
      x <- na.omit(x)
      if(length(x)<2) return(NA)
      mean(x) - 1.96*sd(x)/sqrt(length(x))
    }),
    CI_high = apply(mat, 2, function(x){
      x <- na.omit(x)
      if(length(x)<2) return(NA)
      mean(x) + 1.96*sd(x)/sqrt(length(x))
    })
  )
  
  # ---- Summarize mu(t) ----
  # Use last point of each mu_hat curve (corresponding to time t)
  mu_last_vec <- sapply(mu_test, function(x) if(is.null(x)) NA else tail(x,1))
  
  stats_mu <- data.frame(
    Parameter = "μ̂(t)",
    Point_Estimate = tail(na.omit(mu_last_vec),1),
    Mean = mean(mu_last_vec, na.rm=TRUE),
    SD   = sd(mu_last_vec, na.rm=TRUE),
    CI_low = {
      x <- na.omit(mu_last_vec)
      if(length(x)<2) NA else mean(x)-1.96*sd(x)/sqrt(length(x))
    },
    CI_high = {
      x <- na.omit(mu_last_vec)
      if(length(x)<2) NA else mean(x)+1.96*sd(x)/sqrt(length(x))
    }
  )
  
  stats_all <- rbind(stats_param, stats_mu)
  
  cat("=== Parameter Report (Test Sample Only) ===\n")
  print(stats_all, row.names=FALSE, digits=4)
  
  invisible(stats_all)
}

#---- Data simulation For Trending OU process ----
simulate_trending_ou <- function(n,
                                 theta,
                                 sigma2,
                                 omega2,
                                 seed) {
  set.seed(seed)
  h <- 1
  t_vec <- 1:n
  mu_t <- 5 + 0.001*t_vec + 0.5*sin(2*pi*t_vec/200)
  
  theta <- max(theta, 1e-6)
  sigma2 <- max(sigma2, 1e-6)
  omega2 <- max(omega2, 1e-6)
  
  Q <- numeric(n)
  Z <- numeric(n)
  Q[1] <- mu_t[1]
  sigma <- sqrt(sigma2)
  omega <- sqrt(omega2)
  Z[1]<-Q[1]+rnorm(1, 0, omega)
  
  for (t in 2:n) {
    dW <- rnorm(1, 0, sqrt(h))
    Q[t] <- Q[t-1] + theta*(mu_t[t-1] - Q[t-1])*h + sigma*dW
    E <- rnorm(1, 0, omega)
    Z[t] <- Q[t] + E
  }
  
  list(Z = Z, Q = Q, mu_t = mu_t,
       true_params = list(theta=theta, sigma2=sigma2, omega2=omega2))
}

# ---- Uploading the gambling data for 2020 and 2024 U.S. Election. Since the data have already done the cleaning, we do not need to do the data cleaning ----.
# ---- 2020 data as training data ----
# This is full dataset,
c <- read.csv("2020USElection.csv",header = T)

# After we exclude the covid time
d <- read.csv("2020USElectionwithoutcovid.csv",header=T)

# 2024 data as test data,
b <- read.csv("2024USElection.csv",header = T)

# Y is sum of implied probabilities of two candidates.
# For 2020 data,
# Y <- d$Total.

# For 2024 data, 
Y <- b$Prob

set.seed(239539)
n <- length(Y)

# Set the training size: 10000, based on the length of the whole data.
train_size <- 10000

# Run the output.
monitor_result <- bootstrap_forecast_calibrated(Y, train_size=10000)

# Try to find the last changeable time.
find_first_changes <- function(time_vec, Z_vec) {
  if (length(Z_vec) == 0) return(data.frame(time = numeric(0), Z = numeric(0)))
  change_idx <- c(1, which(diff(Z_vec) != 0) + 1)
  data.frame(time = time_vec[change_idx], Z = Z_vec[change_idx])
}
l <- monitor_result$results
l_true <- l[l$Viol_cal == TRUE, ]

# Find the signal trading points filtered.
change_points <- find_first_changes(l_true$time, l_true$Z)
filtered <- change_points$time
filtered <- na.omit(filtered)

#For smoothing part of mu.
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
      mu_next[t] <- as.numeric(predict(fit_t, newdata = data.frame(tt = t + 1)))
    }
    
    mu_history[[t]] <- mu_hat_t
    
    if (verbose && t %% 50 == 0) {
      message("Online GAM update at t = ", t, "/", n)
    }
  }
  
  list(mu_history = mu_history, mu_next = mu_next)
}
res <- fit_mu_online(Y, k = 15, verbose = TRUE)
r <- unlist(lapply(res$mu_history, function(x) tail(x, 1)))

# Gam for all training data.
fit_mu_all <- function(Z, k = 15) {
  df <- data.frame(Z=Z, t=seq_along(Z))
  gam(Z ~ s(t, bs="cs", k=min(k, max(3, nrow(df)-1))), data=df)
}
t <- as.numeric(predict(fit_mu_all(Y[1:10000])))
r[1:10000] <- t
r <- c(r,NA) # exclude the training size.

# Define the English time zone.
Sys.setlocale("LC_TIME", "English")
time <- as.POSIXct(b$timestampLON,
                   format="%Y/%m/%d %H:%M")

# If we deal with 2020 dataset, we should exclude the Covid period.
# time <- as.POSIXct(c$timestampLON,format="%Y/%m/%d %H:%M")
# Y <-  append(Y,rep(NA,1079),after=11593)
# r <- append(r,rep(NA,1079),after=11593)
# monitor_result$upper<-append(monitor_result$upper,rep(NA,1079),after=11593)
# Ul <- l$Upper_cal
# Ul <- append(Ul,rep(NA,1079),after=11593)
# Ul_break <- Ul
# filtered <- ifelse(filtered < 11594, filtered, filtered + 1079)

# Define the 2020 dataset
# df_Yp <- data.frame(time = time, Y = Y)
# df_covid <- data.frame(time = time, covid = r)
# df_Ul <- data.frame(time = time, Ul = Ul_break)

# df_points <- data.frame(time = time[filtered],Y = Y[filtered])


# q <- ggplot() +
#  geom_line(data = df_Yp, aes(x = time, y = Y), color = "gray60") +
#  geom_line(data = df_covid, aes(x = time, y = r), color = "green3") +
  
#  geom_line(data = df_Ul, aes(x = time, y = Ul), color = "red") +
  
#  geom_point(
#    data = df_points,
#    aes(x = time, y = Y),
#    shape = 4, color = "blue", size = 1
#  ) +
#  geom_vline(
#    xintercept = time[train_size],
#   color = "magenta",
#    linetype = "dashed"
#  ) +
#  geom_vline(xintercept = time[11593], color = "orange", linetype = "dashed") +
#  geom_vline(xintercept = time[12673], color = "orange", linetype = "dashed") +
  
#  theme_bw() +
#  labs(x = "Date", y = "Sum Of Probability")

# Define the 2024 dataset,
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
#Draw the plot
p <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = Y), color = "gray60") +
  geom_line(aes(y = Upper), color = "red", linewidth = 0.7) +
  geom_line(aes(y = r), color = "green3", linewidth = 0.7) +
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

# We find the outside prediction interval points.
filtered <- na.omit(change_points$time)

# Calculate the prediction rate.
# change_count1 <- length(filtered)
# change_count2 <- sum(na.omit(diff(Y[-(1:10000)])) != 0)
# change_count1/change_count2

# Which Bet to place, Using 2020 data to build the Bradley–Terry Model, we have two candidates Trump and Biden from two parties.
Trump <- as.numeric(na.omit(c$probTrump))
Biden <- as.numeric(na.omit(c$probBiden))

# For 2024 test data, we need to check the model.
Trump1 <- as.numeric(b$Donald.Trump.Prob)
Harris1 <- as.numeric(b$Kamala.Harris.Prob)

# According to the outside points find the nearest last changeable points.
get_last_change_index <- function(x, targets) {
  sapply(targets, function(i) {
    if (i <= 1) return(NA)
    prev_idx <- which(x[1:(i - 1)] != x[i])
    if (length(prev_idx) == 0) return(NA)
    return(max(prev_idx))
  })
}

# We need to find the trading points filtered1.
find_next_change <- function(x, pos) {
  start_val <- x[pos]
  idx <- which(x[(pos + 1):length(x)] != start_val)[1]
  if (is.na(idx)) {
    return(NA)   
  } else {
    return(pos + idx)  
  }
}

# Filtered is the signal points and fitered1 is the trading points. 
filtered1 <- sapply(filtered, function(p) find_next_change(Y, p))

# We use the same methods for 2020 dataset which exclude the Covid time(points: 11593 to 12673).
# For Trump we have three points, outside points (predicted points), one point ahead outside points and the nearest changeable points.
last_change_indices <- get_last_change_index(Trump, filtered)
T1 <- Trump[filtered1]
T2 <- Trump[filtered]
T3 <- Trump[last_change_indices]

#For Biden we have the same three points.
last_change_indices <- get_last_change_index(Biden, filtered)
B1 <- Biden[filtered1]
B2 <- Biden[filtered]
B3 <- Biden[last_change_indices]

# First step, build pairwise dataset with feature vectors.
# Feature vectors of Delta Trump and Biden
DeltaTrump1 <- T2-T3
DeltaBiden1 <- B2-B3

feats <- data.frame(
  DeltaTrump1 = DeltaTrump1,
  T2 = T2,
  DeltaBiden1 = DeltaBiden1,
  B2 = B2
)

# Give the lable whether to place bet on Trump or Biden.
utility_1 <- 1/T1-1/T2
utility_2 <- 1/B1-1/B2
feats$label <- as.integer(utility_1>utility_2)

build_pairwise_data <- function(feats) {
  feats$Delta1 <- feats$DeltaTrump1 - feats$DeltaBiden1
  feats$Delta2 <- feats$T2 - feats$B2
  feats
}

pair_df <- build_pairwise_data(feats)
bt_data <- data.frame(
  Trump = feats$label,
  Biden = 1 - feats$label
)

bt_data$diff_Delta <- feats$DeltaTrump1 - feats$DeltaBiden1
bt_data$diff_T2 <- feats$T2 - feats$B2

# Fit Bradley–Terry logistic model.
bt_model <- glm(Trump ~ diff_Delta + diff_T2, data = bt_data, family = binomial())
summary(bt_model)

# Extract fitted probabilities from the Bradley–Terry model.
probs <- predict(bt_model, type = "response")

# Extract the observed place indicator.
#    (1 = Place Trump, 0 = Place Biden)
obs <- bt_data$Trump

# Hosmer–Lemeshow goodness-of-fit test.
hl_test <- hoslem.test(obs, probs)

# Then we use 2024 data to do the prediction, Trump data is T and Harris data is H. Before this we need to use bootstrap_forecast_calibrated(Y, train_size=10000) to find the signal points for 2024.
last_change_indices <- get_last_change_index(Trump1, filtered)
T1 <- Trump1[filtered1]
T2 <- Trump1[filtered]
T3 <- Trump1[last_change_indices]

#For Harris we have the same three points.
last_change_indices <- get_last_change_index(Harris1, filtered)
H1 <- Harris1[filtered1]
H2 <- Harris1[filtered]
H3 <- Harris1[last_change_indices]

# Build the covariates
DeltaTrump1 <- T2-T3
DeltaHarris1 <- H2-H3

newdata <- data.frame(
  diff_Delta = DeltaTrump1-DeltaHarris1,
  diff_T2 = T2-H2
)

bt_pred_prob <- predict(bt_model, newdata = newdata, type = "response")
bt_pred <- as.integer(bt_pred_prob > 0.5)

# We need to check the expected return and sharpe ratio.
selected_value <- ifelse(bt_pred > 0, 1/T1, 1/H1)
selected_baseline <- ifelse(bt_pred > 0, 1/T2, 1/H2)
profit <- selected_value - selected_baseline

# Return_metrics.R,
calc_returns_metrics <- function(o_buy, o_sell, method = c("probability","odds")) {
  # o_buy, o_sell : numeric vectors of buy/sell odds (same length)
  # method : "probability" uses R = 1/o_sell - 1/o_buy (recommended)
  #          "odds"        uses R = (o_sell - o_buy)/o_buy
  method <- match.arg(method)
  if(length(o_buy) != length(o_sell)) stop("Lengths must match")
  n <- length(o_buy)
  
  if(method == "probability") {
    R <- 1 / o_sell - 1 / o_buy   # Absolute change in implied probability
  } else {
    R <- (o_sell - o_buy) / o_buy # Relative change in odds price
  }
  
  # Sample statistics
  ER <- mean(R, na.rm = TRUE)         # Sample mean (expected return)
  sdR <- sd(R, na.rm = TRUE)          # Sample std dev
  SR <- ifelse(sdR > 0, ER / sdR, NA) # Sharpe ratio (risk-free rate assumed zero)
  
  list(
    method = method,
    returns = R,
    mean_return = ER,
    sd_return = sdR,
    sharpe = SR,
    n = n
  )
}

o_buy  <- selected_baseline
o_sell <- selected_value
res_odds <- calc_returns_metrics(o_buy, o_sell, method = "odds")

# After that we build the reliable strategy to place the bet for US election.


rm(list = ls())

# =========================
# Load scripts
# =========================
source("functions.R")

source("data_2020.R")
source("data_2024.R")

source("forecasting_ou.R")
source("modeling_bt.R")
source("plotting.R")

# =========================================================
# PART I : 2020 DATA
# Bradley–Terry Training
# =========================================================

cat("=====================================\n")
cat("2020 Bradley–Terry Training\n")
cat("=====================================\n")

# ---------------------------------------------------------
# Build 2020 signal points
# ---------------------------------------------------------

monitor_result_2020 <- run_ou_forecast(
  Y=Y_2020,
  train_size = 10000,
  B = 200
)

filtered_2020 <- get_signal_points(
  monitor_result = monitor_result_2020
)

# ---------------------------------------------------------
# Plot 2020 results
# ---------------------------------------------------------

plot_2020_results(
  Y = Y_2020,
  monitor_result = monitor_result_2020,
  filtered = filtered_2020,
  time = time_2020,
  train_size = 10000
)

# ---------------------------------------------------------
# Build BT training dataset
# ---------------------------------------------------------

bt_train <- build_bt_dataset(
  candidate1 = Trump_2020,
  candidate2 = Biden_2020,
  filtered = filtered_2020,
  Y = Y_2020
)

# ---------------------------------------------------------
# Fit Bradley–Terry model
# ---------------------------------------------------------

bt_fit <- fit_bt_model(bt_train)

summary(bt_fit$model)

# ---------------------------------------------------------
# Hosmer–Lemeshow Test
# ---------------------------------------------------------

train_probs <- predict(
  bt_fit$model,
  type = "response"
)

train_obs <- bt_train$Candidate1

hl_test <- hoslem.test(
  train_obs,
  train_probs
)

print(hl_test)

# =========================================================
# PART II : 2024 DATA
# Forecasting + Prediction
# =========================================================

cat("=====================================\n")
cat("2024 Forecasting and Prediction\n")
cat("=====================================\n")

# ---------------------------------------------------------
# Forecasting for 2024
# ---------------------------------------------------------

monitor_result_2024 <- run_ou_forecast(
  Y = Y_2024,
  train_size = 10000,
  B = 200
)

# ---------------------------------------------------------
# Find Signal Points
# ---------------------------------------------------------

filtered_2024 <- get_signal_points(
  monitor_result = monitor_result_2024
)

# ---------------------------------------------------------
# Plot 2024 Results
# ---------------------------------------------------------

plot_2024_results(
  Y = Y_2024,
  monitor_result = monitor_result_2024,
  filtered = filtered_2024,
  time = time_2024,
  train_size = 10000
)

# ---------------------------------------------------------
# Build 2024 Test Dataset
# ---------------------------------------------------------

bt_test <- build_bt_dataset(
  candidate1 = Trump_2024,
  candidate2 = Harris_2024,
  filtered = filtered_2024,
  Y = Y_2024
)

# ---------------------------------------------------------
# Predict betting direction
# ---------------------------------------------------------

test_probs <- predict(
  bt_fit$model,
  newdata = bt_test,
  type = "response"
)

bt_pred <- as.integer(
  test_probs > 0.5
)

# =========================================================
# PART III : Strategy Evaluation
# =========================================================

cat("=====================================\n")
cat("Strategy Evaluation\n")
cat("=====================================\n")

strategy_result <- evaluate_strategy(
  bt_pred = bt_pred,
  bt_test = bt_test
)

# ---------------------------------------------------------
# Print results
# ---------------------------------------------------------

cat("Total Profit:\n")
print(strategy_result$profit_sum)

cat("Performance Metrics:\n")
print(strategy_result$metrics)

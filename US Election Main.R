rm(list = ls())

# ==================================================
# Load Scripts
# ==================================================

source("functions.R")

source("data_2020.R")
source("data_2024.R")

source("forecasting_ou.R")
source("modeling_bt.R")
source("plotting.R")

output_dir <- "results1"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!exists("time_index_2020")) {
  time_index_2020 <- convert_2020_index_to_full(
    seq_along(Y_2020)
  )
}

if (!exists("time_index_2024")) {
  time_index_2024 <- seq_along(Y_2024)
}

p_candidate_2020 <- plot_2020_candidate_probabilities(
  time = time_2020,
  trump = Trump_2020,
  biden = Biden_2020
)

p_candidate_2024 <- plot_2024_candidate_probabilities(
  time = time_2024,
  trump = Trump_2024,
  harris = Harris_2024
)

print_signal_diagnostics <- function(monitor_result,
                                     label,
                                     n_top = 10) {
  
  l <- monitor_result$results
  
  cat(label, "coverage:", monitor_result$coverage, "\n")
  cat(label, "signals:", sum(l$Signal, na.rm = TRUE), "\n")
  
  if ("ProbDecrease" %in% names(l)) {
    cat(label, "ProbDecrease summary:\n")
    print(summary(l$ProbDecrease))
  }
  
  score <- if ("ProbDecrease" %in% names(l)) {
    l$ProbDecrease
  } else {
    l$Z - l$Upper
  }
  
  keep <- is.finite(score)
  
  if (sum(keep) > 0) {
    ord <- order(score[keep], decreasing = TRUE)
    diagnostic_df <- l[keep, ]
    cols <- intersect(
      c(
        "effective_time",
        "original_index",
        "time_index",
        "h",
        "Z",
        "Upper",
        "ProbDecrease",
        "PredMean",
        "PredVar",
        "theta",
        "sigma2",
        "omega2"
      ),
      names(diagnostic_df)
    )
    
    cat(label, "top signal diagnostics:\n")
    print(
      head(
        diagnostic_df[ord, cols],
        n_top
      )
    )
  }
}

# ==================================================
# PART I
# 2020 TRAINING
# ==================================================

cat("=====================================\n")
cat("2020 Training Period\n")
cat("=====================================\n")

# --------------------------------------------------
# Forecasting
# --------------------------------------------------

monitor_result_2020 <- run_ou_forecast(
  Y = Y_2020,
  train_size = 10000,
  B = 200,
  skip_flat = TRUE,
  time_index = time_index_2020
)

ci_2020 <- get_bootstrap_param_ci(
  monitor_result_2020,
  method = "percentile"
)

write.csv(
  monitor_result_2020$results,
  file = file.path(output_dir, "monitor_result_2020.csv"),
  row.names = FALSE
)

saveRDS(
  monitor_result_2020,
  file = file.path(output_dir, "monitor_result_2020.rds")
)

last_param_2020 <- tail(
  na.omit(
    monitor_result_2020$results[, c("theta", "sigma2", "omega2")]
  ),
  1
)

cat("Final 2020 parameter estimates:\n")
print(last_param_2020)

print_signal_diagnostics(
  monitor_result = monitor_result_2020,
  label = "2020"
)

# --------------------------------------------------
# Signal Points
# --------------------------------------------------

filtered_2020 <- get_signal_points(
  monitor_result_2020
)

cat(
  "Number of test signals:",
  length(filtered_2020),
  "\n"
)

if (length(filtered_2020) == 0) {
  stop(
    paste(
      "No 2020 signal points were found.",
      "Bradley-Terry training cannot proceed.",
      "Check the printed ProbDecrease diagnostics and the COVID-aware h values."
    )
  )
}

# --------------------------------------------------
# Plot 2020 results
# --------------------------------------------------

p_2020 <- plot_2020_results(
  Y = Y_2020,
  monitor_result = monitor_result_2020,
  filtered = filtered_2020,
  time = time_2020,
  train_size = 10000
)

ggsave(
  filename = file.path(output_dir, "US_Election_2020_forecast.eps"),
  plot = p_2020,
  device = "eps",
  width = 8,
  height = 5
)


# --------------------------------------------------
# Build BT Dataset
# --------------------------------------------------

bt_train <- build_bt_dataset(
  candidate1 = Trump_2020,
  candidate2 = Biden_2020,
  filtered = filtered_2020,
  Y = Y_2020,
  index_map = convert_2020_index_to_full
)

if (nrow(bt_train) == 0) {
  stop("The 2020 Bradley-Terry training data are empty after filtering.")
}

# --------------------------------------------------
# Bradley-Terry
# --------------------------------------------------

bt_fit <- fit_bt_model(
  bt_train
)

summary(
  bt_fit$model
)

# --------------------------------------------------
# Hosmer-Lemeshow
# --------------------------------------------------

train_prob <- predict(
  bt_fit$model,
  type = "response"
)

hl_test <- hoslem.test(
  bt_train$Candidate1,
  train_prob
)

print(hl_test)

# ==================================================
# PART II
# 2024 TESTING
# ==================================================

cat("=====================================\n")
cat("2024 Testing Period\n")
cat("=====================================\n")

# --------------------------------------------------
# Forecasting
# --------------------------------------------------

monitor_result_2024 <- run_ou_forecast(
  Y = Y_2024,
  train_size = 10000,
  B = 200,
  skip_flat = TRUE,
  time_index = time_index_2024
)

ci_2024 <- get_bootstrap_param_ci(
  monitor_result_2024,
  method = "percentile"
)

write.csv(
  monitor_result_2024$results,
  file = file.path(output_dir, "monitor_result_2024.csv"),
  row.names = FALSE
)

saveRDS(
  monitor_result_2024,
  file = file.path(output_dir, "monitor_result_2024.rds")
)

last_param_2024 <- tail(
  na.omit(
    monitor_result_2024$results[, c("theta", "sigma2", "omega2")]
  ),
  1
)

cat("Final 2024 parameter estimates:\n")
print(last_param_2024)

print_signal_diagnostics(
  monitor_result = monitor_result_2024,
  label = "2024"
)

# --------------------------------------------------
# Signal Points
# --------------------------------------------------

filtered_2024 <- get_signal_points(
  monitor_result_2024
)

cat(
  "Number of testing signals:",
  length(filtered_2024),
  "\n"
)

# --------------------------------------------------
# Plot 2024 results
# --------------------------------------------------

p_2024 <- plot_2024_results(
  Y = Y_2024,
  monitor_result = monitor_result_2024,
  filtered = filtered_2024,
  time = time_2024,
  train_size = 10000
)

ggsave(
  filename = file.path(output_dir, "US_Election_2024_forecast.eps"),
  plot = p_2024,
  device = "eps",
  width = 8,
  height = 5
)


# --------------------------------------------------
# Build Test Dataset
# --------------------------------------------------

bt_test <- build_bt_dataset(
  candidate1 = Trump_2024,
  candidate2 = Harris_2024,
  filtered = filtered_2024,
  Y = Y_2024,
  index_map = NULL
)

cat(
  "Testing observations:",
  nrow(bt_test),
  "\n"
)

if (nrow(bt_test) == 0) {
  stop("The 2024 Bradley-Terry test data are empty after filtering.")
}

# --------------------------------------------------
# Predict Direction
# --------------------------------------------------

test_prob <- predict(
  bt_fit$model,
  newdata = bt_test,
  type = "response"
)

bt_pred <- as.integer(
  test_prob > 0.5
)

# ==================================================
# PART III
# Strategy Evaluation
# ==================================================

cat("=====================================\n")
cat("Strategy Evaluation\n")
cat("=====================================\n")

strategy_result <- evaluate_strategy(
  bt_pred = bt_pred,
  bt_test = bt_test
)

# --------------------------------------------------
# Results
# --------------------------------------------------

cat(
  "Total Profit:\n"
)

print(
  strategy_result$profit_sum
)

cat(
  "Performance Metrics:\n"
)

print(
  strategy_result$metrics
)

# --------------------------------------------------
# Classification Metrics
# --------------------------------------------------

accuracy <- mean(
  bt_pred ==
    bt_test$Candidate1
)

cat(
  "Prediction Accuracy:",
  round(accuracy,4),
  "\n"
)

# After that we build the reliable strategy to place the bet for US election.

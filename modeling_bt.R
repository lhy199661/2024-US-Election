# =========================
# Build BT Dataset
# =========================

build_bt_dataset <- function(candidate1,
                             candidate2,
                             filtered,
                             Y) {
  
  filtered1 <- sapply(filtered, function(p) {
    find_next_change(Y, p)
  })
  
  last_change_indices1 <- get_last_change_index(
    candidate1,
    filtered
  )
  
  last_change_indices2 <- get_last_change_index(
    candidate2,
    filtered
  )
  
  C1_trade  <- candidate1[filtered1]
  C1_signal <- candidate1[filtered]
  C1_prev   <- candidate1[last_change_indices1]
  
  C2_trade  <- candidate2[filtered1]
  C2_signal <- candidate2[filtered]
  C2_prev   <- candidate2[last_change_indices2]
  
  Delta1 <- C1_signal - C1_prev
  Delta2 <- C2_signal - C2_prev
  
  utility_1 <- 1 / C1_trade - 1 / C1_signal
  utility_2 <- 1 / C2_trade - 1 / C2_signal
  
  label <- as.integer(utility_1 > utility_2)
  
  data.frame(
    Candidate1 = label,
    Candidate2 = 1 - label,
    diff_Delta = Delta1 - Delta2,
    diff_signal = C1_signal - C2_signal,
    C1_trade = C1_trade,
    C1_signal = C1_signal,
    C2_trade = C2_trade,
    C2_signal = C2_signal
  )
}

# =========================
# Fit BT Model
# =========================

fit_bt_model <- function(bt_data) {
  
  model <- glm(
    Candidate1 ~ diff_Delta + diff_signal,
    data = bt_data,
    family = binomial()
  )
  
  probs <- predict(
    model,
    type = "response"
  )
  
  list(
    model = model,
    probabilities = probs
  )
}

# =========================
# Strategy Evaluation
# =========================

evaluate_strategy <- function(bt_pred,
                              bt_test) {
  
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
    method = "odds"
  )
  
  list(
    profit = profit,
    profit_sum = sum(profit),
    metrics = metrics
  )
}
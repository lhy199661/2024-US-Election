# =========================================================
# modeling_bt.R
# Bradley-Terry training, prediction, and strategy evaluation
# =========================================================

build_bt_training_data_2020 <- function(Trump_2020,
                                        Biden_2020,
                                        filtered_2020,
                                        Y_2020) {
  
  build_bt_dataset(
    candidate1 = Trump_2020,
    candidate2 = Biden_2020,
    filtered = filtered_2020,
    Y = Y_2020,
    index_map = convert_2020_index_to_full
  )
}

build_bt_test_data_2024 <- function(Trump_2024,
                                    Harris_2024,
                                    filtered_2024,
                                    Y_2024) {
  
  build_bt_dataset(
    candidate1 = Trump_2024,
    candidate2 = Harris_2024,
    filtered = filtered_2024,
    Y = Y_2024,
    index_map = NULL
  )
}

fit_bt_training_model <- function(bt_train) {
  
  if (nrow(bt_train) == 0) {
    stop("bt_train is empty. No Bradley-Terry model can be fitted.")
  }
  
  fit_bt_model(
    bt_data = bt_train
  )
}

predict_bt_direction <- function(bt_fit,
                                 bt_test,
                                 cutoff = 0.5) {
  
  if (nrow(bt_test) == 0) {
    stop("bt_test is empty. No prediction can be made.")
  }
  
  probs <- predict(
    bt_fit$model,
    newdata = bt_test,
    type = "response"
  )
  
  pred <- as.integer(
    probs > cutoff
  )
  
  return(list(
    probabilities = probs,
    prediction = pred
  ))
}

evaluate_bt_strategy <- function(bt_pred,
                                 bt_test,
                                 return_method = "absolute") {
  
  evaluate_strategy(
    bt_pred = bt_pred,
    bt_test = bt_test,
    return_method = return_method
  )
}

run_hosmer_lemeshow <- function(bt_fit,
                                bt_train) {
  
  probs <- predict(
    bt_fit$model,
    type = "response"
  )
  
  obs <- bt_train$Candidate1
  
  hoslem.test(
    obs,
    probs
  )
}

run_bt_test_metrics <- function(bt_test,
                                test_probs,
                                bt_pred) {
  
  accuracy <- mean(
    bt_pred == bt_test$Candidate1,
    na.rm = TRUE
  )
  
  if (length(unique(bt_test$Candidate1)) < 2) {
    return(list(
      accuracy = accuracy,
      auc = NA_real_,
      roc = NULL
    ))
  }
  
  roc_obj <- pROC::roc(
    response = bt_test$Candidate1,
    predictor = test_probs,
    quiet = TRUE
  )
  
  auc_value <- as.numeric(
    pROC::auc(roc_obj)
  )
  
  return(list(
    accuracy = accuracy,
    auc = auc_value,
    roc = roc_obj
  ))
}
library(caret)
library(randomForest)
library(boot)
library(dplyr)

# Define RMSE function
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
}

# Function to average results across folds (adapted for regression)
library(ggplot2)

# Function to average results across folds (adapted for regression)
library(ggplot2)
library(dplyr)

# Function to average results across folds (adapted for regression)
average_rf_regression = function(train_rmse_scores, test_rmse_scores,
                                 all_labels_train, all_labels_test,
                                 feature_importance_values,
                                 final_models = NULL){  # Pass list of fitted models if partial dependence desired
  
  # Bootstrap RMSE scores
  boot_rmse_train = boot(data = train_rmse_scores, statistic = function(d, i) mean(d[i]), R = 1000)
  boot_rmse_test = boot(data = test_rmse_scores, statistic = function(d, i) mean(d[i]), R = 1000)
  
  avg_rmse_train = boot_rmse_train$t0
  avg_rmse_test = boot_rmse_test$t0
  
  # Confidence intervals
  ci_train = boot.ci(boot_rmse_train, type = "perc")$percent[4:5]
  ci_test = boot.ci(boot_rmse_test, type = "perc")$percent[4:5]
  
  if(is.null(ci_train)) ci_train = c(1,1)
  
  # Combine test predictions
  test_labels = bind_rows(all_labels_test)
  
  # Combine train predictions (averaged since each sample appears in multiple folds)
  train_labels = bind_rows(all_labels_train) %>%
    group_by(row, true_labels) %>%
    summarise(predicted_values = mean(predicted_values), .groups = 'drop')
  
  # Average feature importances
  mean_feature_importance = Reduce("+", feature_importance_values) / length(feature_importance_values)
  
  # Create importance dataframe using %IncMSE
  importance_df = data.frame(
    Feature = rownames(mean_feature_importance),
    Importance = mean_feature_importance[, "%IncMSE"]
  ) %>%
    arrange(desc(Importance))
  
  rownames(importance_df) = NULL
  
  # Performance metrics
  correlation <- cor(test_labels$true_labels, test_labels$predicted_values)
  r_squared <- correlation^2
  
  sd_y <- sd(test_labels$true_labels)
  rmse <- sqrt(mean((test_labels$true_labels - test_labels$predicted_values)^2))
  mae <- mean(abs(test_labels$true_labels - test_labels$predicted_values))
  
  performance_metrics <- tibble::tribble(
    ~Metric, ~Value,
    "R-squared", round(r_squared, 3),
    "RMSE", round(rmse, 3),
    "MAE", round(mae, 3),
    "SD of True Outcome", round(sd_y, 3)
  )
  
  results = list(
    rmse_train = avg_rmse_train,
    rmse_test = avg_rmse_test,
    rmse_train_ci = ci_train,
    rmse_test_ci = ci_test,
    test_performance_metrics = performance_metrics,
    test_labels = test_labels,
    train_labels = train_labels,
    importance = importance_df
  )
  
  # === VISUALIZATIONS === #
  
  # Feature Importance Bar Chart
  p1 = ggplot(importance_df[1:10, ], aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_col() +
    coord_flip() +
    labs(title = "Top 10 Most Important Features", x = "Features", y = "% Increase in MSE")
  # print(p1)
  
  # Residuals Plot (Test Data Only)
  test_residuals_plot = ggplot(test_labels, aes(x = true_labels, y = true_labels - predicted_values)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(title = "Residuals vs True Values (Test Set)",
         x = "True Outcome",
         y = "Residuals")
  # print(test_residuals_plot)
  
  # Predicted vs Actual Scatter Plot (Test Data Only)
  pred_vs_actual_plot = ggplot(test_labels, aes(x = true_labels, y = predicted_values)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
    labs(title = "Predicted vs Actual Outcomes (Test Set)",
         x = "Actual Outcome",
         y = "Predicted Outcome")
  # print(pred_vs_actual_plot)
  
  # QQ Plot of Residuals (Base R version captured)
  residuals = test_labels$true_labels - test_labels$predicted_values
  qqnorm(residuals, main = "QQ Plot of Residuals")
  qqline(residuals, col = "red")
  
  # Capture the current plot
  qq_plot = recordPlot()
  
  # Training vs Testing RMSE per fold
  rmse_comparison_df = data.frame(
    Fold = seq_along(train_rmse_scores),
    Train_RMSE = train_rmse_scores,
    Test_RMSE = test_rmse_scores
  )
  
  rmse_long = tidyr::pivot_longer(rmse_comparison_df, cols = c(Train_RMSE, Test_RMSE), names_to = "Type", values_to = "RMSE")
  
  rmse_plot = ggplot(rmse_long, aes(x = factor(Fold), y = RMSE, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Per-Fold RMSE Comparison",
         x = "Fold Index",
         y = "Root Mean Square Error")
  # print(rmse_plot)
  
  # Partial Dependence Plots (requires passing final_models)
  if (!is.null(final_models) && length(final_models) > 0) {
    top_features = head(importance_df$Feature, 3)
    
    par(mfrow = c(1, length(top_features)))
    for (feat in top_features) {
      partialPlot(final_models[[1]], pred.data = test_labels, x.var = feat,
                  main = paste("Partial Dependence:", feat))
    }
    par(mfrow = c(1,1))
  }
  
  return(list(results=results,plots=list(
    feature_importance_plot = p1,
    residuals_plot = test_residuals_plot,
    pred_vs_actual_plot = pred_vs_actual_plot,
    residuals_qq_plot = qq_plot,
    rmse_comparison_plot = rmse_plot
  )))
}



# Main function to run Random Forest with cross-validation and RMSE tracking
run_rf = function(X, y, fold_list,
                  hyper, rngseed = 421,
                  kfold = TRUE) {
  
  # X = predictors; y = outcome
  # fold_list = folds; hyper = tune_grid;
  # rngseed=421; kfold=T
  
  # Ensure y is numeric
  stopifnot(is.numeric(y))
  
  number_of_folds = length(fold_list)
  
  train_rmse_scores = c()
  test_rmse_scores = c()
  
  all_labels_train = list()
  all_labels_test = list()
  
  feature_importance_values = list()
  
  for (f in seq_along(fold_list)) {
    fold = fold_list[[f]]

    # Split data into train/test
    X_train_fold = X[-fold, ]
    y_train_fold = y[-fold]
    X_test_fold = X[fold, ]
    y_test_fold = y[fold]
    
    # Train control setup
    train_control = trainControl(
      method = "cv",
      number = number_of_folds,
      savePredictions = "final"
    )
    
    # Hyperparameter tuning using caret/ranger
    set.seed(rngseed)
    rf_model = suppressWarnings(train(
      x = X_train_fold,
      y = y_train_fold,
      method = "ranger",
      trControl = train_control,
      tuneGrid = hyper,
      metric = "RMSE"
    ))
    
    # Final model fit with best params
    set.seed(rngseed)
    final_model = randomForest(
      x = X_train_fold,
      y = y_train_fold,
      mtry = rf_model$bestTune$mtry,
      splitrule = rf_model$bestTune$splitrule,
      min.node.size = rf_model$bestTune$min.node.size,
      importance = TRUE
    )
    
    # Predictions
    train_pred = predict(final_model, X_train_fold)
    test_pred = predict(final_model, X_test_fold)
    
    # Compute RMSEs
    train_rmse = rmse(y_train_fold, train_pred)
    test_rmse = rmse(y_test_fold, test_pred)
    
    train_rmse_scores = c(train_rmse_scores, train_rmse)
    test_rmse_scores = c(test_rmse_scores, test_rmse)
    
    # Store predictions for later aggregation
    temp_train = tibble(
      row = c(1:nrow(X))[-fold],
      true_labels = y_train_fold,
      predicted_values = train_pred
    )
    all_labels_train[[length(all_labels_train) + 1]] = temp_train
    
    temp_test = tibble(
      row = c(1:nrow(X))[fold],
      true_labels = y_test_fold,
      predicted_values = test_pred
    )
    all_labels_test[[length(all_labels_test) + 1]] = temp_test
    
    # Save feature importance
    feature_importance_values[[length(feature_importance_values) + 1]] = final_model$importance
  }
  
  # Aggregate results
  avg_result = average_rf_regression(
    train_rmse_scores,
    test_rmse_scores,
    all_labels_train,
    all_labels_test,
    feature_importance_values
  )
  
  return(avg_result)
}

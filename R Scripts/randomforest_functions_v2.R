# Run random forest on each fold

# Average all the folds
average_rf = function(train_auc_scores,test_auc_scores,
                      all_labels_train,all_labels_test,
                      feature_importance_values){
  
  # Combine the AUC values using a bootstrap approach
  boot_auc_train = boot(data = train_auc_scores, statistic = function(d, i) mean(d[i]), R = 1000)
  boot_auc_test = boot(data = test_auc_scores, statistic = function(d, i) mean(d[i]), R = 1000)
  
  avg_auc_train = boot_auc_train$t0
  avg_auc_test = boot_auc_test$t0
  
  # Calculate AUC confidence intervals
  ci_train = boot.ci(boot_auc_train, type = "perc")$percent[4:5] 
  ci_test = boot.ci(boot_auc_test, type = "perc")$percent[4:5]
  
  # ci_train: if all AUC values are equal to 1 in the training dataset, the above fails. Interval is c(1,1).
  # We'll assign this manually instead:
  if(is.null(ci_train)) ci_train = c(1,1)
  
  # Combine the test data labels and predictions (each fold is currently in a separate df)
  # Each sample is only included in a test dataset once, so no need to average.
  test_labels = bind_rows(all_labels_test)
  
  # Combine the train data labels and predictions (each fold is currently in a separate df)
  # Each sample is included 9 times, so need to average.
  train_labels = bind_rows(all_labels_train) %>% 
    group_by(row,true_labels) %>% 
    summarize(predicted_probabilities = mean(predicted_probabilities)) %>% 
    ungroup()
  
  # Calculate the average importance values for each variable.
  # We will use Reduce to add together all values that are at the same row/column coordinate across all datasets in feature_importance_values, then divide by the number of datasets to get the mean.
  mean_feature_importance = Reduce("+", feature_importance_values) / length(feature_importance_values)
  
  # Convert to a data frame and add a Feature column
  importance_df = data.frame(Feature = rownames(mean_feature_importance), mean_feature_importance)
  # Remove row names
  rownames(importance_df) = NULL
  # Sort the results from most to least important. We will use MeanDecreaseGini.
  importance_df = importance_df %>% arrange(-MeanDecreaseGini)
  
  results = list(auc_train = avg_auc_train,
                 auc_test = avg_auc_test,
                 auc_train_ci = ci_train,
                 auc_test_ci = ci_test,
                 test_labels = test_labels,
                 train_labels = train_labels,
                 importance = importance_df)
  return(results)
}


run_rf = function(X, y, fold_list,
                  hyper, rngseed = 421,
                  kfold=T) {
  
  # X = predictors; y = outcome
  # fold_list = folds; hyper = tune_grid;
  # rngseed=421; kfold=T
  
  # Calculate number of total folds
  number_of_folds = length(fold_list)
  
  # Create series of empty vectors. We'll store model outputs here.
  # AUC scores of the ROC curves
  train_auc_scores = c()
  test_auc_scores = c()
  
  # These will contain data frames of the following:
  # For every sample, we will record a) its actual outcome (PD or Control) and 
  # b) the outcome predicted by the random forest model.
  # This will be used to calculate how accurate our model is.
  all_labels_train = list()
  all_labels_test = list()
  
  # Importance values
  feature_importance_values = list()
  
  for (fold in fold_list) {
    if (kfold == F){
      fold = fold_list[[1]]
    }
    
    # Create train and test datasets.
    X_train_fold = X[-fold, ]
    y_train_fold = y[-fold]
    X_test_fold = X[fold, ]
    y_test_fold = y[fold]
    
    # This will tell the RF command how to perform the RF.
    train_control = trainControl(method = "cv", # K-fold cross validation
                                 number = number_of_folds, # 10 folds
                                 classProbs = TRUE, # Predicted class probabilities 
                                 # are returned instead of just class labels.
                                 summaryFunction = twoClassSummary) # compute AUC
    
    # Use hyperparameter tuning to optimize each parameter.
    # Note that optimal settings are chosen based on ROC/AUC - prone to overfitting!
    set.seed(rngseed) # Reproducible randomness
    rf_model = suppressWarnings(train(X_train_fold, y_train_fold, # training dataset
                                      method = "ranger",
                                      trControl = train_control, # Perform tuning
                                      tuneGrid = hyper,
                                      metric = "ROC"
    ))
    
    # Finally, run random forest using the optimal settings
    set.seed(rngseed) # Reproducible randomness
    final_model = randomForest(X_train_fold, y_train_fold, 
                               mtry = rf_model$bestTune$mtry,
                               splitrule = rf_model$bestTune$splitrule,
                               min.node.size = rf_model$bestTune$min.node.size,
                               importance = TRUE)
    
    # Calculate and save model statistics (TRAINING DATA)
    train_pred_proba = predict(final_model, type = "prob")[, 2]
    train_auc = auc(roc(y_train_fold, train_pred_proba))
    train_auc_scores = c(train_auc_scores, train_auc)
    # Save predictions
    temp = tibble(row = c(1:nrow(X))[-fold], # Row from the original training dataset
                  true_labels = y_train_fold, # Actual outcomes
                  predicted_probabilities = train_pred_proba) # Predicted outcomes
    all_labels_train[[length(all_labels_train) + 1]] = temp
    
    # Calculate and save model statistics (TESTING DATA)
    test_pred_proba = predict(final_model, X_test_fold, type = "prob")[, 2]
    test_auc = auc(roc(y_test_fold, test_pred_proba))
    test_auc_scores = c(test_auc_scores, test_auc)
    # Save predictions
    temp = tibble(row = c(1:nrow(X))[fold], # Row from the original testing dataset
                  true_labels = y_test_fold, # Actual outcomes
                  predicted_probabilities = test_pred_proba) # Predicted outcomes
    all_labels_test[[length(all_labels_test) + 1]] = temp
    
    # Save feature importance values
    # Ctrl, PD: higher number = higher in that group
    # MeanDecrease: how does the model quality decrease if the variable is removed? Higher=more important
    feature_importance_values[[length(feature_importance_values) + 1]] = final_model$importance
  }
  
  # Average all the folds
  avg_result = average_rf(train_auc_scores,test_auc_scores,
                          all_labels_train,all_labels_test,
                          feature_importance_values)
  
  return(avg_result)
}
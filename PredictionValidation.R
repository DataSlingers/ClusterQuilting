library(randomForest)
library(mclust)

pred_val <- function(dataset, clust_ests, test_frac = 0.2) {

  if (!is.data.frame(dataset)) {
    dataset <- as.data.frame(dataset)
  }
  
  nn <- nrow(dataset)
  n_test <- max(1, floor(test_frac * nn))
  test_idx <- sample(1:nn, n_test, replace = FALSE)
  train_idx <- setdiff(seq_len(n), test_idx)
  
  x_train <- dataset[train_idx, ]
  x_test  <- dataset[test_idx, ]
  
  y_train <- clust_ests[train_idx]
  y_test  <- clust_ests[test_idx]
  
  # Ensure classification (clusters) not regression
  y_train <- as.factor(y_train)
  
  # Train RF classifier on training set
  rf_fit <- randomForest(x = x_train, y = y_train)
  
  # Predict on test set (class labels)
  y_pred <- predict(rf_fit, newdata = x_test, type = "response")
  
  ari <- adjustedRandIndex(as.numeric(y_pred), 
                           as.numeric(y_test))
  
  return(ari)
}

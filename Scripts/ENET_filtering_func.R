enet_stability_block <- function(Matrix, Definition, alpha_grid, n_repeats, nfolds, threshold) {
  
  X <- as.matrix(Matrix)                                  # Matrix containing samples (rows) and features (columns)
  Y <- as.factor(Definition)                              # Definition (Sample groups, must be equal to nrows)
  
  feature_names <- colnames(X)                            # Feature names
  n_features    <- ncol(X)                                # Total features in block
  selection_matrix <- matrix(0,                           # Create matrix with just 0's
    nrow = n_repeats,                                     # The number of rows is equal to the times the selection will run
    ncol = n_features,                                    # The number of columns will be equal to the number of features
    dimnames = list(NULL, feature_names))                 # Assign names to matrix dims
  
  alpha_record <- numeric(n_repeats)                      # Will set the number of iterations (repeats)
  for(i in 1:n_repeats){
    set.seed(100 + i)                                     # Make every run different by setting different seed
    best_alpha <- NULL                                    # Make starting alpha
    best_error <- Inf                                     # Make starting error
    best_model <- NULL                                    # Make starting model
    
    for (a in alpha_grid){                                # Iterate over every alpha value in the selected range
      cv_fit <- cv.glmnet(                                # k-fold cross-validation for glmnet
        X, Y,                                             # Matrix and definition
        family = "binomial",                              # Since our response is binary
        alpha = a,                                        # Alpha value for iteration
        nfolds = nfolds,                                  # Number of subsets the dataset is divided to train and validate the model
        type.measure = "class")                           # Since is its bionomial
      current_error <- min(cv_fit$cvm)                    # Assign the current mean cross-validated error
      
      if(current_error < best_error){                     # In case the current error is actually better than the previous one
        best_error <- current_error                       # Set new error if it is better
        best_alpha <- a                                   # Save which alpha value yielded the best error
        best_model <- cv_fit                              # Save the best model
      }
    }
    alpha_record[i] <- best_alpha                         # Signal which was the best alpha
    
    beta <- coef(best_model, s = "lambda.1se")            # Extract model coefficients from best model
    selected_idx <- which(beta != 0)[-1]                  # Remove intercept
    if(length(selected_idx) > 0){
      selection_matrix[i, selected_idx] <- 1
    }
  }

  selection_frequency <- colMeans(selection_matrix)       # Mean times the feature was selected
  stable_features <- names(selection_frequency[           # Select features that comply with the threshold
    selection_frequency >= threshold])
  X_filtered <- X[, stable_features, drop = FALSE]        # Keep only stable features
  
  return(list(                                            # Return object
    X_filtered = X_filtered,
    stable_features = stable_features,
    selection_frequency = selection_frequency,
    alpha_distribution = table(alpha_record)
  ))
}
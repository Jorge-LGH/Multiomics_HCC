# Fit possible penalties
penalty_fit <- function(omic_list, penalty_vector){
  library(mixOmics)
  
  # Make SGCCA model
  model_structure <- wrapper.sgcca(                              
    omic_list,                                                 # List of omics
    penalty = penalty_vector,                                  # Penalties evaluated (for each omic)
    scheme = "centroid",                                       # Maximize 
    scale = T)                                                 # Already standardized 
  
  # Extract average variance explained (AVE)
  model_ave <- do.call(rbind, model_structure$AVE$AVE_X)       # AVE retrieved from model
  
  # Count selected features per run
  n_features <- sapply(features, function(x) sum(x != 0))      # Total features selected from each omic block
  
  # Understandable summary
  model_summary <- data.frame(
    omic = names(n_features),                                  # Feature names
    AVE = model_ave[,1],                                       # AVE by model
    n_features = n_features,                                   # Total selected features
    penalty = penalty_vector                                   # Evaluated penalties for all omics
  )
  return(model_summary)
}

# Fit penalties on subsets
subset_sgcca <- function(omic_list, penalty_vector, n_repeats, subset_size){
  # Get penalties for data subsets
  results <- list()                                            # Make empty list
  for(i in 1:n_repeats){                                       # Iterate over n number of times
    set.seed(100+i)                                            # Make sure every iteration is different
    n_samples <- nrow(omic_list[[1]])                          # Get the number of samples
    idx <- idx <- sample(seq_len(n_samples),                   # Sample the features based on selected size
                         size = floor(n_samples * subset_size))
    data_subset <- lapply(omic_list, function(X) X[idx, ])     # Extract data subset
    results[[i]] <- penalty_fit(data_subset, penalty_vector)   # Evaluate penalty fits based on subsets
  }
  
  # Save results
  results_df <- do.call(rbind, results)
  return(results_df)
}

# Plot possible penalties on a heatmap
penalty_heatmap <- function(results_df){
  library(ggplot2)
  # Heatmap
  ggplot(results_df,
         aes(x = penalty, y = omic, fill = AVE)) +
    geom_tile() +
    facet_wrap(~omic) +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(
      title = "SGCCA penalty performance",
      x = "Penalty",
      y = "Omic block",
      fill = "AVE"
    )
}


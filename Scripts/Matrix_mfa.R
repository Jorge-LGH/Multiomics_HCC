# Scaling, centering, and normalize matrices
mfa_normalize <- function(omic_mat){                              # Take omic matrix as input (features x samples)
  t_omic_mat <- t(omic_mat)                                       # Tranpose matrix so it is samples x features
  t_omic_scaled <- scale(t_omic_mat, center = TRUE, scale = TRUE) # Center and Scale (z-score for features) 
  svd_res <- svd(t_omic_scaled)                                   # SVD (computes eigenvalues)
  w <- svd_res$d[1]                                               # First eigen value
  mat_norm <- omic_mat/w                                          # Normalization
  return(mat_norm)
}

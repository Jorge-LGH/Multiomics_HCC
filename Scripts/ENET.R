# This script is intended as a variable-filter for the selected omic data types. In this case, they are
# miRNA, mRNA, and CpG sites. This filtering of data types is based on the ENET and MFA methods.

#-----------------------Load Function------------------------
source("Scripts/Matrix_mfa.R")
source("Scripts/ENET_filtering_func.R")

#-----------------------Load data----------------------------
# Load each omic's block data
exp_data <- read.table("Data/norm_exp_data.tsv",sep=',',row.names=T)     # mRNA expression data
mir_data <- read.table("Data/norm_mir_data.tsv",sep=',',row.names=T)     # miRNA expression data
cpg_data <- read.table("Data/norm_met_data.tsv",sep=',',row.names=T)     # CpG methylation data

#--------------------Centering and scaling----------------
# Based on doi: 10.1093/bib/bbx060 to make the impact of components' variable comparable between data levels independently
# from the number of variables in each block. In this case, it was decided to divide the data blocks by the square root
# of the first eigenvalue of the transpose block multiplied by the original block and divided by the number of samples (MFA)
exp_eigen <- mfa_normalize(exp_data)
mir_eigen <- mfa_normalize(mir_data)
met_eigen <- mfa_normalize(cpg_data)

#-----------------------Elastic Net--------------------------
# Parameters used by my function
definition <- samples_data$sample_type                                   # Define groups (Tumor and control)
grid_def <- seq(0.1, 0.9, by = 0.1)                                      # Range for shrinkage parameter
reps <- 1000                                                             # Number of iterations
n_folds <- 5                                                             # Number of subsets the dataset is divided to train and validate the model
thresh <- 0.8                                                            # Select features that comply with the threshold (percentage of being kept)

# Execute the ENET selection on each omic block
exp_enet <- enet_stability_block(exp_eigen, definition, grid_def, reps, n_folds, thresh)
mir_enet <- enet_stability_block(mir_eigen, definition, grid_def, reps, n_folds, thresh)
met_enet <- enet_stability_block(met_eigen, definition, grid_def, reps, n_folds, thresh)

# Extract filtered matrices
exp_filtered <- exp_enet$X_filtered
mir_filtered <- mir_enet$X_filtered
met_filtered <- met_enet$X_filtered

# Create list with every filtered omic
omic_list <- list(
  CpGs = met_filtered,
  transcripts = exp_filtered,
  miRNAs = mir_filtered)
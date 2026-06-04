# This script is intended as the base for penalty fit and penlaty selection for each omic based on the SGCCA method.
# In this case, the functions made and used for these omics will have different values as to have harsher punishments
# on large omics like CpG and more relaxed conditions for smaller omics like miRNA.

# Load functions
source("Scripts/Penalty_functions.R")

# Load data
omic_list <- readRDS("Data/omic_list.rds")                         # All omics storted as a list of data frames

# Define parameters
penalty_vector <- expand.grid(CpG = seq(0.1, 0.6, by = 0.1),       # All possible combinations of penalties for each omic 
                              mRNA = seq(0.2, 0.8, by = 0.1), 
                              miRNA = seq(0.3, 0.9, by = 0.1))
reps <- 50                                                         # Number of repetitions
sub_size <- 0.7                                                    # Data subset size

# Execute penalty fitting
p_fits <- subset_sgcca(omic_list, penalty_vector, reps, sub_size)  # Every possible penalty combination will be tried on subsets 

# Visualize penalty fits
penalty_heatmap(p_fits)                                            # Heatmap visualization
penalty_scatter(p_fits)                                            # Scatter plot visualization


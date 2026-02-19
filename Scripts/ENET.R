# This script is intended as a variable-filter for the selected omic data types. In this case, they are
# miRNA, mRNA, and CpG sites. This filtering of data types is based on the ENET and MFA methods.

#-----------------------Load Function------------------------
source("Scripts/Matrix_mfa.R)

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


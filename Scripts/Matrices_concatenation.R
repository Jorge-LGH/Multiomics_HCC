# This script is intended as a concatenation for the selected omic data types. In this case, they are
# miRNA, mRNA, and CpG sites. This concatenation of data types is based on the SGCCA method.

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
# of the first eigenvalue of the transpose block multiplied by the original block and divided by the number of samples
exp_eigen <- mfa_normalize(exp_data)
mir_eigen <- mfa_normalize(mir_data)
met_eigen <- mfa_normalize(cpg_data)

#-----------------------Concatenate--------------------------
# Check if matrices' column names are equally sorted in every information block. While it is most likely
# that the data is sorted in the same way due to the pre-processing for each omic, making sure doesn't hurt
col_names <- colnames(exp_eigen)
if((col_names == colnames(mir_eigen)) & (col_names == colnames(met_eigen))){
  print("Everything matches, concatenating...")
  concatenated <- rbind(exp_data, mir_data, cpg_data) # Concatenate as to have ALL samples as columns and ALL features as rows
}else{
  print("The column names across the data are not consistent, please check it.")
}

#--------------------Save concatenated matrix-------------
write.table(concatenated,"3_Data/concatenated_matrix.tsv",sep=',',row.names=T) 

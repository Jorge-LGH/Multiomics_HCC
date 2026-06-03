# This script is meant as the methylation values pre-processing guideline

#--------------------Load libraries-----------------------
library(TCGAbiolinks)                                                # Version: 2.36.0
library(SummarizedExperiment)                                        # Version: 1.48.1
library(tidyverse)                                                   # Version: 2.0.0
library(methyLImp2)                                                  # Version: 1.2.0
library(BiocParallel)                                                # Version: 1.40.2
library(sesame)                                                      # Version: 1.24.0
library(NOISeq)                                                      # Version: 2.52.0
library(sva)                                                         # Version: 3.54.0
library(ChAMPdata)                                                   # Version: 2.40.0
library(minfi)                                                       # Version: 1.52.1
library(jaffelab)                                                    # Version: 0.99.34

#--------------------Load object--------------------------
# To understand where this object came from, check the 1_get_data.R script
samples_data <- read.table("Data/samples_data.tsv", header = T, sep='\t')

#--------------------Prepare data-------------------------
# Get miRNA expression values
met_query <-  GDCquery(project = "TCGA-LIHC",                        # Liver hepatocellular carcinoma project
                       data.category = "DNA Methylation",            # DNA methylation data
                       platform = "Illumina Human Methylation 450",  # CpG detection platform
                       data.type = "Methylation Beta Value",         # Data type
                       barcode = samples_data$Barcode)               # Barcode

# The resulting object already has the methylation beta values by sample and probe as row names
met_data <- GDCprepare(met_query,                                    # Query object
                       directory = "Data/GDCdata/",                  # Directory where files are stored
                       summarizedExperiment = T)                     # Create summarized experiment (includes little metadata)

#--------------------QC-----------------------------------
# Starting object qualities
dim(met_data)                                                        # 485,577 probes across 407 samples
sum(is.na(assay(met_data)))                                          # 31,589,946 NA's (need to impute values after some filtering)
sum(colSums(assay(met_data), na.rm = T) == 0)                        # No empty columns (Should worry if it wasn't the case)
sum(rowSums(assay(met_data), na.rm = T) == 0)                        # There are 64,213 rows with only 0's if we drop the NA's

# Removing MASKED probes (Check your object for the MASK status)
## This will be removing probes with SNPs and with low uniqueness, meaning they map for multiple sites
masked_cols <- grep("MASK", names(mcols(rowRanges(met_data))),       # Get the columns that have the MASK information
                    value = TRUE)

for(mask in masked_cols){                                            # Remove every single probe that has a True value in any MASK criteria
  met_data <- met_data[which(                                        # 296,864 probes remain (29-04-2026)
    mcols(rowRanges(met_data))[[mask]] == F),]
}

# Remove probes that only appear on the Y chromosome due to sex unbalance or ambiguous mapping
met_data <- met_data[which(                                          # 296,664 probes remain (29-04-2026)
  as.character(seqnames(rowRanges(met_data))) != "chrY"), ]

met_data <- met_data[which(                                          # No ambiguously mapped probes (29-04-2026)
  as.character(seqnames(rowRanges(met_data))) != "*"), ]

# Filter data with just 0's aside from NAs
met_data <- met_data[rowSums(met_data@assays@data@listData[[1]],     # 296,624 probes remain (29-04-2026)
                             na.rm=T) != 0, ,drop = FALSE]

# Filter probes with missing values in 20% or more of the samples
met_data <- met_data[which(!rowSums(                                 # 281,239 probes remain (29-04-2026)
  is.na(met_data@assays@data@listData[[1]])) >
    (ncol(met_data@assays@data@listData[[1]])*0.2)),]

#--------------------Methylation data imputation----------
# Impute missing values with regression based method. See: https://doi.org/10.1186/s12859-020-03592-5
# Separating complete probes from incomplete ones will make the process more efficient and less time-consuming
complete_probes <- met_data[which(                                   # 224,656 complete probes (29-04-2026)
  rowSums(is.na(assay(met_data))) == 0),]

incomplete_probes <- met_data[which(                                 # 56,583 incomplete probes (29-04-2026)
  rowSums(is.na(assay(met_data))) != 0),]
  
# Actual imputation
met_imputed <- methyLImp2(incomplete_probes,
                          type    = "450K",
                          groups  = samples_data$sample_type,
                          BPPARAM = MulticoreParam(workers = 2))

# Merge with complete probes
met_complete <- rbind(complete_probes, met_imputed)                  # Same 281,239 probes with imputations (29-04-2026)

write.table(assay(met_complete), "Data/met_imputed.tsv",
            sep='/t',row.names=T)                                    # Save object with the combined imputations

# Plot beta values for tumor samples and adjacent tissues
## By each sample
png("Figures/CpG/cpg_beta_distribution.png", width=1000) 
densityPlot(assay(met_complete), sampGroups = met_complete@colData$sample_type)
dev.off()

## By sample type
tumor_beta <- assay(met_complete[,                                   # Beta values for tumor samples
  which(met_complete@colData$sample_type == "Primary Tumor")])

control_beta <- assay(met_complete[,                                 # Beta values for adjacent tissues
  which(met_complete@colData$sample_type != "Primary Tumor")])

beta_matrix <- data.frame(betas = c(as.vector(tumor_beta),           # Set beta values as matrix
                                    as.vector(control_beta)), 
                          type = c(rep("Tumor", ncol(tumor_beta)), 
                                   rep("Control", 
                                       ncol(control_beta))))

png("Figures/CpG/cpg_beta_distribution_type.png", width=1000)
ggplot(beta_matrix, aes(x = betas, color = type)) +                  # Actual plot
  geom_density() + 
  labs(x = "Beta values", y = "Density", colour = "Sample type")
dev.off()

#--------------------B-values to M-values-----------------
# See https://doi.org/10.1186/1471-2105-11-587 for decision basis
## Cap beta values to avoid ±Inf in M-value conversion
beta_mat <- assay(met_complete)
beta_mat[beta_mat == 0] <- 1e-6                                      # Make it so it is basically 0 but not really
beta_mat[beta_mat == 1] <- 1-(1e-6)                                  # Make it so it is basically 1 but not really

## Actually converting to M-values
m_values <- BetaValueToMValue(beta_mat)

#--------------------Check batch effect-------------------
noiseqData <- readData(data = m_values, factor = samples_data)       # Create noiseq object
myPCA <- NOISeq::dat(noiseqData, type ="PCA",norm =T, logtransf =T)  # Perform PCA to watch for batch effects, already logtransformed so T

# The samples partially aggregate by cancerous and control
# PC1 explains 28% and PC2 explains 8%
png("Figures/CpG/cpg_batch_effect.png", width = 1000)                # Plot batch effect with PCA
explo.plot(myPCA, factor = "sample_type")
dev.off()

#--------------------Remove btach effect------------------
# I'll be using SVA (Surrogate Variable Analysis) to remove unknown batch effects
com_mod <- model.matrix(~sample_type, data = samples_data)           # Complete model (adjustment and variable of interest cancer or control), it tells which variation we want to keep
nul_mod <- model.matrix(~1, data = samples_data)                     # Null model, only include intercept. It explains nothing biologically and checks for technical variation

## The number of calculated SV's by the method (both "leek" and "be") were to high and risked overfitting.
## I used an elbow plot to check for the "optimal" number of SV's compared to the ones returned by the function
H <- com_mod %*% solve(t(com_mod) %*% com_mod) %*% t(com_mod)        # Compute residual matrix (meaning only technical variation)
res <- m_values - t(H %*% t(m_values))                               # Remove the technical variation
svd_res <- svd(res)                                                  # Singular value decomposition
variance_explained <- svd_res$d^2 / sum(svd_res$d^2) * 100           # How much variance per component
num_la_f <- num.sv(m_values, com_mod, method = "be")                 # Identify number of latent factors to estimate (SVs) based on the "be" method

x_limit <- max(50, num_la_f + 5)                                     # Set the matrix to plot and check for variance and the calculated SV's 
scree_df <- data.frame(SV=1:min(x_limit,length(variance_explained)),
  Var=variance_explained[1:min(x_limit,length(variance_explained))])

## Actual plot
png("Figures/CpG/cpg_scree_plot.png")
ggplot(scree_df, aes(x = SV, y = Var)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = num_la_f,
             linetype = "dashed",
             color = "red") +
  annotate("text",
           x = num_la_f + 1,
           y = max(scree_df$Var) * 0.95,
           label = paste0("BE estimate\n(n = ", num_la_f, ")"),
           color = "red",
           hjust = 0,
           size = 3.5) +
  labs(x = "Surrogate Variable",
       y = "Variance explained (%)",
       title = "Scree plot — residual M-value matrix") +
  theme_classic()
dev.off()

# The plot shows that the variance explained by the SV's starts stabilizing around 7-10 SV's. Therefore, to make the analysis as
# robust as possible, instead of using the 40 SV's recommended by the 'num.sv' function, we'll be using 10 SV's. 
num_la_f <- 10

## Remove batch effect 
sva_obj <- sva(m_values, com_mod, nul_mod, n.sv = num_la_f)          # Estimate surrogate variables' values
clean_m <- cleaningY(y = m_values,                                   # Adjust the m values by regressing the data with the SV's                      
                     mod = cbind(com_mod, sva_obj$sv),
                     P   = 2)                

write.table(clean_m, "Data/clean_m.tsv", sep='/t', row.names=T)      # Save object with clean M-values

noiseqData_clean <- readData(data = clean_m, factor = samples_data)  # Plot again after batch effect correction
myPCA_clean <- NOISeq::dat(noiseqData_clean,type="PCA",
                           norm=T,logtransf =T)
png("Figures/CpG/cpg_batch_effect_after_correction.png")
explo.plot(myPCA_clean, factor = "sample_type")
dev.off()

#--------------------Differential methylation-------------
# First, it is important to take the adjacent tissue as the control for the differential anaylsis
samples_data$sample_type <- relevel(factor(samples_data$sample_type),
                                    ref = "Solid Tissue Normal")

# Now we have to build a design matrix
des_matrix <- model.matrix(~sample_type, data = samples_data)

# We fit the linear model to the M-values
fit <- lmFit(clean_m, des_matrix)

# Empirical Bayes moderation
fit <- eBayes(fit)

# Check which change between conditions
results <- decideTests(fit,
                       coef           = "sample_typePrimary Tumor",
                       adjust.method  = "BH",
                       method         = "separate",
                       lfc            = 1)   
summary(results)
#        (Intercept) sample_typePrimary Tumor
# Down        130798                    30508
# NotSig       41427                   245819
# Up          109014                     4912

# Extract results' info
diff_meth <- topTable(fit,
                      coef = "sample_typePrimary Tumor",
                      number = Inf,                                  # Return all probes
                      adjust.method = "BH")                          # Benjamini-Hochberg FDR

# Change the names for better interpretation
diff_meth$Methylation <- factor(
  results[rownames(diff_meth), "sample_typePrimary Tumor"],
  levels = c(-1, 0, 1),
  labels = c("Hypomethylated", "Unchanged", "Hypermethylated"))

# Convert m clean values to beta values just for biological interpretation
beta_corrected <- MValueToBetaValue(clean_m) 

# Get samples types
tumor_idx <- which(samples_data$sample_type == "Primary Tumor")
control_idx <- which(samples_data$sample_type == "Solid Tissue Normal")

# Get beta value means for each probe on each condition
mean_beta_tumor <- rowMeans(beta_corrected[, tumor_idx], na.rm = T)
mean_beta_control <- rowMeans(beta_corrected[, control_idx], na.rm = T)

# Match probe order to diff_meth before subtracting
diff_meth$delta_beta <- mean_beta_tumor[rownames(diff_meth)] -
                        mean_beta_control[rownames(diff_meth)]
# Volcano plot
png("Figures/CpG/cpg_diff_meth_volcano.png", width = 1000)
ggplot(diff_meth, aes(x = delta_beta, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Methylation), size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c("Hypermethylated" = "firebrick3",
                                "Hypomethylated"  = "dodgerblue3",
                                "Unchanged"       = "gray50")) +
  geom_hline(yintercept = -log10(0.05),
             linetype   = "dashed",
             color      = "black") +
  geom_vline(xintercept = c(-0.2, 0.2),                              # Delta beta reference lines
             linetype   = "dashed",                                  # for biological context only
             color      = "black") +                                 # not used for classification
  labs(x     = "Delta Beta (Tumor - Normal)",
       y     = "-log10(adj. p-value)",
       title = "Differential Methylation") +
  theme_bw()
dev.off()


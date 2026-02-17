# This script is meant as the miRNA pre-processing guideline. Normalization, feature selection, and other
# processing techniques will be applied to ensure the miRNA expression values  are valid across every sample.

#--------------------Load libraries-----------------------
library(TCGAbiolinks)              # Version: 2.36.0
library(SummarizedExperiment)      # Version: 1.48.1
library(tidyverse)                 # Version: 2.0.0
library(biomaRt)                   # Version: 2.64.0
library(NOISeq)                    # Version: 2.52.0
library(edgeR)                     # Version: 4.6.2
library(EDASeq)                    # Version: 2.42.0
library(cqn)                       # Version: 1.45.0
library(DESeq2)                    # Version: 1.48.1
library(ggrepel)                   # Version: 0.9.6

#--------------------Load object--------------------------
# To understand where this object came from, check the 1_get_data.R script
samples_data <- read.table("Data/samples_data.tsv", header = T, sep='\t')

#--------------------Prepare data-------------------------
# Get miRNA expression values
mir_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "miRNA Expression Quantification", # Gene expression 
                      barcode = samples_data$Barcode)                # Barcode
mir_data <- GDCprepare(mir_query,                                    # Query object
                       directory = "Data/GDCdata/",                  # Directory where files are stored
                       summarizedExperiment = F)                     # Do not create a summarized experiment

# Set miRNA-ID as row names
rownames(mir_data) <- mir_data$miRNA_ID

# Keep only raw counts for each sample
keep_cols <- colnames(mir_data)[sapply(colnames(mir_data), function(col){
  strsplit(col, "_")[[1]][1] == "read"})]
mir_data <- mir_data %>% dplyr::select(keep_cols)

# Change column names to barcode
colnames(mir_data) <- samples_data$Barcode

# Remove genes with no expression among every sample
dim(mir_data)                                                        # 1,881 miRNA genes before filtering (11-08-2025)
mir_data <- mir_data[rowSums(mir_data) != 0, , drop = FALSE]         # 1,526 genes remain (11-08-2025)

#--------------------Annotation data----------------------
# Get annotation data
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 114)  # Data set to use
ann_data <- getBM(attributes = c("mirbase_id",                                   # mirbase ID
                                 "percentage_gene_gc_content",                   # miRNA GC content
                                 "start_position",
                                 "chromosome_name",
                                 "end_position"),
                  mart=mart)
ann_data$length <- ann_data$end_position - ann_data$start_position          # Add miRNA length
ann_data <- ann_data[which(ann_data$mirbase_id %in% rownames(mir_data)),]   # Keep annotation of present miRNAs only (1,609 genes 12-08-2025)
ann_data <- ann_data[which(!duplicated(ann_data$mirbase_id)),]              # Remove duplicated miRNAs (1,499 miRNAs remain 12-06-2025)
mir_data <- mir_data[which(rownames(mir_data) %in% ann_data$mirbase_id),]   # Keep only miRNAs with annotation (Same 1,499 miRNAs)

#--------------------Check for biases---------------------
# Create noiseq object for it to be compatible with selected workflow
noiseqData <- NOISeq::readData(mir_data,
                               factors = samples_data,
                               gc = ann_data[, c("mirbase_id", "percentage_gene_gc_content")],
                               length = ann_data[, c("mirbase_id", "length")])

counts_data <- dat(noiseqData, type = "countsbio", factor = "sample_type")

# Expression values for cancer samples and controls
png("Figures/miRNA/mir_vals_before_norm.png",width=1000)
explo.plot(counts_data, plottype = "boxplot")
dev.off()

# Expression values in CPM for cancer samples and controls
png("Figures/miRNA/mir_vals_bar_before_norm.png",width=1000)
explo.plot(counts_data, plottype = "barplot")
dev.off()

# Visualize low CPM values
cpm_hist <- ggplot(mir_data, aes(x = rowMeans(cpm(mir_data, log = T)))) + 
  geom_histogram(colour = "blue", fill = "lightblue") + xlab("CPM") + 
  ylab("Genes") + 
  theme_classic() + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red")
cpm_hist
ggsave("Figures/miRNA/mir_cpm_before_norm.png", plot = cpm_hist)           # Save plot 
sum(rowMeans(cpm(mir_data, log = T))>0)/nrow(mir_data)*100                 # ~27% genes have CPM>0

# Check for transcript composition bias
cd_data <- dat(noiseqData, type = "cd", norm = F)                          # Reference sample is: TCGA-ZS-A9CG-01A-11R-A37K-07
table(cd_data@dat$DiagnosticTest[, "Diagnostic Test"])                     # 405 Failed and 1 Passed (12-08-2025)

png("Figures/miRNA/mir_transcript_bias_before_norm.png",width=1000)
explo.plot(cd_data, samples = sample(1:ncol(mir_data),10))                 # The result clearly shows composition bias
dev.off()

# Check for GC bias
# The results show a large effect on expression values based on the GC content.The fit for the cancer samples
# is 90.77% with a p value of 0.046. The fit for the controls is of 81.85% and a p value of 0.12
gc_content <- dat(noiseqData, type = "GCbias", k = 0, factor = "sample_type")
png("Figures/miRNA/mir_gc_bias_before_norm.png",width=1000)
explo.plot(gc_content)  
dev.off()

# Check for length bias
# The results show a little effect on expression values based on the length. The fit for the cancer samples is 47.67% 
# with a p value of 0.53. The fit for the controls is of 25% and a p value of 0.8
len_bias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "sample_type")
png("Figures/miRNA/mir_length_bias_before_norm.png",width=1000)
explo.plot(len_bias)
dev.off()

# Check for batch effect
# The samples do not really aggregate by cancerous and control samples
# PC1 explains 15% and PC2 explains 14%
myPCA = dat(noiseqData, type = "PCA", norm = F, logtransf = F)
png("Figures/miRNA/mir_batch_effect_before_norm.png",width=1000)
explo.plot(myPCA, factor = "sample_type")
dev.off()

#--------------------Solve biases-------------------------
# Filter genes with low counts (CPM) < 0
mir_data <- filtered.data(mir_data,                   # Data to filter    
                          factor = "sample_type",     # Samples' conditions
                          norm = F,                   # If normalized or not
                          method = 1,                 # Filtering method (CPM)
                          cpm = 0,                    # CPM threshold
                          p.adj = "fdr")              # Correction method
dim(mir_data)                                         # 218 genes remain

# Filter annotation data
ann_data <- ann_data[which(ann_data$mirbase_id %in% rownames(mir_data)),]

# Solve GC and length bias using cqn (conditional quantile normalization)
counts <- as.matrix(mir_data)                         # Raw counts 
lengths <- ann_data$length                            # Gene lengths 
gc <- ann_data$percentage_gene_gc_content             # GC content 

# Apply CQN 
cqn_res <- cqn(counts = counts,                       # Expression counts
               lengths = lengths,                     # Gene length
               x = gc,                                # Covariate to remove
               sizeFactors = colSums(counts),         # Library sizes
               verbose = T)

# Get normalized expression values
mir_data <- cqn_res$y + cqn_res$offset

# cqn NOISeq object to perform batch correction
cqn_noiseq <- NOISeq::readData(data = mir_data,
                               factors = samples_data,
                               gc = ann_data[, c("mirbase_id", "percentage_gene_gc_content")],
                               length = ann_data[, c("mirbase_id", "length")])

# Batch effect correction
cqn_arsyn <- ARSyNseq(cqn_noiseq, 
                      factor = "sample_type", 
                      batch = F,  
                      norm = "n", 
                      logtransf = T)

# Check batch effect removal visually
myPCA <- dat(cqn_arsyn, type = "PCA", norm = T, logtransf = T)       # Perform PCA
png("Figures/miRNA/mir_batch_effect_after_norm.png",width=1000)
explo.plot(myPCA, factor = "sample_type")                            # PCA1 = 16%, PCA2 = 2%
dev.off()

# New NOISeq object to check for GC and length bias
new_noiseq <- NOISeq::readData(cqn_arsyn,
                               factors = samples_data,
                               gc = ann_data[, c("mirbase_id", "percentage_gene_gc_content")],
                               length = ann_data[, c("mirbase_id", "length")])

new_counts_data <- dat(new_noiseq, type = "countsbio",               # Will check expression values
                       factor = "sample_type", norm = T)
png("Figures/miRNA/mir_val_bar_after_norm.png",width=1000)
explo.plot(new_counts_data, plottype = "boxplot")                    # Expression values
dev.off()

# Check for GC and length bias
# I had to check it manually since the EDASeq functions need more features. Since there are only 218 in this case, I can do it this way
expr_mean <- rowMeans(cqn_arsyn@assayData$exprs)

bias_df <- data.frame(expr = expr_mean,
                      GC = ann_data$percentage_gene_gc_content,
                      length = ann_data$length,
                      mirbase_id = ann_data$mirbase_id)

gc_after <- ggplot(bias_df, aes(x = GC, y = expr)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  theme_classic() +
  labs(title = "GC content bias after CQN",
       x = "GC% (miRNA)",
       y = "Normalized expression")

gc_after

len_after <- ggplot(bias_df, aes(x = length, y = expr)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  theme_classic() +
  labs(title = "Length bias after CQN",
       x = "miRNA length (nt)",
       y = "Normalized expression")

len_after

ggsave("Figures/miRNA/mir_gc_bias_after_norm.png", plot = gc_after) 
ggsave("Figures/miRNA/mir_length_bias_after_norm.png", plot = len_after) 

#--------------------Save expression matrices-------------
write.table(new_noiseq@assayData$exprs,"Data/norm_mir_data.tsv",sep=',',row.names=T) # Normalized miRNA expression values

#--------------------Differential expression analysis-----
# Differential expression analysis performed with DESeq2
# Create DESeq object
diff_des <- DESeqDataSetFromMatrix(countData = counts,            # Un-normalized counts matrix
                                   colData = samples_data,        # Meta data
                                   design = ~ sample_type)        # Variable of interest
diff_des$sample_type <- relevel(diff_des$sample_type,             # Set the control samples as the reference
                                ref = "Solid Tissue Normal")

# Perform differential expression analysis
diff_anl <- DESeq(diff_des)

# Extract results
diff_res <- results(diff_anl)

# Make plot
diff_res <- as.data.frame(diff_res)

## Assign labels if unchanged, up-regulated or down-regulated 
diff_res <- diff_res %>% mutate(Expression = case_when(log2FoldChange >= 1 & pvalue <= 0.05 ~ "Up-regulated",
                                                       log2FoldChange <= -1 & pvalue <= 0.05 ~ "Down-regulated",
                                                       TRUE ~ "Unchanged"))
# Check for differential expression results
table(diff_res$Expression)
# Down-regulated      Unchanged   Up-regulated 
# 7                   181         30 

## Extract ID from 50 top and 50 bottom differentially expressed
bot_50_genes <- diff_res[which(
  diff_res$Expression == "Down-regulated"),][order(diff_res[which(
    diff_res$Expression == "Down-regulated"),]$log2FoldChange, decreasing = F), ] %>% head(.,50)

top_50_genes <- diff_res[which(
  diff_res$Expression == "Up-regulated"),][order(diff_res[which(
    diff_res$Expression == "Up-regulated"),]$log2FoldChange, decreasing = T), ] %>% head(.,50)

## Plot itself
png("Figures/miRNA/mir_diff_plot.png",width=1000)
ggplot(diff_res, aes(log2FoldChange, -log(pvalue,10))) +
  geom_point(aes(color = Expression), size = 1.5, alpha=.7) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  labs(x = "log_2(FC)", y = "-log10(p-Value)", title = "Differential Expression") + 
  geom_hline(yintercept=0.30103, linetype="dashed", color = "black") +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color = "black") +
  theme_bw() +
  geom_label_repel(data=subset(diff_res,rownames(diff_res) %in% 
                                 rownames(top_50_genes[1:10,]) | 
                                 rownames(diff_res) %in% rownames(bot_50_genes[1:10,])),
                   aes(log2FoldChange, -log(pvalue,10), 
                       label = rownames(subset(diff_res,rownames(diff_res) %in%
                                                 rownames(top_50_genes[1:10,]) |
                                                 rownames(diff_res) %in% rownames(bot_50_genes[1:10,])))),
                   arrow = arrow(length = unit(0.02, "npc")),
                   nudge_x = 0.2)
dev.off()

#--------------------Save objects-------------------------
write.table(diff_res,"Data/mir_diff.tsv",sep=',',row.names=T)

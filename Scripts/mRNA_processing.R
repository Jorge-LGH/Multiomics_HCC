# This script is meant as the mRNA pre-processing guideline. Normalization, feature selection, and other
# processing techniques will be applied to ensure the mRNA expression values across every sample.

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
# Get mRNA expression values
exp_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "Gene Expression Quantification",  # Gene expression 
                      workflow.type = "STAR - Counts",               # How reads are set as counts per gene
                      barcode = samples_data$Barcode)                # Barcode
exp_data <- GDCprepare(exp_query,                                    # Query object
                       directory = "Data/GDCdata/",                  # Directory where files are stored
                       summarizedExperiment = F)                     # Do not create a summarized experiment

# Keep only protein coding RNAs
exp_data <- exp_data[which(exp_data$gene_type == "protein_coding"),] # 19,962 genes identified

# Remove rows with no transcript reads 
exp_data <- exp_data[rowSums(exp_data[,-c("gene_id",                 # 19,507 genes remain
                                          "gene_name", 
                                          "gene_type")]) != 0,]

# Set ensembl gene id as row names
rownames(exp_data) <- exp_data$gene_id

# Remove ensembl id column as well as the gene name and gene type
exp_data <- dplyr::select(exp_data, -c("gene_id", "gene_name", "gene_type"))

# Keep only unstranded reads
keep_cols <- colnames(exp_data)[sapply(colnames(exp_data), function(col){
  strsplit(col, "_")[[1]][1] == "unstranded"})]
exp_data <- exp_data %>% dplyr::select(keep_cols)

# Remove the "unstranded_" part of the columns' names
colnames(exp_data) <- unlist(strsplit(colnames(exp_data), "_"))[
  unlist(strsplit(colnames(exp_data), "_")) != "unstranded"]

# Remove the transcript version and only keep ensembl gene id
rownames(exp_data) <- sapply(strsplit(rownames(exp_data),".",fixed=T),
                             function(x) x[1])

# Remove rows with no transcript reads again
exp_data <- as.data.frame(exp_data,row.names = rownames(exp_data))
exp_data <- exp_data[rowSums(exp_data) != 0, , drop = FALSE]

#--------------------Annotation data----------------------
# Get annotation data
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", version = 114)
ann_data <- getBM(attributes = c("ensembl_gene_id", 
                                 "percentage_gene_gc_content", 
                                 "gene_biotype",
                                 "start_position",
                                 "end_position",
                                 "hgnc_id",
                                 "chromosome_name",
                                 "hgnc_symbol"),
                  filters = "ensembl_gene_id", 
                  values=rownames(exp_data), 
                  mart=mart)
ann_data <- ann_data[which(!duplicated(ann_data$ensembl_gene_id)),]            # Remove duplicated annotations
ann_data$length <- abs(ann_data$end_position - ann_data$start_position)        # Add length

# Remove non protein coding genes from annotation
ann_data <- ann_data[which(ann_data$gene_biotype == "protein_coding"),]        # 19,350 genes remain

# Remove transcripts with no annotation data
exp_data <- exp_data[which(rownames(exp_data) %in% ann_data$ensembl_gene_id),] # Still 19,350

#--------------------Check for biases---------------------
# Create noiseq object for it to be compatible with selected workflow
noiseqData <- NOISeq::readData(exp_data,
                               factors = samples_data,
                               gc = ann_data[, c("ensembl_gene_id", "percentage_gene_gc_content")],
                               biotype = ann_data[, c("ensembl_gene_id", "gene_biotype")],
                               length = ann_data[, c("ensembl_gene_id", "length")])

counts_data <- dat(noiseqData, type = "countsbio", factor = "sample_type")

# Expression values for cancer samples and controls
png("Figures/mRNA/exp_vals_before_norm.png",width=1000)
explo.plot(counts_data, plottype = "boxplot")
dev.off()

# Expression values in CPM for cancer samples and controls
png("Figures/mRNA/exp_vals_bar_before_norm.png",width=1000)
explo.plot(counts_data, plottype = "barplot")
dev.off()

# Visualize low CPM values
cpm_hist <- ggplot(exp_data, aes(x = rowMeans(cpm(exp_data, log = T)))) + 
  geom_histogram(colour = "blue", fill = "lightblue") + xlab("CPM") + 
  ylab("Genes") + 
  theme_classic() + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red")
cpm_hist
ggsave("Figures/mRNA/exp_cpm_before_norm.png", plot = cpm_hist)            # Save plot 
sum(rowMeans(cpm(exp_data, log = T))>0)/nrow(exp_data)*100                 # ~64% genes have CPM>0

# Check for transcript composition bias
cd_data <- dat(noiseqData, type = "cd", norm = F)                          # Reference sample is: TCGA-ZS-A9CG-01A-11R-A37K-07
table(cd_data@dat$DiagnosticTest[, "Diagnostic Test"])                     # 374 Failed and 32 Passed (01-08-2025)

png("Figures/mRNA/exp_transcript_bias_before_norm.png",width=1000)
explo.plot(cd_data, samples = sample(1:ncol(exp_data),10))                 # The result clearly shows composition bias
dev.off()

# Check for GC bias
# The results show a little effect on expression values based on the GC content.The fit for the cancer samples
# is 51.06% with a p value of 7.5e-09. The fit for the controls is of 37.71% and a p value of 4e-05
gc_content <- dat(noiseqData, type = "GCbias", k = 0, factor = "sample_type")
png("Figures/mRNA/exp_gc_bias_before_norm.png",width=1000)
explo.plot(gc_content)  
dev.off()

# Check for length bias
# The results show a little effect on expression values based on the length. The fit for the cancer samples is 53.07% 
# with a p value of 5.2e-10. The fit for the controls is of 52.58% and a p value of 9.7e-07
len_bias <- dat(noiseqData, k = 0, type = "lengthbias", factor = "sample_type")
png("Figures/mRNA/exp_length_bias_before_norm.png",width=1000)
explo.plot(len_bias)
dev.off()

# Check for batch effect
# The samples do aggregate mostly by cancerous and control samples
# PC1 explains 20% and PC2 explains 8%
myPCA = dat(noiseqData, type = "PCA", norm = F, logtransf = F)
png("Figures/mRNA/exp_batch_effect_before_norm.png",width=1000)
explo.plot(myPCA, factor = "sample_type")
dev.off()

#--------------------Solve biases-------------------------
# Filter genes with low counts (CPM) < 0
exp_data <- filtered.data(exp_data,                   # Data to filter    
                          factor = "sample_type",     # Samples' conditions
                          norm = F,                   # If normalized or not
                          method = 1,                 # Filtering method (CPM)
                          cpm = 0,                    # CPM threshold
                          p.adj = "fdr")              # Correction method
dim(exp_data)                                         # 9,801 genes remain

# Filter annotation data
ann_data <- ann_data[which(ann_data$ensembl_gene_id %in% rownames(exp_data)),]

# Column names must match
colnames(exp_data) <- samples_data$Barcode

# Solve GC and length bias using cqn (conditional quantile normalization)
counts <- as.matrix(exp_data)                         # Raw counts 
lengths <- ann_data$length                            # Gene lengths 
gc <- ann_data$percentage_gene_gc_content             # GC content 

# Apply CQN 
cqn_res <- cqn(counts = counts,                       # Expression counts
               lengths = lengths,                     # Gene length
               x = gc,                                # Covariate to remove
               sizeFactors = colSums(counts),         # Library sizes
               verbose = T)

# Get normalized expression values
exp_data <- cqn_res$y + cqn_res$offset

# cqn NOISeq object to perform batch correction
cqn_noiseq <- NOISeq::readData(data = exp_data,
                               factors = samples_data,
                               gc = ann_data[, c("ensembl_gene_id", "percentage_gene_gc_content")],
                               length = ann_data[, c("ensembl_gene_id", "length")])

# Batch effect correction
cqn_arsyn <- ARSyNseq(cqn_noiseq, 
                      factor = "sample_type", 
                      batch = F,  
                      norm = "n", 
                      logtransf = T)

# Check batch effect removal visually
myPCA <- dat(cqn_arsyn, type = "PCA", norm = T, logtransf = T)       # Perform PCA
png("Figures/mRNA/exp_batch_effect_after_norm.png",width=1000)
explo.plot(myPCA, factor = "sample_type")                            # PCA1 = 13%, PCA2 = 1%
dev.off()

# New NOISeq object to check for GC and length bias
new_noiseq <- NOISeq::readData(cqn_arsyn,
                               factors = samples_data,
                               gc = ann_data[, c("ensembl_gene_id", "percentage_gene_gc_content")],
                               biotype = ann_data[, c("ensembl_gene_id", "gene_biotype")],
                               length = ann_data[, c("ensembl_gene_id", "length")])

new_counts_data <- dat(new_noiseq, type = "countsbio",               # Will check expression values
                       factor = "sample_type", norm = T)
png("Figures/mRNA/exp_val_bar_after_norm.png",width=1000)
explo.plot(new_counts_data, plottype = "boxplot")                    # Expression values
dev.off()

# Check for GC content bias
new_gc_content <- dat(new_noiseq, type = "GCbias", k = 0, factor = "sample_type", norm = T)
png("Figures/mRNA/exp_gc_bias_after_norm.png",width=1000)
explo.plot(new_gc_content)  
dev.off()

# Check for length bias
new_len_bias <- dat(new_noiseq, k = 0, type = "lengthbias", factor = "sample_type", norm = T)
png("Figures/mRNA/exp_length_bias_after_norm.png",width=1000)
explo.plot(new_len_bias)
dev.off()

#--------------------Save expression matrices-------------
write.table(new_noiseq@assayData$exprs,"Data/norm_exp_data.tsv",sep=',',row.names=T) # Normalized mRNA expression values

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

## Extract ID from 50 top and 50 bottom differentially expressed
bot_50_genes <- diff_res[which(
  diff_res$Expression == "Down-regulated"),][order(diff_res[which(
    diff_res$Expression == "Down-regulated"),]$log2FoldChange, decreasing = F), ] %>% head(.,50)

top_50_genes <- diff_res[which(
  diff_res$Expression == "Up-regulated"),][order(diff_res[which(
    diff_res$Expression == "Up-regulated"),]$log2FoldChange, decreasing = T), ] %>% head(.,50)

## Plot itself
png("Figures/mRNA/exp_diff_plot.png",width=1000)
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
                   nudge_x = 0.5)
dev.off()

#--------------------Save objects-------------------------
write.table(diff_res,"Data/exp_diff_exp.tsv",sep=',',row.names=T)

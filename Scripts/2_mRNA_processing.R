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
## Get mRNA expression values
exp_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "Gene Expression Quantification",  # Gene expression 
                      workflow.type = "STAR - Counts",               # How reads are set as counts per gene
                      barcode = samples_data$barcode)                # Barcode
exp_data <- GDCprepare(exp_query,                                    # Query object
                       directory = "Data/GDCdata/",                  # Directory where files are stored
                       summarizedExperiment = F)                     # Do not create a summarized experiment

#--------------------Initial pre-processing--------------------
## Keep only protein coding RNAs
exp_data <- exp_data[which(exp_data$gene_type == "protein_coding"),] # 19,962 genes identified

## Remove rows with no transcript reads 
exp_data <- exp_data[rowSums(exp_data[,-c("gene_id",                 # 19,498 genes remain
                                          "gene_name", 
                                          "gene_type")]) != 0,]

## Set ensembl gene id as row names
rownames(exp_data) <- exp_data$gene_id

## Remove ensembl id column as well as the gene name and gene type
exp_data <- dplyr::select(exp_data, -c("gene_id", "gene_name", "gene_type"))

## Keep only unstranded reads
keep_cols <- colnames(exp_data)[sapply(colnames(exp_data), function(col){
  strsplit(col, "_")[[1]][1] == "unstranded"})]
exp_data <- exp_data %>% dplyr::select(keep_cols)

## Remove the "unstranded_" part of the columns' names
colnames(exp_data) <- unlist(strsplit(colnames(exp_data), "_"))[
  unlist(strsplit(colnames(exp_data), "_")) != "unstranded"]

## Set column names just as the patient_id
colnames(exp_data) <- substr(colnames(exp_data), 1, 19)

## Remove the transcript version and only keep ensembl gene id
rownames(exp_data) <- sapply(strsplit(rownames(exp_data),".",fixed=T),
                             function(x) x[1])

## Remove samples that are not in the samples' patient_id
exp_data <- as.data.frame(exp_data, row.names = row.names(exp_data))
exp_data <- exp_data[,which((colnames(exp_data) %in% samples_data$patient_id))]

## Remove rows with no transcript reads again
exp_data <- exp_data[rowSums(exp_data) != 0, , drop = FALSE]         # 19,448 genes remain

#--------------------Annotation data----------------------
# Get annotation data
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
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
ann_data$length <- abs(ann_data$end_position - ann_data$start_position)        # Add length

## Remove non protein coding genes from annotation
ann_data <- ann_data[which(ann_data$gene_biotype == "protein_coding"),]        # 19,340 genes remain

## Remove transcripts with no annotation data
exp_data <- exp_data[which(rownames(exp_data) %in% ann_data$ensembl_gene_id),] # Kept the 19,340 from the annotation data

#--------------------Check for biases---------------------
## Create noiseq object for it to be compatible with selected workflow
noiseqData <- NOISeq::readData(exp_data,
                               factors = samples_data,
                               gc = ann_data[, c("ensembl_gene_id", "percentage_gene_gc_content")],
                               biotype = ann_data[, c("ensembl_gene_id", "gene_biotype")],
                               length = ann_data[, c("ensembl_gene_id", "length")])

counts_data <- dat(noiseqData, type = "countsbio", factor = "sample_type")

## Expression values for cancer samples and controls
png("Figures/mRNA/exp_vals_before_norm.png",width=1000)
explo.plot(counts_data, plottype = "boxplot")
dev.off()

## Expression values in CPM for cancer samples and controls
png("Figures/mRNA/exp_vals_bar_before_norm.png",width=1000)
explo.plot(counts_data, plottype = "barplot")
dev.off()

## Visualize low CPM values
cpm_hist <- ggplot(exp_data, aes(x = rowMeans(cpm(exp_data, log = T)))) + 
  geom_histogram(colour = "blue", fill = "lightblue") + xlab("CPM") + 
  ylab("Genes") + 
  theme_classic() + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red")
cpm_hist
ggsave("Figures/mRNA/exp_cpm_before_norm.png", plot = cpm_hist)            # Save plot 
sum(rowMeans(cpm(exp_data, log = T))>0)/nrow(exp_data)*100                 # ~64% genes have CPM>0

## Check for transcript composition bias
cd_data <- dat(noiseqData, type = "cd", norm = F)                          # Reference sample is: TCGA-2V-A95S-01A-11
table(cd_data@dat$DiagnosticTest[, "Diagnostic Test"])                     # 332 Failed and 21 Passed (01-08-2025)

png("Figures/mRNA/exp_transcript_bias_before_norm.png",width=1000)
explo.plot(cd_data, samples = sample(1:ncol(exp_data),10))                 # The result clearly shows composition bias
dev.off()

## Check for GC bias
# The results show the effect on expression values based on the GC content.The fit for the cancer samples
# is 58.77% with a p value of 1.2e-11. The fit for the controls is of 51.27% and a p value of 6.4e-09
gc_content <- dat(noiseqData, type = "GCbias", k = 0, factor = "sample_type")
png("Figures/mRNA/exp_gc_bias_before_norm.png",width=1000)
explo.plot(gc_content)  
dev.off()

# Check for length bias
# The results show the effect on expression values based on the length. The fit for the cancer samples is 58.75% 
# with a p value of 3.5e-12. The fit for the controls is of 56.51% and a p value of 2.8e-11
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

#--------------------Descriptive plot-------------------------
## Descriptive data per sample
### Create empty data frame
samples_summary <- data.frame(matrix(NA, 
  nrow = ncol(exp_data), 
  ncol= 6),
  row.names = colnames(exp_data))

### Add summary informatio per sample
for(i in 1:ncol(exp_data)){
  samples_summary[i,] <- summary(exp_data[,i])
}

### Rename columns
colnames(samples_summary) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")

## Violin plot
log(samples_summary + 1) %>%
  pivot_longer(cols = colnames(samples_summary)) %>%
  ggplot(., aes(x = reorder(name, value), y= value, fill= name)) + 
  geom_boxplot() + ylim(0,17) +
  ylab("log reads per sample") + xlab(" ") + 
  labs(dictionary = c(name = "Measures")) +
  theme_classic()

#--------------------Solve biases-------------------------
# Filter genes with low counts (CPM) < 0
exp_data <- filtered.data(exp_data,                   # Data to filter    
                          factor = "sample_type",     # Samples' conditions
                          norm = F,                   # If normalized or not
                          method = 1,                 # Filtering method (CPM)
                          cpm = 0,                    # CPM threshold
                          p.adj = "fdr")              # Correction method
dim(exp_data)                                         # 9,840 genes remain

# Filter annotation data
ann_data <- ann_data[which(ann_data$ensembl_gene_id %in% rownames(exp_data)),]

# Column names must match
#colnames(exp_data) <- samples_data$barcode

# Solve GC and length bias using cqn (conditional quantile normalization)
counts <- as.matrix(exp_data)                         # Raw counts 
lengths <- ann_data$length                            # Gene lengths 
gc_content <- ann_data$percentage_gene_gc_content     # GC content 

# Apply CQN 
cqn_res <- cqn(counts = counts,                       # Expression counts
               lengths = lengths,                     # Gene length
               x = gc_content,                        # Covariate to remove
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
cqn_arsyn <- ARSyNseq(cqn_noiseq,                                    # Object
                      factor = "sample_type",                        # Factor containing batch information
                      batch = F,                                     # F 'cause the sample type is not a bacth in itself
                      norm = "n",                                    # n because it has already been normalized
                      logtransf = T)                                 # T because we don't want a log-transform

# Check batch effect removal visually
myPCA <- dat(cqn_arsyn, type = "PCA", norm = T, logtransf = T)       # Perform PCA
png("Figures/mRNA/exp_batch_effect_after_norm.png",width=1000)
explo.plot(myPCA, factor = "sample_type")                            # PCA1 = 6%, PCA2 = 1%
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
# Now it is R2=32.75% and p=0.016 for tumor and R2=28.71% and p=0.4 
png("Figures/mRNA/exp_gc_bias_after_norm.png",width=1000)
explo.plot(new_gc_content)  
dev.off()

# Check for length bias
new_len_bias <- dat(new_noiseq, k = 0, type = "lengthbias", factor = "sample_type", norm = T)
# Now it is R2=43.08% and p=0.00096 for tumor and R2=48.39% and p=0.00017
png("Figures/mRNA/exp_length_bias_after_norm.png",width=1000)
explo.plot(new_len_bias)
dev.off()

#--------------------Save expression matrices-------------
write.table(new_noiseq@assayData$exprs,"Data/norm_exp_data.tsv",sep='\t',row.names=T) # Normalized mRNA expression values

#--------------------Differential expression analysis-----
# Differential expression analysis performed with DESeq2
## Create DESeq object
diff_des <- DESeqDataSetFromMatrix(countData = counts,            # Un-normalized counts matrix
                                   colData = samples_data,        # Meta data
                                   design = ~ sample_type)        # Variable of interest
diff_des$sample_type <- relevel(diff_des$sample_type,             # Set the control samples as the reference
                                ref = "Solid Tissue Normal")

## Perform differential expression analysis
diff_anl <- DESeq(diff_des)

## Extract results
diff_res <- results(diff_anl)

## Make plot
diff_res <- as.data.frame(diff_res)

### Assign labels if unchanged, up-regulated or down-regulated 
diff_res <- diff_res %>% mutate(Expression = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "Up-regulated",
                                                       log2FoldChange <= -1 & padj <= 0.05 ~ "Down-regulated",
                                                       TRUE ~ "Unchanged"))

#### Check for differential expression results
table(diff_res$Expression)
# Down-regulated      Unchanged   Up-regulated 
#            108           9677             55 

### Extract ID from 50 top and 50 bottom differentially expressed
bot_50_genes <- diff_res[which(
  diff_res$Expression == "Down-regulated"),][order(diff_res[which(
    diff_res$Expression == "Down-regulated"),]$log2FoldChange, decreasing = F), ] %>% head(.,50)

top_50_genes <- diff_res[which(
  diff_res$Expression == "Up-regulated"),][order(diff_res[which(
    diff_res$Expression == "Up-regulated"),]$log2FoldChange, decreasing = T), ] %>% head(.,50)

## Plot itself
png("Figures/mRNA/exp_diff_plot.png",width=1000)
ggplot(diff_res, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 1.5, alpha=.7) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  labs(x = "log_2(FC)", y = "-log10(adj. p-Value)", title = "Differential Expression") + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "black") +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color = "black") +
  theme_bw() +
  geom_label_repel(data=subset(diff_res,rownames(diff_res) %in% 
                                 rownames(top_50_genes[1:10,]) | 
                                 rownames(diff_res) %in% rownames(bot_50_genes[1:10,])),
                   aes(log2FoldChange, -log(padj,10), 
                       label = rownames(subset(diff_res,rownames(diff_res) %in%
                                                 rownames(top_50_genes[1:10,]) |
                                                 rownames(diff_res) %in% rownames(bot_50_genes[1:10,])))),
                   arrow = arrow(length = unit(0.02, "npc")),
                   nudge_x = 0.5)
dev.off()

## More relaxed differential expression anaylisis
diff_res <- results(diff_anl)
diff_res <- as.data.frame(diff_res)
### Assign labels if unchanged, up-regulated or down-regulated based on lfc  of 0.5 instead of 1
diff_res_relaxed <- diff_res %>% 
  mutate(Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Up-regulated",
                                log2FoldChange <= -0.5 & padj <= 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))
table(diff_res_relaxed$Expression)
# Down-regulated      Unchanged   Up-regulated 
#            792           8147            901 
bot_50_genes_relaxed <- diff_res_relaxed[which(
  diff_res_relaxed$Expression == "Down-regulated"),][order(diff_res_relaxed[which(
    diff_res_relaxed$Expression == "Down-regulated"),]$log2FoldChange, decreasing = F), ] %>% head(.,50)

top_50_genes_relaxed <- diff_res_relaxed[which(
  diff_res_relaxed$Expression == "Up-regulated"),][order(diff_res_relaxed[which(
    diff_res_relaxed$Expression == "Up-regulated"),]$log2FoldChange, decreasing = T), ] %>% head(.,50)

png("Figures/mRNA/exp_diff_relaxed_plot.png",width=1000)
ggplot(diff_res_relaxed, aes(log2FoldChange, -log(padj,10))) +
  geom_point(aes(color = Expression), size = 1.5, alpha=.7) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  labs(x = "log_2(FC)", y = "-log10(adj. p-Value)", title = "Differential Expression") + 
  geom_hline(yintercept=0.05, linetype="dashed", color = "black") +
  geom_vline(xintercept=c(-0.5,0.5), linetype="dashed", color = "black") +
  theme_bw() +
  geom_label_repel(data=subset(diff_res_relaxed,rownames(diff_res_relaxed) %in% 
                                 rownames(top_50_genes[1:10,]) | 
                                 rownames(diff_res_relaxed) %in% rownames(bot_50_genes[1:10,])),
                   aes(log2FoldChange, -log(padj,10), 
                       label = rownames(subset(diff_res_relaxed,rownames(diff_res_relaxed) %in%
                                                 rownames(top_50_genes[1:10,]) |
                                                 rownames(diff_res_relaxed) %in% rownames(bot_50_genes[1:10,])))),
                   arrow = arrow(length = unit(0.02, "npc")),
                   nudge_x = 0.5)
dev.off()

#--------------------Save objects-------------------------
write.table(diff_res,"Data/exp_diff.tsv",sep='\t',row.names=T)
write.table(diff_res_relaxed,"Data/exp_diff_relaxed.tsv",sep='\t',row.names=T)


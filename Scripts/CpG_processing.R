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
library(minfi)

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

# Filter probes with missing values in 20% or more of the samples
met_data <- met_data[which(!rowSums(                                 # 396,413 probes remain
  is.na(met_data@assays@data@listData[[1]])) >
    (ncol(met_data@assays@data@listData[[1]])*0.2)),]

# Filter data with just 0's aside from NAs
met_data <- met_data[rowSums(met_data@assays@data@listData[[1]],     # 396,413 probes remain
                             na.rm=T) != 0, ,drop = FALSE]

# Removing MASKED probes (Check your object for the MASK status)
## Remove recommended masked probes
met_data <- met_data[which(met_data@rowRanges@                       # 340,838 probes remain         
                             elementMetadata@listData[[
                               "MASK_rmsk15"]] == F),]

## Remove probes overlapping with common SNP's
met_data <- met_data[which(met_data@rowRanges@                       # 305,328 probes remain           
                             elementMetadata@listData[[
                               "MASK_snp5_common"]] == F),]

## Check for probes with ambiguous chromosome mapping
seqnames(rowRanges(met_data)) %>% table()                            # No ambiguous mapping apparently (25-11-2025)

#--------------------Methylation data imputation----------
# Impute missing values with regression based method. See: https://doi.org/10.1186/s12859-020-03592-5
# I will be separating the data into chunks and doing parallel processing for each chunk since I've had
# some problems running it. If you are more comfortable and knowledgeable, run it in one single chunk 

# Before actually imputing, I will separate the probes that do not have any missing (NA) value to make this faster
complete_probes <- met_data[which(                                   # 243,900 complete probes (22-04-2026)
  rowSums(is.na(assay(met_data))) == 0),]

incomplete_probes <- met_data[which(                                 # 61, 428 incomplete probes (22-04-2026)
  rowSums(is.na(assay(met_data))) != 0),]
  
# Set the basic for running by chunks
workers <- MulticoreParam(workers = 2)                               # Number of workers
chunk_size <- 13000                                                  # Chunk size
outdir <- "Data/methy_chunks"                                        # Directory for chunk output
dir.create(outdir)                                                   # Create directory

# Creating chunks
groups <- samples_data$sample_type                                   # Separate data by tumor and control
total_cpgs <- nrow(incomplete_probes)                                # Total CpG islands
chunk_start <- seq.int(1, total_cpgs, chunk_size)                    # Create chunks' ranges/sizes
chunk_end <- pmin(chunk_start + chunk_size - 1, total_cpgs)          # Create the end of the chunks

# Annotation since there was a problem when subsetting internally by chromosome
data("hm450.manifest.hg19", package = "ChAMPdata")
ann <- data.frame(cpg = rownames(hm450.manifest.hg19),
                  chr = hm450.manifest.hg19$CpG_chrm)
removed_probes_list <- list()                                        # Will store removed probes just to have some info if needed
removed_counts <- integer(length(chunk_start))

# Running chunks
  # Set iterations
for(i in seq_along(chunk_start)){                                    # Iterate over each chunk
  outfile <- file.path(outdir, sprintf("chunk_%02d.rds", i))         # Create a new file with the chunk's imputed data into the out directory
  if(file.exists(outfile)){                                          # Skip if the chunk has already been imputed, since it may crash
    message("Skipping chunk ", i, " (already exists)")
    next
  }
  
  # Set chunks' ids
  idx <- chunk_start[i]:chunk_end[i]                                 # Create chunk id by range
  met_chunk <- incomplete_probes[idx, , drop = F]                    # Create chunk
  message("Running chunk ", i," (", chunk_start[i], ":",             # Tell which chunk is running
  chunk_end[i], ")")
  
  # Extract only valid probes
  met_mat <- assay(met_chunk)                                        # Extract only data matrix from selected chunk
  chunk_anno <- ann[match(rownames(met_mat), ann$cpg), ]             # Match probe annotation data with chunk's probes
  valid_match <- !is.na(chunk_anno$chr)                              # Make sure probe is matched to a chromosome
  chr_counts <- table(chunk_anno$chr[valid_match])                   # Tell how many probes there are per chromosome per chunk
  valid_chr <- names(chr_counts[chr_counts > 1])                     # Make sure there are at least 2 robes to actually run
  keep <- valid_match & (chunk_anno$chr %in% valid_chr)              # Select only fully annotated CpGs

  # Removed probes info
  removed_probes <- rownames(met_mat)[!keep]
  removed_probes_list[[i]] <- removed_probes
  removed_counts[i] <- length(removed_probes)

  # Actual imputation over each chunk
  met_chunk <- met_chunk[keep, , drop = FALSE]
  if (nrow(met_chunk) == 0) {
    message("Skipping chunk ", i, " (no valid probes after filtering)")
    next
  }
  res <- methyLImp2(met_chunk,
                    type    = "450K",
                    groups  = groups,
                    BPPARAM = workers)
  saveRDS(res, outfile)                                              # Save output
  rm(res, met_chunk)                                                 # Remove object from environment as to save space
  invisible(gc())                                                    # Free memory and avoid printing output
}

total_removed <- sum(removed_counts, na.rm = T)                      # Just check the removed probes
table(removed_counts)                                                # Only one probe was removed

# Load chunks into environment and concatenate everything
coldata_ref <- colData(complete_probes)                              # Set metadata reference
path_name <- "Data/methy_chunks/"                                    # Path of chunks

# Read all chunks into a list
chunk_files <- list.files(path_name, full.names = TRUE)              # Get every chunk's name
chunk_list <- lapply(chunk_files, function(f){                       # Read imputed chunks and set reference metadata as main
  x <- readRDS(f)
  colData(x) <- coldata_ref
  x
})

# Combine all imputed chunks (row-wise = probes)
imputed_combined <- do.call(rbind, chunk_list)

# Merge with complete probes
met_imputed <- rbind(complete_probes, imputed_combined)

write.table(assay(met_imputed), "Data/met_imputed.tsv",
            sep=',',row.names=T)                                     # Save object with the combined imputations

# Plot beta values for tumor samples and adjacent tissues
## By each sample
png("Figures/CpG/cpg_beta_distribution.png", width=1000) 
densityPlot(assay(met_imputed), sampGroups = met_imputed@colData$sample_type)
dev.off()

## By sample type
tumor_beta <- assay(met_imputed[,                                    # Beta values for tumor samples
  which(met_imputed@colData$sample_type == "Primary Tumor")])

control_beta <- assay(met_imputed[,                                  # Beta values for adjacent tissues
  which(met_imputed@colData$sample_type != "Primary Tumor")])

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
m_values <- BetaValueToMValue(assay(met_imputed))
write.table(m_values,"Data/m_values.tsv",sep=',',row.names=T)        # Save object with M-values

#--------------------Check batch effect-------------------
noiseqData <- readData(data = m_values, factor = samples_data)       # Create noiseq object
myPCA <- dat(noiseqData, type = "PCA", norm = T, logtransf = T)      # Perform PCA to watch for batch effects

# The samples partially aggregate by cancerous and control
# PC1 explains 29% and PC2 explains 28%
png("Figures/CpG/cpg_batch_effect.png")                              # Plot batch effect with PCA
explo.plot(myPCA, factor = "sample_type")
dev.off()

#--------------------Remove btach effect------------------
com_mod <- model.matrix(~sample_type, data = samples_data)           # Complete model (adjustment and variable of interest cancer or control)
nul_mod <- model.matrix(~1, data = samples_data)                     # Null model, only include intercept
num_la_f <- num.sv(m_values, com_mod, method = "leek")               # Identify number of latent factors to estimate


#--------------------Differential methylation-------------
#TCGAanalyze_DMC(met_data,
#                groupCol = "definition",
#                p.cut = 0.5,
#                title = "Differential methylation",
#                save.directory = "3_Data/")

# This script is meant as the methylation values pre-processing guideline

#--------------------Load libraries-----------------------
library(TCGAbiolinks)            # Version: 2.36.0
library(SummarizedExperiment)    # Version: 1.48.1
library(tidyverse)               # Version: 2.0.0
library(methyLImp2)              # Version: 1.2.0
library(BiocParallel)            # Version: 1.40.2

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
# some problems running it. If you are more comfortable and knowledegable, run it in one single chunk 

# Set the basic for running by chunks
workers <- MulticoreParam(workers = 2)                               # Number of workers
chunk_size <- 50000                                                  # Chunk size
outdir <- "Data/methy_chunks"                                        # Directory for chunk output
dir.create(outdir)                                                   # Create directory

# Creating chunks
groups <- met_data$definition                                        # Separate data by tumor and control
total_cpgs <- nrow(met_data)                                         # Total CpG islands
chunk_start <- seq.int(1, total_cpgs, chunk_size)                    # Create chunks' ranges/sizes
chunk_end <- pmin(starts + chunk_size - 1, total_cpgs)               # Create the end of the chunks

# Running chunks
for (i in seq_along(chunk_start)) {                                  # Iterate over each chunk
  outfile <- file.path(outdir, sprintf("chunk_%02d.rds", i))         # Create a new file with the chunk's imputed data into the out directory
  if (file.exists(outfile)) {                                        # Skip if the chunk has already been imputed, since it may crash
    message("Skipping chunk ", i, " (already exists)")
    next
  }
  idx <- chunk_start[i]:chunk_end[i]                                 # Create chunk id by range
  met_chunk <- met_data[idx, ]                                       # Create chunk
  message("Running chunk ", i," (", starts[i], ":", ends[i], ")")    # Tell which chunk is running
  res <- methyLImp2(                                                 # Impute data for chunk
    met_chunk,
    type    = "450K",
    groups  = groups,
    BPPARAM = workers)
  saveRDS(res, outfile)                                              # Save output
  rm(res, met_chunk)                                                 # Remove object from environment as to save space
  invisible(gc())                                                    # Free memory and avoid printing output
}

#--------------------B-values to M-values-----------------
# See https://doi.org/10.1186/1471-2105-11-587 for decision basis


#--------------------Differential methylation-------------
#TCGAanalyze_DMC(met_data,
#                groupCol = "definition",
#                p.cut = 0.5,
#                title = "Differential methylation",
#                save.directory = "3_Data/")

# This script is meant for downloading all the necessary data fro the project. Keep in mind that, while this script
# is designed to work with the current (03-07-2025) data available in the TCGA platform, it will most likely be
# outdated in some aspects in the future. Please keep this in mind and modify it if necessary.

#--------------------Load libraries-----------------------
library(TCGAbiolinks)            # Version: 2.36.0
library(SummarizedExperiment)    # Version: 1.48.1
library(tidyverse)               # Version: 2.0.0
library(RColorBrewer)            # Version: 1.1.3
library(VennDiagram)             # Version: 1.8.2

#--------------------Query preparation--------------------
# This section is in charge of making the queries to the TCGA database
## mRNA
exp_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "Gene Expression Quantification",  # Gene expression 
                      workflow.type = "STAR - Counts")               # How reads are set as counts per gene

## miRNA
mir_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "miRNA Expression Quantification") # miRNA gen expression 

## Methylation data
met_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "DNA Methylation",             # DNA methylation data
                      platform="Illumina Human Methylation 450")     # CpG detection platform

#--------------------Getting query results----------------
# These functions get the results for the query in a manner that can be later parsed and used to download data
exp_res <- getResults(exp_query)
mir_res <- getResults(mir_query)
met_res <- getResults(met_query)

#--------------------Data pre-selection-------------------
# Identify the cases that have available data for mRNA, CpG islands, and miRNA
# The current TCGA cases' ids are composed with various characters (Eg."TCGA-G3-A3CJ-01A-11R-A213-07"). 
# The way to make sure the same samples are selected despite their data type (RNA or CpG) is cropping 
# the id and selecting the common part. Right now, the first 19 characters are enough to provide this information.
# Additionally, the work will only focus on primary tumors and any available controls, so filtering is required.
## Check for sample types
table(exp_res$sample_type) # (371 primary, 3 recurrent, and 50 normal) (08-07-2025)

## Keep only primary and normal
exp_res <- filter(exp_res, sample_type != "Recurrent Tumor")

## Extract initial samples' names (the first 19 characters are shared amongst all three data types)
cases <- substr(exp_res$cases, 1, 19)

## Check how many are of the cases are shared among the three data types
cases <- cases[cases %in% substr(met_res$cases, 1, 19) & 
                   cases %in% substr(mir_res$cases, 1 ,19)]
length(cases)                                                              # 407 samples shared across all thre data types  

## Plotting which samples are shared amongst all three data types
venn.diagram(x = list(cases, substr(met_res$cases, 1, 19), substr(mir_res$cases, 1 ,19)),
             category.names = c("mRNA", "CpG", "miRNA"),
             filename = "Figures/samples_venn.png",
             fill = brewer.pal(3, "Set1"),
             disable.logging = T,
             output = T,
                      lwd = 3,
                      col = brewer.pal(3, "Set1"))

#--------------------Clinical data------------------------
# Available clinical data is also important as it can provide more insight into each sample. Not every
# characteristic will be used, but be sure to check which columns you wish to keep for further analysis.
## Get clinical data
cli_data <- GDCquery_clinic(project = "TCGA-LIHC",                   # Liver hepatocellular carcinoma project 
                            type = "clinical")                       # Acquire clinical data
cli_data <- select(cli_data,!where(~ all(is.na(.x))))                # Remove empty columns to parse only through available clinical data
for(column in colnames(cli_data)){                                   # Remove columns with only "NA" as string
  if(all(cli_data[, column] == "NA", na.rm = T)){
    cli_data <- select(cli_data, -column)
  }
}

## Select clinical features based on availability and interest
cli_data <- cli_data %>% select(c("bcr_patient_barcode",
                                  "primary_diagnosis",
                                  "synchronous_malignancy",
                                  "prior_malignancy",
                                  "age_at_diagnosis",
                                  "child_pugh_classification",
                                  "prior_treatment",
                                  "follow_ups_disease_response",
                                  "tumor_grade",                     # Numeric value to express the degree of abnormality of cancer cells, a measure of differentiation and aggressiveness.
                                  "ajcc_pathologic_stage",           # The extent of a cancer, especially whether the disease has spread from the original site to other parts of the body based on AJCC staging criteria.
                                  "ajcc_staging_system_edition",     # AJCC version
                                  "race",
                                  "age_at_index",
                                  "sex_at_birth",
                                  "vital_status"))

## Format race's data
cli_data[which(cli_data$race == "not reported"),]$race <- "Unknown"                                # Change "not reported" to "Unknown"

#--------------------Exclusion criteria--------------------
## First remove from clinical data
cli_data <- filter(cli_data, primary_diagnosis == "Hepatocellular carcinoma, NOS")                 # Only keep NOS subtype
cli_data <- cli_data[which(!cli_data$prior_treatment == "Yes"),]                                   # Remove samples that had treatment previous to being samapled
cli_data <- cli_data[which(cli_data$prior_malignancy == "no"),]                                    # Only keep patients with certainty they haven't had aprior malignancies
cli_data <- cli_data[which(cli_data$synchronous_malignancy == "No"),]                              # Remove patients with synchronous malignancies

## Move on to samples, both tumor and control tissues
cases <- cases[substr(cases, 1, 12) %in% cli_data$bcr_patient_barcode]                             # 354 total samples remain

## Create samples' data frame
samples_data <- data.frame(cbind(substr(cases, 1, 12),                                             # Barcode, patient id, sample type
                                cases,
                                exp_res[substr(exp_res$cases, 1, 19) %in% cases,]$sample_type))
colnames(samples_data) <- c("barcode", "patient_id", "sample_type")                                # Rename columns

## Merge with clinical data
samples_data <- merge(samples_data, cli_data, by.x="barcode", by.y="bcr_patient_barcode")          # 354 total samples
table(samples_data$sample_type)                                                                    # 315 tumor samples and 39 solid tissue normals

#--------------------Save object-------------------------
write.table(samples_data, "Data/samples_data.tsv", sep='\t', quote=F, row.names=F)

#--------------------Data download------------------------
# This section is in charge of downloading the queries to the TCGA database. The barcode tells the TCGA
# platform which specific samples I want to request.
## mRNA
exp_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "Gene Expression Quantification",  # Gene expression 
                      workflow.type = "STAR - Counts",               # How reads are set as counts per gene
                      barcode = samples_data$barcode)                # Barcode
GDCdownload(exp_query, directory = "Data/GDCdata")                   # Downloading files

## miRNA
mir_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "Transcriptome Profiling",     # Refers to RNA
                      data.type = "miRNA Expression Quantification", # miRNA gen expression
                      barcode = samples_data$barcode)                # Barcode
GDCdownload(mir_query, directory = "Data/GDCdata")                   # Downloading files

## Methylation data
met_query <- GDCquery(project = "TCGA-LIHC",                         # Liver hepatocellular carcinoma project
                      data.category = "DNA Methylation",             # DNA methylation data
                      platform = "Illumina Human Methylation 450",   # CpG detection platform
                      data.type = "Methylation Beta Value",          # Data type
                      barcode = samples_data$barcode)                # Barcode
GDCdownload(met_query, directory = "Data/GDCdata",                   # Downloading files. I had to select a "files
            files.per.chunk = 50)                                    # per chunk" argument due to the data's size
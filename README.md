# Inferencia de agentes reguladores en el carcinoma hepatocelular mediante la integración de datos multiómicos

This repository will contain all the work pertaining the Hepatocellular Carcinoma (HCC) multiomic integration project.

## Pipeline Overview

**1. Data Acquisition:** Download mRNA, miRNA, and methylation values for samples of HCC available in the TCGA databases.

- [Get_data](./Scripts/Get_data.R): Download data of selected cases.

**2. Data pre-processing:** Filter data, reduce present biases, normalize data, perform differential expression analyses, and store pre-processed data for all three data types.

- [mRNA_processing.R](./Scripts/mRNA_processing.R): Pre-process mRNA data.
- [miRNA_processing.R](./Scripts/miRNA_processing.R): Pre-process miRNA data.
- [CpG_processing.R](./Scripts/CpG_processing.R): Pre-process CpG data.
- [ENET.R](./Scripts/ENET.R): Scale, center, normalize, and perform ENET for each omic data block.

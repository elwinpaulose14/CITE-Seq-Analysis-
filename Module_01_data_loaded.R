library(Seurat)
library(SeuratDisk)
library(BoneMarrowMap)
library(AUCell)
library(tidyverse)
library(BiocManager)
library(symphony)
library(ggplot2)
# Define data paths
getwd()
setwd("//10.100.75.55/NASShare/ELWIN/BDRhapsody")
data_dir <- "//10.100.75.55/NASShare/ELWIN/BDRhapsody"
rds_file <- file.path(data_dir, "CITE-Seq-Lane3_Seurat.rds")
csv_file <- file.path(data_dir, "CITE-Seq-Lane3_Sample_Tag_Calls.csv")

# Verify files exist
if (!file.exists(rds_file)) {
  stop("ERROR: RDS file not found at ", rds_file)
}

if (!file.exists(csv_file)) {
  stop("ERROR: CSV file not found at ", csv_file)
}

cat("Loading CITE-seq data from RDS file...\n")
RNA_data <- readRDS(rds_file)

cat("Initial data dimensions:\n")
cat("  Cells:", ncol(RNA_data), "\n")
cat("  Genes:", nrow(RNA_data), "\n")

# Check available assays
cat("\nAvailable assays:\n")
print(Assays(RNA_data))

# Load sample metadata
cat("\nLoading sample tag metadata...\n")
smk <- read.table(csv_file, sep = ",", header = TRUE, row.names = 1)

cat("Sample metadata dimensions:\n")
cat("  Samples:", nrow(smk), "\n")
cat("  Columns:", ncol(smk), "\n")

cat("\nSample metadata columns:\n")
print(colnames(smk))

# Add metadata to Seurat object
RNA_data <- AddMetaData(RNA_data, metadata = smk)

cat("\nMetadata successfully added to Seurat object\n")
cat("\nCurrent metadata columns:\n")
print(colnames(RNA_data@meta.data))

# Save intermediate object
cat("\nSaving intermediate object...\n")
saveRDS(RNA_data, file = "01_data_loaded.rds")

cat("\n✓ Step 1 completed successfully\n")

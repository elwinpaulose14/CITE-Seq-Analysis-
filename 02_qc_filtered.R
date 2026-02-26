# Load data from previous step
RNA_data <- readRDS("01_data_loaded.rds")

cat("Initial cell count:", ncol(RNA_data), "\n")


# 2.1: Calculate Mitochondrial Percentage

cat("\nCalculating mitochondrial gene percentage...\n")
mito_genes <- grep(pattern = "^MT-", x = rownames(RNA_data), value = TRUE)
cat("  Mitochondrial genes detected:", length(mito_genes), "\n")

RNA_data[["percent.mt"]] <- PercentageFeatureSet(RNA_data, pattern = "^MT-")

# Generate QC plots
cat("\nGenerating QC violin plots...\n")
qc_plot <- VlnPlot(
  RNA_data,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
)
ggsave(
  filename = file.path(data_dir, "01_qc_before_filtering.png"),
  plot = qc_plot,
  width = 12,
  height = 4,
  dpi = 300
)

cat("\n--- Calculating Dynamic QC Thresholds ---\n")


# 2.3: Apply Filtering Steps

lb <- quantile(RNA_data[["nFeature_RNA"]]$nFeature_RNA, probs = 0.01)
ub <- quantile(RNA_data[["nFeature_RNA"]]$nFeature_RNA, probs = 0.99)

RNA_data <- RNA_data[, RNA_data[["nFeature_RNA"]] > lb & RNA_data[["nFeature_RNA"]] < ub]
cat("  Cells after quantile filtering:", ncol(RNA_data), "\n")

# Step 2.3.2: Apply feature count thresholds

cat("\nStep 2.3.2: Applying feature count thresholds ...\n")
RNA_data <- subset(RNA_data, subset = nFeature_RNA > 150 & nFeature_RNA < 2700 & percent.mt < 45)
cat("  Cells after feature filtering:", ncol(RNA_data), "\n")

# Step 2.3.3: Remove multiplets and undetermined cells

cat("\nStep 2.3.3: Removing multiplets and undetermined cells...\n")
print(table(RNA_data$Sample_Name))

RNA_data <- subset(RNA_data, Sample_Name != "Undetermined" & Sample_Name != "Multiplet")
cat("  Cells after multiplet removal:", ncol(RNA_data), "\n")

# Step 2.3.4: Select specific sample tag (SampleTag11_hs)
cat("\nStep 2.3.4: Selecting SampleTag11_hs...\n")
if ("Sample_Tag" %in% colnames(RNA_data@meta.data)) {
  RNA_data <- subset(RNA_data, Sample_Tag %in% c("SampleTag11_hs"))
  cat("  Sample tag distribution after filtering:\n")
  print(table(RNA_data$Sample_Tag))
} else {
  cat("  Warning: Sample_Tag column not found, skipping sample selection\n")
}

# 2.4: Post-Filtering Summary

cat("\nFinal cell count:", ncol(RNA_data), "\n")

# Rename for pipeline continuity
allcombine <- RNA_data
rm(RNA_data)

# Save filtered object
cat("\nSaving filtered object...\n")
saveRDS(allcombine, file = "02_qc_filtered.rds")

cat("\n✓ Step 2 completed successfully\n")

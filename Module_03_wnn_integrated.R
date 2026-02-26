 Multimodal Integration (RNA + ADT)

# Load filtered data
allcombine <- readRDS("02_qc_filtered.rds")

cat("Working with", ncol(allcombine), "cells\n")


# 3.1: RNA Normalization and Dimensionality Reduction

cat("Normalizing RNA data...\n")
DefaultAssay(allcombine) <- "RNA"
allcombine <- NormalizeData(allcombine)

cat("Identifying highly variable features...\n")

allcombine <- FindVariableFeatures(allcombine, nfeatures = 3000)

cat("Number of variable features:", length(VariableFeatures(allcombine)), "\n")

cat("Scaling RNA data...\n")
allcombine <- ScaleData(allcombine)

cat("Running PCA on RNA data...\n")

allcombine <- RunPCA(allcombine, npcs = 50, reduction.name = "rpca")

cat("  RNA PCA dimensions:", ncol(Embeddings(allcombine, "rpca")), "\n")

# Save elbow plot for dimension selection
png("03_rna_elbow_plot.png", width = 800, height = 600, res = 150)
ElbowPlot(allcombine, ndims = 50, reduction = "rpca")
dev.off()
cat("  Elbow plot saved\n")


# 3.2: ADT (Protein) Normalization and Dimensionality Reduction

DefaultAssay(allcombine) <- "ADT"

cat("Available ADT features:\n")
print(rownames(allcombine[["ADT"]]))

cat("\nNormalizing ADT data using CLR method...\n")
# margin = 2: Normalize across cells (recommended for CITE-seq)
allcombine <- NormalizeData(
  allcombine,
  normalization.method = "CLR",
  margin = 2
)

cat("Scaling ADT data...\n")
allcombine <- ScaleData(allcombine)

cat("Running PCA on ADT data...\n")
allcombine <- RunPCA(
  allcombine,
  features = rownames(allcombine[["ADT"]]),
  npcs = 9,
  reduction.name = "apca"
)

cat("  ADT PCA dimensions:", ncol(Embeddings(allcombine, "apca")), "\n")

# Save ADT elbow plot
png("04_adt_elbow_plot.png", width = 800, height = 600, res = 150)
ElbowPlot(allcombine, ndims = 9, reduction = "apca")
dev.off()
cat("ADT elbow plot saved\n")

# 3.3: Weighted Nearest Neighbor (WNN) Analysis

cat("Computing multimodal nearest neighbors...\n")
cat("  RNA dimensions: 1-30\n")
cat("  ADT dimensions: 1-9\n")

allcombine <- FindMultiModalNeighbors(
  allcombine,
  reduction.list = list("rpca", "apca"),
  dims.list = list(1:30, 1:9),
  modality.weight.name = "RNA.weight"
)

# Add metadata for tracking
allcombine$integration_method <- "RNA_ADT_WNN"
allcombine$object_version <- "allcombine_WNN_v1"

# Verify modality weights
cat("\n=== RNA Modality Weights ===\n")
print(summary(allcombine$RNA.weight))

cat("\n=== ADT Modality Weights ===\n")
print(summary(allcombine$ADT.weight))

cat("\nMetadata columns:\n")
print(colnames(allcombine@meta.data))

# 3.4: Visualization of Modality Weights

cat("\nGenerating modality weight histograms...\n")

png("05_modality_weights.png", width = 1200, height = 500, res = 150)
par(mfrow = c(1, 2))

hist(allcombine$RNA.weight,
     breaks = 50,
     col = "steelblue",
     main = "RNA Modality Weights",
     xlab = "RNA Weight",
     xlim = c(0, 1))
abline(v = median(allcombine$RNA.weight), col = "red", lwd = 2, lty = 2)

hist(allcombine$ADT.weight,
     breaks = 50,
     col = "coral",
     main = "ADT Modality Weights",
     xlab = "ADT Weight",
     xlim = c(0, 1))
abline(v = median(allcombine$ADT.weight), col = "red", lwd = 2, lty = 2)

dev.off()
cat("  Modality weight plots saved\n")

# Save integrated object
cat("\nSaving WNN-integrated object...\n")
saveRDS(allcombine, file = "03_wnn_integrated.rds")

cat("\n✓ Step 3 completed successfully\n")
cat("\nNext: Project onto bone marrow reference map\n")

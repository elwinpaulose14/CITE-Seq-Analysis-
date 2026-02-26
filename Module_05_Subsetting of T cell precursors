##Subsetting of Tcell precursors and T/M MPAL primed cells

print_step_header("IDENTIFICATION OF IMMATURE T-LIKE CELLS", 5)

# Load reference-mapped data (or non-reference data if mapping was skipped)
if (file.exists("04_reference_mapped.rds")) {
  query <- readRDS("04_reference_mapped.rds")
} else if (file.exists("/app/r_analysis/results/04_query_no_reference.rds")) {
  query <- readRDS("/app/r_analysis/results/04_query_no_reference.rds")
} else {
  query <- readRDS("/app/r_analysis/results/03_wnn_integrated.rds")
}

cat("Working with", ncol(query), "cells\n")


# 5.1: Define Gene Signatures for Cell Type Identification

T_prog_genes <- c(
  "CD3D",   # T-cell receptor delta
  "CD3E",   # T-cell receptor epsilon
  "CD7",    # Early T-cell marker
  "TRAC",   # T-cell receptor alpha constant
  "TRBC1",  # T-cell receptor beta constant 1
  "IL7R",   # IL-7 receptor (CD127)
  "LCK",    # Lymphocyte-specific protein tyrosine kinase
  "BCL11B", # T-cell commitment factor
  "TCF7",   # T-cell factor 7 (TCF-1)
  "LEF1"    # Lymphoid enhancer-binding factor 1
)

cat("\nT-cell progenitor genes (n=", length(T_prog_genes), "):\n", sep = "")
cat(paste(T_prog_genes, collapse = ", "), "\n")

# Stem/progenitor signature

Stem_prog_genes <- c(
  "CD34",  # Hematopoietic progenitor marker
  "HOPX",  # Hematopoietic stem cell marker
  "SOX4",  # Lymphoid progenitor transcription factor
  "MEF2C", # Myeloid and lymphoid development
  "LMO2"   # Stem cell factor
)

cat("\nStem/progenitor genes (n=", length(Stem_prog_genes), "):\n", sep = "")
cat(paste(Stem_prog_genes, collapse = ", "), "\n")

# Myeloid progenitor signature

Myeloid_prog_genes <- c(
  "CEBPA", # CCAAT/enhancer binding protein alpha
  "MPO",   # Myeloperoxidase
  "LYZ",   # Lysozyme
  "CTSD",  # Cathepsin D
  "FCN1",  # Ficolin 1
  "CD33",  # Myeloid marker
  "ANPEP", # CD13
  "KIT",   # CD117
  "FUT4"   # CD15
)

cat("\nMyeloid progenitor genes (n=", length(Myeloid_prog_genes), "):\n", sep = "")
cat(paste(Myeloid_prog_genes, collapse = ", "), "\n")


# 5.2: Calculate Module Scores


DefaultAssay(query_seurat) <- "RNA"
query_seurat <- NormalizeData(query_seurat)
cat("\nCalculating T-cell progenitor score...\n")
query_seurat <- AddModuleScore(
  query_seurat,
  features = list(T_prog_genes),
  name = "T_prog",
  ctrl = 20
)

cat("Calculating stem/progenitor score...\n")

query_seurat <- AddModuleScore(
  query_seurat,
  features = list(Stem_prog_genes),
  name = "Stem_prog",
  ctrl = 20
)

cat("Calculating myeloid progenitor score...\n")

query_seurat <- AddModuleScore(
  query_seurat,
  features = list(Myeloid_prog_genes),
  name = "Myeloid_prog",
  ctrl = 20
)

cat("\nModule scores added to metadata\n")
cat("New columns:", grep("prog", colnames(query_seurat@meta.data), value = TRUE), "\n")

# Visualize module score distributions
cat("\nGenerating module score distributions...\n")

png("11_module_score_distributions.png", 
    width = 1400, height = 500, res = 150)
par(mfrow = c(1, 3))

hist(query_seurat$T_prog1, breaks = 50, col = "lightblue",
     main = "T-cell Progenitor Score",
     xlab = "Module Score")
abline(v = 0.3, col = "red", lwd = 2, lty = 2)
abline(v = 0.4, col = "darkred", lwd = 2, lty = 2)

hist(query_seurat$Stem_prog1, breaks = 50, col = "lightgreen",
     main = "Stem/Progenitor Score",
     xlab = "Module Score")
abline(v = 0.1, col = "red", lwd = 2, lty = 2)
abline(v = 0.2, col = "darkred", lwd = 2, lty = 2)

hist(query_seurat$Myeloid_prog1, breaks = 50, col = "lightyellow",
     main = "Myeloid Progenitor Score",
     xlab = "Module Score")
abline(v = 0, col = "red", lwd = 2, lty = 2)

dev.off()

cat("  Module score histograms saved\n")

# 5.3: Subset T Cell Prescussors

# Initialize classification
query_seurat$Immature_T_like <- "Other"

cat("\n--- Classification Criterion 1: T cell precussors ---\n")
cat("Criteria: T_prog > 0.3 AND Stem_prog > 0.1\n")

immature_t_cells <- query_seurat$T_prog1 > 0.3 & query_seurat$Stem_prog1 > 0.1
query_seurat$Immature_T_like[immature_t_cells] <- "Immature_T_like"

cat("Cells classified as Immature_T_like:", sum(immature_t_cells), "\n")

cat("\n--- Classification Criterion 2: Immature T-MPAL-like ---\n")
cat("Criteria: T_prog > 0.4 AND Stem_prog > 0.2 AND Myeloid_prog > 0\n")

mpal_like_cells <- query_seurat$T_prog1 > 0.4 & query_seurat$Stem_prog1 > 0.2 & query_seurat$Myeloid_prog1 > 0
query_seurat$Immature_T_like[mpal_like_cells] <- "Immature_T_MPAL_like"

cat("Cells classified as Immature_T_MPAL_like:", sum(mpal_like_cells), "\n")

# Summary
cat("\n=== Final Classification Summary ===\n")
print(table(query_seurat$Immature_T_like))


# 5.4: Visualize Classification

# Determine which UMAP to use
if ("umap_projected" %in% Reductions(query_seurat)) {
  umap_reduction <- "umap_projected"
} else if ("umap" %in% Reductions(query_seurat)) {
  umap_reduction <- "umap"
} else {
  cat("WARNING: No UMAP reduction found. Skipping visualization.\n")
  umap_reduction <- NULL
}

if (!is.null(umap_reduction)) {
  cat("\nGenerating classification UMAP...\n")
  
  p_classification <- DimPlot(
    query,
    reduction = umap_reduction,
    group.by = "Immature_T_like",
    cols = c(
      "Other" = "grey85",
      "Immature_T_like" = "steelblue",
      "Immature_T_MPAL_like" = "firebrick"
    ),
    pt.size = 0.5
  ) +
    ggtitle("Immature T-like Cell Classification") +
    theme_minimal()
  
  ggsave("/app/r_analysis/figures/12_immature_t_classification.png",
         p_classification, width = 10, height = 7, dpi = 300)
  
  # Plot module scores on UMAP
  cat("\nGenerating module score feature plots...\n")
  
  p_modules <- FeaturePlot(
    query,
    reduction = umap_reduction,
    features = c("T_prog1", "Stem_prog1", "Myeloid_prog1"),
    ncol = 3,
    order = TRUE
  ) & theme_minimal()
  
  ggsave("13_module_scores_umap.png",
         p_modules, width = 18, height = 5, dpi = 300)
}


# 5.5: Export Classification Results

cat("\nExporting classification results...\n")

classification_results <- data.frame(
  Cell_Barcode = colnames(query_seurat),
  T_prog_score = query_seurat$T_prog1,
  Stem_prog_score = query_seurat$Stem_prog1,
  Myeloid_prog_score = query_seurat$Myeloid_prog1,
  Classification = query_seurat$Immature_T_like
)

write.csv(
  classification_results,
  file = "05_immature_t_classification.csv",
  row.names = FALSE
)

cat("  Classification results exported\n")

# Save updated query object
cat("\nSaving classified query object...\n")
saveRDS(query_seurat, file = "05_immature_t_classified.rds")

cat("\n✓ Step 5 completed successfully\n")
cat("\nNext: Cluster and annotate immature T-like cells\n")

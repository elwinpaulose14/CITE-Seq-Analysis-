###Clustering and Annotation

# Load classified data
query <- readRDS("05_immature_t_classified.rds")

cat("Total cells:", ncol(query), "\n")
cat("Immature T-like classification:\n")
print(table(query$Immature_T_like))

# 6.1: Subset Immature T-like Cells

print_step_header("Subsetting Immature T-like Cells")

cat("Subsetting to immature T-like populations...\n")
immT <- subset(query, subset = Immature_T_like != "Other")

cat("Immature T-like cells:", ncol(immT), "\n")

if (ncol(immT) < 50) {
  cat("\nWARNING: Insufficient immature T-like cells for clustering (<50)\n")
  cat("Consider adjusting classification thresholds or analyzing full dataset\n")
  saveRDS(immT, file = "/app/r_analysis/results/06_immature_t_subset.rds")
  quit(save = "no", status = 0)
}

# Remove old reductions to force recalculation
cat("\nRemoving old dimensional reductions...\n")
if ("umap_projected" %in% Reductions(immT)) immT@reductions$umap_projected <- NULL
if ("pca_projected" %in% Reductions(immT)) immT@reductions$pca_projected <- NULL
if ("harmony_projected" %in% Reductions(immT)) immT@reductions$harmony_projected <- NULL
if ("umap" %in% Reductions(immT)) immT@reductions$umap <- NULL
if ("pca" %in% Reductions(immT)) immT@reductions$pca <- NULL


# 6.2: Re-process for High-Resolution Clustering

print_step_header("Re-processing RNA Data")

DefaultAssay(immT) <- "RNA"

cat("Finding variable features (n=2000)...\n")
# Reduced from 3000 to 2000 due to smaller cell number
immT <- FindVariableFeatures(immT, nfeatures = 2000)

cat("  Variable features identified:", length(VariableFeatures(immT)), "\n")

cat("Scaling data...\n")
immT <- ScaleData(immT)

cat("Running PCA (npcs=20)...\n")
# Reduced PCs for a smaller dataset
immT <- RunPCA(immT, npcs = 20)

# Visualize PCA variance
cat("\nGenerating elbow plot for PC selection...\n")
png("14_immature_t_elbow_plot.png", 
    width = 800, height = 600, res = 150)
ElbowPlot(immT, ndims = 20)
dev.off()

cat("  Elbow plot suggests using 15 PCs\n")


# 6.3: Graph-Based Clustering

print_step_header("Graph-Based Clustering")

cat("Building nearest neighbor graph (dims 1:15)...\n")
immT <- FindNeighbors(immT, dims = 1:15)

cat("Finding clusters (resolution=1.1)...\n")
# Resolution justification:

immT <- FindClusters(immT, resolution = 1.1)

cat("\nClusters identified:\n")
print(table(immT$seurat_clusters))


# 6.4: UMAP Visualization

print_step_header("UMAP Visualization")

cat("Running UMAP (dims 1:15)...\n")
immT <- RunUMAP(immT, dims = 1:15, reduction.name = "umap_clean")

cat("Generating cluster UMAP...\n")
p_clusters <- DimPlot(
  immT,
  reduction = "umap_clean",
  label = TRUE,
  label.size = 6,
  pt.size = 1.5
) +
  ggtitle("Immature T-like Cells - Unsupervised Clustering") +
  theme_minimal()

ggsave("15_immature_t_clusters.png",
       p_clusters, width = 9, height = 7, dpi = 300)

# 6.5: Marker Gene Expression Analysis
print_step_header("Marker Gene Expression")

# Define stage-specific markers
stage_markers <- c(
  "TCF7",    # Early T-cell progenitor
  "BCL11B",  # T-lineage commitment
  "IL7R",    # Common lymphoid progenitor
  "CD3D",    # Mature T-cell
  "CD4",     # Helper T-cell
  "CD8A"     # Cytotoxic T-cell
)

cat("\nKey developmental markers:\n")
cat(paste(stage_markers, collapse = ", "), "\n")

# Filter available markers
available_markers <- stage_markers[stage_markers %in% rownames(immT)]

if (length(available_markers) > 0) {
  cat("\nGenerating marker feature plots...\n")
  
  p_markers <- FeaturePlot(
    immT,
    reduction = "umap_clean",
    features = available_markers,
    ncol = 3,
    order = TRUE
  ) & theme_minimal()
  
  ggsave("16_developmental_markers.png",
         p_markers, 
         width = 15, 
         height = ceiling(length(available_markers)/3) * 4, 
         dpi = 300)
  
  # Calculate average expression per cluster
  cat("\nCalculating average expression per cluster...\n")
  
  avg_exp <- AverageExpression(
    immT,
    features = available_markers,
    return.seurat = FALSE
  )
  
  cat("\nAverage expression matrix:\n")
  print(round(avg_exp$RNA, 2))
  
  # Export to CSV
  write.csv(
    avg_exp$RNA,
    file = "06_cluster_marker_expression.csv"
  )
}


# 6.6: Developmental Stage Annotation

print_step_header("Developmental Stage Annotation")

cat("\nAnnotating clusters based on marker expression...\n\n")

# Initialize stage annotation
immT$Stage <- "Unassigned"

# Annotation logic (customize based on your average expression matrix)
# These are example annotations - adjust based on actual cluster characteristics

cat("Annotation Criteria:\n")
cat("-------------------\n")

cat("Cluster 0 -> DP (Double Positive)\n")

immT$Stage[immT$seurat_clusters == 0] <- "DP"

immT$Stage[immT$seurat_clusters == 1] <- "CLP"

immT$Stage[immT$seurat_clusters == 2] <- "DN"

immT$Stage[immT$seurat_clusters == 3] <- "ETP"

# Handle additional clusters if present
if (length(unique(immT$seurat_clusters)) > 4)

cat("\n=== Final Stage Annotation Summary ===\n")
print(table(immT$Stage))

# Visualize annotation
cat("\nGenerating annotated UMAP...\n")

p_annotated <- DimPlot(
  immT,
  reduction = "umap_clean",
  group.by = "Stage",
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  pt.size = 1.5,
  cols = c(
    "ETP" = "#E64B35FF",
    "DN" = "#4DBBD5FF",
    "DP" = "#00A087FF",
    "CLP" = "#3C5488FF",
    "Unassigned" = "grey70"
  )
) +
  ggtitle("T-cell Developmental Stages") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("17_developmental_stages_annotated.png",
       p_annotated, width = 10, height = 7, dpi = 300)


# 6.7: Protein (ADT) Expression Validation

print_step_header("ADT Validation")

if ("ADT" %in% Assays(immT)) {
  DefaultAssay(immT) <- "ADT"
  
  cat("\nAvailable ADT features:\n")
  print(rownames(immT[["ADT"]]))
  
  # Key surface markers
  adt_validation_markers <- c(
    "CD3:UCHT1-CD3E-AHS0231-pAbO",
    "CD4:SK3-CD4-AHS0032-pAbO",
    "CD8:RPA-T8-CD8A-AHS0027-pAbO",
    "CD34:581-CD34-AHS0061-pAbO"
  )
  
  adt_available <- adt_validation_markers[adt_validation_markers %in% rownames(immT[["ADT"]])]
  
  if (length(adt_available) > 0) {
    cat("\nGenerating ADT validation plots...\n")
    
    p_adt <- FeaturePlot(
      immT,
      reduction = "umap_clean",
      features = adt_available,
      min.cutoff = "q10",
      max.cutoff = "q90",
      ncol = 2
    ) & theme_minimal()
    
    ggsave("18_adt_validation.png",
           p_adt, 
           width = 12, 
           height = ceiling(length(adt_available)/2) * 5, 
           dpi = 300)
  }
  
  # Reset to RNA
  DefaultAssay(immT) <- "RNA"
}


# 6.8: Export Results

print_step_header("Exporting Results")

cat("\nExporting cell annotations...\n")

annotation_results <- data.frame(
  Cell_Barcode = colnames(immT),
  Cluster = immT$seurat_clusters,
  Stage = immT$Stage,
  Classification = immT$Immature_T_like,
  T_prog_score = immT$T_prog1,
  Stem_prog_score = immT$Stem_prog1,
  nFeature_RNA = immT$nFeature_RNA,
  nCount_RNA = immT$nCount_RNA,
  percent_mt = immT$percent.mt
)

write.csv(
  annotation_results,
  file = "06_final_annotations.csv",
  row.names = FALSE
)

cat("  Annotations exported\n")

# Save final object
cat("\nSaving final annotated object...\n")
saveRDS(immT, file = "06_immature_t_final.rds")

cat("✓ ANALYSIS PIPELINE COMPLETED SUCCESSFULLY\n")

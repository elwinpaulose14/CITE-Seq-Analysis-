Reference Mapping to Bone Marrow Atlas

# Load integrated data
allcombine <- readRDS("03_wnn_integrated.rds")

cat("Working with", ncol(allcombine), "cells\n")

# 4.1: Load Bone Marrow Reference Atlas

projection_path <- "//10.100.75.55/NASShare/ELWIN/BDRhapsody"
ref_file <- file.path(projection_path, "BoneMarrow_RefMap_SymphonyRef.rds")

if (!file.exists(ref_file)) {
  cat("WARNING: Reference file not found at:", ref_file, "\n")
  cat("Proceeding with analysis pipeline without reference mapping\n")
  cat("NOTE: For publication, ensure bone marrow reference is available\n")
  cat("      Reference can be obtained from BoneMarrowMap package\n")
  
  # Save query object and continue
  query <- allcombine
  saveRDS(query, file = "04_query_no_reference.rds")
  
  cat("\n✓ Step 4 completed (without reference mapping)\n")
  quit(save = "no", status = 0)
}

cat("Loading bone marrow reference atlas...\n")
ref <- readRDS(ref_file)

# Set UMAP model path for Symphony projection
ref$save_uwot_path <- file.path(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot')

cat("Creating reference object for mapping...\n")

ReferenceSeuratObj <- create_ReferenceObject(ref)

cat("\nReference atlas cell types:\n")
if ("CellType_Annotation_formatted" %in% colnames(ReferenceSeuratObj@meta.data)) {
  print(table(ReferenceSeuratObj$CellType_Annotation_formatted))
}

# Visualize reference atlas
cat("\nGenerating reference atlas UMAP...\n")
p_ref <- DimPlot(
  ReferenceSeuratObj,
  reduction = 'umap',
  group.by = 'CellType_Annotation_formatted',
  raster = FALSE,
  label = TRUE,
  label.size = 3
) +
  ggtitle("Bone Marrow Reference Atlas") +
  theme_minimal()

ggsave("06_reference_atlas_celltypes.png",
       p_ref, width = 12, height = 8, dpi = 300)

# Plot cell cycle and pseudotime
if ("CyclePhase" %in% colnames(ReferenceSeuratObj@meta.data) &&
    "Pseudotime" %in% colnames(ReferenceSeuratObj@meta.data)) {
  
  p1 <- DimPlot(ReferenceSeuratObj,
                reduction = 'umap',
                group.by = 'CyclePhase',
                raster = FALSE) +
    ggtitle("Cell Cycle Phase")
  
  p2 <- FeaturePlot(ReferenceSeuratObj,
                    reduction = 'umap',
                    features = 'Pseudotime',
                    raster = FALSE) +
    ggtitle("Developmental Pseudotime")
  
  p_combined <- p1 + p2
  
  ggsave("07_reference_atlas_metadata.png",
         p_combined, width = 16, height = 6, dpi = 300)
}


# 4.2: Map Query Data to Reference

cat("Extracting query RNA counts...\n")
query_counts <- GetAssayData(allcombine, assay = "RNA", layer = "counts")

cat("  Query dimensions:", dim(query_counts), "\n")

batchvar <- "Sample_Name"
if (!batchvar %in% colnames(allcombine@meta.data)) {
  cat("WARNING: Batch variable", batchvar, "not found. Using default.\n")
  batchvar <- NULL
}

cat("\nProjecting query onto reference...\n")


query <- mapQuery(
  exp_query = query_counts,
  metadata_query = allcombine@meta.data,
  ref_obj = ref,
  vars = "Sample_Name"
)

# 4.3 Converting Symphony output into Seuratobject
query_seurat <- CreateSeuratObject(
  counts = query_counts,
  meta.data = query$meta_data
)

query_seurat[["ref.umap"]] <- CreateDimReducObject(
  embeddings = query$umap,
  key = "refUMAP_",
  assay = "RNA"
)
query_seurat[["ADT"]] <- CreateAssayObject(counts = adt_counts)

cat("  Projection complete\n")

cat("\nAvailable reductions in query:\n")
Reductions(query_seurat)

# 4.3: Transfer ADT Data to Query Object

cat("Extracting ADT counts from original object...\n")
adt_counts <- GetAssayData(allcombine, assay = "ADT", layer = "counts")

cat("  ADT features:", nrow(adt_counts), "\n")
cat("  ADT cells:", ncol(adt_counts), "\n")

# Verify cell barcode alignment
if (!all(colnames(adt_counts) == colnames(query))) {
  stop("ERROR: Cell barcodes do not match between ADT and query object")
}

cat("  Cell barcode verification: PASS\n")

cat("Creating ADT assay in query object...\n")
query[["ADT"]] <- CreateAssayObject(counts = adt_counts)

cat("\nQuery object assays:\n")


# 4.4: Normalize ADT in Query Object

cat("\nNormalizing ADT data in query object...\n")
DefaultAssay(query_seurat) <- "ADT"

query <- NormalizeData(
  query_seurat,
  normalization.method = "CLR",
  margin = 2
)

query <- ScaleData(query)

cat("ADT normalization complete\n")


# 4.5: Visualize Query Projection


# Switch to RNA for visualization
DefaultAssay(query) <- "RNA"

cat("Generating projected UMAP plots...\n")

# Plot by predicted cell types (if available)
if ("predicted.celltype" %in% colnames(query@meta.data)) {
  p_projected <- DimPlot(
    query,
    reduction = "umap_projected",
    group.by = "predicted.celltype",
    label = TRUE,
    label.size = 3,
    repel = TRUE
  ) +
    ggtitle("Query Cells - Predicted Cell Types") +
    theme_minimal()
  
  ggsave("08_query_projected_celltypes.png",
         p_projected, width = 12, height = 8, dpi = 300)
}

# Plot selected ADT markers
DefaultAssay(query) <- "ADT"

cat("\nGenerating ADT feature plots on projected UMAP...\n")

adt_markers <- c(
  "CD34:581-CD34-AHS0061-pAbO",
  "CD7-CD7-AHS0043-pAbO",
  "CD3:UCHT1-CD3E-AHS0231-pAbO"
)

# Filter markers that exist
adt_markers_available <- adt_markers[adt_markers %in% rownames(query_seurat[["ADT"]])]

if (length(adt_markers_available) > 0) {
  p_adt <- FeaturePlot(
    query_seurat,
    features = adt_markers_available,
    reduction = "ref.umap",
    ncol = 3
  )
  
  ggsave("09_query_adt_markers.png",
         p_adt, width = 15, height = 5, dpi = 300)
}

# Compare RNA and ADT for key marker
DefaultAssay(query_seurat) <- "RNA"

if ("CEBPA" %in% rownames(query)) {
  p1 <- FeaturePlot(
    query_seurat,
    features = "CEBPA",
    reduction = "umap_projected"
  ) + ggtitle("RNA: CEBPA")
  
  DefaultAssay(query) <- "ADT"
  
  if ("CD34:581-CD34-AHS0061-pAbO" %in% rownames(query[["ADT"]])) {
    p2 <- FeaturePlot(
      query,
      features = "CD34:581-CD34-AHS0061-pAbO",
      reduction = "umap_projected"
    ) + ggtitle("Protein: CD34")
    
    p_compare <- p1 | p2
    
    ggsave("10_rna_vs_adt_comparison.png",
           p_compare, width = 14, height = 6, dpi = 300)
  }
}

# Save query object
cat("\nSaving reference-mapped query object...\n")
saveRDS(query, file = "04_reference_mapped.rds")

cat("\n✓ Step 4 completed successfully\n")
cat("\nNext: Identify T cell precussors\n")

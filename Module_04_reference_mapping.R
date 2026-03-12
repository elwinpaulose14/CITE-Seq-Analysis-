#Reference Mapping to Bone Marrow Atlas

# Load integrated data
allcombine <- readRDS("03_wnn_integrated.rds")

cat("Working with", ncol(allcombine), "cells\n")

# 4.1: Load Bone Marrow Reference Atlas
projection_path <- "\\\\10.100.75.55/NASShare/ELWIN/BDRhapsody/"
ref <- readRDS(file.path(projection_path, "BoneMarrow_RefMap_SymphonyRef.rds"))
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot')
ReferenceSeuratObj <- create_ReferenceObject(ref)                      
DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', raster = FALSE, label = TRUE, label.size = 4)
p1 <- DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CyclePhase', raster = FALSE)
p2 <- FeaturePlot(ReferenceSeuratObj, reduction = 'umap', features = 'Pseudotime', raster = FALSE)
p1 + p2

# 4.2: Map Query Data to Reference

cat("Extracting query RNA counts...\n")
batchvar <- 'Sample_Name'
query <- map_Query(
  exp_query = GetAssayData(allcombine, assay = "RNA", slot = "counts"),
  metadata_query = allcombine@meta.data,
  ref_obj = ref,
  vars = "Sample_Name"
)
Reductions(query)
query[["ADT"]] <- allcombine[["ADT"]]
class(query)
Assays(query)

adt_counts <- GetAssayData(
  allcombine,
  assay = "ADT",
  slot = "counts"
)
all(colnames(adt_counts) == colnames(query))
# Should be TRUE

query[["ADT"]] <- CreateAssayObject(
  counts = adt_counts
)

names(query@assays)
# should show: RNA ADT

# 4.3: Normalize ADT in Query Object

cat("\nNormalizing ADT data in query object...\n")
DefaultAssay(query) <- "ADT"

query <- NormalizeData(
  query,
  normalization.method = "CLR",
  margin = 2
)

query <- ScaleData(query)
FeaturePlot(
  query,
  features = "CD3:UCHT1-CD3E-AHS0231-pAbO",
  reduction = "umap_projected"
)
FeaturePlot(
  query,
  features = c(
    "CD34:581-CD34-AHS0061-pAbO",
    "CD7-CD7-AHS0043-pAbO",
    "CD3:UCHT1-CD3E-AHS0231-pAbO"
  ),
  reduction = "umap_projected",
  ncol = 3
)

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


# Compare RNA and ADT for key marker
DefaultAssay(query_seurat) <- "RNA"

if ("CEBPA" %in% rownames(query)) {
  p1 <- FeaturePlot(
    query,
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

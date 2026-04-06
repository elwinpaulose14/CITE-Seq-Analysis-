# ===============================
# CLEAN MODULE 4: Mapping + Visualization + Precursor Detection
# ===============================

library(Seurat)
library(patchwork)
library(FNN)

# -------------------------------
# LOAD DATA
# -------------------------------
allcombine <- readRDS("03_wnn_integrated.rds")

projection_path <- "\\\\10.100.75.55/NASShare/ELWIN/BDRhapsody/"
ref <- readRDS(file.path(projection_path, "BoneMarrow_RefMap_SymphonyRef.rds"))
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot')

ReferenceSeuratObj <- create_ReferenceObject(ref)

# -------------------------------
# MAP QUERY
# -------------------------------
batchvar <- "Sample_Name"

query <- map_Query(
  exp_query = GetAssayData(allcombine, assay = "RNA", slot = "counts"),
  metadata_query = allcombine@meta.data,
  ref_obj = ref,
  vars = batchvar
)

# -------------------------------
# ADD ADT
# -------------------------------
query[["ADT"]] <- allcombine[["ADT"]]

DefaultAssay(query) <- "ADT"
query <- NormalizeData(query, normalization.method = "CLR", margin = 2)
query <- ScaleData(query)

# -------------------------------
# TRANSFER REFERENCE LABELS (KNN)
# -------------------------------
ref_embed <- t(ref$Z_corr)
query_embed <- Embeddings(query, "pca_projected")

# Ensure same dimensions
query_embed <- query_embed[, 1:ncol(ref_embed)]

nn <- get.knnx(ref_embed, query_embed, k = 5)

ref_labels <- ref$meta_data$CellType_Annotation_formatted

query$predicted_CellType <- apply(nn$nn.index, 1, function(idx) {
  names(sort(table(ref_labels[idx]), decreasing = TRUE))[1]
})

# -------------------------------
# OUTPUT FOLDER
# -------------------------------
outdir <- "Step4_Visualizations"
if (!dir.exists(outdir)) dir.create(outdir)

# -------------------------------
# RNA FEATURES
# -------------------------------
DefaultAssay(query) <- "RNA"

rna_features <- c(
  "TCF7", "LEF1", "IL7R", "CCR7",
  "GATA3", "BCL11B", "SELL",
  "MPO", "LYZ", "CD33", "CEBPA"
)

rna_features <- rna_features[rna_features %in% rownames(query)]

p_rna <- FeaturePlot(
  query,
  features = rna_features,
  reduction = "umap_projected",
  min.cutoff = 0,
  max.cutoff = 3,
  cols = c("grey90", "red"),
  ncol = 4
)

ggsave(file.path(outdir, "RNA_expression.png"), p_rna, width = 14, height = 10)

# -------------------------------
# ADT FEATURES (SAME SCALE 0–3)
# -------------------------------
DefaultAssay(query) <- "ADT"

adt_features <- c(
  "CD34:581-CD34-AHS0061-pAbO",
  "CD3:UCHT1-CD3E-AHS0231-pAbO",
  "CD7-CD7-AHS0043-pAbO",
  "CD33:WM53-CD33-AHS0044-pAbO",
  "CD13-ANPEP-AHS0050-pAbO",
  "CD15-FUT4-AHS0196-pAbO"
)

adt_features <- adt_features[adt_features %in% rownames(query[["ADT"]])]

# Clamp ADT to 0–3
adt_mat <- GetAssayData(query, slot = "scale.data")
adt_mat[adt_mat < 0] <- 0
adt_mat[adt_mat > 3] <- 3
query[["ADT"]]@scale.data <- adt_mat

p_adt <- FeaturePlot(
  query,
  features = adt_features,
  reduction = "umap_projected",
  slot = "scale.data",
  cols = c("grey90", "blue"),
  ncol = 3
)

ggsave(file.path(outdir, "ADT_expression.png"), p_adt, width = 14, height = 10)

# -------------------------------
# ANNOTATED UMAP
# -------------------------------
DefaultAssay(query) <- "RNA"

p_annot <- DimPlot(
  query,
  reduction = "umap_projected",
  group.by = "predicted_CellType",
  label = TRUE,
  repel = TRUE
)

ggsave(file.path(outdir, "Annotated_UMAP.png"), p_annot, width = 10, height = 8)

# -------------------------------
# PRECURSOR DETECTION
# -------------------------------

# RNA
rna_data <- FetchData(query, vars = c("TCF7", "LEF1", "IL7R", "MPO", "LYZ", "CD33"))

# ADT
DefaultAssay(query) <- "ADT"
adt_data <- FetchData(query, vars = "CD34:581-CD34-AHS0061-pAbO", slot = "scale.data")

df <- cbind(rna_data, CD34 = adt_data[,1])

df$CD34_pos <- df$CD34 > 1
df$Tcell_program <- rowMeans(df[, c("TCF7","LEF1","IL7R")], na.rm = TRUE) > 0.5
df$Myeloid_negative <- rowMeans(df[, c("MPO","LYZ","CD33")], na.rm = TRUE) < 0.5

df$precursor <- "Other"
df$precursor[df$CD34_pos & df$Tcell_program & df$Myeloid_negative] <- "Precursor"

query$precursor_status <- df$precursor

# -------------------------------
# SAVE PRECURSOR PLOT
# -------------------------------
p_precursor <- DimPlot(
  query,
  reduction = "umap_projected",
  group.by = "precursor_status",
  cols = c("grey", "red")
)

ggsave(file.path(outdir, "Precursor_cells.png"), p_precursor, width = 10, height = 8)

# -------------------------------
# COUNT PRECURSORS
# -------------------------------
cat("\nPrecursor counts:\n")
print(table(query$precursor_status))

cat("\nPercentage:\n")
print(prop.table(table(query$precursor_status))*100)

# -------------------------------
# SAVE OBJECT
# -------------------------------
saveRDS(query, file = "04_reference_mapped_clean.rds")

cat("\n✓ CLEAN MODULE 4 COMPLETED\n")

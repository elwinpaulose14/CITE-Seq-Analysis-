# CITE-Seq-Analysis

```text
# CITE-Seq Analysis Pipeline

> Multimodal single-cell analysis of RNA and protein co-expression for T-cell progenitor identification and developmental staging in bone marrow

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.2-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5-green.svg)](https://seuratproject.org/)

---

## Overview

This pipeline performs an end-to-end single-cell multimodal analysis of CITE-seq data from cryopreserved bone marrow samples, integrating simultaneous RNA and surface protein (ADT) expression. The primary goal is to identify **T-cell progenitors** and **T/Myeloid MPAL-like cells**, map , and resolve their developmental trajectory into ETP → DN → DP stages.

### Key Features

- Simultaneous RNA + protein (ADT) processing via Weighted Nearest Neighbor (WNN) integration
- Symphony-based reference mapping onto a curated bone marrow atlas
- Gene module scoring for T-cell, stem, and myeloid progenitor programs
- High-resolution clustering and developmental stage annotation (ETP, DN, DP, CLP)


---

## Workflow

```text
Raw CITE-seq Data
       │
       ▼

┌─────────────────────┐
│  Module 01          │
│  Data Loading       │
└────────┬────────────┘
         ├─ Load RDS + Sample Tag Metadata
         └─ Filter: Undetermined / Multiplets

         ▼

┌─────────────────────┐
│  Module 02          │
│  Quality Control    │
└────────┬────────────┘
         ├─ Mitochondrial % calculation
         ├─ Dynamic quantile thresholds (1st–99th)
         └─ nFeature: 150–2700 | MT < 45%

         ▼

┌─────────────────────┐
│  Module 03          │
│  Multimodal         │
│  Integration (WNN)  │
└────────┬────────────┘
         ├─ RNA: NormalizeData → PCA (50 PCs)
         ├─ ADT: CLR normalization → PCA (9 PCs)
         └─ WNN graph: RNA dims 1:30, ADT dims 1:9

         ▼

┌─────────────────────┐
│  Module 04          │
│  Reference Mapping  │
└────────┬────────────┘
         ├─ Symphony reference mapping
         ├─ KNN label transfer (k = 5)
         └─ Precursor detection: CD34+ / TCF7+ / MPO−

         ▼

┌─────────────────────┐
│  Module 05          │
│  Progenitor         │
│  Identification     │
└────────┬────────────┘
         ├─ AddModuleScore: T_prog
         ├─ AddModuleScore: Stem_prog
         ├─ AddModuleScore: Myeloid_prog
         └─ Classify: Immature_T_like / Immature_T_MPAL_like

         ▼

┌─────────────────────┐
│  Module 06          │
│  Clustering &       │
│  Annotation         │
└────────┬────────────┘
         ├─ PCA (20 PCs)
         ├─ FindNeighbors → FindClusters
         ├─ UMAP visualization
         ├─ Stage labels: ETP / DN / DP / CLP
         └─ ADT validation + export
```

---

## Directory Structure

```text
r_analysis/
├── data/
│   ├── CITE-Seq-Lane3_Seurat.rds
│   └── CITE-Seq-Lane3_Sample_Tag_Calls.csv
│
├── scripts/
│   ├── Module_01_data_loaded.R
│   ├── Module_02_qc_filtered.R
│   ├── Module_03_wnn_integrated.R
│   ├── Module_04_reference_mapping.R
│   ├── Module_05_progenitor_identification.R
│   └── Module_06_clustering_annotation.R
│
├── Main_Pipeline_Target.R
│
├── results/
│   ├── 01_data_loaded.rds
│   ├── 02_qc_filtered.rds
│   ├── 03_wnn_integrated.rds
│   ├── 04_reference_mapped_clean.rds
│   ├── 05_immature_t_classified.rds
│   ├── 06_immature_t_final.rds
│   ├── 05_immature_t_classification.csv
│   ├── 06_cluster_marker_expression.csv
│   └── 06_final_annotations.csv
│
└── figures/
    ├── 01_qc_before_filtering.png
    ├── 03_rna_elbow_plot.png
    ├── 04_adt_elbow_plot.png
    ├── 05_modality_weights.png
    ├── 11_module_score_distributions.png
    ├── 12_immature_t_classification.png
    ├── 13_module_scores_umap.png
    ├── 14_immature_t_elbow_plot.png
    ├── 15_immature_t_clusters.png
    ├── 16_developmental_markers.png
    ├── 17_developmental_stages_annotated.png
    └── 18_adt_validation.png
```

---

## Dependencies

### R Packages

| Package | Version | Purpose |
|---------|---------|---------|
| `Seurat` | ≥ 5.0 | Core single-cell analysis framework |
| `SeuratDisk` | ≥ 0.0.5 | H5Seurat / RDS I/O |
| `BoneMarrowMap` | latest | Reference atlas and Symphony mapping utilities |
| `symphony` | ≥ 0.1.0 | Reference-based cell embedding |
| `AUCell` | ≥ 1.20 | Gene set activity scoring |
| `FNN` | ≥ 1.1 | Fast k-nearest neighbor search |
| `patchwork` | ≥ 1.1 | Plot composition |
| `ggplot2` | ≥ 3.4 | Visualization |
| `tidyverse` | ≥ 2.0 | Data wrangling and transformation |
| `BiocManager` | ≥ 1.30 | Bioconductor package management |
| `targets` | ≥ 1.0 | Pipeline orchestration |
| `tarchetypes` | ≥ 0.7 | Targets helper functions |

### Installation

```r
# CRAN packages
install.packages(c(
  "Seurat", "patchwork", "FNN", "ggplot2",
  "tidyverse", "targets", "tarchetypes"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("AUCell"))

# GitHub packages
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("immunogenomics/symphony")

#For Reference Bone Marrow Map
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_SymphonyRef.rds', destfile = paste0(projection_path, 'BoneMarrow_RefMap_SymphonyRef.rds'))
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot', destfile = paste0(projection_path, 'BoneMarrow_RefMap_uwot_model.uwot'))
```

---

## Configuration

Before running, update the network path to your NAS share in `Module_01_data_loaded.R` and `Module_04_reference_mapping.R`:

```r
# Module_01_data_loaded.R
data_dir <- "//YOUR_NAS/Share/BDRhapsody"

# Module_04_reference_mapping.R
projection_path <- "//YOUR_NAS/Share/BDRhapsody/"
```

Ensure the following files are present on the NAS:

- `CITE-Seq-Lane3_Seurat.rds`
- `CITE-Seq-Lane3_Sample_Tag_Calls.csv`
- `BoneMarrow_RefMap_SymphonyRef.rds`
- `BoneMarrow_RefMap_uwot_model.uwot`

---

## Usage

### Run the Full Pipeline

```r
source("Main_Pipeline_Target.R")
```

This will execute all six modules in sequence using the `{targets}` framework, with automatic dependency tracking and caching. If a module's inputs haven't changed, it will be skipped on re-runs.

### Run Individual Modules

Each module is self-contained and reads/writes intermediate `.rds` files. Run them sequentially:

```r
source("scripts/Module_01_data_loaded.R")       # → 01_data_loaded.rds
source("scripts/Module_02_qc_filtered.R")        # → 02_qc_filtered.rds
source("scripts/Module_03_wnn_integrated.R")     # → 03_wnn_integrated.rds
source("scripts/Module_04_reference_mapping.R")  # → 04_reference_mapped_clean.rds
source("scripts/Module_05_Subsetting of T cell precursors.R")  # → 05_immature_t_classified.rds
source("scripts/Module_06_Clustering_Annotation.R")            # → 06_immature_t_final.rds
```

### Skip Reference Mapping (Module 04)

If a Symphony reference is unavailable, Module 05 falls back to the WNN-integrated object:

```r
# In Module_05, the pipeline automatically attempts:
query <- readRDS("04_reference_mapped_clean.rds")
# Falls back to:
query <- readRDS("03_wnn_integrated.rds")
```

---

## Pipeline Steps

### Module 01 — Data Loading

Reads the raw Seurat RDS object and merges sample demultiplexing metadata from the BD Rhapsody Sample Tag CSV. Validates file paths before loading.

**Output:** `01_data_loaded.rds`

---

### Module 02 — Quality Control

Computes per-cell mitochondrial percentage and applies a multi-step filtering strategy:

- Dynamic quantile filtering (1st–99th percentile on `nFeature_RNA`)
- Hard thresholds: `nFeature_RNA` > 150 and < 2700; `percent.mt` < 45%
- Removal of `Undetermined` and `Multiplet` sample tags
- Retention of `SampleTag11_hs` only

**Output:** `02_qc_filtered.rds`, `01_qc_before_filtering.png`

---

### Module 03 — Multimodal Integration (WNN)

Processes both RNA and ADT modalities independently before joint integration:

- **RNA:** NormalizeData → FindVariableFeatures (n=3000) → ScaleData → PCA (50 PCs)
- **ADT:** CLR normalization (margin=2) → ScaleData → PCA (9 PCs)
- **WNN:** `FindMultiModalNeighbors` combining RNA dims 1:30 and ADT dims 1:9

**Output:** `03_wnn_integrated.rds`, elbow plots, modality weight histograms

---

### Module 04 — Reference Mapping

Projects query cells onto a curated bone marrow Symphony reference atlas:

- `map_Query()` with batch correction by `Sample_Name`
- KNN label transfer (k=5) from reference `CellType_Annotation_formatted`
- CD34+ / TCF7+ / MPO− precursor detection logic
- Feature plots for RNA markers and ADT surface proteins

**Output:** `04_reference_mapped_clean.rds`, UMAP plots, `Precursor_cells.png`

---

### Module 05 — Progenitor Identification

Scores each cell against three gene programs using `AddModuleScore`:

| Program | Key Genes | Threshold |
|---------|-----------|-----------|
| T-cell progenitor | TCF7, LEF1, IL7R, BCL11B, LCK | > 0.3 |
| Stem/progenitor | CD34, HOPX, SOX4, MEF2C, LMO2 | > 0.1 |
| Myeloid progenitor | MPO, LYZ, CEBPA, CD33, KIT | > 0 |

Cell classifications:
- **Tcell_prog:** T_prog > 0.3 AND Stem_prog > 0.1
- **T_M_MPAL_like:** T_prog > 0.4 AND Stem_prog > 0.2 AND Myeloid_prog > 0.1

**Output:** `05_immature_t_classified.rds`, `05_immature_t_classification.csv`

---

### Module 06 — Clustering and Annotation

High-resolution re-analysis of the immature T-like subset:

- PCA (20 PCs) → FindNeighbors (dims 1:15) → FindClusters (resolution=1.1)
- UMAP embedding (`umap_clean`)
- Developmental stage annotation: ETP, DN, DP, CLP
- ADT surface marker validation (CD3, CD4, CD8, CD34)
- Per-cluster average expression matrix export

**Output:** `06_immature_t_final.rds`, `06_final_annotations.csv`, annotated UMAPs

---

## Expected Outputs

### Key Result Files

| File | Description |
|------|-------------|
| `06_immature_t_final.rds` | Final annotated Seurat object |
| `06_final_annotations.csv` | Per-cell barcodes, cluster, stage, module scores |
| `06_cluster_marker_expression.csv` | Average marker expression per cluster |
| `05_immature_t_classification.csv` | Module scores and classification labels |

---

## Troubleshooting

**NAS path not found**
Ensure the UNC path (`\\server\share`) is accessible from your machine. On Linux/macOS, mount the share with CIFS/SMB before running.

**`map_Query` fails**
Verify the Symphony reference `.rds` and `.uwot` model files are present at `projection_path`. The uwot path must be set correctly before calling `create_ReferenceObject`.

**Too few cells after QC**
The thresholds (`nFeature_RNA` 150–2700, `percent.mt` < 45%) are tuned for bone marrow CITE-seq. Adjust in Module 02 if your dataset has different characteristics.

**No cells classified as Immature_T_MPAL_like**
The MPAL classification uses strict thresholds. Inspect the module score distributions in `11_module_score_distributions.png` and adjust cutoffs in Module 05 accordingly.

---

**Reference:**
>Zeng AGX, Iacobucci I, Shah S, Mitchell A, Wong G, Bansal S, Chen D, Gao Q, Kim H, Kennedy JA, Arruda A, Minden MD, Haferlach T, Mullighan CG, Dick JE. Single-cell transcriptional mapping reveals genetic and non-genetic determinants of aberrant differentiation in AML. Nature Genetics. 2022. https://doi.org/10.1038/s41588-022-01043-0


**Core tools:**
> Hao, Y., et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573–3587. https://doi.org/10.1016/j.cell.2021.04.048

> Kang, J.B., et al. (2021). Efficient and precise single-cell reference atlas mapping with Symphony. *Nature Communications*, 12, 5890. https://doi.org/10.1038/s41467-021-25957-x

> Stoeckius, M., et al. (2017). Simultaneous epitope and transcriptome measurement in single cells. *Nature Methods*, 14, 865–868. https://doi.org/10.1038/nmeth.4380

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

*Developed by Elwin Paulose · 2026*
```

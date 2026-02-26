# CITE-Seq-Analysis

## Project Structure

```text
r_analysis/
├── data/
│   ├── CITE-Seq-Lane3_Seurat.rds
│   └── CITE-Seq-Lane3_Sample_Tag_Calls.csv
│
├── scripts/
│   ├── helper_functions.R
│   ├── 01_data_loading.R
│   ├── 02_quality_control.R
│   ├── 03_multimodal_integration.R
│   ├── 04_reference_mapping.R
│   ├── 05_immature_t_identification.R
│   └── 06_clustering_annotation.R
│
├── main_pipeline.R
│
├── results/
│   ├── 01_data_loaded.rds
│   ├── 02_qc_filtered.rds
│   ├── 03_wnn_integrated.rds
│   ├── 04_reference_mapped.rds
│   ├── 05_immature_t_classified.rds
│   ├── 06_immature_t_final.rds
│   ├── 06_final_annotations.csv
│   └── session_info.txt
│
└── figures/
    ├── 01_qc_before_filtering.png
    ├── 02_qc_after_filtering.png
    ├── 03_rna_elbow_plot.png
    ├── 15_immature_t_clusters.png
    └── 17_developmental_stages_annotated.png
```

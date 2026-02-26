# CITE-Seq-Analysis
r_analysis/
├── data/                                   # Input data files
│   ├── CITE-Seq-Lane3_Seurat.rds           # Seurat object with RNA + ADT
│   └── CITE-Seq-Lane3_Sample_Tag_Calls.csv # Sample metadata
│
├── scripts/                                 # Analysis scripts
│   ├── helper_functions.R                  # Utility functions
│   ├── 01_data_loading.R                   # Data import
│   ├── 02_quality_control.R                # QC filtering
│   ├── 03_multimodal_integration.R         # RNA + ADT WNN
│   ├── 04_reference_mapping.R              # BoneMarrowMap projection
│   ├── 05_immature_t_identification.R      # Gene signature classification
│   └── 06_clustering_annotation.R          # Developmental stage annotation
│
├── main_pipeline.R                         
│
├── results/                                # Output files (.rds, .csv)
│   ├── 01_data_loaded.rds
│   ├── 02_qc_filtered.rds
│   ├── 03_wnn_integrated.rds
│   ├── 04_reference_mapped.rds
│   ├── 05_immature_t_classified.rds
│   ├── 06_immature_t_final.rds
│   ├── 06_final_annotations.csv
│   └── session_info.txt
│
├── figures/                                
│   ├── 01_qc_before_filtering.png
│   ├── 02_qc_after_filtering.png
│   ├── 03_rna_elbow_plot.png
│   ├── 15_immature_t_clusters.png
│   ├── 17_developmental_stages_annotated.png
│   └── ... (18 figures total)

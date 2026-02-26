library(targets)
library(tarchetypes)
source("Module_01_data_loaded.R")
source("Module_02_qc_filtered.R")
source("Module_03_wnn_integrated.R")
source("Module_04_reference_mapping.R")
source("Module_05_Subsetting of T cell precursors.R")
source("Module_06_Clustering_Annotation.R")

list(
  tar_target(raw_data, load_data()),
  tar_target(qc_obj, qc_filter(raw_data)),
  tar_target(norm_obj, normalize_data(qc_obj)),
  tar_target(cluster_obj, cluster_cells(norm_obj)),
  tar_target(scored_obj, calculate_module_scores(cluster_obj)),
  tar_target(final_obj, classify_cells(scored_obj))
)

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")

cat("\n\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("✓ PIPELINE COMPLETED SUCCESSFULLY\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("Analysis completed at:", as.character(end_time), "\n")
cat("Total runtime:", round(runtime, 2), "minutes\n\n")

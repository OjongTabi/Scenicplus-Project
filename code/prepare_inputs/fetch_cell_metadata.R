#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(Signac)

# Path to input and output files
input_rds <- "/n/sci/SCI-004375-NYUDATA/Kenzie/Multiome-Merged.rds"
output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/cell_metadata/"

multiome_merged <- readRDS(input_rds)

for (i in 1:42) {
  
  metacluster <- paste0("M", i)
  
  M <- subset(multiome_merged, metacluster.idents == metacluster)
  
  # Save the subsetted Seurat object
  output_path <- paste0(output_dir, metacluster, "_cell_data.tsv")
  
  write.table(M@meta.data, file = output_path, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  
  cat("Completed processing for", metacluster, "\n")
}

cat("All metaclusters processed successfully.\n")


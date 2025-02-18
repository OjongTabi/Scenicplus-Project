#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(Signac)

# Path to input and output files
input_rds <- "/n/sci/SCI-004375-NYUDATA/Kenzie/Multiome-Merged.rds"
input_fragments <- "/n/sci/SCI-004375-NYUDATA/Kenzie/atac_fragments.tsv.gz"
output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/Fragments/"

multiome_merged <- readRDS(input_rds)

for (i in 1:42) {
  
  metacluster <- paste0("M", i)
  
  M <- subset(multiome_merged, metacluster.idents == metacluster)
  
  # Set the active assay to 'ATAC'
  DefaultAssay(M) <- "ATAC"
  
  # Subset fragments file
  output_file <- paste0(output_dir, "atac_fragments_", metacluster, ".tsv.gz")
  FilterCells(
    fragments = input_fragments,
    cells = Cells(M),
    outfile = output_file
  )
  
  cat("Completed processing for", metacluster, "\n")
}

cat("All metaclusters processed successfully.\n")









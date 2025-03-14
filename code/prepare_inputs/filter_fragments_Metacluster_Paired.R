#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(Signac)

# Path to input and output files
input_rds <- "/n/sci/SCI-004375-NYUDATA/Data/Multiome/Multiome-Merged.rds"
input_fragments <- "/n/sci/SCI-004375-NYUDATA/Kenzie/atac_fragments.tsv.gz"
output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/Fragments/"

multiome_merged <- readRDS(input_rds)

selected_metaclusters <- list(
  c("M28", "M13"),
  c("M35", "M27"),
  c("M25", "M29", "M30"),
  c("M34", "M36"),
  c("M19", "M20"),
  c("M15", "M17"),
  c("M31", "M32", "M33"),
  c("M7", "M8", "M22"),
  c("M23", "M24", "M26")
)


for (group in selected_metaclusters) {
  
  metacluster_name <- paste(group, collapse = "_")  # Create a combined name for the output file
  
  
  M <- subset(multiome_merged, metacluster.idents %in% group)
  
  # Set the active assay to 'ATAC'
  DefaultAssay(M) <- "ATAC"
  
  # Subset fragments file
  output_file <- paste0(output_dir, "atac_fragments_", metacluster_name, ".tsv.gz")
  FilterCells(
    fragments = input_fragments,
    cells = Cells(M),
    outfile = output_file
  )
  
  cat("Completed processing for", metacluster_name, "\n")
}

cat("All metaclusters processed successfully.\n")









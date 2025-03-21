#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(Signac)

# Path to input and output files
input_rds <- ("/n/sci/SCI-004375-NYUDATA/Data/Multiome/P0/Dm/P0-Dm-M24.rds")
input_fragments <- "/n/sci/SCI-004375-NYUDATA/Kenzie/atac_fragments.tsv.gz"
output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/Fragments/"

M <- readRDS(input_rds)

  
DefaultAssay(M) <- "ATAC"
  

output_file <- paste0(output_dir, "atac_fragments_", "P0_Dm_M24", ".tsv.gz")

FilterCells(
  fragments = input_fragments,
  cells = Cells(M),
  outfile = output_file
)









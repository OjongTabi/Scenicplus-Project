#!/usr/bin/env Rscript

rm(list = ls())


library(Seurat)
library(Signac)
library(SeuratDisk)

# Path to input and output files
input_rds <- "/n/sci/SCI-004375-NYUDATA/Kenzie/Multiome-Merged.rds"
output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/GEX_h5ad/"

multiome_merged <- readRDS(input_rds)

for (i in 1:42) {
  
  metacluster <- paste0("M", i)
  
  M <- subset(multiome_merged, metacluster.idents == metacluster)
  
  
  # New Seurat object with only the RNA assay
  rna_seurat_object <- CreateSeuratObject(counts = M[["RNA"]]@counts)
  
  rna_seurat_object@meta.data <- M@meta.data[rownames(rna_seurat_object@meta.data), ]
  
  rna_seurat_object$AnnotatedJul24 <- as.character(rna_seurat_object$AnnotatedJul24)
  
  rna_seurat_object <- FindVariableFeatures(rna_seurat_object)# Needed for conversion in next step

  # Convert "Assay5" to "Assay" for SeuratDisk compatibility
  rna_seurat_object[["RNA"]] <- as(object = rna_seurat_object[["RNA"]], Class = "Assay")
  
  # Save the subsetted Seurat object
  output_path <- paste0(output_dir, metacluster, "_RNA.h5Seurat")
  
  SaveH5Seurat(rna_seurat_object, filename = output_path, overwrite = TRUE)
  
  # Convert the H5Seurat file to H5AD
  Convert(output_path, dest = "h5ad", overwrite = TRUE)
  
  cat("Completed processing for", metacluster, "\n")
}

cat("All metaclusters processed successfully.\n")


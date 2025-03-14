#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(Signac)
library(SeuratDisk)

# Path to input and output files
input_rds <- "/n/sci/SCI-004375-NYUDATA/Data/Multiome/Multiome-Merged.rds"
output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/GEX_h5ad/"


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
  
  # New Seurat object with only the RNA assay
  rna_seurat_object <- CreateSeuratObject(counts = M[["RNA"]]@counts)
  
 
  rna_seurat_object@meta.data <- M@meta.data[rownames(rna_seurat_object@meta.data), ]
  
  
  rna_seurat_object$AnnotatedJul24 <- as.character(rna_seurat_object$AnnotatedJul24)
  
  
  rna_seurat_object <- FindVariableFeatures(rna_seurat_object)# Needed for conversion in next step
  
  # Convert "Assay5" to "Assay" for SeuratDisk compatibility
  rna_seurat_object[["RNA"]] <- as(object = rna_seurat_object[["RNA"]], Class = "Assay")
  
  
  output_path <- paste0(output_dir, metacluster_name, "_RNA.h5Seurat")
  
  # Save the subsetted Seurat object
  SaveH5Seurat(rna_seurat_object, filename = output_path, overwrite = TRUE)
  
  # Convert H5Seurat to H5AD format
  Convert(output_path, dest = "h5ad", overwrite = TRUE)
  
  #Save the metadata
  write.table(rna_seurat_object@meta.data, file = paste0("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/cell_metadata/", metacluster_name, "_cell_data.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

  cat("Completed processing for", metacluster_name, "\n")
}

cat("All selected metaclusters processed successfully.\n")

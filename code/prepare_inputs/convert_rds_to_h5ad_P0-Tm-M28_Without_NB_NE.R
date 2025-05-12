#!/usr/bin/env Rscript

rm(list = ls())


library(Seurat)
library(Signac)
library(SeuratDisk)


input_rds <- ("/n/sci/SCI-004375-NYUDATA/Data/Multiome/P0/Tm/P0-Tm-M28.rds")

output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/GEX_h5ad/"

M <- readRDS(input_rds)
M <- subset(M, idents = setdiff(levels(M), c("NE-OPC", "NB-OPC")))
M$Idents <- as.character(Idents(M))

write.table(M@meta.data, file = paste0("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/cell_metadata/", "P0-Tm-M28_Without_NB_NE", "_cell_data.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


rna_seurat_object <- CreateSeuratObject(counts = M[["RNA"]]@counts)

rna_seurat_object@meta.data <- M@meta.data[rownames(rna_seurat_object@meta.data), ]

rna_seurat_object$AnnotatedJul24 <- as.character(rna_seurat_object$AnnotatedJul24)

rna_seurat_object <- FindVariableFeatures(rna_seurat_object)# Needed for conversion in next step

# Convert "Assay5" to "Assay" for SeuratDisk compatibility
rna_seurat_object[["RNA"]] <- as(object = rna_seurat_object[["RNA"]], Class = "Assay")

# Save the subsetted Seurat object
output_path <- paste0(output_dir, "P0-Tm-M28_Without_NB_NE", "_RNA.h5Seurat")

SaveH5Seurat(rna_seurat_object, filename = output_path, overwrite = TRUE)

# Convert the H5Seurat file to H5AD
Convert(output_path, dest = "h5ad", overwrite = TRUE)



  
  
  




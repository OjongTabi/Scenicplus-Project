#!/usr/bin/env Rscript

rm(list = ls())


library(Seurat)
library(Signac)
library(SeuratDisk)


input_rds <- ("/n/sci/SCI-004375-NYUDATA/Data/Multiome/P0/Dm/P0-Dm-M24.rds")

output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/GEX_h5ad/Comprehensive/P0/"

M <- readRDS(input_rds)
M <- subset(M, idents = setdiff(levels(M), c("NE-OPC", "NB-OPC")))
M$Idents <- as.character(Idents(M))

P15_Dm = readRDS("/n/sci/SCI-004375-NYUDATA/Data/Ozel21/P15_Dm.rds")
P15_Dm$AnnotatedJul24 <- as.character(Idents(P15_Dm))
P15_Dm$Idents <- as.character(Idents(P15_Dm))

# Merge Multiome data with the external RNA-seq data
combined_data <- merge(M, y = list(P15_Dm), 
                         add.cell.ids = c("Multiome", "Ozel21"), project = "Comprehensive_RNA")

rna_seurat_object <- CreateSeuratObject(counts = combined_data[["RNA"]]@counts)


rna_seurat_object@meta.data <- combined_data@meta.data[rownames(rna_seurat_object@meta.data), ]

rna_seurat_object$AnnotatedJul24 <- as.character(rna_seurat_object$AnnotatedJul24)

rna_seurat_object <- FindVariableFeatures(rna_seurat_object)# Needed for conversion in next step

# Convert "Assay5" to "Assay" for SeuratDisk compatibility
rna_seurat_object[["RNA"]] <- as(object = rna_seurat_object[["RNA"]], Class = "Assay")

# Save the subsetted Seurat object
output_path <- paste0(output_dir, "P0_Dm_M24_Without_NB_NE", "_RNA.h5Seurat")

SaveH5Seurat(rna_seurat_object, filename = output_path, overwrite = TRUE)

# Convert the H5Seurat file to H5AD
Convert(output_path, dest = "h5ad", overwrite = TRUE)
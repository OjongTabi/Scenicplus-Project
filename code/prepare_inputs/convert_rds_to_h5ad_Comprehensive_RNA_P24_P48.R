#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(Signac)
library(SeuratDisk)
library(readxl)
library(dplyr)

# Path to input and output files
input_rds <- "/n/sci/SCI-004375-NYUDATA/Kenzie/Multiome-Merged.rds"
output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/GEX_h5ad/Comprehensive/P24_P48/"
file_path_cell_type <- "/n/sci/SCI-004375-NYUDATA/Data/Multiome/metacluster-identities-final-5.23.24.xlsx" # External RNA cell type annotation

multiome_merged <- readRDS(input_rds)

data <- read_excel(file_path_cell_type) %>% 
  slice(-38) # Remove the second header in file

for (i in setdiff(seq_len(nrow(data)), 13)) { ## Skip M13 as there are merging issues
  
  metacluster <- data$Metacluster[i]
  # Extract cell IDs
  cell_ids <- as.numeric(unlist(strsplit(data$`Cluster Identities (Clusters)`[i], ",\\s*")))
  
  multiome_data <- subset(multiome_merged, subset = metacluster.idents == metacluster)
  multiome_data <- subset(multiome_data, subset = stage %in% c("P24","P48"))
  DefaultAssay(multiome_data) <- "RNA"
  
  # Load Ozel21 objects
  rds_directory <- "/n/sci/SCI-004375-NYUDATA/Data/Ozel21"
  
  specific_files <- c("P30.rds", "P40.rds", "P50.rds")
  
  for (file in specific_files) {
    
    full_path <- file.path(rds_directory, file)
    
    base_name <- gsub("\\.rds$", "", file)
    base_name <- paste0("Ozel21_", base_name)
    
    seurat_obj <- readRDS(full_path)
    
    # Subset Ozel21 objects
    seurat_obj <- subset(seurat_obj, subset = MultiomeNN %in% cell_ids)
     
    assign(base_name, seurat_obj)
    
  }
  
  # Load Zipursky objects
  rds_directory <- "/n/sci/SCI-004375-NYUDATA/Data/Zipursky"
  
  specific_files <- c(
    "Zipursky_P24.rds",
    "Zipursky_P36.rds",
    "Zipursky_P48.rds")
  
  for (file in specific_files) {
    
    full_path <- file.path(rds_directory, file)
    
    base_name <- gsub("\\.rds$", "", file)
    
    seurat_obj <- readRDS(full_path)
    
    # Subset Zipursky objects
    seurat_obj <- subset(seurat_obj, subset = MultiomeNN %in% cell_ids)
    
    assign(base_name, seurat_obj)
    
  }
  
  ozel_data <- merge(Ozel21_P30, y = list( Ozel21_P40, Ozel21_P50), 
                   add.cell.ids = c("P30","P40", "P50"), project = "Ozel21")
  
  ozel_data$AnnotatedJul24 <- ozel_data$MultiomeAnnotated # Consistent cell type annotation column name
  
  zipursky_data <- merge(Zipursky_P24, y = list(Zipursky_P36, Zipursky_P48), 
                       add.cell.ids = c("P24", "P36", "P48"), 
                       project = "Zipursky")
  
  zipursky_data$AnnotatedJul24 <- zipursky_data$MultiomeAnnotated # Consistent cell type annotation column name
  
  # Merge Multiome data with the external RNA-seq data
  combined_data <- merge(multiome_data, y = list(ozel_data, zipursky_data), 
                         add.cell.ids = c("Multiome", "Ozel21", "Zipursky"), project = "Comprehensive_RNA")
  
  # Keep only the RNA assay
  rna_assay <- combined_data[["RNA"]]
  
  # Create a new Seurat object with only the RNA assay
  rna_seurat_object <- CreateSeuratObject(counts = rna_assay@counts)
  rna_seurat_object@meta.data <- combined_data@meta.data[rownames(rna_seurat_object@meta.data), ]
  
  rna_seurat_object$AnnotatedJul24 <- as.character(rna_seurat_object$AnnotatedJul24)
  
  rna_seurat_object <- FindVariableFeatures(rna_seurat_object) # Needed for conversion in next step
  
  # Convert "Assay5" to "Assay" for SeuratDisk compatibility
  rna_seurat_object[["RNA"]] <- as(object = rna_seurat_object[["RNA"]], Class = "Assay")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  output_path <- paste0(output_dir, metacluster, "_Comprehensive_RNA.h5Seurat")
  
  SaveH5Seurat(rna_seurat_object, filename = output_path, overwrite = TRUE)
  
  # Convert the H5Seurat file to H5AD
  Convert(output_path, dest = "h5ad", overwrite = TRUE)
  
}

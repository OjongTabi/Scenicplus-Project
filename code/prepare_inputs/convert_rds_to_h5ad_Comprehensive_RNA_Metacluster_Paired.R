#!/usr/bin/env Rscript

rm(list = ls())


library(Seurat)
library(Signac)
library(SeuratDisk)
library(readxl)
library(dplyr)

# Path to input and output files
input_rds <- "/n/sci/SCI-004375-NYUDATA/Data/Multiome/Multiome-Merged.rds"
output_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/GEX_h5ad/Comprehensive/P24_P48_Adult/"

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
  
  metacluster <- paste(group, collapse = "_")  # Create a combined name for the output file

  
  
  multiome_data <- subset(multiome_merged, metacluster.idents %in% group)
  DefaultAssay(multiome_data) <- "RNA"
  
  # Extract cell cluster IDs
  cell_ids <- unique(multiome_data$Clusters)
  
  # Load Ozel21 objects
  rds_directory <- "/n/sci/SCI-004375-NYUDATA/Data/Ozel21"
  
  specific_files <- c("Adult.rds", "P30.rds", "P40.rds", "P50.rds", "P70.rds")
  
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
    "Zipursky_P48.rds",
    "Zipursky_P60.rds",
    "Zipursky_P72.rds",
    "Zipursky_P84.rds",
    "Zipursky_P96.rds"
  )
  
  for (file in specific_files) {
    
    full_path <- file.path(rds_directory, file)
    
    base_name <- gsub("\\.rds$", "", file)
    
    seurat_obj <- readRDS(full_path)
    
    # Subset Zipursky objects
    seurat_obj <- subset(seurat_obj, subset = MultiomeNN %in% cell_ids)
    
    assign(base_name, seurat_obj)
    
  }
  
  ozel_data <- merge(Ozel21_Adult, y = list(Ozel21_P30, Ozel21_P40, Ozel21_P50, Ozel21_P70), 
                     add.cell.ids = c("Adult", "P30", "P40", "P50", "P70"), project = "Ozel21")
  
  ozel_data$AnnotatedJul24 <- ozel_data$MultiomeAnnotated # Consistent cell type annotation column name
  
  zipursky_data <- merge(Zipursky_P24, y = list(Zipursky_P36, Zipursky_P48, Zipursky_P60, 
                                                Zipursky_P72, Zipursky_P84, Zipursky_P96), 
                         add.cell.ids = c("P24", "P36", "P48", "P60", "P72", "P84", "P96"), 
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
  
  output_path <- paste0(output_dir,metacluster, "_Comprehensive_RNA.h5Seurat")
  
  SaveH5Seurat(rna_seurat_object, filename = output_path, overwrite = TRUE)
  
  # Convert the H5Seurat file to H5AD
  Convert(output_path, dest = "h5ad", overwrite = TRUE)
  
}

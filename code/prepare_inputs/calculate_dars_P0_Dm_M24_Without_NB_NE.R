#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(Signac)
library(GenomicRanges)
library(rtracklayer)


input_file <- ("/n/sci/SCI-004375-NYUDATA/Data/Multiome/P0/Dm/P0-Dm-M24.rds")
output_dir <-("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/dars/P24_P48_Adult/P0_Dm_M24_Without_NB_NE")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

atac_data <- readRDS(input_file)

DefaultAssay(atac_data) <- "ATAC"
atac_data <- subset(atac_data, idents = setdiff(levels(atac_data), c("NE-OPC", "NB-OPC")))


# Find Differentially Accessible Regions (DARs)
DARs <- FindAllMarkers(
  atac_data,
  assay = "ATAC",
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments',
  logfc.threshold = 0.5, 
  min.pct = 0.05
)


DARs_filtered <- DARs[DARs$p_val_adj < 0.05, ]

# Export all markers as BED
all.markers <- unique(DARs_filtered$gene)
consensus <- StringToGRanges(all.markers)
consensus <- sort(consensus)
export.bed(consensus, con = paste0(output_dir, "P24_P48_Adult_peaks.bed")) #Name consistency was maintained, and SCENIC+ run scripts were written to use the 'P24_P48_Adult_peaks.bed' file name  

# Export BED files for each unique cluster
unique_clusters <- unique(DARs_filtered$cluster)
cluster_dir <- paste0(output_dir, "/DARs_cell_type/")
dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

for (cluster in unique_clusters) {
  cluster_data <- DARs_filtered[DARs_filtered$cluster == cluster, ]
  peaks <- StringToGRanges(cluster_data$gene)
  peaks <- sort(peaks)
  export.bed(peaks, con = paste0(cluster_dir, cluster, ".bed"))}
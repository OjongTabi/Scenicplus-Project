#!/usr/bin/env Rscript

rm(list = ls())

library(Seurat)
library(Signac)
library(GenomicRanges)
library(rtracklayer)


metaclusters <- paste0("M", 1:42)

for (metacluster in metaclusters) {
  
  
  input_file <- paste0("/n/sci/SCI-004375-NYUDATA/Kenzie/PeaksCounts/Recounts/", metacluster, "-ATACrecounts.rds")
  output_dir <- paste0("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/dars/P24_P48_Adult/", metacluster, "/")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  atac_data <- readRDS(input_file)
  
  # Subset to retain only stages P24,P48,Adult
  atac_data <- subset(atac_data, stage %in% c("P24", "P48","Adult"))
  atac_data <- subset(atac_data, metacluster.idents == metacluster)
  
  
  # Combine the 'stage' and 'AnnotatedJul24' columns
  atac_data$stage_AnnotatedJul24 <- paste(atac_data$stage, atac_data$AnnotatedJul24, sep = "_")
  
  
  
  DefaultAssay(atac_data) <- "ATAC"
  Idents(atac_data) <- "stage_AnnotatedJul24"
  
  
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
  export.bed(consensus, con = paste0(output_dir, "P24_P48_Adult_peaks.bed"))
  
  # Export BED files for each unique cluster
  unique_clusters <- unique(DARs_filtered$cluster)
  cluster_dir <- paste0(output_dir, "/DARs_cell_type/")
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (cluster in unique_clusters) {
    cluster_data <- DARs_filtered[DARs_filtered$cluster == cluster, ]
    peaks <- StringToGRanges(cluster_data$gene)
    peaks <- sort(peaks)
    export.bed(peaks, con = paste0(cluster_dir, cluster, ".bed"))
  }
  
  cat("Finished processing", metacluster, "\n")
}

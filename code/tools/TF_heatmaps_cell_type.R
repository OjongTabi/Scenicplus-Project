library(ComplexHeatmap)
library(circlize)

input_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M28/P24_P48_Adult/scplus_pipeline/Snakemake"

files <- list.files(input_dir, pattern = "^RNA_Counts_tfs_.*\\.tsv$", full.names = TRUE)

for (file in files) {
  
  category <- gsub("RNA_Counts_tfs_|\\.tsv", "", basename(file))
  
  df <- read.delim(file, row.names = 1, check.names = FALSE)
  
  
  # Z-score normalization per gene
  df_scaled <- scale(as.matrix(df))
  df_scaled[is.na(df_scaled)] <- 0
  
  hm <- Heatmap(
    df_scaled,
    name = category,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = FALSE,
    column_title = category,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  )
  
  # Save to PNG
  png_filename <- file.path(input_dir, paste0("Heatmap_", category, ".png"))
  png(png_filename, width = 1000, height = 800)
  draw(hm, heatmap_legend_side = "right")
  dev.off()
  
}

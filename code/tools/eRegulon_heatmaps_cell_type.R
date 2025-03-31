library(ComplexHeatmap)
library(circlize)

input_dir <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results/M28/P24_P48_Adult/scplus_pipeline/Snakemake"

files <- list.files(input_dir, pattern = "^eRegulon_gene_auc_.*\\.tsv$", full.names = TRUE)

for (file in files) {
  
  category <- gsub("eRegulon_gene_auc_|\\.tsv", "", basename(file))
  
  df <- read.delim(file, row.names = 1, check.names = FALSE)
  
  
  
  
  hm <- Heatmap(
    as.matrix(df),
    name = category,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = FALSE,
    column_title = category,
    col = colorRamp2(c(0, 1), c( "white", "red"))
  )
  
  # Save to PNG
  png_filename <- file.path(input_dir, paste0("Heatmap_eRegulon_", category, ".png"))
  png(png_filename, width = 1200, height = 900)
  draw(hm, heatmap_legend_side = "right")
  dev.off()
  
}

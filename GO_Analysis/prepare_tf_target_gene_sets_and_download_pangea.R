rm(list = ls())
library(tidyverse)
library(readr)
library(dplyr)


base_path <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results"


file_paths <- paste0(base_path, "/M", 1:42, "/P24_P48_Adult_Comprehensive/eRegulon_direct.tsv")


unique_tf_gene_df <- file_paths %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE) %>% dplyr::select(TF, Gene)) %>% # Read and combine the TF and Gene columns from all files
  distinct() %>%
  arrange(TF)



ConversionTable <- read_rds("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/ConversionTable-new.rds")

# Create named vector for lookup
conversion_vector <- setNames(ConversionTable$converted_id, ConversionTable$reference_symbol)

# Map TF and Gene to their converted names
unique_tf_gene_df$TF_Converted <- conversion_vector[unique_tf_gene_df$TF]
unique_tf_gene_df$Gene_Converted <- conversion_vector[unique_tf_gene_df$Gene]


background_df <- data.frame(Gene_Converted = unique(unique_tf_gene_df$Gene_Converted))
background_genes <- background_df$Gene_Converted

sheet_list <- list()
sheet_list[["background"]] <- background_df


unique_tfs <- unique(unique_tf_gene_df$TF)

for (tf in unique_tfs) {
  tf_genes <- unique_tf_gene_df %>%
    filter(TF == tf) %>%
    dplyr::select(Gene_Converted) %>%
    distinct()
  
  if (nrow(tf_genes) >= 50) {  #Minimum 50 targets
    sheet_list[[tf]] <- tf_genes
  }
}




output_file <- file("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/TF_Targets.txt", open = "wt")

for (i in 2:length(sheet_list)) {
  
  group_id <- names(sheet_list)[i]
  
  genes <- sheet_list[[i]]
  
  for (gene in genes) {
    writeLines(paste0(group_id, ", ", gene), output_file)
  }
}


close(output_file)




library(rvest)
library(httr)


main_url <- "https://www.flyrnai.org/tools/pangea/web/multiple_enrichment_display/multiple_20250428_122409s6"

page <- read_html(main_url)
str(page)

csv_links <- page %>%
  html_elements("a") %>%     #select all <a> tags
  html_attr("href") %>%      #get href attributes
  (\(x) x[startsWith(x, "/tools/pangea/web/show_results_from_file/")])() 


full_csv_links <- paste0("https://www.flyrnai.org", csv_links)


for (i in seq_along(full_csv_links)) {
  file_url <- full_csv_links[i]
  
  
  filename <- paste0("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/","PANGEA_TF_", names(sheet_list)[-1][i], ".csv")
  destfile <- file.path(filename)
  
  
  #download.file(file_url, destfile, mode = "wb")
  
  #cat("Downloaded:", filename, "\n")
}


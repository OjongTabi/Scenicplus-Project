rm(list = ls())

library(tidyverse)
library(writexl)
library(dplyr)

base_path <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results"


file_paths <- paste0(base_path, "/M", 1:42, "/P24_P48_Adult_Comprehensive/eRegulon_direct.tsv")


unique_tf_gene_df <- file_paths %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE) %>% select(TF, Gene)) %>% # Read and combine the TF and Gene columns from all files
  distinct() %>%
  arrange(TF)



ConversionTable <- read_rds("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/data/ConversionTable-new.rds")

# Create named vector for lookup
conversion_vector <- setNames(ConversionTable$converted_id, ConversionTable$reference_symbol)

# Map TF and Gene to their converted names
unique_tf_gene_df$TF_Converted <- conversion_vector[unique_tf_gene_df$TF]
unique_tf_gene_df$Gene_Converted <- conversion_vector[unique_tf_gene_df$Gene]


background_df <- data.frame(Gene_Converted = unique(unique_tf_gene_df$Gene_Converted))


sheet_list <- list()
sheet_list[["background"]] <- background_df


unique_tfs <- unique(unique_tf_gene_df$TF)

for (tf in unique_tfs) {
  tf_genes <- unique_tf_gene_df %>%
    filter(TF == tf) %>%
    select(Gene_Converted) %>%
    distinct()
  sheet_list[[tf]] <- tf_genes
}


write_xlsx(sheet_list, path = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/code/tools/TF_GeneConverted_List.xlsx")




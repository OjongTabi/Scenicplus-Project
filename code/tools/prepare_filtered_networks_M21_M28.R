library(dplyr)

flybase_tfs <- readRDS("/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/FlyBaseTFsConverted.rds")
comprehensive_df <- read.delim("/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/M21_P24_P48_Adult_Comprehensive.tsv", 
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)


comprehensive_df <- comprehensive_df %>%
  filter(Gene %in% flybase_tfs, triplet_rank <= 2000) %>%
  group_by(TF, Gene) %>% # Get top interaction for same interaction
  slice_min(triplet_rank, with_ties = FALSE) %>%
  ungroup()


comprehensive_df$importance_x_rho_TF2G <- comprehensive_df$importance_TF2G * comprehensive_df$rho_TF2G

write.table(comprehensive_df,
            file = "/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/Filtered_Comprehensive_TFs_M21.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


comprehensive_df <- read.delim("/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/M28_P24_P48_Adult_Comprehensive.tsv", 
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)


comprehensive_df <- comprehensive_df %>%
  filter(Gene %in% flybase_tfs, triplet_rank <= 2000) %>%
  group_by(TF, Gene) %>% # Get top interaction for same interaction
  slice_min(triplet_rank, with_ties = FALSE) %>%
  ungroup()

comprehensive_df$importance_x_rho_TF2G <- comprehensive_df$importance_TF2G * comprehensive_df$rho_TF2G

write.table(comprehensive_df,
            file = "/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/Filtered_Comprehensive_TFs_M28.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
################
comprehensive_df <- read.delim("/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/M21_P24_P48_Adult_Comprehensive.tsv", 
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)


comprehensive_df <- comprehensive_df %>%
  filter(Gene %in% flybase_tfs, triplet_rank <= 1000) %>%
  group_by(TF, Gene) %>% # Get top interaction for same interaction
  slice_min(triplet_rank, with_ties = FALSE) %>%
  ungroup()


comprehensive_df$importance_x_rho_TF2G <- comprehensive_df$importance_TF2G * comprehensive_df$rho_TF2G

write.table(comprehensive_df,
            file = "/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/Filtered_Comprehensive_TFs_M21_1000.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


comprehensive_df <- read.delim("/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/M28_P24_P48_Adult_Comprehensive.tsv", 
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE)


comprehensive_df <- comprehensive_df %>%
  filter(Gene %in% flybase_tfs, triplet_rank <= 1000) %>%
  group_by(TF, Gene) %>% # Get top interaction for same interaction
  slice_min(triplet_rank, with_ties = FALSE) %>%
  ungroup()

comprehensive_df$importance_x_rho_TF2G <- comprehensive_df$importance_TF2G * comprehensive_df$rho_TF2G

write.table(comprehensive_df,
            file = "/n/sci/SCI-004375-NYUDATA/Ojong/Scenicplus_v1.0a2/Filtered_Comprehensive_TFs_M28_1000.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)





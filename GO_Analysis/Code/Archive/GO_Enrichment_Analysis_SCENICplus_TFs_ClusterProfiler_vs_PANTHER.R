rm(list = ls())
library(tidyverse)
library(patchwork)
library(readr)
library(writexl)
library(rbioapi)  # For PANTHER enrichment via API
library(clusterProfiler)
library(org.Dm.eg.db)  # Drosophila annotation package
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

# run enrichment using PANTHER via rbioapi
run_panther_enrichment <- function(genes) {
  genes <- na.omit(as.character(genes))
  
  ontology_ids <- list(
    BP = "GO:0008150",
    CC = "GO:0005575",
    MF = "GO:0003674"
  )
  
  results <- list()
  
  for (ont in names(ontology_ids)) {
    enriched <- rba_panther_enrich(
      genes = genes,
      organism = 7227,
      annot_dataset = ontology_ids[[ont]],
      cutoff = 0.05,
      ref_genes = background_genes,
      ref_organism = 7227
    )
    results[[ont]] <- enriched$result
  }
  
  return(results)
}

# Run enrichment for each TF and each ontology
go_results_list_PANTHER <- list()

for (tf in names(sheet_list)[-1]) {  # skip "background"
  genes <- sheet_list[[tf]]$Gene_Converted
  result_list <- run_panther_enrichment(genes)
  
  for (ont in names(result_list)) {
    go_results_list_PANTHER[[paste0(tf, "_", ont)]] <- result_list[[ont]]
  }
}



# Function to perform enrichment using clusterProfiler
run_go_enrichment <- function(genes, background) {
  
  genes <- na.omit(as.character(genes))
  background <- na.omit(as.character(background))
  
  results <- list(
    BP = enrichGO(
      gene          = genes,
      universe      = background,
      OrgDb         = org.Dm.eg.db,
      keyType       = "FLYBASE",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1,
      readable      = TRUE
    ),
    CC = enrichGO(
      gene          = genes,
      universe      = background,
      OrgDb         = org.Dm.eg.db,
      keyType       = "FLYBASE",
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1,
      readable      = TRUE
    ),
    MF = enrichGO(
      gene          = genes,
      universe      = background,
      OrgDb         = org.Dm.eg.db,
      keyType       = "FLYBASE",
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1,
      readable      = TRUE
    )
  )
  
  return(results)
}



go_results_list_clusterProfiler <- list()

for (tf in names(sheet_list)[-1]) {  # skip "background"
  genes <- sheet_list[[tf]]$Gene_Converted
  result_list <- run_go_enrichment(genes, background_genes)
  
  for (ont in names(result_list)) {
    go_results_list_clusterProfiler[[paste0(tf, "_", ont)]] <- as.data.frame(result_list[[ont]])
  }
}



write_xlsx(go_results_list_PANTHER, path = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Results_PANTHER.xlsx")


write_xlsx(go_results_list_clusterProfiler, path = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis//GO_Results_clusterProfiler.xlsx")


library(tidyverse)
library(patchwork)


prepare_go_plot_data <- function(tf, ontology, list_cp, list_panther) {
  
  key_cp <- paste0(tf, "_", ontology)
  key_panther <- paste0(tf, "_", ontology)
  
  
  df_cp <- list_cp[[key_cp]] %>%
    dplyr::slice(1:20) %>%
    mutate(Source = paste0("clusterProfiler_", ontology),
           Term = Description,
           log10_p = -log10(pvalue))
  
  
  df_panther <- list_panther[[key_panther]] %>%
    dplyr::slice(1:20) %>%
    mutate(Source = paste0("PANTHER_", ontology),
           Term = term.label,
           log10_p = -log10(pValue))
  
  return(list(cluster = df_cp, panther = df_panther))
}


plot_top20 <- function(df) {
  ggplot(df, aes(x = reorder(Term, log10_p), y = log10_p)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = unique(df$Source), x = NULL, y = "-log10(p-value)") +
    theme_minimal(base_size = 13)
}



for (tf in names(sheet_list)[-1]) {
  
  ontologies <- c("BP", "MF", "CC")
  
  plots <- list()
  
  for (ont in ontologies) {
    data_pair <- prepare_go_plot_data(tf, ont, go_results_list_clusterProfiler, go_results_list_PANTHER)
    plots[[paste0("cp_", ont)]] <- plot_top20(data_pair$cluster)
    plots[[paste0("panther_", ont)]] <- plot_top20(data_pair$panther)
  }
  
  final_plot <- (plots$cp_BP | plots$panther_BP) /
    (plots$cp_MF | plots$panther_MF) /
    (plots$cp_CC | plots$panther_CC)
  
  # Save plot to file
  ggsave(
    filename = paste0("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Enrichment_Top20_", tf, ".png"),
    plot = final_plot,
    width = 16,
    height = 12
  )
}



######################################### No Background

# run enrichment using PANTHER via rbioapi
run_panther_enrichment <- function(genes) {
  genes <- na.omit(as.character(genes))
  
  ontology_ids <- list(
    BP = "GO:0008150",
    CC = "GO:0005575",
    MF = "GO:0003674"
  )
  
  results <- list()
  
  for (ont in names(ontology_ids)) {
    enriched <- rba_panther_enrich(
      genes = genes,
      organism = 7227,
      annot_dataset = ontology_ids[[ont]],
      cutoff = 0.05,
      #ref_genes = background_genes,
      #ref_organism = 7227
    )
    results[[ont]] <- enriched$result
  }
  
  return(results)
}

# Run enrichment for each TF and each ontology
go_results_list_PANTHER <- list()

for (tf in names(sheet_list)[-1]) {  # skip "background"
  genes <- sheet_list[[tf]]$Gene_Converted
  result_list <- run_panther_enrichment(genes)
  
  for (ont in names(result_list)) {
    go_results_list_PANTHER[[paste0(tf, "_", ont)]] <- result_list[[ont]]
  }
}



# Function to perform enrichment using clusterProfiler
run_go_enrichment <- function(genes, background) {
  
  genes <- na.omit(as.character(genes))
  background <- na.omit(as.character(background))
  
  results <- list(
    BP = enrichGO(
      gene          = genes,
      #universe      = background,
      OrgDb         = org.Dm.eg.db,
      keyType       = "FLYBASE",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1,
      readable      = TRUE
    ),
    CC = enrichGO(
      gene          = genes,
      #universe      = background,
      OrgDb         = org.Dm.eg.db,
      keyType       = "FLYBASE",
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1,
      readable      = TRUE
    ),
    MF = enrichGO(
      gene          = genes,
      #universe      = background,
      OrgDb         = org.Dm.eg.db,
      keyType       = "FLYBASE",
      ont           = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.1,
      readable      = TRUE
    )
  )
  
  return(results)
}



go_results_list_clusterProfiler <- list()

for (tf in names(sheet_list)[-1]) {  # skip "background"
  genes <- sheet_list[[tf]]$Gene_Converted
  result_list <- run_go_enrichment(genes, background_genes)
  
  for (ont in names(result_list)) {
    go_results_list_clusterProfiler[[paste0(tf, "_", ont)]] <- as.data.frame(result_list[[ont]])
  }
}




write_xlsx(go_results_list_PANTHER,
           path = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Results_PANTHER_No_Background.xlsx")

write_xlsx(go_results_list_clusterProfiler,
           path = "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Results_clusterProfiler_No_Background.xlsx")


library(tidyverse)
library(patchwork)


prepare_go_plot_data <- function(tf, ontology, list_cp, list_panther) {
  
  key_cp <- paste0(tf, "_", ontology)
  key_panther <- paste0(tf, "_", ontology)
  
  
  df_cp <- list_cp[[key_cp]] %>%
    dplyr::slice(1:20) %>%
    mutate(Source = paste0("clusterProfiler_", ontology),
           Term = Description,
           log10_p = -log10(pvalue))
  
  
  df_panther <- list_panther[[key_panther]] %>%
    dplyr::slice(1:20) %>%
    mutate(Source = paste0("PANTHER_", ontology),
           Term = term.label,
           log10_p = -log10(pValue))
  
  return(list(cluster = df_cp, panther = df_panther))
}


plot_top20 <- function(df) {
  ggplot(df, aes(x = reorder(Term, log10_p), y = log10_p)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = unique(df$Source), x = NULL, y = "-log10(p-value)") +
    theme_minimal(base_size = 13)
}



for (tf in names(sheet_list)[-1]) {
  
  ontologies <- c("BP", "MF", "CC")
  
  plots <- list()
  
  for (ont in ontologies) {
    data_pair <- prepare_go_plot_data(tf, ont, go_results_list_clusterProfiler, go_results_list_PANTHER)
    plots[[paste0("cp_", ont)]] <- plot_top20(data_pair$cluster)
    plots[[paste0("panther_", ont)]] <- plot_top20(data_pair$panther)
  }
  
  final_plot <- (plots$cp_BP | plots$panther_BP) /
    (plots$cp_MF | plots$panther_MF) /
    (plots$cp_CC | plots$panther_CC)
  
  # Save plot to file
  ggsave(
    filename = paste0("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Enrichment_Top20_No_Background_", tf, ".png"),
    plot = final_plot,
    width = 16,
    height = 12
  )
}


rm(list = ls())
library(tidyverse)
library(patchwork)
library(readr)
library(writexl)
library(rbioapi)  # For PANTHER enrichment via API
library(clusterProfiler)
library(org.Dm.eg.db)  # Drosophila annotation package
library(dplyr)

setwd("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/")
files <- c(
  "Terms2GO_BLIMP_1_No_Background_GO_MF.txt",
  "Terms2GO_BLIMP_1_With_Background_GO_BP.txt",
  "PANGEA_BLIMP_1_No_Background.csv",
  "Terms2GO_BLIMP_1_With_Background_GO_CC.txt",
  "PANGEA_BLIMP_1_With_Background.csv",
  "Terms2GO_BLIMP_1_With_Background_GO_MF.txt",
  "Terms2GO_BLIMP_1_No_Background_GO_BP.txt",
  "Terms2GO_BLIMP_1_No_Background_GO_CC.txt"
)

for (file in files) {
  name <- tools::file_path_sans_ext(basename(file))
  assign(name, if (grepl("\\.csv$", file)) read.csv(file) else read.delim(file))
}


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
  sheet_list[[tf]] <- tf_genes
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



tp_bp1 <- Terms2GO_BLIMP_1_With_Background_GO_BP %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "Terms2GO_BP", Term = Description, log10_p = -log10(pvalue))

tp_bp2 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO BP") %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANGEA_SLIM2 GO BP", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_bp3 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "EXP GO Biological Processes") %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANGEA_EXP GO Biological Processes", Term = Gene.Set.Name, log10_p = -log10(P.value))


tp_bp4 <- go_results_list_clusterProfiler[["Blimp-1_BP"]] %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "clusterProfiler_BP", Term = Description, log10_p = -log10(pvalue))


tp_bp5 <- go_results_list_PANTHER[["Blimp-1_BP"]] %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANTHER_BP", Term = term.label, log10_p = -log10(pValue))


tp_cc1 <- Terms2GO_BLIMP_1_With_Background_GO_CC %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "Terms2GO_CC", Term = Description, log10_p = -log10(pvalue))

tp_cc2 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO CC") %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANGEA_SLIM2 GO CC", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_cc3 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "EXP GO Cellular Component") %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANGEA_EXP GO Cellular Component", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_cc4 <- go_results_list_clusterProfiler[["Blimp-1_CC"]] %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "clusterProfiler_CC", Term = Description, log10_p = -log10(pvalue))


tp_cc5 <- go_results_list_PANTHER[["Blimp-1_CC"]] %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANTHER_CC", Term = term.label, log10_p = -log10(pValue))


tp_mf1 <- Terms2GO_BLIMP_1_With_Background_GO_MF %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "Terms2GO_MF", Term = Description, log10_p = -log10(pvalue))

tp_mf2 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO MF") %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANGEA_SLIM2 GO MF", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_mf3 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "EXP GO Molecular Function") %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANGEA_EXP GO MF", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_mf4 <- go_results_list_clusterProfiler[["Blimp-1_MF"]] %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "clusterProfiler_MF", Term = Description, log10_p = -log10(pvalue))


tp_mf5 <- go_results_list_PANTHER[["Blimp-1_MF"]] %>%
  dplyr::slice(1:20) %>%
  mutate(Source = "PANTHER_MF", Term = term.label, log10_p = -log10(pValue))



plot_top20 <- function(df) {
  ggplot(df, aes(x = reorder(Term, log10_p), y = log10_p)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = unique(df$Source), x = NULL, y = "-log10(p-value)") +
    theme_minimal(base_size = 13)
}

p1 <- plot_top20(tp_bp1)
p2 <- plot_top20(tp_bp2)
p3 <- plot_top20(tp_bp3)
p4 <- plot_top20(tp_bp4)
p5 <- plot_top20(tp_bp5)

p6 <- plot_top20(tp_cc1)
p7 <- plot_top20(tp_cc2)
p8 <- plot_top20(tp_cc3)
p9 <- plot_top20(tp_cc4)
p10 <- plot_top20(tp_cc5)

p11 <- plot_top20(tp_mf1)
p12 <- plot_top20(tp_mf2)
p13 <- plot_top20(tp_mf3)
p14 <- plot_top20(tp_mf4)
p15 <- plot_top20(tp_mf5)

final_plot_BP <- (p1 | p2 | p3) / (p4 | p5) 

final_plot_CC <- ( p6 | p7 | p8) / (p9 | p10)

final_plot_MF <- ( p11 | p12 | p13)/ (p14 | p15)

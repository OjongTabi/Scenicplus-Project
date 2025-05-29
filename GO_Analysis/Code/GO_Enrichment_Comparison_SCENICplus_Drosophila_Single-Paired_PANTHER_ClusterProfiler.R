rm(list = ls())
library(tidyverse)
library(patchwork)
library(readr)
library(writexl)
library(rbioapi)  # For PANTHER enrichment via API
library(clusterProfiler)
library(org.Dm.eg.db)  # Drosophila annotation package
library(dplyr)
library(GO.db)


base_path <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/results"


all_dirs <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)
all_dirs <- grep("^M\\d+", all_dirs, value = TRUE)


paired_dirs <- grep("_", all_dirs, value = TRUE)
single_dirs <- setdiff(all_dirs, paired_dirs)

paired_components <- unique(unlist(strsplit(paired_dirs, "_")))


final_dirs <- setdiff(all_dirs, intersect(single_dirs, paired_components))


file_paths <- file.path(base_path, final_dirs, "P24_P48_Adult_Comprehensive", "eRegulon_direct.tsv")


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





prepare_go_plot_data <- function(tf, ontology, list_cp, list_panther) {
  
  key_cp <- paste0(tf, "_", ontology)
  key_panther <- paste0(tf, "_", ontology)
  
  
  df_cp <- list_cp[[key_cp]] %>%
    filter(p.adjust < 0.05) %>%
    dplyr::slice(1:20) %>%
    mutate(Source = paste0("clusterProfiler_", ontology),
           Term = Description,
           FoldEnrichment = FoldEnrichment)
  
  
  df_panther <- list_panther[[key_panther]] %>%
    filter(fdr < 0.05) %>%
    dplyr::slice(1:20) %>%
    mutate(Source = paste0("PANTHER_", ontology),
           Term = term.label,
           FoldEnrichment = fold_enrichment)
  
  return(list(cluster = df_cp, panther = df_panther))
}


plot_top20 <- function(df) {
  ggplot(df, aes(x = reorder(Term, FoldEnrichment), y = FoldEnrichment)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = unique(df$Source), x = NULL, y = "FoldEnrichment") +
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
    filename = paste0("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Enrichment_Single_Paired_Metaclusters/GO_Enrichment_Top20_", tf, ".png"),
    plot = final_plot,
    width = 16,
    height = 12
  )
}

# Read terminal selector TFs
common_tf_markers <- read.csv(
  "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Enrichment_Single_Paired_Metaclusters/GO_Term_FoldEnrichment_By_TF/CommonTFmarkersMM_Table.csv",
  stringsAsFactors = FALSE
)

marker_tfs <- tolower(common_tf_markers$X) #case insensitive

# Merge all GO results, tag with TF and ontology
merged_df <- imap_dfr(go_results_list_clusterProfiler, function(df, name) {
  tf <- sub("_(BP|MF|CC)$", "", name)
  ont <- sub(".*_", "", name)
  df %>%
    mutate(TF = tf, Ontology = ont)
})

filtered_df <- merged_df %>%
  filter(p.adjust <= 0.05)

fe_df <- filtered_df %>%
  dplyr::select(Description, FoldEnrichment, TF, Ontology)

# Get top 20 GO terms per TF by Fold Enrichment
top20_per_tf <- fe_df %>%
  group_by(TF) %>%
  arrange(desc(FoldEnrichment)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Get unique set of top terms across all TFs
unique_top_terms <- top20_per_tf %>%
  distinct(Description)

# Filter full Fold Enrichment data to only those terms
final_df <- fe_df %>%
  filter(Description %in% unique_top_terms$Description) %>%
  mutate(
    is_terminal_selector = tolower(TF) %in% marker_tfs,
    TF_status = ifelse(is_terminal_selector, "Terminal Selector", "Other TF")
  )

unique_terms <- unique(final_df$Description)

walk(unique_terms, function(term) {
  term_df <- final_df %>%
    filter(Description == term)
  
  ontology_label <- unique(term_df$Ontology)
  ontology_label_text <- paste("GO Ontology:", paste(unique(ontology_label), collapse = ", "))
  
  p <- ggplot(term_df, aes(x = reorder(TF, FoldEnrichment), y = FoldEnrichment, fill = TF_status)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("GO Term:", term),
      subtitle = ontology_label_text,
      y = "Fold Enrichment",
      x = "TF",
      fill = "TF Type"
    ) +
    scale_fill_manual(values = c("Terminal Selector" = "#E41A1C", "Other TF" = "#377EB8")) +
    theme_minimal(base_size = 12) + 
    theme(
      axis.text.y = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Safe filename
  safe_term <- gsub("[^a-zA-Z0-9]", "_", term)
  
  ggsave(
    filename = paste0(
      "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Enrichment_Single_Paired_Metaclusters/GO_Term_FoldEnrichment_By_TF/",
      safe_term, "_FoldEnrichment_GO_Results_clusterProfiler-With_Background.png"
    ),
    plot = p,
    width = 10,
    height = max(6, 0.2 * nrow(term_df)),
    bg = "white"
  )
})


###################################################

library(GO.db)
library(AnnotationDbi)


go_ids <- keys(GO.db, keytype = "GOID")
go_df <- select(GO.db, keys = go_ids, columns = c("GOID", "TERM", "ONTOLOGY"), keytype = "GOID")

# Prepare empty vectors
top_anc_id <- character(length(go_ids))
top_anc_term <- character(length(go_ids))
distance <- integer(length(go_ids))

# Loop through each GO term
for (i in seq_along(go_ids)) {
  id <- go_ids[i]
  ont <- go_df$ONTOLOGY[i]
  
  # Get ancestors
  ancestors <- switch(ont,
                      "BP" = GOBPANCESTOR[[id]],
                      "MF" = GOMFANCESTOR[[id]],
                      "CC" = GOCCANCESTOR[[id]],
                      NULL)
  
  # Only assign if ancestors exist
  if (!is.null(ancestors) && length(ancestors) > 0) {
    top_anc_id[i] <- tail(ancestors, 2)[1]
    top_anc_term[i] <- Term(GOTERM[[top_anc_id[i]]])
    distance[i] <- length(ancestors)
  } else {
    top_anc_id[i] <- NA
    top_anc_term[i] <- NA
    distance[i] <- 0
  }
}



go_df$Top_Ancestor_ID <- top_anc_id
go_df$Top_Ancestor_Term <- top_anc_term

annotated_final_df <- final_df %>%
  right_join(go_df, by = c("GO_ID" = "GOID"))

final_df$GO_ID <- sub("\\.\\.\\..*$", "", rownames(final_df))

final_df_summary <- final_df %>%
  group_by(GO_ID, Description, Ontology) %>%
  summarise(
    TF_Count = n_distinct(TF),  # Count how many unique TFs contributed to this GO term
    .groups = "drop"
  )

annotated_final_df <- final_df_summary %>%
  left_join(go_df, by = c("GO_ID" = "GOID"))

save.image("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Enrichment_Single_Paired_Metaclusters/GO_Enrichment_Single_Paired_Metaclusters.RData")



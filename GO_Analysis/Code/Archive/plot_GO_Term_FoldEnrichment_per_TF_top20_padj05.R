library(tidyverse)
library(readxl)


col_types <- c(
  "text",    # ID
  "text",    # Description
  "text",    # GeneRatio
  "text",    # BgRatio
  "numeric", # pvalue
  "numeric", # p.adjust
  "numeric", # qvalue
  "numeric", # geneID
  "numeric", # Count
  "numeric", # FoldEnrichment
  "text",    # Category
  "numeric"  # GeneCount
)


xlsx_path <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Results_clusterProfiler.xlsx"
sheet_names <- excel_sheets(xlsx_path)

go_results_list <- setNames(
  lapply(sheet_names, function(sheet) {
    read_excel(xlsx_path, sheet = sheet, col_types = col_types)
  }),
  sheet_names
)

#Merge all GO results, tag with TF and ontology
merged_df <- imap_dfr(go_results_list, function(df, name) {
  tf <- sub("_(BP|MF|CC)$", "", name)
  ont <- sub(".*_", "", name)
  df %>%
    mutate(TF = tf, Ontology = ont)
})


filtered_df <- merged_df %>%
  filter(p.adjust <= 0.05)


fe_df <- filtered_df %>%
  select(Description, FoldEnrichment, TF, Ontology)

#Get top 20 GO terms per TF by Fold Enrichment
top20_per_tf <- fe_df %>%
  group_by(TF) %>%
  arrange(desc(FoldEnrichment)) %>%
  slice_head(n = 20) %>%
  ungroup()

#Get unique set of top terms across all TFs
unique_top_terms <- top20_per_tf %>%
  distinct(Description)

#Filter full Fold Enrichment data to only those terms
final_df <- fe_df %>%
  filter(Description %in% unique_top_terms$Description)


unique_terms <- unique(final_df$Description)

walk(unique_terms, function(term) {
  term_df <- final_df %>%
    filter(Description == term)
  
  p <- ggplot(term_df, aes(x = reorder(TF, FoldEnrichment), y = FoldEnrichment, fill = Ontology)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("GO Term:", term),
      y = "Fold Enrichment",
      x = "TF"
    ) +
    theme_minimal(base_size = 12)
  
  #Safe filename by replacing problematic characters
  safe_term <- gsub("[^a-zA-Z0-9]", "_", term)
  
  
  ggsave(
    filename = paste0("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Term_FoldEnrichment_By_TF/", safe_term, "_FoldEnrichment_GO_Results_clusterProfiler-With_Background.png"),
    plot = p,
    width = 10,
    height = 6,
    bg = "white",
  ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
})


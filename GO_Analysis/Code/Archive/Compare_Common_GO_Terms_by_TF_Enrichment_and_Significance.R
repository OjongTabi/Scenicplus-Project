library(tidyverse)
library(patchwork)
library(readxl)



col_types <- c(
  "text",    
  "text",    
  "text",    
  "text",    
  "numeric", 
  "numeric", 
  "numeric", 
  "numeric", 
  "numeric", 
  "numeric", 
  "text",    
  "numeric"  
)


xlsx_path <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Results_clusterProfiler.xlsx"
sheet_names <- excel_sheets(xlsx_path)

go_results_list_clusterProfiler <- setNames(
  lapply(sheet_names, function(sheet) {
    read_excel(xlsx_path, sheet = sheet, col_types = col_types)
  }),
  sheet_names
)



#Merge GO results across transcription factors, tagging each row with its TF and GO category
term_df <- imap_dfr(go_results_list_clusterProfiler, function(df, name) {
  tf <- sub("_(BP|MF|CC)$", "", name)
  ont <- sub(".*_", "", name)
  df %>%
    dplyr::select(Description, p.adjust) %>% 
    mutate(TF = tf, Ontology = ont)
})


term_tf_counts <- term_df %>%
  count(Description, name = "TF_Count")


top_terms <- term_tf_counts %>%
  arrange(desc(TF_Count)) %>%
  dplyr::slice(1:3) %>%
  pull(Description)


plot_list <- map(top_terms, function(term) {
  top_term_df <- term_df %>%
    filter(Description == term) %>%
    mutate(log10_padj = -log10(p.adjust))
  
  ggplot(top_term_df, aes(x = reorder(TF, log10_padj), y = log10_padj, fill = Ontology)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("GO term:", term),
      y = "-log10(padj-value)",
      x = "TF"
    ) +
    theme_minimal(base_size = 13)
})


final_plot <- wrap_plots(plotlist = plot_list, ncol = 3)


print(final_plot)

ggsave("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/Top3_Common_GO_Terms_across_TFs_GO_Results_clusterProfiler.png", plot = final_plot, width = 18, height = 22)


term_df_fe <- imap_dfr(go_results_list_clusterProfiler, function(df, name) {
  tf <- sub("_(BP|MF|CC)$", "", name)
  ont <- sub(".*_", "", name)
  df %>%
    dplyr::select(Description, FoldEnrichment) %>%
    mutate(TF = tf, Ontology = ont)
})


term_tf_counts_fe <- term_df_fe %>%
  count(Description, name = "TF_Count")

# 3. Get top N common GO terms
top_terms_fe <- term_tf_counts_fe %>%
  arrange(desc(TF_Count)) %>%
  dplyr::slice(1:3) %>%
  pull(Description)


plot_list_fe <- map(top_terms_fe, function(term) {
  top_term_df <- term_df_fe %>%
    filter(Description == term)
  
  ggplot(top_term_df, aes(x = reorder(TF, FoldEnrichment), y = FoldEnrichment, fill = Ontology)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("GO term:", term),
      y = "Fold Enrichment",
      x = "TF"
    ) +
    theme_minimal(base_size = 13)
})


final_plot_fe <- wrap_plots(plotlist = plot_list_fe, ncol = 3)

print(final_plot_fe)

ggsave("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/Top3_Common_GO_Terms_FoldEnrichment_GO_Results_clusterProfiler.png",
        plot = final_plot_fe, width = 18, height = 22)


##########


col_types <- c(
  "text",    
  "text",    
  "text",    
  "text",    
  "numeric", 
  "numeric", 
  "numeric", 
  "numeric", 
  "numeric", 
  "numeric", 
  "text",    
  "numeric"  
)


xlsx_path <- "/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/GO_Results_clusterProfiler_No_Background.xlsx"
sheet_names <- excel_sheets(xlsx_path)

go_results_list_clusterProfiler <- setNames(
  lapply(sheet_names, function(sheet) {
    read_excel(xlsx_path, sheet = sheet, col_types = col_types)
  }),
  sheet_names
)



#Merge GO results across transcription factors, tagging each row with its TF and GO category
term_df <- imap_dfr(go_results_list_clusterProfiler, function(df, name) {
  tf <- sub("_(BP|MF|CC)$", "", name)
  ont <- sub(".*_", "", name)
  df %>%
    dplyr::select(Description, p.adjust) %>% 
    mutate(TF = tf, Ontology = ont)
})


term_tf_counts <- term_df %>%
  count(Description, name = "TF_Count")


top_terms <- term_tf_counts %>%
  arrange(desc(TF_Count)) %>%
  dplyr::slice(1:3) %>%
  pull(Description)


plot_list <- map(top_terms, function(term) {
  top_term_df <- term_df %>%
    filter(Description == term) %>%
    mutate(log10_padj = -log10(p.adjust))
  
  ggplot(top_term_df, aes(x = reorder(TF, log10_padj), y = log10_padj, fill = Ontology)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("GO term:", term),
      y = "-log10(padj-value)",
      x = "TF"
    ) +
    theme_minimal(base_size = 13)
})


final_plot <- wrap_plots(plotlist = plot_list, ncol = 3)


print(final_plot)

ggsave("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/Top3_Common_GO_Terms_across_TFs_GO_Results_clusterProfiler_No_Background.png", plot = final_plot, width = 18, height = 22)


term_df_fe <- imap_dfr(go_results_list_clusterProfiler, function(df, name) {
  tf <- sub("_(BP|MF|CC)$", "", name)
  ont <- sub(".*_", "", name)
  df %>%
    dplyr::select(Description, FoldEnrichment) %>%
    mutate(TF = tf, Ontology = ont)
})


term_tf_counts_fe <- term_df_fe %>%
  count(Description, name = "TF_Count")

# 3. Get top N common GO terms
top_terms_fe <- term_tf_counts_fe %>%
  arrange(desc(TF_Count)) %>%
  dplyr::slice(1:3) %>%
  pull(Description)


plot_list_fe <- map(top_terms_fe, function(term) {
  top_term_df <- term_df_fe %>%
    filter(Description == term)
  
  ggplot(top_term_df, aes(x = reorder(TF, FoldEnrichment), y = FoldEnrichment, fill = Ontology)) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste("GO term:", term),
      y = "Fold Enrichment",
      x = "TF"
    ) +
    theme_minimal(base_size = 13)
})


final_plot_fe <- wrap_plots(plotlist = plot_list_fe, ncol = 3)

print(final_plot_fe)

ggsave("/n/sci/SCI-004375-NYUDATA/Ojong/scenicplus_project/GO_Analysis/Top3_Common_GO_Terms_FoldEnrichment_GO_Results_clusterProfiler_No_Background.png",
       plot = final_plot_fe, width = 18, height = 22)

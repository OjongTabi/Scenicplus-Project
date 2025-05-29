rm(list = ls())
library(tidyverse)
library(patchwork)


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


tp_bp1 <- Terms2GO_BLIMP_1_No_Background_GO_BP %>%
  slice(1:10) %>%
  mutate(Source = "Terms2GO_BP", Term = Description, log10_p = -log10(pvalue))

tp_bp2 <- PANGEA_BLIMP_1_No_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO BP") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_SLIM2 GO BP", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_bp3 <- PANGEA_BLIMP_1_No_Background %>%
  filter(Gene.Set.Category == "EXP GO Biological Processes") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_EXP GO Biological Processes", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_cc1 <- Terms2GO_BLIMP_1_No_Background_GO_CC %>%
  slice(1:10) %>%
  mutate(Source = "Terms2GO_CC", Term = Description, log10_p = -log10(pvalue))

tp_cc2 <- PANGEA_BLIMP_1_No_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO CC") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_SLIM2 GO CC", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_cc3 <- PANGEA_BLIMP_1_No_Background %>%
  filter(Gene.Set.Category == "EXP GO Cellular Component") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_EXP GO Cellular Component", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_mf1 <- Terms2GO_BLIMP_1_No_Background_GO_MF %>%
  slice(1:10) %>%
  mutate(Source = "Terms2GO_MF", Term = Description, log10_p = -log10(pvalue))

tp_mf2 <- PANGEA_BLIMP_1_No_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO MF") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_SLIM2 GO MF", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_mf3 <- PANGEA_BLIMP_1_No_Background %>%
  filter(Gene.Set.Category == "EXP GO Molecular Function") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_EXP GO MF", Term = Gene.Set.Name, log10_p = -log10(P.value))


plot_top10 <- function(df) {
  ggplot(df, aes(x = reorder(Term, log10_p), y = log10_p)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = unique(df$Source), x = NULL, y = "-log10(p-value)") +
    theme_minimal(base_size = 13)
}


p1 <- plot_top10(tp_bp1)
p2 <- plot_top10(tp_bp2)
p3 <- plot_top10(tp_bp3)

p4 <- plot_top10(tp_cc1)
p5 <- plot_top10(tp_cc2)
p6 <- plot_top10(tp_cc3)

p7 <- plot_top10(tp_mf1)
p8 <- plot_top10(tp_mf2)
p9 <- plot_top10(tp_mf3)

final_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9)


final_plot


########

tp_bp1 <- Terms2GO_BLIMP_1_With_Background_GO_BP %>%
  slice(1:10) %>%
  mutate(Source = "Terms2GO_BP", Term = Description, log10_p = -log10(pvalue))

tp_bp2 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO BP") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_SLIM2 GO BP", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_bp3 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "EXP GO Biological Processes") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_EXP GO Biological Processes", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_cc1 <- Terms2GO_BLIMP_1_With_Background_GO_CC %>%
  slice(1:10) %>%
  mutate(Source = "Terms2GO_CC", Term = Description, log10_p = -log10(pvalue))

tp_cc2 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO CC") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_SLIM2 GO CC", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_cc3 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "EXP GO Cellular Component") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_EXP GO Cellular Component", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_mf1 <- Terms2GO_BLIMP_1_With_Background_GO_MF %>%
  slice(1:10) %>%
  mutate(Source = "Terms2GO_MF", Term = Description, log10_p = -log10(pvalue))

tp_mf2 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "SLIM2 GO MF") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_SLIM2 GO MF", Term = Gene.Set.Name, log10_p = -log10(P.value))

tp_mf3 <- PANGEA_BLIMP_1_With_Background %>%
  filter(Gene.Set.Category == "EXP GO Molecular Function") %>%
  slice(1:10) %>%
  mutate(Source = "PANGEA_EXP GO MF", Term = Gene.Set.Name, log10_p = -log10(P.value))


p1 <- plot_top10(tp_bp1)
p2 <- plot_top10(tp_bp2)
p3 <- plot_top10(tp_bp3)

p4 <- plot_top10(tp_cc1)
p5 <- plot_top10(tp_cc2)
p6 <- plot_top10(tp_cc3)

p7 <- plot_top10(tp_mf1)
p8 <- plot_top10(tp_mf2)
p9 <- plot_top10(tp_mf3)

final_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9)


final_plot

###############```r
###############
# --- Setup ---
###############
source("src/00_init.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")

library(tidyverse)
library(enrichR)
library(purrr)
library(fgsea)
library(latex2exp)
library(readr)

base    <- "Ag_ScRNA_12_enrichR_nonsigniicant/"
basedir <- dirout(base)

###############
# --- Load & preprocess limma ---
###############
limmaRes <- readRDS(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))

limmaRes <- limmaRes %>%
  mutate(
    genotype  = gsub("ex.vivo|in.vivo|interaction", "", coef),
    condition = case_when(
      grepl("interaction", coef) ~ "interaction",
      grepl("ex.vivo", coef)     ~ "ex.vivo",
      grepl("in.vivo", coef)     ~ "in.vivo",
      TRUE                       ~ "other"
    )
  )

###############
# --- Fast relevance computation (no slow mutate chains) ---
###############
# summary_df <- limmaRes %>%
#   group_by(ensg, genotype, celltype) %>%
#   summarise(
#     max_main_logFC = max(abs(logFC[condition %in% c("ex.vivo","in.vivo")]), na.rm = TRUE),
#     main_sig       = any(group != "n.s" & condition %in% c("ex.vivo","in.vivo")),
#     interaction_logFC = logFC[condition == "interaction"][1],
#     .groups = "drop"
#   ) %>%
#   mutate(
#     relevance = ifelse(
#       main_sig & !is.na(interaction_logFC) &
#         abs(interaction_logFC) < 0.5 * max_main_logFC,
#       "no_int", "int_relevant"
#     )
#   )
summary_df <- limmaRes %>%
  group_by(ensg, genotype, celltype) %>%
  summarise(
    max_main_logFC = max(abs(logFC[condition %in% c("ex.vivo","in.vivo")]),
                         na.rm = TRUE),
    main_sig = any(group != "n.s" & condition %in% c("ex.vivo","in.vivo")),
    interaction_sig = any(group != "n.s" & condition == "interaction"),
    .groups = "drop"
  ) %>%
  mutate(
    relevance = ifelse(
      main_sig & !interaction_sig,
      "no_int",
      "int_relevant"
    )
  )

limmaRes <- limmaRes %>%
  left_join(summary_df, by = c("ensg","genotype","celltype"))

###############
# --- Non-significant interaction genes ---
###############
non_sig_genes <- limmaRes %>%
  filter(condition == "interaction",
         relevance == "no_int",
         group == "n.s") %>%
  group_by(celltype, genotype) %>%
  summarise(non_sig = list(ensg), .groups = "drop")

###############
# --- EnrichR ---
###############
enrichr_dbs <- c(""KEGG_2019_Mouse"")

run_enrichr <- function(gene_list) {
  if (length(gene_list) == 0) return(NULL)
  enrichr(gene_list, enrichr_dbs)
}

non_sig_enrichment <- non_sig_genes %>%
  mutate(enrichr_res = map(non_sig, run_enrichr))

###############
# --- Flatten EnrichR results safely ---
###############
non_sig_enrichment_flat <- non_sig_enrichment %>%
  mutate(enrichr_flat = map(enrichr_res, function(res) {
    
    if (is.null(res)) return(tibble())
    
    bind_rows(lapply(names(res), function(db) {
      df <- res[[db]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$db <- db
      df
    }))
  })) %>%
  select(celltype, genotype, enrichr_flat) %>%
  unnest(enrichr_flat, keep_empty = TRUE)

###############
# --- Save ---
###############
write_rds(non_sig_enrichment_flat, basedir("enrichment_non_sig_genes.rds"))
write.table(non_sig_enrichment_flat,
            basedir("enrichment_non_sig_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

###############
# --- Filter enriched terms ---
###############
filtered_terms <- non_sig_enrichment_flat %>%
  filter(Odds.Ratio > 10, Adjusted.P.value < 0.01) %>%
  pull(Term) %>%
  unique()

###############
# --- Prepare plotting table ---
###############
pDT <- non_sig_enrichment_flat %>%
  filter(db == "MSigDB_Hallmark_2020",
         Term %in% filtered_terms) %>%
  mutate(overlap_num = as.numeric(sub("/.*", "", Overlap))) %>%
  filter(overlap_num > 1) %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 100)

###############
# --- Plot ---
###############
Fig1C <- ggplot(
  pDT,
  aes(
    x = genotype,
    y = Term,
    color = log2(Odds.Ratio),
    size = pmin(5, -log10(Adjusted.P.value))
  )
) +
  geom_point() +
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = TeX("log2(Odds.Ratio)")
  ) +
  scale_size_continuous(
    range = c(0, 1.8),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    x = NULL,
    y = "Pathways",
    title = "Enriched pathways"
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  facet_grid(~celltype)

ggsave(basedir("non_sig_interactions.pdf"), Fig1C, width = 10, height = 6)


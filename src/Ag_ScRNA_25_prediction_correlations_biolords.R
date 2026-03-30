# ============================================================
# Load libraries & custom functions
# ============================================================
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
source("src/Ag_ko_classification.R")

library(tidyverse)
library(enrichR)
library(scales)
library(patchwork)
library(cowplot)
library(latex2exp)
library(ggridges)
library(data.table)

# ============================================================
# Directories
# ============================================================
InDir_pred   <- dirout("Ag_ScRNA_15_celltype_biolord_limma_v2_mye")
InDir_limma  <- dirout("Ag_ScRNA_11_limma_all_ko_ex.vivo_vs_in.vivo_guide")
InDir_cor    <- dirout("Ag_ScRNA_11_limma_all_ko_ex.vivo_vs_in.vivo_correlation")
InDir_meta   <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
out      <- dirout("Ag_ScRNA_25_prediction_Biolord_correlations")



# ============================================================
# Metadata
# ============================================================
meta <- fread(InDir_meta("meta_cleaned.tsv")) %>% as.data.frame()
rownames(meta) <- meta[[1]]
meta <- meta[, -1, drop = FALSE]
colnames(meta) <- gsub("rowname", "sample", colnames(meta))
meta$sample
# ============================================================
# KO validity flags
# ============================================================
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarise(n = n_distinct(meta$sample), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = n, values_fill = 0) %>%
  mutate(valid_ko = in.vivo >= 3 & ex.vivo >= 3) %>%
  group_by(genotype, celltype) %>%
  summarise(valid_ko = any(valid_ko), .groups = "drop")

# ============================================================
# Load limma results
# ============================================================
limma_pred <- read_rds(InDir_pred(
  "limma_ex.vivo_vs_in.vivo_per_CT_all_coef_pred.rds"
))

limma_act <- read_rds(InDir_limma(
  "limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"
)) %>%
  filter(coef %in% limma_pred$coef)

# ============================================================
# In vivo logFC comparison (prediction vs actual)
# ============================================================
act_in <- limma_act %>%
  filter(str_detect(coef, "in.vivo")) %>%
  dplyr::select(logFC, ensg, celltype, coef)

pred_in <- limma_pred %>%
  filter(str_detect(coef, "in.vivo")) %>%
  filter(coef %in% act_in$coef) %>%
  dplyr::rename(logFC_pred = logFC) %>%
  dplyr::select(logFC_pred, ensg, celltype, coef)

combined <- left_join(pred_in, act_in,
                      by = c("coef", "ensg", "celltype")) %>%
  drop_na()

write_rds(combined, out("LogFC_invivo_pred_biolord.rds"))

# ============================================================
# Correlations per KO & celltype
# ============================================================
cor_pred <- combined %>%
  group_by(celltype, coef) %>%
  summarise(correlation = cor(logFC_pred, logFC),
            .groups = "drop") %>%
  mutate(
    genotype = gsub("in.vivo", "", coef),
    data = "prediction"
  ) %>%
  dplyr::select(celltype, genotype, correlation, data)

write_rds(cor_pred, out("Prediction-actual_cor_biolord.rds"))

# ============================================================
# Plot: prediction vs actual correlation
# ============================================================
ggplot(cor_pred, aes(genotype, correlation)) +
  geom_col(fill = "#4C889C") +
  facet_grid(rows = vars(celltype)) +
  optimized_theme_fig()

ggsave(out("percelltype_correlation.pdf"))

# ============================================================
# Ex vivo vs in vivo correlation (biological)
# ============================================================
cor_ex <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds")) %>%
  replace(is.na(.), 0)

cor_ex_long <- cor_ex %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(-celltype,
               names_to = "genotype",
               values_to = "correlation") %>%
  mutate(data = "ex.vivo")

# ============================================================
# DEG counts
# ============================================================
deg_in <- limma_act %>%
  filter(str_detect(coef, "in.vivo"),
         abs(logFC) > 1,
         adj.P.Val < 0.05) %>%
  group_by(coef, celltype) %>%
  summarise(num_degs_act = n(), .groups = "drop") %>%
  mutate(genotype = gsub("in.vivo", "", coef)) %>%
  dplyr::select(celltype, genotype, num_degs_act)

# ============================================================
# Merge all correlation data
# ============================================================
cor_all <- bind_rows(cor_pred, cor_ex_long) %>%
  inner_join(deg_in, by = c("celltype", "genotype")) %>%
  inner_join(ko_flags, by = c("celltype", "genotype")) %>%
  filter(valid_ko)

write_rds(cor_all, out("correlation_to_invivo_exvivo_vs_pred_biolord.rds"))

# ============================================================
# Final correlation plot
# ============================================================
ggplot(cor_all) +
  geom_point(
    aes(
      x = data,
      y = celltype,
      fill = correlation,
      size = pmin(3, log10(num_degs_act))
    ),
    shape = 21,
    color = "black",
    stroke = 0.4
  ) +
  facet_grid(cols = vars(genotype)) +
  scale_fill_gradient2(low = "#4C889C", mid = "white", high = "#D0154E",
                       limits = c(-1, 1)) +
  scale_size_continuous(
    range = c(0, 2),
    limits = c(0, 3),
    name = TeX("$\\log_{10}$(No. of DEGs)")
  ) +
  optimized_theme_fig() +
  theme(strip.text.x = element_text(angle = 90))

ggsave(out("correlation_pred_vs_exvivo_biolord.pdf"), width = 18, height = 5, units = "cm")

# ============================================================
# Example scatter plot
# ============================================================
combined %>%
  filter(coef == "in.vivoSetdb1", celltype == "GMP") %>%
  ggplot(aes(logFC_pred, logFC)) +
  geom_point(color = "#4C889C", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", linetype = "dashed") +
  
  theme_minimal()

ggsave(out("setdb1.pdf"))

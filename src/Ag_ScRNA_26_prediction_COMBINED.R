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
InDir_biolord   <- dirout("Ag_ScRNA_15_celltype_biolord_limma_v2_mye")
InDir_perturb   <- dirout("Ag_ScRNA_23_PerturbNet")
InDir_limma  <- dirout("Ag_ScRNA_11_limma_all_ko_ex.vivo_vs_in.vivo_guide")
InDir_cor    <- dirout("Ag_ScRNA_11_limma_all_ko_ex.vivo_vs_in.vivo_correlation")
InDir_meta   <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
out      <- dirout("Ag_ScRNA_26_prediction_COMBINED")
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
limma_biolord <- read_rds(InDir_biolord(
  "limma_ex.vivo_vs_in.vivo_per_CT_all_coef_pred.rds"
))

limma_perturb <- read_rds(InDir_perturb(
  "limma_ex.vivo_vs_in.vivo_per_CT_all_coef_pred_perturb.rds"
))
limma_act <- read_rds(InDir_limma(
  "limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"
)) 

# ============================================================
# In vivo logFC comparison (prediction vs actual)
# ============================================================
act_in <- limma_act %>%
  filter(str_detect(coef, "in.vivo")) %>%
  dplyr::select(logFC, ensg, celltype, coef)

biolord_in <- limma_biolord %>%
  filter(str_detect(coef, "in.vivo")) %>%
  filter(coef %in% act_in$coef) %>%
  dplyr::rename(logFC_biolord = logFC) %>%
  dplyr::select(logFC_biolord, ensg, celltype, coef)
perturb_in <- limma_perturb %>%
  filter(str_detect(coef, "in.vivo")) %>%
  filter(coef %in% act_in$coef) %>%
  dplyr::rename(logFC_perturb = logFC) %>%
  dplyr::select(logFC_perturb, ensg, celltype, coef)

combined_biolord <- left_join(act_in, biolord_in,
                      by = c("coef", "ensg", "celltype")) 

combined_perturb <- left_join( act_in, perturb_in,
                               by = c("coef", "ensg", "celltype")) 

################################################################



cor_ex <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds")) %>%
  replace(is.na(.), 0)

# Correlation of actual vs biolord
cor_biolord <- combined_biolord %>%
  group_by(celltype, coef) %>%
  summarise(
    correlation = cor(logFC_biolord, logFC, use = "pairwise.complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    genotype = str_remove(coef, "in\\.vivo_?"),  # remove "in.vivo" from coef
    data = "biolord"
  ) %>%
  dplyr::select(celltype, genotype, correlation, data)


# Correlation of actual vs biolord
cor_perturb <- combined_perturb %>%
  group_by(celltype, coef) %>%
  summarise(
    correlation = cor(logFC_perturb, logFC, use = "pairwise.complete.obs"),
    .groups = "drop"
  ) %>%
  mutate(
    genotype = str_remove(coef, "in\\.vivo_?"),  # remove "in.vivo" from coef
    data = "PerturbNet"
  ) %>%
  dplyr::select(celltype, genotype, correlation, data)


cor_ex_long <- cor_ex %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(-celltype,
               names_to = "genotype",
               values_to = "correlation") %>%
  mutate(data = "ex vivo")
##################


# Combine all correlation results into one tidy table
cor_all <- bind_rows(
  cor_biolord,
  cor_perturb,
  cor_ex_long
)

# Check the result

cor_all$celltype_genotype <- paste(cor_all$celltype, cor_all$genotype, sep = "_")

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
cor_all <- cor_all %>%
  inner_join(deg_in, by = c("celltype", "genotype")) %>%
  inner_join(ko_flags, by = c("celltype", "genotype")) %>%
  filter(valid_ko)

cor_all_clean <- cor_all %>% 
  filter(!is.na(correlation)) %>%
  mutate(
    data = fct_relevel(data, "ex vivo", "biolord", "PerturbNet")  # reorder factor levels
  )
cor_all_clean %>% write_rds(out("correlation_to_invivo_exvivo_vs_pred_all.rds"))

# Create boxplot with jitter
ggplot(cor_all_clean, aes(x = data, y = correlation)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5) +  # boxplot
  geom_jitter(aes(color = genotype), width = 0.2, size = 2, alpha = 0.8) +  # jitter points
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # vertical line at 1
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) + # y-axis scale
  labs(
    x = "Dataset",
    y = "Correlation with actual logFC",
    title = "Comparison of predicted vs actual logFC per dataset"
  ) +
  theme_bw() +
  facet_grid(cols = vars(celltype), scales = "free")+
  optimized_theme_fig()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
###############################

combined_all <- act_in %>%
  left_join(biolord_in,  by = c("coef", "ensg", "celltype")) %>%
  left_join(perturb_in, by = c("coef", "ensg", "celltype")) %>%
  pivot_longer(
    cols = c(logFC_biolord, logFC_perturb),
    names_to = "method",
    values_to = "logFC_pred"
  ) %>%
  mutate(
    method = recode(
      method,
      logFC_biolord = "Biolord",
      logFC_perturb = "PerturbNet"
    )
  ) %>%
  drop_na(logFC_pred, logFC)
combined_all %>% write_rds(out("combined_logFCs.rds"))
selected_coef <- c("in.vivoBrd9","in.vivoKmt2d","Setd1b")
selected_celltype <- "GMP"

plot_df <- combined_all %>%
  filter(
    coef %in% selected_coef,
    celltype == selected_celltype
  )%>%
  mutate(coef = gsub("in.vivo","", coef))

corr_df <- plot_df %>%
  group_by(method,coef) %>%
  summarise(
    correlation = cor(logFC_pred, logFC, use = "complete.obs")
  )
ggplot(plot_df, aes(x = logFC_pred, y = logFC)) +
  geom_point(alpha = 0.6, size = 2, color = "#4C889C") +
  geom_smooth(method = "lm", se = FALSE,
              linetype = "dashed", color = "black") +
  facet_grid(cols = vars(method), rows = vars(coef) ) +
  geom_text(
    data = corr_df,
    aes(
      x = -Inf, y = Inf,
      label = paste0("r = ", round(correlation, 2))
    ),
    hjust = -0.1, vjust = 1.2,
    inherit.aes = FALSE,
    size = 3
  ) +
  labs(
    title = paste("In vivo KO prediction vs observed:", selected_coef, "-", selected_celltype),
    x = "Predicted logFC",
    y = "Observed logFC"
  ) +
  optimized_theme_fig()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0))

ggsave(out("correlation_examples.pdf"), w = 8, h = 8, units = "cm")

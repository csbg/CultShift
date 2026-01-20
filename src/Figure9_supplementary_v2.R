source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
library(tidyverse)
library(enrichR)
library(purrr)
library("scales")
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library(readr)

#directories ------
base <-"Figure3"
basedir <- dirout("Figure3")
## InDir4 <- dirout("Figure2_Mye")
InDir6 <- dirout("Ag_ScRNA_12_fgsea_overlap")
#
###############
correlation_matrix_data <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(InDir_cor("DEGs_per_tissue.rds")) %>%
  dplyr::select(-correlation)

# Clean correlation matrix
correlation_matrix_no_na <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)

# Hierarchical clustering for column order
hc_cols <- hclust(dist(t(correlation_matrix_no_na)), method = "ward.D2")

# Reshape correlation data
correlation <- correlation_matrix_data %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(
    cols = colnames(correlation_matrix_data),
    names_to = "genotype",
    values_to = "correlation",
    values_drop_na = FALSE
  )
correlation$dataset <- str_wrap("Perturb-seq(Hematopoietic)", width = 11)
correlation_perturb <- correlation
# Prepare DEG data: select max num_degs per genotype/celltype
deg_plot_data <- deg_plot_data %>%
  na.omit() %>%
  group_by(genotype, celltype) %>%
  filter(num_degs == max(num_degs, na.rm = TRUE)) %>%
  dplyr::select(-condition) %>%
  unique()

# Merge DEG data with correlation data
correlation_deg <- inner_join(deg_plot_data, correlation,
                              by = c("celltype", "genotype"))

# Merge with KO flags to include valid KO status
correlation_deg_flagged <- correlation_deg %>%
  inner_join(ko_flags, by = c("genotype", "celltype")) %>%
  filter(valid_ko) %>%
  na.omit()

# Update genotype factor levels based on hierarchical clustering and mean correlation
ko_order <- correlation_deg_flagged %>%
  group_by(genotype) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_corr)) %>%
  pull(genotype)

correlation_deg_flagged <- correlation_deg_flagged %>%
  mutate(genotype = factor(genotype, levels = rev(ko_order)))
mean_corr <- correlation_deg_flagged %>%
  group_by(coef) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(coef = factor(coef, levels = rev(ko_order)))


#enrichment plot----------

# Load limma results
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds")) %>%
  mutate(coef = gsub("interaction", "", coef))

# Step 1: Filter for significant genes
limmaRes_significant <- limmaRes %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# summary_df described in KO classification
data <- summary_df %>%
  filter(coef %in% koi) %>%
  filter(Count >= 10) %>%
  inner_join(ko_flags, by = c("celltype", "coef")) %>%
  filter(valid_ko)

# Load enrichment results
enrichment_ntc_in.vivo <- read_rds(InDir6("enrichment_to_NTC_genes.rds"))

# Combine KO/celltype with enrichment results
combined <- data %>%
  dplyr::select(coef, celltype) %>%
  distinct() %>%
  left_join(enrichment_ntc_in.vivo, by = c("coef", "celltype"))

# Significance stars
combined <- combined %>%
  mutate(significance_en = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))

combined <- combined %>%
  mutate(log2.odds.ratio = log2(odds.ratio))

# Special case
combined <- combined %>%
  mutate(log2.odds.ratio = case_when(
    coef == "Smc3" ~ 0,
    TRUE ~ log2.odds.ratio
  ))

# Filter: keep only if overlap > 5
combined <- combined %>%
  mutate(log2.or.filtered = ifelse(overlap > 5, pmin(log2.odds.ratio, 7), NA))

# Cap -log10(p) at 10
combined <- combined %>%
  mutate(minuslog10p = pmin(-log10(p.value), 10))

# === NEW: Apply identical KO ordering to BOTH datasets ===
ko_levels <- rev(ko_order)

combined$coef  <- factor(combined$coef,  levels = ko_levels)
mean_corr$coef <- factor(mean_corr$coef, levels = ko_levels)

# optional: drop unused KOs in heatmap
mean_corr$coef <- droplevels(mean_corr$coef)

# --- Enrichment dot plot ---
p <- ggplot(
  combined,
  aes(
    x = coef,
    y = celltype,
    fill = log2.or.filtered,
    size = minuslog10p
  )
) +
  geom_point(shape = 21, colour = "black", stroke = 0.2) +
  scale_fill_gradient2(
    low  = "#CCE5FF",
    mid  = "white",
    high = "#FF9933",
    midpoint = 0,
    na.value = "grey80",
    name = "log2(OR)"
  ) +
  scale_size_continuous(range = c(0, 2), name = "-log10(Padj)") +
  labs(
    title = str_wrap("Overlap enrichment (NTC in vivo) per KO and cell type", width = 90),
    x = "KO",
    y = "Cell type"
  ) +
  optimized_theme_fig() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ko_present <- combined %>%
  droplevels() %>%        # remove unused factor levels
  pull(coef) %>%          # extract factor
  levels()                # get only levels that occur in the data

ko_present

heatmap_corr <- ggplot(mean_corr %>%
                         filter(coef %in% ko_present)
                       , aes(x = coef, y = 1, fill = mean_corr)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    midpoint = 0,
    na.value = "grey80",
    name = "Mean correlation",
    limits = c(-1, 1)
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


final_plot <- heatmap_corr / p + plot_layout(heights = c(0.2, 1))
final_plot

ggsave(basedir("Overlap_enrichment_mean_cor.pdf"),
       plot = final_plot,
       width = 10, height = 6, units = "cm")



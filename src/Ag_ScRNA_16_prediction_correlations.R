#load libraries and functions--------------------------------------------------
#
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
library(tidyverse)
library(enrichR)
library(purrr)
library("scales")
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library(ggridges)
source("src/Ag_ko_classification.R")

###############
InDir1 <- dirout("Ag_ScRNA_15_celltype_biolord_limma/")
#InDir2 <- dirout("Ag_ScRNA_11_limma_all_ko_ex.vivo_vs_in.vivo_guide/")
#Indir3 <- dirout("Ag_ScRNA_11_limma_all_ko_ex.vivo_vs_in.vivo_correlation/")
#InDir5 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
basedir <- dirout("Ag_ScRNA_16_prediction_correlations_new/")
source("src/Ag_Optimized_theme_fig.R")
#
# meta <- fread(InDir5("meta_cleaned.tsv")) # Read data
# meta <- as.data.frame(meta)               # Convert to dataframe (optional)
# rownames(meta) <- meta[[1]]   
# 
# meta <- meta[, -1, drop = FALSE] 
# colnames(meta) <- gsub("rowname","sample1", colnames(meta))
# # Check if there are at least 2 distinct samples per tissue for each genotype and celltype
# ko_flags <- meta %>%
#   group_by(genotype, celltype, tissue) %>%
#   summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
#   pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
#   mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
#   group_by(genotype, celltype) %>%
#   summarize(valid_ko = any(valid_ko), .groups = "drop")%>%
#   mutate(coef = genotype)

# replicates_per_ko <- meta %>%
#   group_by(genotype, celltype, tissue) %>%
#   summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
#   pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
#   mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
#   group_by(genotype, celltype) %>%
#   summarize(
#     valid_ko = any(valid_ko),
#     total_in_vivo = sum(in.vivo, na.rm = TRUE),
#     total_ex_vivo = sum(ex.vivo, na.rm = TRUE)
#     , .groups = "drop") %>%
#   mutate(coef = genotype)
limmaRes_pred <- read_rds(InDir1("limma_ex.vivo_vs_in.vivo_per_CT_all_coef_pred.rds"))
limmaRes_act <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))# %>%
  #filter(coef %in% unique(limmaRes_pred$coef))
unique(limmaRes_pred$coef)
unique(limmaRes_act$coef)
#processing
limmaRes_act_in.vivo <- limmaRes_act %>%
  filter(coef %in% grep("in.vivo",coef,value = T)) %>%
  dplyr::select(logFC, celltype,coef, ensg)

unique(limmaRes_act_in.vivo$coef)
limmaRes_pred$coef <- gsub("genotype","in.vivo",limmaRes_pred$coef)
limmaRes_pred_in.vivo <- limmaRes_pred %>%
  filter(coef %in% grep("in.vivo",coef,value = T)) %>%
  filter(coef %in% limmaRes_act_in.vivo$coef) %>%
  dplyr::rename(logFC_pred = logFC) %>%
  dplyr::select(logFC_pred, celltype,coef, ensg)
combined_genotype <- limmaRes_pred_in.vivo %>%
  left_join(limmaRes_act_in.vivo, by = c("coef","ensg","celltype")) %>%
  na.omit()
combined_genotype %>% write_rds(basedir("LogFC_invivo_pred.rds"))
combined_genotype_with_corr <- combined_genotype %>%
  group_by(celltype, coef) %>%
  mutate(correlation = cor(logFC_pred, logFC, use = "complete.obs")) %>%
  ungroup() %>%
  dplyr::select(celltype, coef, correlation) %>%
  distinct()

combined_genotype_with_corr$coef <- gsub("in.vivo","",combined_genotype_with_corr$coef )
combined_genotype_with_corr %>% write_rds(basedir("Prediction-actual_cor.rds"))
ggplot(combined_genotype_with_corr,
       aes(x = coef, y = correlation))+
  geom_bar(stat = "identity", fill = "#4C889C") +
  facet_grid(rows = vars(celltype))+optimized_theme_fig()
ggsave(basedir("percelltype_correlation.pdf"))
###########
#exvivo pred vs biolord pred ----------------------
###########

correlation_matrix_data <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds"))
deg_in <- limmaRes_act %>%
  filter(coef %in% grep("in.vivo",coef,value = T)) %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(coef, celltype) %>%
  summarise(num_degs_act = n())
deg_in$genotype <- gsub("in.vivo","", deg_in$coef)

# Prepare correlation data
correlation_matrix_no_na <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)

# Perform hierarchical clustering
hc_cols <- hclust(dist(t(correlation_matrix_no_na)), method = "ward.D2")
column_order <- colnames(correlation_matrix_data)[hc_cols$order]

# Reshape correlation data
correlation <- correlation_matrix_data %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(cols = colnames(correlation_matrix_data),
               names_to = "genotype",
               values_to = "correlation",
               values_drop_na = FALSE)%>%
  mutate(data = "ex.vivo")

combined_genotype_with_corr <- combined_genotype_with_corr %>%
  dplyr::rename(genotype = coef) 
combined_genotype_with_corr$data <- "prediction"
pred_vs_ex_corr <- rbind(combined_genotype_with_corr,correlation)
# Merge DEG data with correlation data
correlation_deg <- inner_join(deg_in[c("celltype","num_degs_act","genotype")], pred_vs_ex_corr,
                              by = c("celltype", "genotype"))

# Merge with KO flags to include valid KO status
correlation_deg_flagged <- correlation_deg %>%
  inner_join(ko_flags, by = c("genotype", "celltype")) %>%
  filter(valid_ko)%>%
  na.omit()%>%
  filter(valid_ko) 
correlation_deg_flagged %>% write_rds(basedir("correlation_to_invivo_exvivo_vs_pred.rds"))
ggplot(correlation_deg_flagged) +
  geom_point(aes(
    x = data,
    y = celltype,
    size = pmin(3,log10(num_degs_act)),
    fill = correlation  # Set transparency based on KO validity
  ),
  shape = 21,           # Use shape 21 to enable fill and color
  color = "black",       # Black outline
  stroke = 0.5 
  ) +
  facet_grid(cols = vars(genotype))+
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E"
  ) +
  scale_size_continuous(
    range = c(0,2),
    limits = c(0,3),
    breaks = c(1,2,3),
    name =TeX("$\\log_{10}\\; (\\No.\\; of \\;DEGs)$"))+
 labs(x = "KOs",
       y = "Celltype") +
  optimized_theme_fig()+
  theme(
    legend.justification = "right",
    strip.text.x = element_text(angle = 90, hjust = 0)
  )


ggsave(basedir("correlation_pred_vs_exvivo.pdf"), w = 18, h = 5, units = "cm")
##########################################
# Filter data for the selected coef and celltype
selected_coef <- "in.vivoSetdb1"
selected_celltype <- "GMP"
filtered_data <- combined_genotype %>%
  filter(coef == selected_coef, celltype == selected_celltype)

# Scatter plot with linear regression line
ggplot(filtered_data, aes(x = logFC_pred, y = logFC)) +
  geom_point(color = "#4C889C", alpha = 0.7, size = 3) +  # Scatter points
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +  # Trendline
  labs(
    title = paste("Correlation for", selected_coef, "in", selected_celltype),
    x = "Predicted logFC",
    y = "Actual logFC"
  ) +
  theme_minimal()
ggsave(basedir("setdb1.pdf"))  

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

basedir <- dirout("Figure4_Interactions_v2")
## InDir4 <- dirout("Figure2_Mye")
InDir6 <- dirout("Ag_ScRNA_12_fgsea_overlap")

limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))

################################################################################
#gene_based analysis

# Step 1: Filter for significant genes (adj.P.Val < 0.05 and abs(logFC) > 1)
limmaRes_significant <- limmaRes %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)  # Only significantly altered genes


data <- summary_df %>%
  filter(coef %in% koi) %>%
  filter(Count >= 10) %>%
  inner_join(ko_flags, by = c("celltype", "coef")) %>%
  filter(valid_ko)
data %>% write_rds(basedir("kos_fig3.rds"))




median_degs <- summary_df %>%
  filter(coef %in% koi)%>%
  filter(Count >= 10) %>%
  
  # 1. Sum Up + Down per celltype and KO
  group_by(celltype, coef) %>%
  summarise(total_degs = sum(Count), .groups = "drop") %>%
  
  # 2. Compute median across cell types per KO
  group_by(coef) %>%
  summarise(median_total_degs = median(total_degs), .groups = "drop")
length(unique(limmaRes_int$coef))
length(unique(median_degs$coef))
median(median(median_degs$median_total_degs))
#Fig3B----------------------------




#Fig3CA---------------------------(Correllation)----------
limmaRes_int <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))%>%
  mutate(genes = ensg)
#########
limmaRes_int %>%
  group_by(celltype,coef) %>%
  summar
#########
limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

correlation_results <- merged_data %>%
  inner_join(ko_flags, by = c("coef","celltype")) %>%
  filter(valid_ko) %>%
  group_by(coef, celltype) %>%
  filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results <- correlation_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(coef %in% koi)

# Step 3: Add significance labels based on p-value thresholds
correlation_results <- correlation_results %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )

correlation_results$sign_level <- factor(
  correlation_results$significance,
  levels = c("***","**","*","")
)

Fig3Ca_heatmap <- ggplot(correlation_results,
                         aes(x = coef, y = celltype, fill = cor_abs,
                             size = pmin(-log10(p_adj),5))) +
  
  geom_point(shape = 21) +
  #geom_text(aes(label = significance), size = 1.5) +
  
  scale_fill_gradient2(
    
    low = "purple",
    mid = "white",
    high = "#00A600",
    limits =c(-1,1),
    name = str_wrap(
      "Correlation",
      width = 20
    )
  )+
  
  scale_size_continuous(
    range = c(0,2),
    name = TeX("$\\log_{10}\\; (Padj)$")
  )+
  labs(title = str_wrap("Culture effects drive discrepancies of KO-effects in ex vivo models", 
                        width = 80),
       x = "KO",
       y = "Cell type") +
  
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
Fig3Ca_heatmap
ggsave(basedir("Fig3Ca_heatmap_perurb_small.pdf"),plot=Fig3Ca_heatmap,
       w=10.5,h = 3.8, units = "cm")
######
#Fig4C--------density-----------
#perturb_seq---------
limmaRes_int <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))%>%
  mutate(genes = ensg) 
deg_counts <- limmaRes_int %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.01) %>%
  group_by(coef, celltype) %>%
  tally(name = "n_DEGs") %>%
  mutate(valid = n_DEGs >= 10)  # TRUE if ≥10 DEGs, else FALSE

limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
merged_data <- limmaRes_int %>%
  inner_join(deg_counts %>% filter(valid), 
             by = c("coef", "celltype"))%>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

correlation_results_perturb <- merged_data %>%
  inner_join(ko_flags, by = c("coef","celltype")) %>%
  filter(valid_ko) %>%
  group_by(coef, celltype) %>%
  filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results_perturb <- correlation_results_perturb %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(coef %in% koi)

# Step 3: Add significance labels based on p-value thresholds
correlation_results_perturb <- correlation_results_perturb %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )

#Spleen-------------
out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"
#*
InDir1 <- dirout("Figure4")

limmaRes <- read_rds(InDir1("combined_jakstat_diff_exp.rds"))

limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))

# Modify the 'coef' column for any KO
limmaRes <- limmaRes %>%
  mutate(coef = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(cell_type = gsub("M", "Macrophage", cell_type)) %>%
  mutate(cell_type = gsub("T8", "T-cells", cell_type)) 
limmaRes$coef <- gsub("treatmentex_vivo","WT",
                      limmaRes$coef)

limmaRes_int <- limmaRes %>%
  filter(grepl("Interaction", limmaRes$coef))%>%
  filter(cell_type == "T-cells")
deg_counts <- limmaRes_int %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.01) %>%
  group_by(coef, cell_type) %>%
  tally(name = "n_DEGs") %>%
  mutate(valid = n_DEGs >= 10)

limmaRes_WT <- limmaRes %>%
  filter(coef == "WT")%>%
  filter(cell_type == "T-cells")%>%
  dplyr::select(-coef)
merged_data <- limmaRes_int %>%
  inner_join(deg_counts %>% filter(valid), 
             by = c("coef", "cell_type"))%>%
  inner_join(limmaRes_WT, by = c("gene","cell_type"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)%>%
  mutate(celltype = cell_type)

unique(merged_data$coef)
unique(merged_data$celltype)
correlation_results_spleen <- merged_data %>%
  group_by(coef, celltype) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results_spleen <- correlation_results_spleen %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Step 3: Add significance labels based on p-value thresholds
correlation_results_spleen <- correlation_results_spleen %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )

Fig3Ca_heatmap_spleen <- ggplot(correlation_results_spleen,
                         aes(x = gsub("Interaction_","",coef), 
                                      y = celltype,
                                      fill = cor_abs,
                             size = pmin(-log10(p_adj),5))) +
  
  geom_point(shape = 21) +
  #geom_text(aes(label = significance), size = 1.5) +
  
  scale_fill_gradient2(
    
    low = "purple",
    mid = "white",
    high = "#00A600",
    limits =c(-1,1),
    name = str_wrap(
      "Correlation",
      width = 20
    )
  )+
  
  scale_size_continuous(
    range = c(0,2),
    limits = c(0,5),
    name = TeX("$\\log_{10}\\; (Padj)$")
  )+
  labs(title = str_wrap("Culture effects drive discrepancies of KO-effects in ex vivo models", 
                        width = 50),
       x = "KO",
       y = "Cell type") +
  
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(basedir("Spleen_interaction_NTC_cor.pdf"), w = 3, h = 4,units = "cm")
#
#Glioblastoma---------------------
InDir2 <- dirout("Glioblastoma_limmaRes.3way.mod")

limmaRes <- read_rds(InDir2("limmaRes_threeway.rds"))
unique(limmaRes$coef)


  



limmaRes_NTC_noRT <- limmaRes%>%
  filter(coef == "tissueex.vivo")%>%
  dplyr::select(-coef)
limmaRes_int_glio <- limmaRes%>% 
  filter(grepl("tissueex.vivo.genotype[A-Za-z0-9]*$",coef))%>%
  mutate(coef = gsub("tissueex.vivo.genotype","", coef))

limmaRes_clean <- limmaRes %>%
  filter(
    coef %in% c(
      grep("tissueex.vivo.genotype[A-Za-z0-9]*$", coef, value = TRUE),
      "tissueex.vivo",
      "tissueex.vivo.RT_statusRT"
    )
  ) %>%
  mutate(
    coef = case_when(
      coef == "tissueex.vivo" ~ "culture_effect",
      coef == "tissueex.vivo.RT_statusRT" ~ "interaction_RT",
      str_detect(coef, "^tissueex\\.vivo\\.genotype") ~ str_replace(coef, "^tissueex\\.vivo\\.genotype", "Interaction_genotype_"),
      str_detect()
      TRUE ~ coef
    )
  )

unique(limmaRes$coef)
library(dplyr)
library(stringr)

# Extract only in.vivo_noRT and ex.vivo_noRT genotype coefficients
genotype_noRT <- limmaRes %>%
  filter(str_detect(coef, "_noRT$")) %>%                  # keep only noRT rows
  filter(str_detect(coef, "in\\.vivo|ex\\.vivo")) %>%    # keep in.vivo or ex.vivo
  filter(!str_detect(coef, "tissueex\\.vivo|RT_status")) # remove tissueex or other unwanted coeffs
unique(genotype_noRT$coef)

limmaRes_clean <- rbind(limmaRes_clean, genotype_noRT)
write_rds(limmaRes_clean,basedir("limmaRes_glio.rds"))

deg_counts <- limmaRes_int_glio %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.01) %>%
  group_by(coef) %>%
  tally(name = "n_DEGs") %>%
  mutate(valid = n_DEGs >= 10) 
merged_data <- limmaRes_int_glio %>%
  inner_join(deg_counts %>% filter(valid), 
             by = c("coef"))%>%
  inner_join(limmaRes_NTC_noRT, by = c("ensg"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)%>%
  mutate(celltype = "glio")

correlation_results_glio <- merged_data %>%
  group_by(coef, celltype) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results_glio <- correlation_results_glio %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Step 3: Add significance labels based on p-value thresholds
correlation_results_glio <- correlation_results_glio %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )
Fig3Ca_heatmap_glio <- ggplot(correlation_results_glio,
                                aes(x = gsub("Interaction","",coef), 
                                    y = celltype,
                                    fill = cor_abs,
                                    size = pmin(-log10(p_adj),5))) +
  
  geom_point(shape = 21) +
  #geom_text(aes(label = significance), size = 1.5) +
  
  scale_fill_gradient2(
    
    low = "purple",
    mid = "white",
    high = "#00A600",
    limits =c(-1,1),
    name = str_wrap(
      "Correlation",
      width = 20
    )
  )+
  
  scale_size_continuous(
    range = c(0,2),
    limits = c(0,5),
    name = TeX("$\\log_{10}\\; (Padj)$")
  )+
  labs(title = str_wrap("Culture effects drive discrepancies of KO-effects in ex vivo models", 
                        width = 70),
       x = "KO",
       y = "Cell type") +
  
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(basedir("Fig3Ca_heatmap_glio.pdf"),plot=Fig3Ca_heatmap_glio,
       w=9.5,h=4.5, units = "cm")
#Combined plot --------------------
correlation_results_glio$dataset <- "Glioblastoma"
correlation_results_spleen$dataset <- "Spleen"
correlation_results_perturb$dataset <- str_wrap("Perturb-seq(Hematiopoietic)", width = 11)
correlation_all <- rbind(correlation_results_perturb,
                         correlation_results_spleen, 
                         correlation_results_glio)

condition_colors <- c("#4E92D0" ,
                      "#478E5C",
                      "#AD1F69"
)
ggplot(correlation_all, aes(x = cor_abs, color = dataset)) +
  geom_density(alpha = 0.4, size = 0.2) +        # semi-transparent fill
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype ="dotted", color = "red", size =0.2)+
  #scale_fill_manual(values = condition_colors)+
  scale_color_manual(values = condition_colors)+
  labs(
    title = str_wrap("Correlation between culture effects and interaction effects",
                     width = 40),
    x = "Correlation",
    y = "Density",
    fill = "Dataset",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "top"
  )+optimized_theme_fig()+theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank())
ggsave(basedir("Correlation_density_interaction.NTC.only_Tcells.pdf"), w = 5.9, h = 3.5, units = "cm")
##############















#alternative
### ------------------------------------------------------------
### LOAD & PREPARE INPUT DATA
### ------------------------------------------------------------



merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype")) %>%
  mutate(
    logFC_KO = logFC.x,
    logFC_NTC = logFC.y,
    adj.P.Val_KO = adj.P.Val.x,
    adj.P.Val_NTC = adj.P.Val.y
  )

### Keep only real genes tested
all_genes <- unique(merged_data$genes)

### ------------------------------------------------------------
### SIGNIFICANT GENE SETS FOR VENN + FISHER TEST
### ------------------------------------------------------------

sig_sets <- merged_data %>%
  mutate(
    sig_NTC = (abs(logFC_NTC) > 1 & adj.P.Val_NTC < 0.05),
    sig_KO  = (abs(logFC_KO)  > 1 & adj.P.Val_KO  < 0.05)
  ) %>%
  group_by(coef, celltype) %>%
  summarise(
    set_NTC = list(genes[sig_NTC]),
    set_KO  = list(genes[sig_KO]),
    .groups = "drop"
  )

### ------------------------------------------------------------
### FISHER'S EXACT TEST FOR OVERLAP
### ------------------------------------------------------------

fisher_results <- sig_sets %>%
  rowwise() %>%
  mutate(
    a = length(intersect(set_NTC, set_KO)),
    b = length(set_NTC) - a,
    c = length(set_KO) - a,
    d = length(all_genes) - (a + b + c),
    fisher_p = fisher.test(matrix(c(a,b,c,d), nrow = 2))$p.value,
    OR = unname(fisher.test(matrix(c(a,b,c,d), nrow = 2))$estimate)
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(fisher_p, method = "BH"),
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE           ~ ""
    )
  )

### ------------------------------------------------------------
### CORRELATION ANALYSIS
### ------------------------------------------------------------

correlation_results <- merged_data %>%
  inner_join(ko_flags, by = c("coef","celltype")) %>%
  filter(valid_ko) %>%
  group_by(coef, celltype) %>%
  filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>%
  summarise(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO)),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO))$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

### ------------------------------------------------------------
### MERGE CORRELATION + FISHER DATA
### ------------------------------------------------------------

cor_fisher <- correlation_results %>%
  inner_join(
    fisher_results %>% dplyr::select(coef, celltype, OR, fisher_p, p_adj, significance),
    by = c("coef","celltype")
  )

### ------------------------------------------------------------
### FIGURE 3C — WITH FISHER SIGNIFICANCE LABELS
### ------------------------------------------------------------

Fig3Ca <- ggplot(cor_fisher, aes(x = coef, y = cor_abs)) +
  geom_col(color = "grey30", fill = NA, width = 0.6) +
  
  # fisher significance (*, **, ***)
  geom_text(
    aes(label = significance),
    y = 0.75,
    color = "black",
    size = 2.4
  ) +
  
  labs(
    title = "Correlation of interaction effect to culture effect",
    x = NULL,
    y = expression(atop("Correlation of interaction effect","to culture effects"))
  ) +
  facet_grid(cols = vars(celltype), scales = "free_x", space = "free_x") +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

Fig3Ca

### ------------------------------------------------------------
### (OPTIONAL) SAVE TABLES
### ------------------------------------------------------------

write_tsv(fisher_results, basedir("Fisher_overlap_results.tsv"))
write_tsv(cor_fisher, basedir("correlation_plus_fisher.tsv"))

#Fig3Cb DEG interaction logFC ---------------------------------------------
ggsave(basedir("Fig3Ca.pdf"),plot=Fig3Ca,
       w=12,h=4, units = "cm")
summary_total <- data %>%
  filter(valid_ko)%>%
  group_by(celltype, coef, genotype, valid_ko) %>%
  summarise(Total_Regulated = sum(Count), .groups = "drop")

Fig3Cb <- ggplot(summary_total,aes(
  x = coef,
  y = log10(Total_Regulated)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  facet_grid(cols = vars(celltype), space = "free", scales = "free") +
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
#strip.text.x = element_blank(),
#axis.ticks.x = element_blank())


# Display the plot
Fig3Cb

ggsave(basedir("Fig3Cb.pdf"),plot=Fig3Cb,
       w=12,h=4, units = "cm")

# combine with enrichment to ntc
#Fig3Cc-------------

enrichment_ntc_in.vivo <- read_rds(InDir6("enrichment_to_NTC_genes.rds"))
combined <- data %>%
  dplyr::select(coef,celltype)%>%
  distinct()%>%
  left_join(enrichment_ntc_in.vivo, by = c("coef","celltype"))

# Step 1: Combine the correlation data and enrichment data
combined <- combined %>%
  mutate(significance_en = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))
combined <- combined %>%
  mutate(log2.odds.ratio =log2(odds.ratio))

combined <- combined %>%
  mutate(log2.odds.ratio = case_when(
    coef == "Smc3" ~ 0,
    TRUE ~ odds.ratio
  ))

combined <- combined %>%
  mutate(
    log2.or.filtered = ifelse(overlap > 5, pmin(log2.odds.ratio, 7), NA)  # keep value only if overlap>5
  )

Fig3Cc <- ggplot(combined, aes(x = coef, y = log2.or.filtered)) +
  geom_col(color = "darkgrey", fill =NA,  width = 0.6) +
  
  # significance labels (for overlap > 5 only, already handled by filtering)
  geom_text(
    data = combined %>% filter(overlap > 5),
    aes(label = significance_en),
    y = 6.5,
    color = "black",
    size = 1.5
  ) +
  
  # add "NA" label for overlap <= 5
  geom_text(
    data = combined %>% filter(overlap <= 5),
    aes(label = "NA"),
    y = 1.5,   # adjust position
    color = "black",
    size = 1.5
  ) +
  
  facet_grid(cols = vars(celltype), scales = "free", space = "free_x") +
  labs(
    x = "KOs",
    y = TeX("$\\log_{2}\\; (Odds ratio)$"),
    title = "KOs with large overlap to culture effect genes show more discordant effects between ex vivo and in vivo models"
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),  # 👈 center facet labels
    panel.spacing = unit(0.02, "cm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


Fig3Cc

ggsave(basedir("Fig3Cc.pdf"),plot=Fig3Cc,
       w=14.8,h=4, units = "cm")


#Fig3C-----------
Fig3C <-  Fig3Ca / Fig3Cb /  Fig3Cc
ggsave(basedir("Fig3C.pdf"),plot=Fig3C,
       w = 16.5,h=12, units = "cm")

#combined-----------
row1 <- (plot_spacer() | Fig3B) +
  plot_layout(widths = c(1, 1.5))

ggsave(basedir("row1.pdf"), plot = row1, w = 18, h = 8, units = "cm" )
ggsave(basedir("row1_title.pdf"), plot = row1, w = 18, h = 11, units = "cm" )



# Create a truly small blank plot
blank_plot <- ggplot() + theme_void()

# Combine with Fig3C
row2 <- plot_grid(
  blank_plot,
  Fig3C,
  rel_widths = c(1, 3),  # Use relative widths
  nrow = 1
)

#ggsave(basedir("row2_height.pdf"), plot = row2, width = 18, height = 11, units = "cm")


ggsave(basedir("row2.pdf"), row2, w = 18, h = 12, units = "cm" )
#ggsave(basedir("row2_height1.pdf"), row2, w = 18, h = 11, units = "cm" )
##################


####################################################################
#***Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01***#
####################################################################
#limma with all 
#data
#
#Aarathy
###############
###############
source("src/00_init.R")
source("src/Ag_ScRNA_11_invivo_exvivo_KO_limma_function.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)

#####################################################################
inDir <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")

InDir_biolord   <- dirout("Ag_ScRNA_15_celltype_biolord_limma_v2_mye")
InDir_perturb   <- dirout("Ag_ScRNA_23_PerturbNet")
InDir_all <- dirout("Ag_ScRNA_26_prediction_COMBINED_dataVoom")

basedir <- dirout("Ag_ScRNA_24_predicted__pergene_plot")
                  

############################
## 0. METADATA
############################
meta_obs  <- meta
meta_biolord <- read_rds(InDir_biolord("meta.rds"))
meta_perturb <- read_rds(InDir_perturb("meta_obs_pred.rds"))
meta_all <- read_rds(InDir_all("meta_all.rds"))
meta_all$tissue <- gsub("invivo","in.vivo",meta_all$tissue)
############################
## 1. PREP LIMMA TABLES
############################
limmaRes_perturb <- read_rds(InDir_perturb("limma_ex.vivo_vs_in.vivo_per_CT_all_coef_pred_perturb.rds"))
limmaRes_biolord <- read_rds(InDir_biolord("limma_ex.vivo_vs_in.vivo_per_CT_all_coef_pred.rds"))
limmaRes_all <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))

# dataVoom_GMP_biolord <- read_rds(InDir_biolord("GMP_dataVoom.rds"))
# dataVoom_GMP_perturb <- read_rds(InDir_perturb("GMP_dataVoom.rds"))
# dataVoom_GMP <-  read_rds(InDir_int("GMP_dataVoom.rds"))
dataVoom_all <- read_rds(InDir_all("dataVoom.rds"))
## ---- PREDICTED ----
limmaRes_biolord$gene <- limmaRes_biolord$ensg
limmaRes_perturb$gene <- limmaRes_perturb$ensg
limmaRes_all$gene <- limmaRes_all$ensg


# Add a column indicating the source of predictions
limmaRes_perturb <- limmaRes_perturb %>%
  mutate(
    model = "PerturbNet",
    tissue = ifelse(
      startsWith(coef, "interaction"), "interaction",
      ifelse(startsWith(coef, "ex.vivo"), "ex.vivo",
             ifelse(startsWith(coef, "in.vivo"), "in.vivo", NA))
    ),
    comparison = sub("^(interaction|ex\\.vivo|in\\.vivo)", "", coef)
  )

unique(limmaRes_perturb$coef)


limmaRes_biolord <- limmaRes_biolord %>%
  mutate(
    model = "biolord",
    tissue = ifelse(
      startsWith(coef, "interaction"), "interaction",
      ifelse(startsWith(coef, "ex.vivo"), "ex.vivo",
             ifelse(startsWith(coef, "in.vivo"), "in.vivo", NA))
    ),
    comparison = sub("^(interaction|ex\\.vivo|in\\.vivo)", "", coef)
  )





# Keep only relevant columns for consistency
limmaRes_biolord <- limmaRes_biolord %>%
  dplyr::select(gene, coef, logFC, adj.P.Val, tissue, comparison, celltype, model)

limmaRes_perturb <- limmaRes_perturb %>%
  dplyr::select(gene, coef, logFC, adj.P.Val, tissue, comparison, celltype, model)

# Combine both
limmaRes_pred <- bind_rows(limmaRes_biolord, limmaRes_perturb)



## ---- OBSERVED ----


limmaRes_all <- limmaRes_all %>%
  mutate(
    model = "Observed",
    tissue = ifelse(
      startsWith(coef, "interaction"), "interaction",
      ifelse(startsWith(coef, "ex.vivo"), "ex.vivo",
             ifelse(startsWith(coef, "in.vivo"), "in.vivo", NA))
    ),
    comparison = sub("^(interaction|ex\\.vivo|in\\.vivo)", "", coef)
  )


############################
## PARAMETERS
############################
KO  <- c("Brd9")
ct  <- "GMP"
goi <- c("Ifit1")

## ---- OBSERVED ----

# goi = vector of genes
# KO  = vector of KO genotypes

dat_all <- meta_all[names(dataVoom_all$E[rownames(dataVoom_all$E) %in% goi, ]),
] %>%
  mutate(
    E = dataVoom_all$E[goi, ],
    gene = goi,
    comparison = KO
  ) %>%
  filter(genotype %in% c(KO, "NTC")) %>%
  mutate(genotype = factor(genotype, levels = c("NTC", KO)))
## 3. MERGE WITH LIMMA
############################

limmaRes_pred_invivo <- limmaRes_pred%>%
  filter(tissue == "in.vivo")
unique(limmaRes_pred_invivo$coef)
unique(limmaRes_pred_invivo$tissue)
colnames(limmaRes_pred)


limmaRes_all <- limmaRes_all[,colnames(limmaRes_pred_invivo)]
limmaRes <- rbind(limmaRes_pred_invivo,limmaRes_all)
unique(dat_all$comparison)
unique(limmaRes$comparison)

merged_all <- dat_all %>%
  left_join(
    limmaRes,
    by = c("gene", "celltype", "tissue","model","comparison")
  )

############################
## 4. SIGNIFICANCE + Y POS
############################

prep_plot_df <- function(df) {
  df %>%
    group_by(tissue,model) %>%
    mutate(
      y_pos = max(E, na.rm = TRUE) * 1.1,
      significance = case_when(
        adj.P.Val < 0.001 ~ "***",
        adj.P.Val < 0.01  ~ "**",
        adj.P.Val < 0.05  ~ "*",
        TRUE ~ "ns") )%>%
    ungroup()
}


plot_all  <- prep_plot_df(merged_all)
plot_all$tissue <- gsub("invivo","in.vivo",plot_all$tissue)
############################
## 5. COMBINE
############################
plot_all <- plot_all %>%
  mutate(
    group = paste0(model, ":", genotype),
    group = factor(
      group,
      levels = c(
        "Observed:NTC", 
        paste0("Observed:", KO),
        paste0("PerturbNet:", KO),
        paste0("biolord:", KO)
      )
    ),
    fill_group = case_when(
      genotype == "NTC" ~ "Observed",
      genotype == "Brd9" & model %in% c("PerturbNet", "biolord") ~ "Predicted",
      genotype == "Brd9" & model == "Observed" ~ "Observed",
      TRUE ~ NA_character_  # catch other cases
    )
  )
                                     
plot_df2 <- plot_all %>%
  filter(celltype == "GMP") 
unique(plot_df2$group)
plot_df2 <- plot_df2 %>%
  mutate(color_group = ifelse(tissue == "ex.vivo", "ex.vivo", fill_group))
plot_df2 %>% write_rds(basedir("Ifit1_Brd9.rds"))
############################
## 6. PLOT
############################


ggsave(basedir("Ifit1_Brd9_in.vivo.pdf"))
###################

############################
## PARAMETERS
############################
KO  <- "Setdb1"
ct  <- "GMP"
goi <- "Myc"

## ---- OBSERVED ----

# goi = vector of genes
# KO  = vector of KO genotypes

dat_all <- meta_all[names(dataVoom_all$E[rownames(dataVoom_all$E) %in% goi, ]),
] %>%
  mutate(
    E = dataVoom_all$E[goi, ],
    gene = goi,
    comparison = KO
  ) %>%
  filter(genotype %in% c(KO, "NTC")) %>%
  mutate(genotype = factor(genotype, levels = c("NTC", KO)))
## 3. MERGE WITH LIMMA
############################

limmaRes_pred_invivo <- limmaRes_pred%>%
  filter(tissue == "in.vivo")
unique(limmaRes_pred_invivo$coef)
unique(limmaRes_pred_invivo$tissue)
colnames(limmaRes_pred)


limmaRes_all <- limmaRes_all[,colnames(limmaRes_pred_invivo)]
limmaRes <- rbind(limmaRes_pred_invivo,limmaRes_all)
unique(dat_all$comparison)
unique(limmaRes$comparison)

merged_all <- dat_all %>%
  left_join(
    limmaRes,
    by = c("gene", "celltype", "tissue","model","comparison")
  )

############################
## 4. SIGNIFICANCE + Y POS
############################

prep_plot_df <- function(df) {
  df %>%
    group_by(tissue,model) %>%
    mutate(
      y_pos = max(E, na.rm = TRUE) * 1.1,
      significance = case_when(
        adj.P.Val < 0.001 ~ "***",
        adj.P.Val < 0.01  ~ "**",
        adj.P.Val < 0.05  ~ "*",
        TRUE ~ "ns") )%>%
    ungroup()
}


plot_all  <- prep_plot_df(merged_all)
plot_all$tissue <- gsub("invivo","in.vivo",plot_all$tissue)
############################
## 5. COMBINE
############################
plot_all <- plot_all %>%
  mutate(
    group = paste0(model, ":", genotype),
    group = factor(
      group,
      levels = c(
        "Observed:NTC", 
        paste0("Observed:", KO),
        paste0("PerturbNet:", KO),
        paste0("biolord:", KO)
      )
    ),
    fill_group = case_when(
      genotype == "NTC" ~ "Observed",
      genotype == KO & model %in% c("PerturbNet", "biolord") ~ "Predicted",
      genotype == KO & model == "Observed" ~ "Observed",
      TRUE ~ NA_character_  # catch other cases
    )
  )

plot_df2 <- plot_all %>%
  filter(celltype == "GMP") 
unique(plot_df2$group)
plot_df2 <- plot_df2 %>%
  mutate(color_group = ifelse(tissue == "ex.vivo", "ex.vivo", fill_group))
plot_df2 %>% write_rds(basedir("Myc_Setdb1.rds"))
#

ggplot(
  plot_df2,
  aes(x = group, y = E, color = color_group)
) +
  geom_boxplot(
    outlier.shape = NA,
    size = 0.2,
    width = 0.6,
    position = position_dodge2(preserve = "single")
  ) +
  scale_color_manual(
    values = c(
      Observed  = "#d38d5fff",
      Predicted = "#46beea",
      "ex.vivo" = "#6a3d9aff"
    )
  ) +
  labs(
    title = paste0(goi, " expression in ", ct, " (", KO, " KO)"),
    y = "Expression",
    x = "Genotype"
  ) +
  facet_grid(cols = vars(tissue), scales = "free", space = "free") +
  optimized_theme_fig() +
  geom_text(
    data = plot_df2 %>% distinct(model, significance),
    aes(x = 1.5, label = significance),
    inherit.aes = FALSE,
    y = 3.5,
    size = 2.5
  ) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

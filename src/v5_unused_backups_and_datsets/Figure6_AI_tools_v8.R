source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
source("src/Ag_ko_classification.R")
library(tidyverse)
library(enrichR)
library(purrr)
library(scales)
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library(ggridges)
library(ggsci)  # Optional for clean palettes like "npg"
#directories ------
#
base <- "Figure6"
basedir <- dirout("Figure6")


###############
InDir1 <- dirout("Ag_ScRNA_25_prediction_Biolord_correlations/")
InDir2 <- dirout("Ag_ScRNA_25_prediction_Perturb_correlations/")
InDir3 <- dirout("Ag_ScRNA_26_prediction_COMBINED")
InDir_all <-  dirout("Ag_ScRNA_24_predicted__pergene_plot")
#Fig5.1
correlation_deg_flagged_all <- read_rds(InDir3("correlation_to_invivo_exvivo_vs_pred_all.rds"))
#setting colours--------
# All genotypes that ever appear anywhere

all_genotypes <- sort(unique(c(
  correlation_deg_flagged_all$genotype,
  correlation_deg_flagged_all$genotype,
  meta$genotype
)))

# Fixed color palette
genotype_colors <- setNames(
  hue_pal()(length(all_genotypes)),
  all_genotypes
)

# Fixed shapes (ggplot supports ~25 clearly)
# Choose visually distinct shapes
distinct_shapes <- c(21, 22, 23, 24, 25, 3, 4, 8)  # filled + open shapes

# Assign shapes cyclically if there are more genotypes than shapes
genotype_shapes <- setNames(
  rep(distinct_shapes, length.out = length(all_genotypes)),
  all_genotypes
)


#Option 2 Fig6A
ggplot(correlation_deg_flagged_all,
       aes(x = data, y = correlation)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.4,
    position = position_dodge2(
      width = 0.6,
      preserve = "single"   # ✅ THIS enforces identical widths
    )
  ) +
  geom_jitter(
    aes(color = genotype, shape = genotype,size = pmin(1,log10(num_degs_act))),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.6),
   
    alpha = 0.7
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_y_continuous(limits = c(-0.7, 1)) +
  scale_x_discrete(drop = TRUE) +
  scale_shape_manual(values = genotype_shapes) +
  scale_size_continuous(
    range = c(0,1),
    limits = c(0,3),
    breaks = c(1,2,3),
    name = expression(atop("No. of genes", log[10](n))))+
  scale_color_manual(values = genotype_colors, drop = TRUE) +
  facet_grid(
    cols = vars(celltype),
    scales = "free_x",   # ✅ now allowed
    space  = "free_x"
  ) +
  labs(
    x = NULL,
    y = "Pearson correlation",
    title = "Correlation to observed in vivo KO effects shows suboptimal prediction by existing AI tools"
  ) +
  optimized_theme_fig()+
  
  theme(
    legend.position = "none")
  #   legend.direction = "vertical",
  #   legend.box = "horizontal",
  #   legend.box.just = "center",
  #   legend.key.height = unit(0.3, "cm"),
  #   legend.key.width  = unit(0.4, "cm"),
  #   legend.spacing.x  = unit(0.2, "cm"),
  #   legend.margin     = margin(t = 2, b = 2)
  # )
  
ggsave(basedir("Fig6B_Overall_prediction_correlation_B.pdf"), w=18 , h = 4.5, units = "cm")
#summary for text
unique(correlation_deg_flagged_all$data)
# 1. Extract ex vivo reference
exvivo_ref <- correlation_deg_flagged_all %>%
  filter(data == "ex vivo") %>%
  dplyr::select(celltype, genotype, exvivo_cor = correlation)

# 2. Join and compare
summary_table <- correlation_deg_flagged_all %>%
  left_join(exvivo_ref, by = c("celltype", "genotype")) %>%
  filter(data != "ex vivo") %>%  # remove ex vivo itself
  group_by(data) %>%
  summarise(
    total_cases = n(),
    better_than_exvivo = sum(correlation > exvivo_cor, na.rm = TRUE),
    better_and_above_0_5 = sum(correlation > exvivo_cor & correlation > 0.5, na.rm = TRUE),
    better_and_above_0_8 = sum(correlation > exvivo_cor & correlation > 0.8, na.rm = TRUE),
    better_but_below_0_5 = sum(correlation > exvivo_cor & correlation < 0.5, na.rm = TRUE),
    frac_above_0_5_given_better = better_and_above_0_5 / better_than_exvivo,
    frac_total = better_and_above_0_5 / total_cases
  ) %>%
  arrange(desc(better_and_above_0_5))

summary_table
##################################

#Fig6B



correlation_long <- correlation_deg_flagged_all %>%
  pivot_wider(
    names_from = data,
    values_from = correlation
  ) %>%
  dplyr::rename(
    cor_exvivo  = `ex vivo`,
    Biolord = biolord,
    PerturbNet = PerturbNet
  ) %>%
  pivot_longer(
    cols = c(Biolord, PerturbNet),
    names_to = "method",
    values_to = "cor_pred"
  ) %>%
  filter(!is.na(cor_exvivo), !is.na(cor_pred))
#
correlation_test <- correlation_long %>%
  na.omit()
# drop unused factor levels
correlation_test <- correlation_test %>%
  droplevels()

correlation_test <- correlation_test %>%
  mutate(
    celltype = as.character(celltype),   # ensure it is character
    celltype = fct_relevel(celltype, c(
      "GMP", "Gran.P", "HSC", "MkP", "Mono", "Eo.Ba", "Gran.", "MEP.early"
    ))
  )
# Set the desired order of methods
correlation_test$method <- factor(
  correlation_test$method,
  levels = c("Observed","PerturbNet", "Biolord")  # adjust as needed
)
# Now plot row-wise (method as rows, celltype as columns)
ggplot(correlation_test, aes(x = cor_exvivo, y = cor_pred)) +
  geom_point(aes(
    color = genotype,
    shape = genotype,
    size  = pmin(2, log10(num_degs_act))
  ),
  alpha = 0.8,
  stroke = 0.4
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "red") +
  facet_grid(
    method ~ celltype,
    scales = "fixed",  # or "free" if you want axis per celltype
    drop = TRUE         # ✅ important
  ) +
  scale_color_manual(values = genotype_colors, drop = FALSE) +
  scale_shape_manual(values = genotype_shapes, drop = FALSE) +
  scale_size_continuous(
    range = c(0,1),
    limits = c(0,3),
    breaks = c(1,2,3),
    name = expression(atop("No. of genes", log[10](n))))+
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    x = "Correlation to Ex Vivo",
    y = "Correlation to Prediction",
    color = "Genotype",
    shape = "Genotype"
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "right",
    strip.text = element_text(size = 5, face = "bold")
    )+
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0.),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.4, "cm"),
    legend.spacing.x  = unit(0.1, "cm"),
    legend.margin     = margin(t = 2, b = 2)
  )
  


ggsave(basedir("Fig.6Cprediction_compared_to_ex.vivo.pdf"), w= 17,
       h = 5, units = "cm")
#######################
#Sup.Fig. 17.B---------------
plot_df2 <- read_rds(InDir_all("Ifit1_Brd9.rds"))
goi <- "Ifit1"
ct <-"GMP"
KO <- "Brd9"
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

ggsave(basedir("Ifit1_Brd9_in.vivo.pdf"), w =8, h= 5, units = "cm")
##########

#######################
#Sup.Fig. 17.B---------------
plot_df2 <- read_rds(InDir_all("Myc_Setdb1.rds"))
goi <- "Myc"
ct <-"GMP"
KO <- "Setdb1"




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

ggsave(basedir("Myc_Setdb1_in.vivo.pdf"), w =8, h= 5, units = "cm")


###########
#summaries----
better_counts <- correlation_long %>%
  
  # keep only valid comparisons
  filter(!is.na(cor_exvivo), !is.na(cor_pred)) %>%
  
  mutate(better_than_exvivo = cor_pred > cor_exvivo) %>%
  
  group_by(method) %>%
  summarise(
    n_better = sum(better_than_exvivo),
    total    = n(),
    frac_better = n_better / total,
    .groups = "drop"
  )
better_counts

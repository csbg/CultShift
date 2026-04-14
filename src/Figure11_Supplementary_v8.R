source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")


require(tidyverse)

basedir <- dirout("Figure11_Supplementary/")
InDir <- dirout("Figure3_correlation_ex_in_v8")
InDir1 <- dirout("Ag_ScRNA_16_Linear_model_predictionNF")
InDir2 <- dirout("Ag_ScRNA_26_prediction_COMBINED")
InDir3 <- dirout("Ag_ScRNA_24_predicted__pergene_plot")
# load, format and merge data ---------------------------------------------------



#Sup.Fig.17A

pred_act <- read_csv(InDir1("pearson correlations_pred_act.csv"))

deg <- read_rds(InDir_cor(paste0("DEGs_per_tissue.rds")))  

deg <- deg %>%
  filter(condition == "In Vivo")%>%
  dplyr::select("celltype", "genotype","num_degs")%>%
  dplyr::rename(knockout = genotype)
pred_act <- merge(pred_act, deg,by = c("celltype", "knockout"))

combine <- pred_act
unique(combine$model)

model_names = c(
  "exVivoOnly" = "0 + logFC_ex.vivo",
  "interactionAndIntercept" = "logFC_ex.vivo * logFC_NTC"
  #"ex vivo" = "baseline"
)
unique(combine$model)

ko_flag <- ko_flags%>%
  dplyr::rename(knockout = genotype)
combine <- combine %>%
  inner_join(ko_flag, by = c("knockout", "celltype"))

combine <- combine %>% filter(valid_ko)
combine <- combine %>%
  filter(model %in% c(
    "exVivoOnly",
    "interactionAndIntercept"  
  ))
column_order <- read_rds(InDir("column_order.rds"))
combine$knockout <- factor(combine$knockout,
                           levels = column_order)
#geom_point

Sup.Fig.11A <- combine %>%
  filter(valid_ko)%>%
  filter(ct_all != "individual cts") %>%
  ggplot(aes(
    x = knockout,
    y = celltype,
    size = pmin(0,log10(sig_in)),
    fill = cor  #
  )) +
  geom_point(aes(
    size = pmin(3,log10(num_degs)),
    fill = cor  # Set transparency based on KO validity
  ),
  shape = 21,           # Use shape 21 to enable fill and color
  color = "black",       # Black outline
  stroke = 0.5          # Adjust the width of the outline
  ) +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    limits = c(-1,1),
    name = expression("Pearson's correlation")
  ) +
  scale_size_continuous(
    range = c(0,2.5),
    #limits = c(0,2.5),
    breaks = c(1,2,3),
    name = expression(atop("No. of genes", log[10](n)))
  )+
  facet_grid(ko_all ~ model,
             space = "free",
             scales = "free",
             labeller = labeller(model = model_names))+
  labs(x = "KOs",
       y = "Cell type",
       title =  "Correlation of KO-effects (actual versus predicted)") +
  optimized_theme_fig()+
  theme(
    
    legend.position  = "bottom",
    panel.spacing = unit(0.1,"cm")
  )
ggsave(basedir("Sup.Fig.11A_dot.pdf"),plot = Sup.Fig.11A, w = 18, h = 10,units = "cm")
###############################
# #Sup.Fig. 17.B---------------
# plot_df2 <- read_rds(InDir3("Ifit1_Brd9.rds"))
# goi <- "Ifit1"
# ct <-"GMP"
# KO <- "Brd9"
# ggplot(
#   plot_df2,
#   aes(x = group, y = E, color = fill_group)
# ) +
#   geom_boxplot(
#     outlier.shape = NA,
#     size = 0.2
#   ) +
#   
#   scale_color_manual(
#     values = c(
#       Observed = "#d38d5fff",
#       Predicted= "#46beea"
#     )
#   ) +
#   labs(
#     title = paste0(goi, " expression in ", ct, " (", KO, " KO)"),
#     y = "Expression",
#     x = NULL
#   ) +
#   theme(legend.position = "none") +
#   optimized_theme_fig() +
#   theme(
#     panel.grid = element_blank(),
#     strip.background = element_blank()
#   ) +
#   geom_text(
#     data = plot_df2 %>%
#       distinct(model, significance),
#     aes(x = 1.5, label = significance),
#     inherit.aes = FALSE,
#     y = 3.5,
#     size = 2.5
#   ) +
#   optimized_theme_fig()+
#   
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# ggsave(basedir("Ifit1_Brd9_in.vivo.pdf"), w = 7, h = 4.5, units = "cm")
##################

#Sup.Fig. 17 C---------
combined_all <- read_rds(InDir2("combined_logFCs.rds"))

cor_all <- combined_all %>%
  group_by(method,coef) %>%
  summarise(
    correlation = cor(logFC_pred, logFC, use = "complete.obs")
  )
selected_coef <- c("in.vivoBrd9","in.vivoSetdb1")
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
  
  geom_hex(bins = 50) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") + 
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.2, color ="#e41a1c") +
  #scale_fill_viridis_c() +
  #scale_fill_viridis_c() +  # Optional: nice continuous color scale
  
  
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

ggsave(basedir("Fig.17B.correlation_examples.pdf"), w = 7, h = 5, units = "cm")

#####################


# Compute correlations
cor_all <- read_rds(InDir2("correlation_to_invivo_exvivo_vs_pred_all.rds"))

# Step 1: pivot wider
cor_wide <- cor_all %>%
  filter(!is.na(correlation)) %>%
  mutate(data = as.character(data)) %>%
  dplyr::select(celltype, coef, data, correlation) %>%
  pivot_wider(names_from = data, values_from = correlation)

# Step 2: rename columns to safe names
cor_wide <- cor_wide %>%
  dplyr::rename(
    biolord = `biolord`,
    PerturbNet = `PerturbNet`,
    ex_vivo = `ex vivo`   # note the backticks to reference the space
  )

# Step 3: determine which method has the max correlation
wins <- cor_wide %>%
  mutate(
    winner = case_when(
      (!is.na(biolord) & biolord >= coalesce(PerturbNet, -Inf)) &
        (!is.na(biolord) & biolord >= coalesce(ex_vivo, -Inf)) ~ "biolord",
      
      (!is.na(PerturbNet) & PerturbNet >= coalesce(biolord, -Inf)) &
        (!is.na(PerturbNet) & PerturbNet >= coalesce(ex_vivo, -Inf)) ~ "PerturbNet",
      
      (!is.na(ex_vivo) & ex_vivo >= coalesce(biolord, -Inf)) &
        (!is.na(ex_vivo) & ex_vivo >= coalesce(PerturbNet, -Inf)) ~ "ex_vivo",
      
      TRUE ~ NA_character_  # tie or missing
    )
  )


# Count wins overall
win_counts <- wins %>%
  select(winner) %>%          # drop all other columns
  filter(!is.na(winner)) %>%
  group_by(winner) %>%
  summarise(n_wins = n()) %>%
  arrange(desc(n_wins))

win_counts

## Count wins per celltype
win_counts_celltype <- wins %>%
  select(celltype, winner) %>%   # drop other columns
  filter(!is.na(winner)) %>%
  group_by(celltype, winner) %>%
  summarise(n_wins = n()) %>%
  arrange(celltype, desc(n_wins))

win_counts_celltype
library(ggplot2)

ggplot(win_counts_celltype, aes(x = celltype, y = n_wins, fill = winner)) +
  geom_col(position = "dodge") +              # side-by-side bars per cell type
  labs(
    x = "Cell Type",
    y = "Number of Wins",
    fill = "Method",
    title = "Number of KO × Celltype Cases Won by Each Method"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # rotate x labels
    legend.position = "top"
  )+
  optimized_theme_fig()
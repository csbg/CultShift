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
base <- "Figure5"
basedir <- dirout("Figure5")


###############
InDir1 <- dirout("Ag_ScRNA_23_prediction_correlations/")
#Fig5.1
correlation_deg_flagged <- read_rds(InDir1("correlation_to_invivo_exvivo_vs_pred.rds"))
correlation_deg_flagged <- correlation_deg_flagged %>%
  filter(!(celltype %in% c("Mono", "MEP.early")))
  
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
    high = "#D0154E",
    name = expression("Pearson's\ncorrelation")
  ) +
  scale_size_continuous(
    range = c(0,2),
    limits = c(0,3),
    breaks = c(1,2,3),
    name = expression(atop("No. of genes", log[10](n))))+
  labs(x = "KOs",
       y = "Cell type",
       title = "Correlation of predicted versus ex vivo KO effect\nto actual (in vivo) KO effect") +
  optimized_theme_fig()+
  theme(
    legend.position = "right",
    legend.justification = "right",
    panel.spacing = unit(0.1,"cm"),
    strip.text.x = element_text(angle = 90, hjust = 0)
  )
ggsave(basedir("Fig5A.pdf"), w = 18, h = 5, units = "cm")

#Fig5B

corr_pred <- read_rds(InDir1("Prediction-actual_cor.rds"))
logFC <- read_rds(InDir1("LogFC_invivo_pred.rds"))

# ggplot(corr_pred,
#        aes(x = coef, y = correlation))+
#   geom_bar(stat = "identity", fill = "#4C889C") +
#   facet_grid(rows = vars(celltype))+optimized_theme_fig()
# 
# 
# plot_logFC_comparison <- function(selected_celltype, logFC) {
#   selected_coefs <- c("in.vivoKmt2d", "in.vivoBrd9")
#   
#   filtered_data <- logFC %>%
#     filter(coef %in% selected_coefs, celltype == selected_celltype) %>%
#     mutate(coef = gsub("in.vivo", "", coef))  # Clean up labels
#   
#   p <- ggplot(filtered_data, aes(x = logFC_pred, y = logFC)) +
#     geom_hex(bins = 40, aes(fill = after_stat(count))) +
#     scale_fill_gradient(low = "#d0e1f2", high = "#08306b") +  # Blue gradient
#     geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c") +
#    
#     facet_wrap(~coef) +
#     labs(
#       title = paste("Genotype Effect in", selected_celltype),
#       x = "Predicted logFC",
#       y = "Actual logFC",
#       fill = "Cell Count",
#       color = "Genotype"
#     ) +
#     optimized_theme_fig()
#   
#   ggsave(basedir(paste0("Fig5B.pdf")), 
#          plot = p, width = 14, height = 5, units = "cm")
# }
# 
# plot_logFC_comparison("Eo.Ba", logFC)
# 

correlation_wide <- correlation_deg_flagged %>%
  pivot_wider(
    names_from = data,
    values_from = correlation
  ) %>%
  dplyr::rename(
    cor_exvivo = `ex.vivo`,
    cor_pred = prediction
  )
correlation_wide <- correlation_wide %>%
  filter(!celltype %in% c("Mono", "MEP.early"))



# Generate a distinct color palette using base R tools (no system deps)
genotypes <- unique(correlation_wide$genotype)
color_list <- hue_pal()(length(genotypes))
names(color_list) <- genotypes

# Create a distinct shape list (max ~25 standard ggplot shapes)
shape_list <- seq_along(genotypes)
names(shape_list) <- genotypes

# Final plot
ggplot(correlation_wide, aes(x = cor_exvivo, y = cor_pred)) +
  geom_point(aes(
    size = pmin(3, log10(num_degs_act)),
    color = genotype,
    shape = genotype
  ),
  stroke = 0.5,
  alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "red") +
  scale_size_continuous(
    range = c(0, 2),
    limits = c(0, 3),
    breaks = c(1, 2, 3),
    name = expression(atop("No. of DEGs", log[10](n)))
  ) +
  scale_color_manual(values = color_list,) +
  scale_shape_manual(values = shape_list) +
  facet_wrap(~celltype, nrow = 2) +
  labs(
    x = "Correlation to Ex Vivo",
    y = "Correlation to Predicted",
    color = "Genotype",
    shape = "Genotype"
  ) +
  optimized_theme_fig() +
  theme(
    legend.position =  "right",
    legend.direction = "horizontal",
    legend.spacing.x = unit(0.1,"cm"),
    strip.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(basedir("Fig5B_scatter_labeled.pdf"), width = 13, height = 5.5, units = "cm")
ggsave(basedir("Fig5B_scatter_labeled_no_legend.pdf"), width = 15, height = 3.5, units = "cm")
############Fig5C-------

pred_act <- read_rds(InDir1("LogFC_invivo_pred.rds"))
pred_act$coef <- gsub("in.vivo","",pred_act$coef)
colnames(pred_act)
unique(pred_act$coef)
data <-pred_act%>%
  filter(coef %in% c("Brd9","Kmt2d"))%>%
  filter(celltype == "Eo.Ba")
ggplot(data, aes(y = logFC_pred, x = logFC)) +
  geom_hex(bins = 40, aes(fill = after_stat(count))) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") +  # Blue gradient
  geom_smooth(method = "lm", linetype = "dashed",
              se = FALSE, size = 0.8, color ="#e41a1c")+  # Linear regression line
  
  facet_grid(cols = vars(coef))+
  labs(
    title = "logFC predicted versus actual (In vivo)",
    x = "logFC (predicted)",
    y = "logFC (actual)",
    fill = "Gene Count"
  )+ 
  optimized_theme_fig() +
  theme(
    plot.title = element_text(hjust = 1, vjust = 1)
  )
ggsave(basedir("Fig5C_scatter.pdf"), width = 6.5, height = 4, units = "cm")


#Figure2
#Invivo vs exvivo 
#Aarathy
#24-09-2024
###############################################################################
#load libraries and functions--------------------------------------------------
###############################################################################
source("src/00_init.R")
source("src/Ag_Optimized_theme.R")
library(tidyverse)
library(enrichR)
library(patchwork)
#library(purrr)
#library(gridExtra)
library(pheatmap)
########################
#directories ------
########################
base<-"Figure2"
basedir<-dirout("Figure2")
###############################################################################
#Fig2.1 Correlation between logFC ---------------------------------------------
###############################################################################
Indir1<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/")
correlation_matrix_data <- read_rds(Indir1("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(Indir1("DEGs_per_tissue.rds"))
#head(correlation_matrix_data)
# replacing NA with 0s
# I need to keep the NA and display as NAs. But during the clustering, NA values are a problem
correlation_matrix_data <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)

# Plot the heatmap

##########
#plot2.1.1
##########
# Load required libraries
#pdf(basedir("Fig2.1.pdf"), w = 16, h = 8)
pheatmap(correlation_matrix_data,
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         clustering_method = "ward.D2",
         display_numbers = TRUE, 
         cellwidth = 20,
         cellheight = 20,
         number_format = "%.1f", 
         main = "Correlation of logFC between ex vivo and in vivo by Genotype",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = T,
         na_col = "darkgrey")  # Display NA values in grey

#dev.off()
#alternative
# color_palette <- colorRamp2(breaks = seq(-1, 1, length = 3), 
#                             colors = c("blue", "white", "red"))
# 
# Create the heatmap
# heatmap_obj <- Heatmap(
#   correlation_matrix_data,
#   name = "Correlation",  # Name
#   col = color_palette, 
#   cluster_rows = TRUE,  
#   cluster_columns = TRUE, 
#   clustering_method_rows = "ward.D2",
#   clustering_method_columns = "ward.D2",
#   show_row_names = TRUE, 
#   show_column_names = TRUE, 
#   row_names_side = "left",  
#   column_names_side = "bottom",  
#   na_col = "darkgrey",  # Color for NA values? how to include NA
#   heatmap_legend_param = list(title = "Correlation Coefficient",
#                               at = seq(-1, 1, by = 0.2)),  # Legend parameters
#   row_title = "Genes",  # Title for row labels
#   column_title = "Samples"  # Title for column labels
# )
# 
# # # Draw the heatmap
#  draw(heatmap_obj)
# Draw the heatmap

hc_cols <- hclust(dist(t(correlation_matrix_data)), method = "ward.D2")
column_order <- colnames(correlation_matrix_data)[hc_cols$order]
deg_plot_data$genotype <- factor(deg_plot_data$genotype, levels = column_order)

deg_plot_data$deg_category <- factor(deg_plot_data$deg_category,levels = c("Below 10",
                                                                           "10-50",
                                                                           "50-100",
                                                                           "Above 100",
                                                                           "Above 1000",
                                                                           "NA"))

##########
#plot2.1.2
##########
# Create the plot
deg_dot_plot <- ggplot(deg_plot_data, aes(x = genotype, y = condition)) +
  geom_point(aes(color = deg_category), size = 5) +  # Fix size and vary shape for condition
  scale_color_manual(values = c("Below 10" = "#4B79A5", 
                                "10-50" = "#B7D4E9", 
                                "50-100" = "#FFC2B5", 
                                "Above 100" = "#ED564E",
                                "Above 1000"="#B7203C",
                                "NA" = "grey")) +  # Discrete color gradient+
  facet_grid(rows = vars(celltype))+
  scale_shape_manual(values = c("Ex Vivo" = 16, "In Vivo" = 17)) +  # Different shapes for Ex Vivo/In Vivo
  labs(
    title = "Categorized DEGs by Genotype, Celltype, and Condition",
    x = "Genotype", 
    y = "Celltype",
    color = "Number of DEGs",
    shape = "Condition"
  ) +
  theme_minimal() +
  optimized_theme()
ggsave(basedir("deg_dot_plot.pdf"),plot=deg_dot_plot)
##################

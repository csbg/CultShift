# Load required libraries
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(ComplexHeatmap)

###############################################################################
# Define directories
base <- "Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/"
Indir <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
basedir <- dirout(base)
###############################################################################

# Read in data
limmaRes <- read_rds(Indir("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))

# Prepare ex.vivo data
ex_vivo_data <- limmaRes %>%
  filter(coef != "ex.vivo") %>%
  filter(str_detect(coef, "^ex.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "ex.vivo\\.", "")) %>%
  mutate(genotype = gsub("ex.vivo", "", coef)) %>%
  select("ensg", "logFC", "celltype", "genotype", "adj.P.Val")

# Prepare in.vivo data
in_vivo_data <- limmaRes %>%
  filter(str_detect(coef, "^in.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "in.vivo\\.", "")) %>%
  mutate(genotype = gsub("in.vivo", "", coef)) %>%
  select("ensg", "logFC", "celltype", "genotype", "adj.P.Val")

# Combine datasets
merged_data <- ex_vivo_data %>%
  inner_join(in_vivo_data, by = c("ensg", "celltype", "genotype"), 
             suffix = c("_ex.vivo", "_in.vivo"))

# Calculate the correlation for each genotype within each cell type, excluding NA
correlation_results <- merged_data %>%
  group_by(celltype, genotype) %>%
  summarise(correlation = cor(logFC_ex.vivo, logFC_in.vivo, use = "complete.obs")) %>%
  ungroup()

# Reshape the data for the heatmap
correlation_matrix <- correlation_results %>%
  pivot_wider(names_from = genotype, values_from = correlation)

# Convert to a matrix
correlation_matrix_data <- as.matrix(correlation_matrix[,-1])
rownames(correlation_matrix_data) <- correlation_matrix$celltype
write_rds(correlation_matrix_data, basedir("correlation_ex.vivo_vs_in.vivo.rds"))

# Check which columns contain NA values
cols_with_na <- apply(correlation_matrix_data, 2, function(col) any(is.na(col)))

# Exclude columns with NA values
filtered_correlation_matrix_data <- correlation_matrix_data[, !cols_with_na]

# Compute distance matrix (using Euclidean distance)
dist_matrix <- dist(t(filtered_correlation_matrix_data), method = "euclidean")

# Perform hierarchical clustering
col_clustering <- hclust(dist_matrix, method = "ward.D")

# Filter for differentially expressed genes (DEGs)
deg_ex_vivo <- ex_vivo_data %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype, celltype) %>%
  summarise(num_degs_ex_vivo = n())

deg_in_vivo <- in_vivo_data %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype, celltype) %>%
  summarise(num_degs_in_vivo = n())

# Merge the DEG counts
deg_counts <- full_join(deg_ex_vivo, deg_in_vivo, by = c("genotype", "celltype"))

# Create a dot plot for DEGs below the heatmap
deg_plot_data <- deg_counts %>%
  pivot_longer(cols = starts_with("num_degs"), names_to = "condition", values_to = "num_degs") %>%
  mutate(condition = ifelse(condition == "num_degs_ex_vivo", "Ex Vivo", "In Vivo"))

# Plot the heatmap and DEG counts together
library(ggplot2)
library(patchwork)  # To combine heatmap and dot plots

# Plotting DEGs as dot plots
deg_dot_plot <- ggplot(deg_plot_data, aes(x = genotype, y = celltype)) +
  geom_point(aes(size = num_degs, color = condition)) +
  scale_size_continuous(range = c(2, 8)) +  # Adjust size for better visibility
  labs(size = "Number of DEGs", color = "Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save heatmap and DEG plot as PDF
pdf(basedir("ex.vivo_vs_in.vivo_ko_effects_with_deg_dotplots.pdf"), width = 16, height = 10)

# Draw heatmap and combine with dot plots
pheatmap(
  filtered_correlation_matrix_data,
  cluster_rows = TRUE, 
  cluster_cols = col_clustering,  # Use the correct clustering object
  display_numbers = TRUE, 
  cellwidth = 20,
  cellheight = 20,
  number_format = "%.1f", 
  main = "Correlation of logFC between ex vivo and in vivo by Genotype",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  show_rownames = TRUE
) + deg_dot_plot  # Combine heatmap with dot plots

dev.off()

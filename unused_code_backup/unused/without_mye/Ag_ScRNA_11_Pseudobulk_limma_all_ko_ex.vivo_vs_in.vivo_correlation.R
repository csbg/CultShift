source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
###############################################################################
base <- "Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/"
Indir <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
basedir <- dirout(base)
###############################################################################
limmaRes <- read_rds(Indir("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))
in.vivo_degs <- limmaRes %>%
  filter(coef %in% grep("in.vivo", limmaRes$coef, value = T))%>%
  filter(group != "n.s")

summary_df <-in.vivo_degs %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < 0.05 & logFC > 1),
    Downregulated = sum(adj.P.Val < 0.05 & logFC < -1)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")


# Prepare ex.vivo data
ex_vivo_data <- limmaRes %>%
  filter(coef != "ex.vivo") %>%
  filter(str_detect(coef, "^ex.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "ex.vivo\\.", "")) %>%
  mutate(genotype = gsub("ex.vivo", "", coef)) %>%
  dplyr::select("ensg", "logFC", "celltype", "genotype", "adj.P.Val")

# Prepare in.vivo data
in_vivo_data <- limmaRes %>%
  filter(str_detect(coef, "^in.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "in.vivo\\.", "")) %>%
  mutate(genotype = gsub("in.vivo", "", coef)) %>%
  dplyr::select("ensg", "logFC", "celltype", "genotype", "adj.P.Val")
unique(limmaRes$coef)
# Combine datasets
merged_data <- ex_vivo_data %>%
  inner_join(in_vivo_data, by = c("ensg", "celltype", "genotype"), 
             suffix = c("_ex.vivo", "_in.vivo"))
saveRDS(merged_data,basedir("in.vivo_ex.vivo_logFC.rds"))
# Calculate the correlation for each genotype within each cell type, excluding NA
correlation_results <- merged_data %>%
  group_by(celltype, genotype) %>%
  summarise(correlation = cor(logFC_ex.vivo, logFC_in.vivo, use = "complete.obs")) %>%
  ungroup()

unique(correlation_results$genotype)

# Reshape the data for the heatmap
correlation_matrix <- correlation_results %>%
  pivot_wider(names_from = genotype, values_from = correlation)

# Convert to a matrix
correlation_matrix_data <- as.matrix(correlation_matrix[,-1])
rownames(correlation_matrix_data) <- correlation_matrix$celltype
write_rds(correlation_matrix_data,basedir("correlation_ex.vivo_vs_in.vivo.rds"))
#
library(tidyverse)

# Assume correlation_matrix_data is your wide data.frame or matrix
# First, convert to long format
corr_long <- as.data.frame(correlation_matrix_data) %>%
  rownames_to_column("genotype") %>%
  pivot_longer(-genotype, names_to = "gene", values_to = "correlation")

# Remove NAs (optional)
corr_long <- corr_long %>% drop_na(correlation)

# Order the genes within each genotype based on correlation (ascending or descending)
corr_long <- corr_long %>%
  group_by(genotype) %>%
  arrange(correlation, .by_group = TRUE) %>%
  mutate(gene = factor(gene, levels = unique(gene))) %>%
  ungroup()
write_rds(corr_long,basedir("correlation_ex.vivo_vs_in.vivo_ordered_per_celltype.rds"))
# Plot
ggplot(corr_long, aes(x = gene, y = correlation, fill = correlation)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ genotype, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation of logFC_ex.vivo vs logFC_in.vivo",
       x = "Gene", y = "Correlation")

# Check which columns contain NA values
cols_with_na <- apply(correlation_matrix_data, 2, function(col) any(is.na(col)))

# Exclude columns with NA values
filtered_correlation_matrix_data <- correlation_matrix_data[, !cols_with_na]
# Compute distance matrix (using Euclidean distance)
dist_matrix <- dist(t(filtered_correlation_matrix_data), method = "euclidean")

# Perform hierarchical clustering
col_clustering <- hclust(dist_matrix, method = "ward.D")
# Verify dimensions of the matrix and clustering object
cat("Number of columns in filtered matrix:", ncol(filtered_correlation_matrix_data), "\n")
cat("Length of clustering object:", length(col_clustering$labels), "\n")

# Plot the heatmap using the filtered matrix and the correct clustering object
pdf(basedir("ex.vivo_vs_invivo_ko_effects_filtered_no_na_cols.pdf"), w = 16, h = 8)

pheatmap(filtered_correlation_matrix_data,
         cluster_rows = TRUE, 
         cluster_cols = col_clustering,  # Use the correct clustering object
         display_numbers = TRUE, 
         cellwidth = 20,
         cellheight = 20,
         number_format = "%.1f", 
         main = "Correlation of logFC between ex vivo and in vivo by Genotype",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = TRUE)

dev.off()
# Convert correlation matrix to a tidy format
correlation_long <- as.data.frame(filtered_correlation_matrix_data) %>%
  rownames_to_column("celltype") %>%
  pivot_longer(cols = -celltype, names_to = "genotype", values_to = "correlation")

# Check the structure
head(correlation_long)


# Assign row names directly from the filtered data
################################################################################
#Number of DEGs-----------------------------------------------------------------
################################################################################
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
deg_plot_data <- deg_counts %>%
  pivot_longer(cols = starts_with("num_degs"), names_to = "condition", values_to = "num_degs") %>%
  mutate(condition = ifelse(condition == "num_degs_ex_vivo", "Ex Vivo", "In Vivo")) %>%
  left_join(correlation_long, by = c("genotype", "celltype"))
# Plotting DEGs as dot plots
# Categorize the number of DEGs
deg_plot_data <- deg_plot_data %>%
  mutate(deg_category = case_when(
    num_degs < 10 ~ "Below 10",
    num_degs >= 10 & num_degs <= 50 ~ "10-50",
    num_degs > 50 & num_degs <= 100 ~ "50-100",
    num_degs > 100 & num_degs <= 1000 ~ "100-1000",
    num_degs > 1000 ~ "Above 1000",
    TRUE ~ "NA"
  ))
write_rds(deg_plot_data,basedir("DEGs_per_tissue.rds"))


##########

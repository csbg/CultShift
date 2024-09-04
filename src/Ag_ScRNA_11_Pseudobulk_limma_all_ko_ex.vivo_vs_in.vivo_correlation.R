source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
###############################################################################
base<-"Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/"
Indir<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
out<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation")
###############################################################################
limmaRes<-read_rds(Indir("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))

# Prepare ex.vivo data
ex_vivo_data <- limmaRes %>%
  filter(coef != "ex.vivo")%>%
  filter(str_detect(coef, "^ex.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "ex.vivo\\.", "")) %>%
  mutate(genotype = gsub("ex.vivo","",coef))%>%
  select("ensg","logFC","celltype","genotype","adj.P.Val")


# For in.vivo data, no prefix adjustment is necessary
in_vivo_data <- limmaRes %>%
  filter(str_detect(coef, "^in.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "in.vivo\\.", "")) %>%
  mutate(genotype = gsub("in.vivo","",coef))%>%
  select("ensg","logFC","celltype","genotype","adj.P.Val")

# Combine datasets on cleaned coefficient names

# Step 1: Merge the ex_vivo_data and in_vivo_data
merged_data <- ex_vivo_data %>%
  inner_join(in_vivo_data, by = c("ensg", "celltype", "genotype"), 
             suffix = c("_ex.vivo", "_in.vivo")) 

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
# Remove the first column (celltype)

# Assign row names directly from the filtered data

# Check dimensions

# Replace NA, NaN, and Inf with a specific value (e.g., 0) or remove rows/columns
# Here, we replace them with 0
correlation_matrix_data <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)



# Plot the heatmap
pdf(out("ex.vivo_vs_invivo_ko_effects1.pdf"), w = 16, h = 8)

pheatmap(correlation_matrix_data,
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = TRUE, 
         cellwidth = 20,
         cellheight = 20,
         number_format = "%.1f", 
         main = "Correlation of logFC between ex vivo and in vivo by Genotype",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = T,
         na_col = "darkgrey")  # Display NA values in grey

dev.off()

#############################
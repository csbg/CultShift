########################################################################
source("src/00_init.R")
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(tidyverse)
source("src/Ag_Optimized_theme.R")

# Input/Output directories-----------
InDir <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
InDir4 <- dirout("Figure1")
base <- "Ag_ScRNA_19_invivo_exvivo_zscore/"
basedir <- dirout("Ag_ScRNA_19_invivo_exvivo_zscore/")

########################################################################
#load data and clean metadata
########################################################################
#metadata from in-vivo ex-vivo
meta <- fread(InDir2("meta_cleaned.tsv"))
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   
meta <- meta %>%
  filter(genotype == "NTC")# Assign first column as row names
meta <- meta[, -1, drop = FALSE] 
meta$sample1 <- rownames(meta)
meta <- meta[,c("sample","celltype","tissue","sample1")]
#metadata fromizzo

# filtering steps
celltypes_to_exclude <- c("CLP",  "EBMP", "unclear","T-cell","MEP","Imm. B-cell")
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                      "Hba-a2","Pim1","Fabp5","Fdps","Cd9")

meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]
# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)\\-]", ".", rownames(meta))



meta <- meta %>%
  filter(!(celltype %in% c("B-cell","Ery","Neu","T-cell-Cd3d+","E/B","Imm. B-cell")))
rownames(meta) <- gsub("[\\ \\(\\)\\-]", ".", rownames(meta))
meta <- meta %>% filter(!grepl("NA",rownames(meta)))
meta <- meta[!grepl("T.cell",rownames(meta)),]



dataVoom <- read_rds(InDir2("dataVoom_perCTex.vivovsin.vivo.rds"))



# modify dataVoom
longer_dataVoom <-  dataVoom$E %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample1",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(meta, by ="sample1")%>%
  mutate(tissue_celltype =paste0(tissue,"_",celltype))
##################
#without replicate
ex_vivo_data <- longer_dataVoom %>% filter(tissue == "ex.vivo")
in_vivo_data <- longer_dataVoom %>% filter(tissue == "in.vivo")

# Step 2: Pivot ex vivo and in vivo data to wide format
# Pivot ex vivo data to wide format
ex_vivo_wide <- ex_vivo_data %>%
  select(genes, tissue_celltype, Expression) %>%
  pivot_wider(
    names_from = tissue_celltype, 
    values_from = Expression,
    values_fn = list(Expression = mean)
  )

# Pivot in vivo data to wide format
in_vivo_wide <- in_vivo_data %>%
  select(genes, tissue_celltype, Expression) %>%
  pivot_wider(
    names_from = tissue_celltype, 
    values_from = Expression,
    values_fn = list(Expression = mean)
  )

# Step 3: Merge ex vivo and in vivo data on the genes column
combined_data <- inner_join(ex_vivo_wide, in_vivo_wide, by = "genes")

# Remove the gene column for correlation calculation
combined_data_for_correlation <- combined_data %>% select(-genes)

# Step 4: Calculate the correlation matrix
cor_matrix <- cor(combined_data_for_correlation, method = "pearson")

# Step 5: Get the exact cell type order from the ex_vivo data
# Ex vivo cell type order
# Filter for ex vivo and in vivo samples
ex_vivo_data <- longer_dataVoom %>% filter(tissue == "ex.vivo")
in_vivo_data <- longer_dataVoom %>% filter(tissue == "in.vivo")



# Pivot ex vivo data to wide format (keeping individual replicates)
ex_vivo_wide <- ex_vivo_data %>%
  select(genes, sample1, Expression) %>%
  pivot_wider(
    names_from = sample1, 
    values_from = Expression
  )

# Pivot in vivo data to wide format (keeping individual replicates)
in_vivo_wide <- in_vivo_data %>%
  select(genes, sample1, Expression) %>%
  pivot_wider(
    names_from = sample1, 
    values_from = Expression
  )

# Merge ex vivo and in vivo data on the genes column
combined_data <- inner_join(ex_vivo_wide, in_vivo_wide, by = "genes")

# Remove the gene column for correlation calculation (we'll use it for labeling later)
combined_data_for_correlation <- combined_data %>% select(-genes)
exvivo <- grep("ex.vivo", rownames(cor_matrix), value = T)
invivo <- gsub("ex.vivo","in.vivo",exvivo)
# Compute the correlation matrix
cor_matrix <- cor(combined_data_for_correlation, method = "pearson")
cor_matrix_1 <- cor_matrix[exvivo,invivo]
# Check the correlation matrix (optional)

# Load ComplexHeatmap library
cor_matrix_1 %>% write_rds(basedir("Correlation_ex.vivo_in.vivo.rds"))

pdf(basedir("correlation_celltypes.pdf"), h=)
# Plot the heatmap with ex vivo samples on x-axis and in vivo samples on y-axis
heatmap_plot <- Heatmap(
  cor_matrix_1, 
  name = "Correlation",
  row_title = "Ex Vivo Samples", 
  column_title = "In Vivo Samples",
  show_row_names = TRUE,  # Show row names (in vivo celltypes)
  show_column_names = TRUE,  # Show column names (ex vivo celltypes)
  row_names_gp = gpar(fontsize = 8),  # Adjust font size for row names (in vivo)
  column_names_gp = gpar(fontsize = 8),
  cluster_columns = F,
  cluster_rows = F,
  row_order = exvivo,
  column_order = invivo,# Adjust font size for column names (ex vivo)
  show_row_dend = FALSE,  # No row clustering
  show_column_dend = FALSE,  # No column clustering
  #clustering_distance_rows = "none",  # Disable row clustering
  #clustering_distance_columns = "none",  # Disable column clustering
  color = colorRampPalette(c("blue", "white", "red"))(100)  # Correlation color scale
)
# Draw the heatmap
draw(heatmap_plot)
dev.off()

###################
#with replicate
# Filter for ex vivo and in vivo samples
ex_vivo_data <- longer_dataVoom %>% filter(tissue == "ex.vivo")
in_vivo_data <- longer_dataVoom %>% filter(tissue == "in.vivo")



# Pivot ex vivo data to wide format (keeping individual replicates)
ex_vivo_wide <- ex_vivo_data %>%
  select(genes, sample1, Expression) %>%
  pivot_wider(
    names_from = sample1, 
    values_from = Expression
  )

# Pivot in vivo data to wide format (keeping individual replicates)
in_vivo_wide <- in_vivo_data %>%
  select(genes, sample1, Expression) %>%
  pivot_wider(
    names_from = sample1, 
    values_from = Expression
  )

# Merge ex vivo and in vivo data on the genes column
combined_data <- inner_join(ex_vivo_wide, in_vivo_wide, by = "genes")

# Remove the gene column for correlation calculation (we'll use it for labeling later)
combined_data_for_correlation <- combined_data %>% select(-genes)
exvivo <- grep("ex.vivo", rownames(cor_matrix), value = T)
invivo <- gsub("ex.vivo","in.vivo",exvivo)
# Compute the correlation matrix
cor_matrix <- cor(combined_data_for_correlation, method = "pearson")
cor_matrix_1 <- cor_matrix[exvivo,invivo]
# Check the correlation matrix (optional)
head(cor_matrix)
# Load ComplexHeatmap library


# Plot the heatmap with ex vivo samples on x-axis and in vivo samples on y-axis
heatmap_plot <- Heatmap(
  cor_matrix_1, 
  name = "Correlation",
  row_title = "Ex Vivo Samples", 
  column_title = "In Vivo Samples",
  show_row_names = TRUE,  # Show row names (in vivo celltypes)
  show_column_names = TRUE,  # Show column names (ex vivo celltypes)
  row_names_gp = gpar(fontsize = 8),  # Adjust font size for row names (in vivo)
  column_names_gp = gpar(fontsize = 8),
  cluster_columns = F,
  cluster_rows = F,
  row_order = exvivo,
  column_order = invivo,# Adjust font size for column names (ex vivo)
  show_row_dend = FALSE,  # No row clustering
  show_column_dend = FALSE,  # No column clustering
  #clustering_distance_rows = "none",  # Disable row clustering
  #clustering_distance_columns = "none",  # Disable column clustering
  color = colorRampPalette(c("blue", "white", "red"))(100)  # Correlation color scale
)
# Draw the heatmap
draw(heatmap_plot)


# # wide format
# wide_data <- longer_dataVoom %>%
#   select(genes, tissue_celltype_replicate, Expression) %>%
#   pivot_wider(names_from = tissue_celltype_replicate, values_from = Expression) %>%
#   column_to_rownames("genes")

# small_test <- longer_dataVoom %>%
#   filter(genes %in% sample(unique(genes), 100)) %>%
#   pivot_wider(names_from = tissue_celltype_replicate, values_from = Expression)
# 
# t <- dataVoom$E
# 
# # mean per tissue_celltype per gene
# gene_mean <- longer_dataVoom %>%
#   group_by(tissue_celltype,tissue,celltype, genes) %>%
#   summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

# ex-in-------------
# based on expression
# data to wide format
wide_data <- gene_mean %>%
  select(genes, tissue_celltype, mean_expr) %>%
  pivot_wider(names_from = tissue_celltype, values_from = mean_expr) %>%
  column_to_rownames("genes")


# Select ex.vivo and in.vivo columns
ex_vivo_columns <- grep("^ex.vivo", colnames(wide_data), value = TRUE)
in_vivo_columns <- grep("^in.vivo", colnames(wide_data), value = TRUE)
# in_vivo_columns <- grep("inVivo", colnames(dataVoom$E), value = TRUE)
# ex_vivo_columns <- colnames(dataVoom[,!colnames(dataVoom$E) %in% in_vivo_columns])

# Subset the data for ex.vivo and in.vivo
ex_vivo_data <- wide_data[, ex_vivo_columns]
in_vivo_data <- wide_data[, in_vivo_columns]
# Compute correlation matrix
# Calculate the correlation between ex.vivo and in.vivo data
cor_matrix <- cor(ex_vivo_data, in_vivo_data)
cor_matrix1 <- cor(dataVoom$E)
# Load required libraries
library(ggplot2)
library(reshape2)

# Compute correlation matrix
cor_matrix1 <- cor(dataVoom$E, method = "pearson")

# Convert correlation matrix to long format
cor_long <- melt(cor_matrix1)

# Remove self-correlations and duplicates
cor_long <- cor_long[cor_long$Var1 != cor_long$Var2, ]  # Remove diagonal
cor_long <- cor_long[!duplicated(t(apply(cor_long, 1, sort))), ]  # Remove duplicates

# Plot scatter plot
ggplot(cor_long, aes(x = Var1, y = Var2, color = value)) +
  geom_point(size = 4) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Scatter Plot of Correlation Matrix",
       x = "Tissue-Cell Type 1",
       y = "Tissue-Cell Type 2",
       color = "Correlation")+optimized_theme()


ggsave(basedir("corlong.pdf"))


diag(cor_matrix1) <- NA

pdf(file = basedir("cor_heatmap_iv_ex_dataVoom.pdf")
)
Heatmap(cor_matrix1 )
dev.off()
# Plot the heatmap
Heatmap(cor_matrix
         )   
#
library(ggplot2)

# Ensure that ex.vivo and in.vivo samples are aligned
ex_vivo_long <- wide_data[, ex_vivo_columns] %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  pivot_longer(cols = -genes, names_to = "sample", values_to = "ex_vivo_expr")

in_vivo_long <- wide_data[, in_vivo_columns] %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  pivot_longer(cols = -genes, names_to = "sample", values_to = "in_vivo_expr")

# Merge by genes and sample (to pair ex.vivo and in.vivo for each sample)
merged_data <- inner_join(ex_vivo_long, in_vivo_long, by = c("genes"))

# Scatter plot
ggplot(merged_data, aes(x = ex_vivo_expr, y = in_vivo_expr)) +
  geom_point(alpha = 0.5, color = "blue") +  # Each dot is a sample
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +  # Regression line
  theme_minimal() +
  labs(
    title = "Expression Correlation: ex.vivo vs. in.vivo",
    x = "Ex Vivo Expression",
    y = "In Vivo Expression"
  )+optimized_theme()
ggsave(basedir("corlong.pdf"))
########################
#by taking mean across similar samples
# Step 1: Calculate Mean per Gene per (Sample, Cell Type, Tissue)
gene_mean <- longer_dataVoom %>%
  group_by(tissue_celltype,tissue,celltype, genes) %>%
  summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop")

# ex-in-------------
# based on expression
# data to wide format
wide_data <- gene_mean %>%
  select(genes, tissue_celltype, mean_expr) %>%
  pivot_wider(names_from = tissue_celltype, values_from = mean_expr) %>%
  column_to_rownames("genes")

# Compute correlation matrix
cor_matrix <- cor(wide_data, method = "pearson")

diag(cor_matrix) <- NA
# Convert correlation to distance (1 - correlation)
dist_matrix <- as.dist(1 - cor_matrix)

# Perform hierarchical clustering
row_clust <- hclust(dist_matrix, method = "ward.D2")

pdf(file = basedir("Onlycorheatmap_ex_in_celltype(expression).pdf")
)
Heatmap(cor_matrix,cluster_rows = row_clust, cluster_columns = row_clust)
dev.off()

#based on zscore within celltype_tissue across genes

zscore_tissue_celltype <- gene_mean %>%
  group_by(tissue_celltype) %>%  # Group by tissue_celltype
  mutate(
    mean_tissue = mean(mean_expr, na.rm = TRUE),  # Mean expression for the group
    sd_tissue = sd(mean_expr, na.rm = TRUE),      # Standard deviation for the group
    zscore = (mean_expr - mean_tissue) / sd_tissue # Z-score for each gene
  ) %>%
  ungroup()


wide_data <- zscore_tissue_celltype %>%
  select(genes, tissue_celltype,zscore) %>%
  pivot_wider(names_from = tissue_celltype, values_from = zscore) %>%
  column_to_rownames("genes")

# Compute correlation matrix
cor_matrix <- cor(wide_data, method = "pearson")

diag(cor_matrix) <- NA
# Convert correlation to distance (1 - correlation)
dist_matrix <- as.dist(1 - cor_matrix)

# Perform hierarchical clustering
row_clust <- hclust(dist_matrix, method = "ward.D2")

pdf(file = basedir("corheatmap_ex_in_celltype(zscore_within_tissue_celltype_across_genes).pdf")
)
Heatmap(cor_matrix,cluster_rows = row_clust, cluster_columns = row_clust)
dev.off()
# based on zscore across sample(tissue_celltype)

filtered_data <- gene_mean %>%
  select(genes, tissue_celltype, mean_expr) %>%
  pivot_wider(names_from = tissue_celltype, values_from = mean_expr) %>%
  column_to_rownames("genes") %>%
  as.matrix()

# Normalize using Z-score transformation
filtered_data <- t(scale(t(filtered_data)))

# Compute Spearman correlation
tissue_celltype_correlation <- cor(filtered_data, use = "pairwise.complete.obs", method = "pearson")


diag(tissue_celltype_correlation) <- NA
pdf(file = basedir("Corr_Heatmap_zscore_across_tissue_celltype_In_Vivo_vs_Ex_Vivo_Heatmap.pdf"))
Heatmap(tissue_celltype_correlation)
dev.off()
# Save plot


############################

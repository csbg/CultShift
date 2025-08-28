#limma_NTC
#Comparison:ex.vivo-in.vivo

###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)

##################################################################################

inDir <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
InDir2 <- dirout("Metadata")
InDir3 <- dirout("SCRNA_06_01_Markers")
InDir4 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")

basedir <- dirout("Ag_ScRNA_18_zscore_plots_celltype_marker/")
source("src/Ag_Optimized_theme_fig.R")
##################################################################################
#load data
########################
#metadata
marker_genes <- fread(InDir2("/FIGS_02_DE_Genes.tsv"))

#
meta <- fread(InDir4("meta_cleaned.tsv")) 

# Convert to a dataframe (optional, depending on use case)
meta <- as.data.frame(meta)  

# Set the first column as row names
rownames(meta) <- meta[[1]]  

# Remove the first column from the data
meta <- meta[, -1, drop = FALSE] 

NTC_meta <- meta %>%
  filter(genotype == "NTC")

NTC_meta <- NTC_meta %>% rownames_to_column("Sample")
dataVoom <- read_rds(InDir4("dataVoom_perCTex.vivovsin.vivo.rds"))

# Load marker gene table
marker_genes <- fread(InDir2("/FIGS_02_DE_Genes.tsv"))

colnames(marker_genes) <- c("Gene","CELL")
# Z-score normalize within each sample
#zscore_normalized_counts <- t(scale(t(dataVoom$E)))  # Scale each gene (rows) within each sample

# Convert to data frame
df <- as.data.frame(dataVoom$E)
df$Gene <- rownames(df)

# Merge with metadata to get cell types
long <- df %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expr") %>%
  left_join(NTC_meta, by = "Sample")  # Ensure correct merging


# Aggregate by cell type
celltype_avg <- long %>%
  group_by(Gene, celltype,tissue) %>%
  summarise(Mean_Expr = mean(Expr, na.rm = TRUE), .groups = "drop")
#With B cells and erythrocytes
# Filter for marker genes
marker_results <- celltype_avg %>%
  filter(Gene %in% marker_genes$Gene)%>%
  left_join(marker_genes, by = "Gene")%>%
  filter(!(CELL == "Bcell"))

# Normalize Z-scores across each celltypewithin tissue for a gene
marker_results_scaled <- marker_results_scaled %>%
  group_by(tissue, Gene) %>%
  mutate(Scaled_Zscore = scale(Mean_Expr)) %>%
  ungroup()
# Normalize Z-scores within celltype and tissue across all genes
# marker_results_scaled <- marker_results %>%
#   group_by(celltype,tissue) %>%
#   mutate(Scaled_Expr = scale(Mean_Expr)) %>%
#   ungroup()

marker_results_scaled %>%
  write_rds(basedir("scaled_zscore_marker_genes.rds"))

ggplot(marker_results_scaled %>% filter(Gene %in% Gene[CELL == unique(CELL)]), 
       aes(x = Gene, y = tissue, fill = pmin(2,Scaled_Zscore)))+#, size = Mean_Zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4C889C",
                       mid = "white",
                       high = "#D0154E",
                       name= "scaled acr ct within tis"
  )+#,
  
  facet_grid(cols = vars(CELL), rows = vars(celltype), scales = "free", space = "free") +
  optimized_theme_fig()+
  theme(strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90, hjust = 0))

ggplot(marker_results_scaled %>% filter(Gene %in% Gene[CELL == unique(CELL)]), 
       aes(x = Gene, y = tissue, fill = pmin(2,Scaled_Expr)))+#, size = Mean_Zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4C889C",
                       mid = "white",
                       high = "#D0154E",
                       name= "scaled Expr"
  )+#,
  
  facet_grid(cols = vars(CELL), rows = vars(celltype), scales = "free", space = "free") +
  optimized_theme_fig()+
  theme(strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90, hjust = 0))
ggplot(marker_results_scaled %>% filter(Gene %in% Gene[CELL == unique(CELL)]), 
       aes(x = Gene, y = tissue, fill = pmin(2,Scaled_Expr)))+#, size = Mean_Zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4C889C",
                       mid = "white",
                       high = "#D0154E",
                       name= "zscore((zscore within ct)"
  )+#,
  
  facet_grid(cols = vars(CELL), rows = vars(celltype), scales = "free", space = "free") +
  optimized_theme_fig()+
  theme(strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90, hjust = 0))

ggsave(basedir("Celltype_marker_incl_Ery.pdf"), w=18, h=10,units = "cm")
#without B-cells and Ery
# Filter for marker genes
marker_results <- celltype_avg %>%
  filter(Gene %in% marker_genes$Gene)%>%
  left_join(marker_genes, by = "Gene")%>%
  filter(CELL != "Bcell")# %>%
  #filter(CELL != "Ery")

# Normalize Z-scores within each celltype
marker_results_scaled <- marker_results %>%
  group_by(tissue, Gene) %>%
  mutate(Scaled_Zscore = scale(Mean_Expr)) %>%
  ungroup()




ggplot(marker_results_scaled %>% filter(Gene %in% Gene[CELL == unique(CELL)]), 
       aes(x = Gene, y = tissue, color = pmin(pmax(Scaled_Zscore, -2), 2)))+#, size = Mean_Expr)) +
  geom_point() +
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name= "scaled Expr"
  )+#,
  scale_size_continuous(
    range = c(0, 1.5),
    limits = c(-3, 2.5),
    name=TeX("zscore within"))+
  facet_grid(cols = vars(CELL), rows = vars(celltype), scales = "free", space = "free") +
  optimized_theme_fig()+
  theme(strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90))
        
ggsave(basedir("Celltype_marker_incl_Ery.pdf"), w=18, h=9,units = "cm")
ggsave(basedir("Celltype_marker.pdf"))
#################################################################################
# #from merged markers
# merged_marker_list
# colnames(merged_marker_list) <- c("Gene","CELL","DB")
# 
# # Convert to data frame
# datavoom_df <- as.data.frame(dataVoom$E)
# datavoom_df$Gene <- rownames(datavoom_df)
# 
# #NTC_meta <- NTC_meta %>% rownames_to_column("Sample")
# # Merge with metadata to get cell types
# dataVoom_long <- datavoom_df %>%
#   pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expr") %>%
#   left_join(NTC_meta, by = "Sample")  # Ensure correct merging
# 
# 
# # Aggregate by cell type
# celltype_avg <- dataVoom_long %>%
#   group_by(Gene, celltype,tissue) %>%
#   summarise(Mean_Expr = mean(Expr, na.rm = TRUE), .groups = "drop")
# #With B cells and erythrocytes
# # Filter for marker genes
# marker_results <- celltype_avg %>%
#   filter(Gene %in% merged_marker_list$Gene)%>%
#   left_join(merged_marker_list, by = "Gene")#%>%
# #filter(CELL != "Bcell")
# 
# 
# # Plot 
# ggplot(marker_results %>% filter(Gene %in% Gene[CELL == unique(CELL)]), 
#        aes(x = Gene, y = celltype, color = Mean_Expr)) +
#   geom_tile() +
#   scale_color_gradient2(low = "blue", mid = "white", high = "red") +
#   facet_grid(cols = vars(CELL), rows = vars(tissue), scales = "free_x", space = "free_x") +
#   optimized_theme_fig()+
#   theme(
#     strip.text = element_text(size = 16, face = "bold", angle = 90)
#   )
# ggsave(basedir("Celltype_marker_incl_Ery.pdf"))
# #without B-cells and Ery
# # Filter for marker genes
# marker_results <- celltype_avg %>%
#   filter(Gene %in% marker_genes$Gene)%>%
#   left_join(marker_genes, by = "Gene")%>%
#   filter(CELL != "Bcell")
# #filter(CELL != "Ery")
# 
# marker_results <- marker_results %>%
#   filter(CELL %in% c(
#     "GMP","Early-Mye","Granulocyte",
#     "Megakaryocytes","Hematopoietic stem cells",
#     "Erythroid-like and erythroid precursor cells","Eo",
#     "Basophils",
#     "Monocytes",
#     "Early-Ery",
#     
#   ))
# 
# # Normalize Z-scores within each celltype
# marker_results_scaled <- marker_results %>%
#   group_by(celltype) %>%
#   mutate(Scaled_Expr = scale(Mean_Expr)) %>%
#   ungroup()
# 
# unique(marker_results$celltype)
# # Plot the scaled Z-scores
# ggplot(marker_results_scaled %>% filter(Gene %in% Gene[CELL == unique(CELL)]), 
#        aes(x = Gene, y = celltype, color = Mean_Expr)) +
#   geom_tile() +
#   scale_color_gradient2(low = "blue", mid = "white", high = "red") +
#   facet_grid(cols = vars(CELL), rows = vars(tissue), scales = "free_x", space = "free_x") +
#   optimized_theme_fig()+
#   theme(
#     strip.text = element_text(size = 16, face = "bold", angle = 90)
#   )
# ggsave(basedir("Celltype_marker.pdf"))
# 
# 

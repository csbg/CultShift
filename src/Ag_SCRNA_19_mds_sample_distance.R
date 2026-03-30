########################################################################
source("src/00_init.R")
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(tidyverse)
source("src/Ag_Optimized_theme.R")
library(ggplot2)
library(reshape2)

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
head(longer_dataVoom)
ex_vivo_data <- longer_dataVoom %>% filter(tissue == "ex.vivo")
in_vivo_data <- longer_dataVoom %>% filter(tissue == "in.vivo")
unique(ex_vivo_data$sample1)
# Step 2: Pivot ex vivo and in vivo data to wide format
# Pivot ex vivo data to wide format
ex_vivo_wide <- ex_vivo_data %>%
  mutate(sample_id = paste0(sample,tissue_celltype))%>%
  dplyr::select(genes, Expression, sample_id) %>%
  pivot_wider(
    names_from = sample_id, 
    values_from = Expression,
    values_fn = list(Expression = mean)
  )

# Pivot in vivo data to wide format
in_vivo_wide <- in_vivo_data  %>%
  mutate(sample_id = paste0(sample,tissue_celltype))%>%
  dplyr::select(genes, Expression, sample_id) %>%
  pivot_wider(
    names_from = sample_id, 
    values_from = Expression,
    values_fn = list(Expression = mean)
  )
# Ex vivo matrix
ex_mat <- ex_vivo_wide %>%
  column_to_rownames("genes") %>%
  as.matrix()

# In vivo matrix
in_mat <- in_vivo_wide %>%
  column_to_rownames("genes") %>%
  as.matrix()
colnames(ex_mat)
meta_lookup <- meta %>%
  mutate(sample_id = paste0(sample, tissue, "_",celltype)) %>%
  dplyr::select(sample_id, celltype, tissue)

#############

plot_celltype_dist <- function(ct) {
  
  samples_ct <- meta_lookup %>%
    filter(celltype == ct) %>%
    pull(sample_id)
  colnames(ex_mat)
  # Extract only the samples for this celltype
  ex_ct <- ex_mat[, colnames(ex_mat) %in% samples_ct, drop = FALSE]
  in_ct <- in_mat[, colnames(in_mat) %in% samples_ct, drop = FALSE]
  
  # Combine ex.vivo + in.vivo
  mat_ct <- cbind(ex_ct, in_ct)
  
  if (ncol(mat_ct) < 2) return(NULL)
  
  dist_mat <- dist(t(mat_ct), method = "euclidean")
  dist_df <- as.matrix(dist_mat)
  
  heat_df <- melt(dist_df, varnames = c("Sample1", "Sample2"), 
                  value.name = "Distance")
  
  ggplot(heat_df, aes(gsub("LSK_|in.vivo_|_Lin-|_lin-|_ckit|ex.vivo_[A-Za-z]*","",Sample1),
                      gsub("LSK_|in.vivo|_Lin-|_lin-|_ckit|ex.vivo_[A-Za-z]*","",Sample2), fill = Distance)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6)
    ) +
    labs(
      title = paste("Euclidean Distance –", ct),
      subtitle = "Within-celltype comparison: ex.vivo vs in.vivo",
      x = "Sample", y = "Sample"
    )
}
all_celltypes <- unique(meta_lookup$celltype)

plots <- lapply(all_celltypes, plot_celltype_dist)

plots[[6]]     # show first plot

# Remove NULL plots (if any celltype had <2 samples)
plots_clean <- plots[!sapply(plots, is.null)]

# Only keep plots that are actual ggplot objects
plots_clean <- plots[sapply(plots, function(x) inherits(x, "gg"))]

length(plots_clean)  # check how many plots remain

ggsave("sample_distance_plot.pdf", combined_plot, width = 30, height = 20, units = "cm")

##############
mat_ct <- cbind(ex_mat, in_mat)
if (ncol(mat_ct) < 2) return(NULL)

dist_mat <- dist(t(mat_ct))
mds <- cmdscale(dist_mat)

mds_df <- as.data.frame(mds) %>%
  rownames_to_column("sample_id") %>%
  mutate(
    tissue = ifelse(grepl("ex.vivo", sample_id), "ex.vivo", "in.vivo"),
    celltype = gsub(".*_", "", sample_id),
    celltype_tissue = paste0(celltype, "_", tissue)
  )
write_rds(mds_df ,basedir("mds_df.rds"))
condition_colors <- c("#3D1778" ,
                          "#82498C", "#D2A9DB","hotpink",
                          "#05547F",
                          "#3690C0",
                          "#6B91C6","grey"
  )
ggplot(mds_df, aes(V1, V2, color = celltype, shape = tissue)) +
  geom_point(size = 1.5, alpha = 0.7) +
  theme_minimal() +
  scale_color_manual(values = condition_colors)+
  labs(
    title = paste("MDS (Euclidean Distance)"),
    x = "Dimension 1", 
    y = "Dimension 2",
    color = "Tissue",
    shape = "Celltype"
  ) +
  optimized_theme_fig()
ggsave(basedir("MDS_plot.pdf"), w = 8, h = 7, units = "cm")
#######################################

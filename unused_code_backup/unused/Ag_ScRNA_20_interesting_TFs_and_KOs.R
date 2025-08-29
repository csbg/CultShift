########################################################################
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(tidyverse)
library(latex2exp)

# Input/Output directories-----------
InDir <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")

base<-"Ag_ScRNA_20_interesting_TFs_and_KOs/"
basedir <- dirout("Ag_ScRNA_20_interesting_TFs_and_KOs/")

########################################################################
#load data and clean metadata
########################################################################
#metadata from in-vivo ex-vivo
meta <- read.delim(InDir("metadata_guide.tsv"),row.names=1)%>%
  filter(genotype == "NTC")

meta <- meta[,c("celltype","tissue")]

meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]

# filtering steps
celltypes_to_exclude <- c("CLP",  "EBMP", "unclear","T-cell","MEP","Imm. B-cell")
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                      "Hba-a2","Pim1","Fabp5","Fdps","Cd9")

meta <- meta %>%
  # Clean row names
  tibble::rownames_to_column("sample") %>%
  mutate(
    sample = gsub("[\\ \\(\\)\\-]", ".", sample),
    sample = gsub("Eo/Ba", "Eo.Ba", sample),
    sample = gsub("IMP1|IMP2", "GMP", sample)
  ) %>%
  # Clean all columns
  mutate(across(everything(), ~ gsub("Eo/Ba", "Eo.Ba", .))) %>%
  mutate(across(everything(), ~ gsub("IMP1|IMP2", "GMP", .))) %>%
  # Update celltype
  mutate(celltype = ifelse(celltype %in% c("Eo", "Ba"), "Eo.Ba", celltype)) %>%
  # Remove unwanted cell types
  filter(!celltype %in% c("B-cell", "Ery", "Neu", "T-cell-Cd3d+", "E/B", "Imm. B-cell","B.cell")) %>%
  # Remove rows with "NA" or "T.cell" in sample names
  filter(!grepl("NA|T.cell", sample))
  
# Normalized_data
dataVoom <- read_rds(InDir2("dataVoom_perCTex.vivovsin.vivo.rds"))

#limmaRes from NTC
limmaRes <- read_rds(InDir2("limma_perCTex.vivovsin.vivo.rds"))

limmaRes_significant <- limmaRes %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05)%>%
  pull(genes)%>%
  unique()

#genes of used KOs


KO_genes <- read.delim(InDir("metadata_guide.tsv"),row.names=1)%>%
  group_by(genotype, tissue, celltype) %>%  
  filter(genotype != "NTC")%>% # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample), .groups = 'drop') %>% # Count distinct samples for each group
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>% # Spread tissue to separate columns (in.vivo and ex.vivo)
  group_by(genotype) %>%                                    # Regroup by genotype
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) %>%              # Keep genotypes that have at least one celltype with 3+ samples in both tissues
  pull(genotype) %>% unique()

# modify dataVoom
longer_dataVoom <-  dataVoom$E %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(meta, by ="sample")%>%
  mutate(tissue_celltype =paste0(tissue,"_",celltype))



# Step 1: Calculate Mean per Gene per (Sample, Cell Type, Tissue)
gene_mean <- longer_dataVoom %>%
  group_by(tissue_celltype,tissue,celltype, genes) %>%
  summarise(mean_expr = mean(Expression, na.rm = TRUE), .groups = "drop")
#

generate_heatmap <- function(data, geneset, filename_prefix, control = NULL) {
  filtered_data <- data %>%
    filter(genes %in% geneset) %>%
    group_by(genes) %>%  
    mutate(scaled_expr = scale(mean_expr)) %>% 
    ungroup() %>%  
    select(genes, tissue_celltype, tissue, celltype, mean_expr, scaled_expr)
  
  # Keep "Spi1" at the top or bottom
  filtered_data <- filtered_data %>%
    mutate(genes = factor(genes, levels = c(setdiff(unique(genes), control), control)))
  limmaRes_data <- limmaRes %>%
    filter(genes %in% c(geneset, control))%>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1)
  # Heatmap without replicates
  p1 <- ggplot(filtered_data, aes(x = tissue, y = genes, fill = scaled_expr)) +
    geom_tile() +
    facet_grid(cols = vars(celltype)) +
    labs(y = "TFs", x = "Tissues") +
    scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
    optimized_theme() +
    coord_equal()+
    theme(strip.text.x = element_text(angle = 90))
  ggsave(basedir(paste0(filename_prefix, "_mean_expr_scaled.pdf")), plot = p1)
  p2 <- ggplot(limmaRes_data,aes(x = celltype,
                                  y = genes,
                                  size = pmin(-log10(adj.P.Val),5),
                                  color = logFC)) +
    geom_point() +
    scale_color_gradient2(high = "red",
                          mid = "white",
                          low = "blue")+
    scale_size_continuous(limits = c(1,5))+
    optimized_theme_fig()+theme(
      legend.position = "right"
    )
    coord_fixed()
  
  ggsave(basedir(paste0(filename_prefix, "_logFC.pdf")), plot = p2)
  
  
}

#  TFs ------------------------------------
TF_JAK_STAT <- c("Stat1", "Stat2", "Irf9", "Irf2", "Irf1", "Spi1")
TF_JAK_STAT <- c(TF_JAK_STAT[TF_JAK_STAT %in% limmaRes_significant],"Spi1")
  
TF_JAK_STAT_plot <- generate_heatmap(gene_mean, TF_JAK_STAT, "TFs_JAK_stat", "Spi1")

# JAK_stat and mtorc
TFs <- c("Stat1", "Stat2", "Irf9", "Irf2", "Irf1", "Spi1",
         "Myc", "Atf4", "Hif1a", "Eif4ebp1", "Eif4ebp2",
         "Eif4ebp3", "Tfe3", "Mitf", "Foxo1", "Foxo3", 
         "Foxo4", "Foxo6", "Mlx", "MondoA (Mlxipl)", "Srebf1", "Srebf2"
)
TFs <- c(TFs[TFs %in% limmaRes_significant],"Spi1")
TFs_plot <- generate_heatmap(gene_mean, TFs, "TFs", "Spi1")
# KOs---------------------------------------------------------------------------# Example usage for KOs
KO_plot <- generate_heatmap(gene_mean, KO_genes, "KOs")


# JAK_stat and mtorc combined
TFs_JAK_STAT <- c("Stat1", "Stat2", "Irf9", "Irf2", "Irf1","Irf8")
TF_control <- "Spi1"
TFs_mTORC1 <-c("Myc", "Atf4", "Hif1a", "Eif4ebp1", "Eif4ebp2",
         "Eif4ebp3", "Tfe3", "Mitf", "Foxo1", "Foxo3", 
         "Foxo4", "Foxo6", "Mlx", "MondoA (Mlxipl)", "Srebf1", "Srebf2")

gene_labels_df <- data.frame(
  genes = c(TFs_JAK_STAT, TF_control, TFs_mTORC1),
  label = c(rep("JAK-STAT", length(TFs_JAK_STAT)),
            rep("Control", length(control)),
            rep("mTORC1", length(TFs_mTORC1)))
)

#Expression plot
TF_data <- gene_mean %>%
  filter(genes %in% gene_labels_df$genes) %>%
  group_by(genes) %>%  
  mutate(scaled_expr = scale(mean_expr)) %>% 
  ungroup() %>%  
  select(genes, tissue_celltype, tissue, celltype, mean_expr, scaled_expr)%>%
  inner_join(gene_labels_df, by = "genes")


ggplot(TF_data, aes(x = tissue, y = genes, fill = scaled_expr)) +
  geom_tile() +
  facet_grid(cols = vars(celltype), rows = vars(label), scales = "free", space = "free") +
  labs(y = "TFs", x = "Tissues") +
  scale_fill_gradient2(high = "red", mid = "white", low = "blue") +
  optimized_theme() +
  #coord_equal()+
  theme(strip.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle =  0))
ggsave(basedir("TFs_JAK_stat_mtorc_Expression.pdf"))

#logFC plot
limmaRes_data <- limmaRes %>%
  filter(genes %in% gene_labels_df$genes)%>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  inner_join(gene_labels_df, by = "genes") 

ggplot(limmaRes_data,aes(x = celltype,
                         y = genes,
                         size = pmin(-log10(adj.P.Val),5),
                         color = logFC)) +
  geom_point() +
  facet_grid(rows = vars(label) ,space = "free", scales = "free_y")+
  scale_color_gradient2(high = "red",
                        mid = "white",
                        low = "blue")+
  scale_size_continuous(limits = c(1,5))+
  optimized_theme_fig()+
  theme(
    legend.position = "right"
  )+
  theme(aspect.ratio = 1)
ggsave(basedir("TFs_JAK_stat_mtorc_logFC.pdf"))


#  with replicates
TF_with_rep <- longer_dataVoom %>%
  filter(genes %in% TF_list) %>%
  group_by(genes) %>%  # Scale within each tissue_celltype
  mutate(scaled_expr = scale(Expression)) %>% 
  ungroup() %>%  # Remove grouping after scaling
  select(genes, tissue_celltype, tissue, celltype, Expression, scaled_expr,sample)
  
ggplot(TF_with_rep, aes(x = gsub("_*","",sample), y = genes, fill = scaled_expr))+
  geom_tile()+
  facet_grid(cols  = vars(celltype),scales = "free_x")+
  labs(
    y = "TFs",
    x = "tissues"
  )+
  scale_fill_gradient2(high = "red",
                       mid = "white",
                       low = "blue")+
  optimized_theme()
ggsave(basedir(paste0("TFs_JAK_stat_expr_scaled", ".pdf")),w=30,h=20)


#KO_genes logFCs in NTC-------------------

KO_genes_logFC <- limmaRes %>%
  filter(genes %in% KO_genes)%>%
  filter(adj.P.Val < 0.05) %>%
  select(genes, logFC, adj.P.Val, celltype)



#ggplot KO_genes_logFC
ggplot(KO_genes_logFC,aes(x = celltype,
                          y = genes,
                          size = pmin(-log10(adj.P.Val),5),
                          color = logFC)) +
  geom_point() +
  scale_color_gradient2(high = "red",
                        mid = "white",
                        low = "blue")+
  scale_size_continuous(limits = c(1,5))+
  optimized_theme()+
  coord_fixed()
ggsave(basedir("differential_expr_of_target_KO.pdf"), )  
  

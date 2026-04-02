source("src/00_init.R")
source("src/Ag_ko_classification.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)
library(ggpubr)
library(ggplot2)


Indir2 <- dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/")
#Indir3 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
Indir5 <- dirout("Ag_ScRNA_13_zscore_plots_celltype_marker/")
Indir6 <- dirout("Ag_ScRNA_19_invivo_exvivo_izzo_zscore/")
out <- "Figure1_Supplementary"
outdir <- dirout("Figure1_Supplementary")
source("src/Ag_Optimized_theme_fig.R")
########################
ENRICHR <- dirout(paste0("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/ENRICHR"))
ENRICHR.DBS <-union(ENRICHR.DBS,
                    c("GO_Biological_Process_2021",
                      "TRRUST_Transcription_Factors_2019",
                      "Reactome_2022",
                      "GO_Molecular_Function_2023",
                      "GO_Biological_Process_2023",
                      "CellMarker_2024"))
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
# # save(enr.terms, file=out("Genesets_Human.RData"))
# Convert to mouse --------------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl){
  dbl <- lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})
#######
#Supp_Fig_1a --------------
#########
marker_results_scaled <- read_rds(Indir5("scaled_zscore_marker_genes.rds"))
marker_results_scaled$CELL <- factor(marker_results_scaled$CELL, 
                                     levels = c("Baso","Early Mye","Granulocyte","HSC","Ery","Mega","Mono"))
#Scaled_Zscore here is the zscore calculated across celltype for the same gene
#within the same tissue 
ggplot(marker_results_scaled, #%>% filter(Gene %in% Gene[CELL == unique(CELL)]), 
       aes(x = Gene, y = tissue, fill = pmin(pmax(Scaled_Zscore, -2), 2)))+#, size = Mean_Zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4C889C",
                       mid = "white",
                       high = "#D0154E",
                       name= "Expression"
  )+#,
  
  facet_grid(cols = vars(CELL), rows = vars(celltype), scales = "free", space = "free") +
  optimized_theme_fig()+
  theme(strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 90, hjust = 0))

ggsave(outdir("Supp.Fig.1A.pdf"), w=18, h=9,units = "cm")
###################

########################
#Fig1B----------
limmaRes_NTC_wide <- limmaRes_NTC %>%
  dplyr::select(genes, celltype, logFC) %>%
  pivot_wider(names_from = celltype, values_from = logFC)
# Remove gene names column and compute correlation
cor_matrix <- cor(limmaRes_NTC_wide[,-1], use = "pairwise.complete.obs")
diag(cor_matrix) <- NA
# Optional: define a color function
col_fun <- colorRamp2(c(-1, 0, 1), c("#4C889C", "white", "#D0154E"))

# Plot the correlation matrix
Heatmap(
  cor_matrix,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE
)

cor_long <- melt(cor_matrix)
colnames(cor_long) <- c("CellType1", "CellType2", "Correlation")


ggplot(cor_long, aes(CellType1, CellType2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4C889C", mid = "white", high = "#D0154E", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name = "Correlation") +
  labs(x = "Cell type",
       y = "Cell type",
       title = "Correlation of culture effects")+
  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
  ) +
  optimized_theme_fig()+
  
  coord_fixed()
ggsave(outdir("Sup_Fig.1B.pdf"), w= 8, h=8, units = "cm")



# Filter out NA correlations
ggplot(cor_long, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.1, fill = "#4C889C", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ CellType2, scales = "free_y") +
  labs(x = "Correlation", y = "Count",
       title = "Distribution of Correlation Values by Cell Type") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 10)) +
  optimized_theme_fig()  # Optional

#Fig1C------------------
InDir <- dirout("Ag_ScRNA_19_invivo_exvivo_zscore/")

mds_df <- read_rds(InDir("mds_df.rds"))
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
  optimized_theme_fig()+
  theme(axis.text.x = element_text(angle = 0))
ggsave(outdir("MDS_plot.pdf"), w = 8, h = 7, units = "cm")

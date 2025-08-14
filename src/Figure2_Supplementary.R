source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(latex2exp)
out <- "Figure2_Supplementary"
basedir <- dirout("Figure2_Supplementary")

Indir3 <- dirout("Ag_ScRNA_19_invivo_exvivo_external_zscore/")
Indir1 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
Indir2 <- dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/")
Indir4 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
Indir5 <- dirout("Figure1")
#Supp_Fig1c----------------
##########
limmaRes_NTC <- read_rds(Indir1("limma_perCTex.vivovsin.vivo.rds"))
dataVoom_NTC_in_ex <- read_rds(Indir1("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(Indir1("NTC_meta.rds"))
##########
get_top_genes <- function(pathway_name,
                          pathway_genes,
                          limma_results, 
                          logFC_threshold = 1,
                          pval_threshold = 0.05, 
                          top_n = 5) {
  # Filter limma results for the pathway genes
  filtered_genes <- limma_results %>%
    filter(group != "n.s") %>%
    group_by(celltype) %>%
    filter( adj.P.Val < pval_threshold) %>%
    filter(toupper(genes) %in% toupper(pathway_genes)) %>%
    arrange(adj.P.Val) %>%
    #mutate(abs_logFC = abs(logFC)) %>%
    slice_head(n = top_n) %>%
    pull(genes)%>%unique()
  
  # Extract and return data for the filtered genes
  pathway_plot <- limma_results %>%
    filter(toupper(genes) %in% toupper(filtered_genes))
  #%>%
  #filter(group != "n.s") 
  
  return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
}


gsea.res <- read_rds(Indir2("NTC_fgsea.rds"))
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][, -c("log2err", "size", "pval"), with = F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, function(vec) paste(vec[1:10], collapse = ","))

gsea.res <- read_rds(Indir2("NTC_fgsea.rds"))
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][, -c("log2err", "size", "pval"), with = F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, function(vec) paste(vec[1:10], collapse = ","))

unique(gsea.res$db)
dbx <- "KEGG_2019_Mouse"
# Iterate over each database
#for (dbx in unique(gsea.res$db)) {
dat <- dirout(paste0(out, "FGSEA/", dbx))

# Define file path
output_file <- dat("GSEA_significant_", dbx, ".tsv")

# Ensure the directory exists before writing
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)  # Create directory if it doesn't exist
}
# Convert list columns to character before writing
df <- gsea.res[db == dbx]

df[] <- lapply(df, function(x) if (is.list(x)) sapply(x, toString) else x)

# Now, write the table
write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Get the pathways for this database
pDT <- gsea.res[db == dbx]

pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n = 5), by = c("celltype")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n = 5), by = c("celltype")]$pathway)
pw.display <- unique(c(pw.display.pos, pw.display.neg))

pDT <- pDT[pathway %in% pw.display]

if (nrow(pDT) > 0) {
  # Aggregate NES values across all cell types
  pDT_agg <- pDT %>%
    group_by(pathway) %>%
    summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
    arrange(desc(average_NES)) 
  
  #pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)
  # Ensure pDT follows the same order as pDT_agg
  pDT <- pDT %>%
    mutate(pathway = factor(pathway, levels = pDT_agg$pathway)) %>%
    arrange(factor(pathway, levels = pDT_agg$pathway))  # Explicitly reorder pDT
  # Step 3: Plot with the new pathway order (highest NES first)
  
  pathways <- unique(pDT$pathway)
  # Split the leadingEdge column into separate genes and create a combined list for each pathway
  combined_genes <- pDT %>%
    group_by(pathway) %>%
    summarise(
      Genes = list(unique(unlist(leadingEdge))),  # Combine and take the union of gene lists
      .groups = "drop"
    )
  
  # View the cleaned-up result
  print(combined_genes)
  
  
  
  # Function to extract top genes for each pathway
  get_top_genes <- function(pathway_name, 
                            limma_results,
                            logFC_threshold = 1,
                            pval_threshold = 0.05,
                            top_n = 5) {
    
    # Extract pathway genes
    pathway_genes <- combined_genes %>%
      filter(pathway == pathway_name) %>%
      pull(Genes) %>%
      unlist() %>%
      unique()
    
    # Filter limma results for genes in the pathway
    filtered_genes <- limma_results %>%
      filter(group != "n.s") %>%
      filter(adj.P.Val < pval_threshold) %>%
      filter(genes %in% pathway_genes) %>%
      arrange(adj.P.Val, desc(abs(logFC))) %>%
      group_by(celltype) %>%  # Ensure top genes per cell type & tissue
      slice_head(n = top_n) %>%
      ungroup() %>%
      mutate(pathway = pathway_name)
    
    return(filtered_genes)
  }
  
  # Initialize final combined results table
  all_top_genes <- list()
  # Define the database of interest
  dbx <- "GO_Biological_Process_2023"  # Replace with your actual database name
  
  # Create directory path
  dat <- dirout(paste0(out, "/FGSEA/", dbx))
  
  
  # Extract unique pathways
  pathways <- unique(pDT$pathway)
  
  # Get combined gene data for all pathways
  results <- map_dfr(pathways, function(pathway) {
    get_top_genes(
      pathway_name = pathway,
      limma_results = limmaRes_NTC
    )
  })
  results <- results %>%
    dplyr::select(genes,pathway)%>%
    distinct()
}

top_genes <- limmaRes_NTC %>%
  inner_join(results, by = "genes")
#interaction logFCs
limmaRes <- read_rds(Indir4("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
#exclude fig1 genes
genes_fig1 <- read_rds(Indir5("genes_fig1.rds"))
# select significant genes fron interaction, present in results
top_int_genes <- limmaRes %>%
  mutate(genes = ensg)%>%
  inner_join(results, by = "genes") %>%
  filter(group != "n.s") %>%
  filter(!(genes %in% genes_fig1)) %>%
  pull(genes)%>%
  unique()

top_genes <- limmaRes_NTC %>%
  inner_join(results, by = "genes") %>%
  filter(genes %in% top_int_genes) %>%
  filter(!(genes %in% genes_fig1$genes)) %>%
  group_by(genes) %>%
  filter(n_distinct(celltype[group != "n.s"]) >= 2) %>%  # Keep genes significant in at least 3 celltypes
  ungroup() %>%
  group_by(pathway, genes) %>%  
  slice_max(order_by = abs(logFC), n = 1) %>%  # Keep the gene entry with the highest logFC per pathway
  ungroup() %>%
  group_by(pathway) %>%
  slice_head(n = 50) %>%  # Select top 10 distinct genes per pathway
  ungroup() %>%
  dplyr::select(genes,pathway)





top_genes_NTC <-  limmaRes_NTC %>%
  filter(genes %in% top_genes$genes) %>%
  inner_join(top_genes, by = "genes")
#pathway_list <- unique(top_genes_NTC$pathway)
top_pathway <- pDT %>%
  group_by(pathway) %>%
  arrange(padj)%>%
  slice_head(n=40) %>%
  summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
  pull(pathway)


#genes_to_plot
top <- top_genes_NTC %>%
  filter(pathway %in% top_pathway) %>%
  mutate(geneset = case_when(
    pathway %in% c("Ribosome", "Ribosome biogenesis in eukaryotes") ~ "Ribosome machinery",
    pathway %in% c("DNA replication", 
                   "Cell cycle" 
    ) ~ "Replication/ cell cycle",
    
    pathway %in% c("Oxidative phosphorylation" 
    ) ~ "Oxphos/ Electron transport",
    
    
    pathway %in% c("Cell adhesion molecules (CAMs)" 
    ) ~ "Cell adhesion molecules (CAMs)",
    TRUE ~ as.character(pathway)  # Convert factor to character to avoid issues
  )) %>%
  group_by(genes) %>%
  filter(sum(group != "n.s.") >= 2) %>%  # Exclude genes n.s. in more than 2 cases
  ungroup() %>%
  dplyr::select(-pathway) #%>%

# Only select those pathways
top <- top %>%
  filter(geneset %in% c(
    "Replication/ cell cycle",
    "Oxphos/ Electron transport",
    "Ribosome machinery",
    "Cell adhesion molecules (CAMs)")) %>%
  filter(!(genes %in% c("Cdh1","Cd40","Mag","Pvr","Itgb7","Cox7a1","E2f2")))
unique(top$geneset)
#include myc
# Append Myc to top with gene set info
top <- rbind(
  top,
  limmaRes_NTC %>%
    filter(genes == "Myc") %>%
    mutate(
      group = ifelse(logFC > 0, "up", "down"),
      geneset = "growth/cellcycle"
    )
)
top$geneset <- factor(top$geneset, levels = c("Cell adhesion molecules (CAMs)" ,
                                              "Oxphos/ Electron transport",
                                              "Ribosome machinery",
                                              "Replication/ cell cycle",
                                              "growth/cellcycle"))
Supp_fig_1c <- ggplot(top, aes(
  x = celltype, 
  y = genes,
  color = pmin(2, pmax(-2, logFC)), 
  size = pmin(3, -log10(adj.P.Val))
)) +
  geom_point() +  # Base layer for all points
  # geom_point(
  #   data = subset(top, adj.P.Val < 0.05),  # Add black ring only for significant points
  #   aes(x = celltype, y = genes),
  #   shape = 21,  # Circle with border
  #   color = "black",  # Black outline
  #   fill = NA,  # No fill inside the ring
  #   stroke = 0.3  # Thickness of the black ring
  # ) +
  scale_color_gradient2(low = "#4C889C", mid = "white", high = "#D0154E") +
  scale_size_continuous(
    range = c(0, 1.5),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    title = "Differentially Expressed Genesets",
    y = "Genes",
    x = "Cell Type",
    color = "logFC",
    size = "-log10(padj)"
  ) +
  facet_grid(rows = vars(geneset), scales = "free", space = "free") +
  optimized_theme_fig() +
  theme(strip.text.y = element_text(angle = 0, hjust = 0))

Supp_fig_1c


#save
ggsave(basedir("Supp.Fig.1C.pdf"), plot = Supp_fig_1c, w = 9.5, h = 9.5, units = "cm")

#external  zscore
#Supp_fig_1d 
zscore_external <- read_rds(Indir3("zscore_plot_external.rds"))
unique(zscore_external$tissue_celltype)
ex.vivo_cells <- grep("ex.vivo",unique(zscore_external$tissue_celltype),value = T)
in.vivo_cells <- grep("in.vivo",unique(zscore_external$tissue_celltype),value = T)
izzo_cells <- grep("izzo",unique(zscore_external$tissue_celltype),value = T)
Anna_cells <-grep("Anna",unique(zscore_external$tissue_celltype),value = T)
zscore_external$tissue_celltype <- factor(zscore_external$tissue_celltype, 
                          levels = c(ex.vivo_cells,in.vivo_cells,izzo_cells,Anna_cells))
unique(zscore_external$tissue_celltype)
#zscore_external <- gsub("izzo","Izzo_et_al",zscore_external )
plot <- ggplot(zscore_external, aes(x = sample, y = genes, fill = zscore)) +
  geom_tile(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) +
  facet_grid(
    cols = vars(tissue_celltype),
    scales = "free",
    space = "free",
    labeller = labeller(tissue_celltype = function(x) gsub("^(ex\\.vivo_|in\\.vivo_|izzo_|Anna_et_al_)", "", x))
  ) +
  labs(title = "Z-score of Gene Expression Across Tissues", 
       x = "Genes", 
       y = NULL) +
  scale_fill_gradient2(low = "#4C889C",
                       mid = "white",
                       high = "#D0154E", midpoint = 0) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(angle = 90, hjust = 0),
    panel.spacing = unit(0, "lines"))         # <--- reduces facet spacing
plot1 <- plot + theme(legend.position = "none")
ggsave(basedir(paste0("Z-score_Distribution_", "_per_Tissue_line_without_guide.pdf")), plot = plot1, w=18,
       h=8, units = "cm")
plot2 <- plot + theme(legend.position = "bottom")
ggsave(basedir(paste0("Z-score_Distribution_", "_per_Tissue_line.pdf")), plot = plot2, w=18,
       h=8, units = "cm")
library(ggplot2)

# Boxplot version of your plot
library(ggplot2)

ggplot(zscore_external, aes(x = genes, y = zscore)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +         # boxplot without outliers (since jitter shows points)
  geom_jitter(aes(color = sample), width = 0.2, size = 0.3, alpha = 0.6) +  # jitter per sample
  facet_wrap(~ tissue_celltype, scales = "free_x", nrow = 1,
             labeller = labeller(tissue_celltype = function(x) gsub("^(ex\\.vivo_|in\\.vivo_|izzo_|Anna_et_al_)", "", x))) +
  labs(title = "Z-score of Gene Expression Across Tissues",
       x = "Genes",
       y = "Z-score") +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"  # or keep legend if you want to show samples colors
  )

print(boxplot)
library(dplyr)
library(ggplot2)

# Step 1: Calculate median zscore per sample per tissue_celltype across genes
median_per_sample <- zscore_external %>%
  group_by(tissue_celltype, sample, tissue) %>%
  summarize(median_zscore = median(zscore, na.rm = TRUE)) %>%
  ungroup()
my_color <- c("#004949", "#009292",  "#006DDB", "#B66DFF")
# Step 2: Plot boxplot per tissue_celltype and jitter points per sample
plot3 <- ggplot(median_per_sample, aes(x = tissue_celltype, y = median_zscore)) +
  geom_boxplot(aes(color = tissue), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(color = "darkgrey", width = 0.2, size = 0.8, alpha = 0.7) +
  labs(
    title = "Median Z-score per Sample across ISG-core genes",
    x = "Tissue-Celltype",
    y = "Median Z-score"
  ) +
  scale_color_manual(
    values = my_color,
    labels = c(
      "Ext.data\n(Anna Konturek-Ciesla et al)",
      "Ex vivo",
      "In vivo",
      "Ext.data\n(Izzo et al)"
    )
  ) +
  scale_x_discrete(labels = function(x) gsub("^(ex\\.vivo_|in\\.vivo_|izzo_|Anna_et_al_)", "", x)) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right"
  )

ggsave(basedir(paste0("Sup_Fig2B", "_per_Tissue_line.pdf")), plot = plot3, w = 18, h = 5, units = "cm")


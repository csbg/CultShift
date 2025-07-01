source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)
library(ggpubr)
library(ggplot2)

Indir1 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
Indir2 <- dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/")
Indir3 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
Indir4 <- dirout("Figure1")
Indir5 <- dirout("Ag_ScRNA_18_zscore_plots_celltype_marker/")
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
ggplot(marker_results_scaled %>% filter(Gene %in% Gene[CELL == unique(CELL)]), 
       aes(x = Gene, y = tissue, fill = pmin(pmax(Scaled_Zscore, -2), 2)))+#, size = Mean_Zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4C889C",
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
        strip.text.x = element_text(angle = 90, hjust = 0))

ggsave(outdir("Supp.Fig.1A.pdf"), w=18, h=9,units = "cm")
###################
limmaRes_NTC <- read_rds(Indir1("limma_perCTex.vivovsin.vivo.rds"))
head(limmaRes_NTC)
dataVoom_NTC_in_ex <- read_rds(Indir1("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(Indir1("NTC_meta.rds"))

top_genes <- limmaRes_NTC[limmaRes_NTC$genes %in% c("Idi1","Oas2","Msmo1"),] %>%
  unique()

#fig-----
# Ensure factor levels are correctly assigned
limmaRes_NTC$celltype <- factor(limmaRes_NTC$celltype,  
                                levels = c("HSC", "MEP.early", "MkP",  
                                           "GMP", "Gran.P", "Gran.",  
                                           "Mono", "Eo.Ba"), ordered = TRUE)
limmaRes_NTC <- limmaRes_NTC %>%
  mutate(group = recode(group,
                        down = "downregulation ex vivo",
                        up = "upregulation ex vivo"
  ))

# Debugging: Check if factor levels are correct
print(levels(limmaRes_NTC$celltype))

Fig1.2 <- ggplot() +
  stat_bin_hex(data = filter(limmaRes_NTC, group == "n.s"), 
               aes(x = logFC, y = -log10(adj.P.Val), fill = ..count..), 
               bins = 20, color = NA, alpha = 0.7) +
  scale_fill_gradient(low = "lightgrey", high = "black",
                      limits = c(1, 5000), name = "Gene Count") +
  
  stat_bin_hex(data = filter(limmaRes_NTC, group == "upregulation ex vivo"), 
               aes(x = logFC, y = -log10(adj.P.Val), fill = ..count..), 
               bins = 20, color = NA, fill = "#D0154E", alpha = 0.7) +
  
  stat_bin_hex(data = filter(limmaRes_NTC, group == "downregulation ex vivo"), 
               aes(x = logFC, y = -log10(adj.P.Val), fill = ..count..), 
               bins = 20, color = NA, fill = "#4C889C", alpha = 0.7) +
  
  # Draw black dots for target genes
  geom_point(
    data = top_genes,
    aes(x = logFC, y = -log10(adj.P.Val)),
    color = "black",
    size = 0.5
  ) +
  
  # Label target genes
  geom_text_repel(
    data = top_genes,
    aes(x = logFC, y = -log10(adj.P.Val), label = genes),
    size = 3, 
    color = "black",
    max.overlaps = 100,
    force = 10,
    force_pull = 0.1,
    max.iter = 3000,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.color = "black",
    segment.size = 0.3,
    min.segment.length = 0.02,
    arrow = arrow(length = unit(0.02, "npc"), type = "closed", angle = 25)
  ) +
  

  
  labs(title = "DEGs (Ex-vivo vs in-vivo)",
       x = "logFC",
       y = "-log10(adj.P)") +
  
  facet_wrap(~ factor(celltype, levels = c("HSC", "MEP.early", "MkP",  
                                           "GMP", "Gran.P", "Gran.",  
                                           "Mono", "Eo.Ba")),
             scales = "free", drop = FALSE) +
  
  optimized_theme_fig()+
  theme(
    panel.spacing = unit(0.0001, "lines")  # reduce facet spacing
  )

print(Fig1.2)

ggsave(outdir(paste0("Fig1.2.pdf")), plot=Fig1.2,w=11,h=8,units = "cm")

#

########################
limmaRes_NTC_wide <- limmaRes_NTC %>%
  select(genes, celltype, logFC) %>%
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
       y = "Cell type")+
  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
  ) +
  optimized_theme_fig()+
  
  coord_fixed()
ggsave(outdir("Correlation_culture_effects_across_cell_type.pdf"), w= 8, h=8, units = "cm")



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
limmaRes <- read_rds(Indir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
#exclude fig1 genes
genes_fig1 <- read_rds(Indir4("genes_fig1.rds"))
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
  select(genes,pathway)





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
   select(-pathway) #%>%
 
# Only select those pathways
top <- top %>%
  filter(geneset %in% c(
    "Replication/ cell cycle",
    "Oxphos/ Electron transport",
    "Ribosome machinery",
    "Cell adhesion molecules (CAMs)")) %>%
  filter(!(genes %in% c("Cdh1","Cd40","Mag","Pvr","Itgb7","Cox7a1","E2f2")))
unique(top$geneset)

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
ggsave(outdir("Supp.Fig.1C.pdf"), plot = Supp_fig_1c, w = 9.5, h = 9.5, units = "cm")
#save genes
supp_fig1_genes <- top %>%
  select(genes,geneset)
supp_fig1_genes %>% write_rds(outdir("supp_fig1_genes.rds"))

genes_ntc <- genes_fig1 %>%
  mutate(geneset = pathway)%>%
  select(genes,geneset) %>%
  filter(!(geneset %in% "mTORC1_or_Cholesterol"))%>%
  rbind(supp_fig1_genes)
top_int <- limmaRes %>%
  filter(ensg %in% genes_ntc$genes)%>%
  mutate(genes = ensg) %>%
  inner_join(genes_ntc, by = "genes")



meta <- fread(Indir1("meta_cleaned.tsv")) # Read data
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   

meta <- meta[, -1, drop = FALSE] 
meta$sample1 <- rownames(meta)


# Check if there are at least 3 distinct samples per tissue for each genotype and celltype
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop")%>%
  mutate(coef = genotype)
summary_df <- limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < 0.05 & logFC > 1),
    Downregulated = sum(adj.P.Val < 0.05 & logFC < -1)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")

unique(top_int$geneset)
top_int <- limmaRes %>%
  filter(ensg %in% genes_ntc$genes)%>%
  mutate(genes = ensg) %>%
  inner_join(genes_ntc, by = "genes")%>%
  inner_join(ko_flags, by = c("coef","celltype"))%>%
  filter(valid_ko)%>%
  #filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>% 
  filter(geneset != "ISG core")
  
Fig3.3.2_supplementary <- ggplot(top_int, 
                                 aes(x = coef, y = ensg,
                                     color = pmin(1.5, pmax(-1.5, logFC)),
                                     size = pmin(3, -log10(adj.P.Val))  # Set size by adjusted p-value
                                 )) +
  geom_point() + 
  
  # Add black ring only for significant points with adjusted stroke size
  # geom_point(
  #   data = subset(top_int, adj.P.Val < 0.05),  # Add black ring only for significant points
  #   aes(x = coef, y = ensg,
  #       size = pmin(5, -log10(adj.P.Val))),
  #   shape = 21,  # Circle with border
  #   color = "black",  # Black outline
  #   fill = NA,  # No fill inside the ring
  #   stroke = 0.5  # Set stroke to match the size
  # ) +  # Create points
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0, 1.5),
   # limits = c(0, 3),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(title = "Differentially expressed genesets",
       x = "KOs",
       y = "Genes") +
  facet_grid(rows = vars(geneset), cols = vars(celltype),
             scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig() + 
  theme(
    axis.text = element_text(size = 5),
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.text.x = element_text(angle = 90, hjust = 0)
  )

Fig3.3.2_supplementary

ggsave(outdir("Supplementary.Fig.3.1.pdf"), w=12, h= 17, units = "cm")
############################################################################
#izzo  zscore
#Supp_fig_1d 
zscore_izzo <- read_rds(Indir6("zscore_plot_izzo.rds"))
ggplot(longer_dataVoom_zscore, aes(x = sample, y = genes , fill = zscore)) +
  geom_tile(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) +  # Scatter plot with jitter for better visualization
  facet_grid(cols = vars(tissue_celltype), scales = "free", space = "free") +  # Separate by tissue, free_x ensures that each tissue has its own x-axis range
  labs(title = "Z-score of Gene Expression Across Tissues", 
       x = "Genes", 
       y = "Z-score") +
  scale_fill_gradient2(low = "#4C889C",
                       mid = "white",
                       high = "#D0154E", midpoint = 0) +
  optimized_theme_fig()+
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90,hjust = 0))+
  theme(panel.spacing = unit(0.2, "lines")) # Rotate x-axis labels for readability
# Remove minor gridlines

ggsave(basedir(paste0("Z-score_Distribution_", "_per_Tissue_line.pdf")), w=18,
       h=8, units = "cm")

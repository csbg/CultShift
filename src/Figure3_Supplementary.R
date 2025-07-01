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
Indir6 <- basedir <- dirout("Ag_ScRNA_19_invivo_exvivo_izzo_zscore/")

outdir <- dirout("Figure3_Supplementary")
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
##################################################################################

InDir7  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide")
gsea.res <- read_rds(InDir7("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
gsea.res$coef <- gsub("interaction","",gsea.res$coef )
# Step 2: Summarize to find KOs with at least one valid cell type
valid_ko_summary <- ko_flags %>%
  group_by(genotype) %>%
  summarize(has_valid_celltype = any(valid_ko), .groups = 'drop')


# Step 3: Filter GSEA results based on valid KOs

db = "MSigDB_Hallmark_2020"
create_gsea_plot <- function(db) {
  
  # Filter the data for the given database
  pDT <- gsea.res %>%
    filter(coef %in% koi) %>%
    filter(db == !!db) %>%  # Correctly reference the current database
    left_join(ko_flags, by = c("coef", "celltype")) %>%
    left_join(summary_df, by = c("coef", "celltype")) %>%
    filter(Count > 5) %>%
    filter(valid_ko)%>%
    dplyr::select(-c("Count", "Regulation"))
  # Merge with KO flags
  
  # Step 3 continued: Keep only valid KOs for the specific cell type
  pDT <- pDT %>% filter(valid_ko == TRUE, padj < 0.05)
  
  # Select the pathways for plotting (both positive and negative)
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n = 7), by = c("coef", "celltype", "pathway")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n = 7), by = c("coef", "celltype", "pathway")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(union(pw.display.pos, pw.display.neg))
  # Remove duplicate rows
  pDT <- pDT %>% distinct()
  
  # Filter pDT to include only selected pathways
  # **Convert list-columns to character format**
  dat <- pDT %>%
    mutate(across(where(is.list), ~sapply(., toString)))  # Converts lists to comma-separated strings
  
  # Save the filtered data table
  write.table(dat, file = basedir(paste0("fgsea_", db, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # Create the plot for the current database
  fig <- ggplot(pDT, aes(x = coef, y = pathway, color = NES, size = pmin(5, -log10(padj)))) +
    geom_point() + 
    scale_color_gradient2(low = "#4C889C", mid = "white", high = "#D0154E", name = TeX("NES")) +
    geom_point(data = pDT[padj < 0.05], shape = 1) +
    scale_size_continuous(range = c(0, 2), limits = c(0, 5), name = TeX("$-\\log_{10}(p_{adj})$")) +
    theme_bw() +
    xRot() +
    labs(x = "KOs") +
    facet_grid(cols = vars(celltype), scales = "free", space = "free") +  # Create facets for each cell type
    optimized_theme_fig()+
    theme(strip.text.x = element_text(angle = 90))
  #strip.text.y = element_text(angle = 0)) +
  
  # Save the plot for the current database
  ggsave(basedir("fig2.2_supplementary_fgsea", db, "_3rep.pdf"), fig,
         w =11, h =12, units = "cm")
  
  
  
}

# Example usage: 
# Create the plot for "MSigDB_Hallmark_2020"
create_gsea_plot("MSigDB_Hallmark_2020")
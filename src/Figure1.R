#figure 1
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtext)
library(latex2exp)
library(circlize)
library(ComplexHeatmap)
##################################################################################
base <- "Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/"
basedir <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
out <- "Figure1"
outdir <- dirout("Figure1")
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
#s
########################################################################
#load data
########################################################################
#Fig1.1
InDir1 <- dirout("FIG_02_scRNA_UMAPs_ar/")
pDT.labels <- read_rds(InDir1("pDT.labels.rds"))
merged_data <- read_rds(InDir1("Cross_projected_on_in.vivo.rds"))
merged_data$functional.cluster <- factor(merged_data$functional.cluster, 
                                         levels = names(cluster_colors))

# color coding
cluster_colors <- c(
  "Mono" = "#E69F00",      # Orange
  "Eo/Ba" = "#56B4E9",     # Sky Blue
  "GMP" = "#009E73",       # Green
  "MEP (early)" = "#F0E442", # Yellow
  "MkP" = "#CC79A7",       # Pink/Purple
  "Gran. P" = "#0072B2",   # Blue
  "Gran." = "#D55E00",     # Reddish Orange
  "HSC" = "#A020F0",       # Purple
  "GMP (early)" = "#999999",  # Light Gray (unchanged)
  "CLP" = "#D9D9D9",       # Light gray for CLP
  "unclear" = "#B0B0B0",    # Gray for unclear
  "Imm. B-cell" = "#8DA0CB", # Soft Blue
  "MEP" = "#D3D3D3",        # Lighter gray for MEP
  "Ery" = "#A9A9A9",        # Slightly darker gray for Ery
  "Imm.B.cell" = "gray"     # Other
)
Fig1.1 <- ggplot(merged_data[tissue != "leukemia"], aes(x = UMAP_1, y = UMAP_2)) + 
  
  geom_point(aes(color = functional.cluster), size = 0.00000001 ) + 
  
  geom_text_repel(data = pDT.labels %>%
                    filter(functional.cluster %in% c("Mono", "Eo/Ba", "GMP", "MEP (early)",
                                                     "MkP", "Gran. P", "Gran.", "HSC"
                    )),
                  aes(x = hex.x, y = hex.y, label = functional.cluster),
                  size = 2,                  # Adjust text size
                  box.padding = 0.21,         # Distance from points
                  point.padding = 0.21,       # Distance from label anchor
                  segment.color = "black",   # Line color
                  segment.size = 0.004,        # Line thickness
                  force = 10,                # Repelling force
                  max.overlaps = Inf) +
  facet_grid(cols = vars(tissue)) + 
  # Defining color manual scale for clusters
  # Defining color manual scale for clusters
  scale_color_manual(name = "Celltype",
                     values = cluster_colors,
                     guide = guide_legend(override.aes = list(size = 3))) + 
  #labs(title = "Cell type composition") +
  # Adjusting the alpha scale for fraction
  #scale_alpha_continuous(name = "Fraction", range = c(0, 1)) +  # To use fraction values for transparency
  labs(x =" UMAP 1", y = "UMAP 2")+
  # Adjusting the theme
  optimized_theme_fig() +
  
  # Positioning legends separately
  # Positioning legends separately
  theme(
    legend.position = "right",  # Color legend at the bottom
    legend.box = "horizontal",   # Horizontal alignment of legends
    legend.text = element_text(size = 5),     # Adjust legend text font size
    legend.title = element_blank(), # Remove title for color legend
    #legend.spacing = unit(0.5, "cm"),  # Adjust spacing between legends
    #legend.key.size = unit(4, "lines")  # Adjust size of the legend keys
  ) 
# Correct axis labels (you can define xu and yu separately if needed)
ggsave(outdir("UMAP_InvivoX_NTC.png"),Fig1.1,dpi=300, w=4, h=2, units = "in")

#fig 1.2 ---------------
limmaRes_NTC <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
head(limmaRes_NTC)
dataVoom_NTC_in_ex <- read_rds(basedir("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(basedir("NTC_meta.rds"))

top_genes <- limmaRes_NTC[limmaRes_NTC$genes %in% c("Idi1","Oas2","Msmo1"),] %>%
  unique()


#fig.1.3-----

filtered_data <- limmaRes_NTC %>%
  filter(group != "n.s")%>%
  filter(abs(logFC)>1)%>%
  mutate(Regulation = group)

# Count the number of up and down genes for each cell type
gene_counts <- filtered_data %>%
  group_by(celltype,Regulation) %>%
  summarise(count = n()) %>%
  ungroup()

# Log10-transform the counts, adding 1 to avoid log(0)
gene_counts <- gene_counts %>%
  mutate(log10_count = log10(count + 1))
#fig 1.1 ---------------
limmaRes_NTC <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
dataVoom_NTC_in_ex <- read_rds(basedir("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(basedir("NTC_meta.rds"))

top_genes <- limmaRes_NTC[limmaRes_NTC$genes %in% c("Idi1","Msmo1","Iigp1","Oas2"),] %>%
  unique()
#fig-----
gene_counts$celltype <- factor(gene_counts$celltype,
                                   levels =c("HSC","MEP.early","MkP" ,
                                             "GMP", "Gran.P", "Gran.", "Mono","Eo.Ba" ))
gene_counts <-gene_counts%>%
  mutate(Regulation = recode(Regulation,
                           "down" = "downregulation ex vivo",
                           "up" = "upregulation ex vivo"))
# Fig
Fig1.3 <- ggplot(gene_counts,aes(y = celltype,
                                 x = ifelse(Regulation == "down",
                                            -log10(count), log10(count)), 
                                 fill = Regulation)) + geom_col() +
  scale_fill_manual(values = c("downregulation ex vivo" =  "#4C889C", "upregulation ex vivo" = "#D0154E")) +
  
  labs(y = "Cell type",
       x = TeX("log_{10}(No. of genes+1)"),
       title = paste("No.of DEGs")) +
  
  coord_flip()+
  optimized_theme_fig()+
  theme(legend.position = "right") 
Fig1.3
ggsave(outdir(paste0("Fig1.3.pdf")), plot=Fig1.3,w=6,h=5, units = "cm")
#Figure 1.4
################################################################################
#FGSEA
celltype_order <- c("HSC","MEP.early","MkP" ,"GMP", "Gran.P", "Gran.", "Mono","Eo.Ba" )
InDir1 <- dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/")
gsea.res <- read_rds(InDir1("NTC_fgsea.rds"))
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge,
                                      function(vec) paste(vec[1:10], collapse = ","))
dbx<-"MSigDB_Hallmark_2020"
pDT <- gsea.res[db == dbx]
## Splitting the task to handle both ends of the NES spectrum-positive and negative
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5),
                                                       by=c("celltype")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5),
                                                      by=c("celltype")]$pathway)

# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(c(pw.display.pos, pw.display.neg))
pDT <- pDT[pathway %in% pw.display]

if (nrow(pDT) > 0){
  # Step 1: Aggregate NES values across all celltypes (mean of NES per pathway)
  pDT_agg <- pDT %>%
    group_by(pathway) %>%
    summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
    arrange(desc(average_NES))  # Ordering pathways by the average NES, highest first
  
  # Step 2: Create a factor for pathway that reflects the aggregated NES order
  pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)
  
  
  
  pDT$celltype <- factor(pDT$celltype,
                         levels =c("HSC","MEP.early","MkP" ,
                                   "GMP", "Gran.P", "Gran.", "Mono","Eo.Ba" ))}
  
  # Step 3: Plot with the new pathway order (highest NES first)
Fig1.4 <- ggplot(pDT, aes(y=celltype, x=pathway, color=pmin(2,NES), size=pmin(5, -log10(padj)))) +
         
         scale_color_gradient2(low = "#4C889C",
                               mid = "white",
                               high = "#D0154E"
         )+#,
         #name=TeX("log_{2}(FC)"))+
         geom_point() +
         scale_size_continuous(
           range = c(0, 2),
           limits = c(0, 5),
           name=TeX("$-\\log_{10}(p_{adj})$"))+
         
         
         #xRot() +
         #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
         labs(y = "celltype",
              x = "Pathways",
              title = "Enriched pathways",
              size = "log10(padj)")+
    #coord_flip()+
  optimized_theme_fig()+
  theme(axis.text.x = element_text(angle = 45))
  
  

Fig1.4
ggsave(outdir("Fig1.4fGSEA_plot_",dbx,"._long.pdf"),plot = Fig1.4, w=11,h=6, units = "cm")

################################################################################
row2 <-  Fig1.3 + Fig1.4 + plot_layout(widths = c(1, 4))
#######################################
#ISG and Cholesterol
# Define your pathways of interest
pathways <- list(
  ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
    filter(L1=="ISG_Core")%>%pull(value),
  mTORC1_or_Cholesterol = union(enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
                                enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`),
  ROC = enr.terms$MSigDB_Hallmark_2020$`Reactive Oxygen Species Pathway`,
  
  #Cholesterol = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  Hypoxia = enr.terms$MSigDB_Hallmark_2020$`Hypoxia`,
  Glycolysis = enr.terms$MSigDB_Hallmark_2020$`Glycolysis`,
  Protein_loc = Reduce(union, list(
    enr.terms$MSigDB_Hallmark_2020$`Myc Targets V1`,
    enr.terms$GO_Biological_Process_2023$`Protein Import (GO:0017038)`,
    enr.terms$GO_Biological_Process_2023$`Protein Insertion Into ER Membrane (GO:0045048)`,
    enr.terms$GO_Biological_Process_2023$`Protein Insertion Into Membrane (GO:0051205)`,
    enr.terms$GO_Biological_Process_2021$`RNA biosynthetic process (GO:0032774)`
  ))
)

# Define the function to get top genes for a given pathway
#here top 5 genes
#IFN genes

get_top_genes_up <- function(pathway_name,
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
    arrange(desc(logFC)) %>%
    slice_head(n = top_n) %>%
    pull(genes)%>%unique()
  
  # Extract and return data for the filtered genes
  pathway_plot <- limma_results %>%
    filter(toupper(genes) %in% toupper(filtered_genes))
  #%>%
  #filter(group != "n.s") 
  
  return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
}
get_top_genes_down <- function(pathway_name,
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
    arrange(logFC) %>%
    slice_head(n = top_n) %>%
    pull(genes)%>%unique()
  
  # Extract and return data for the filtered genes
  pathway_plot <- limma_results %>%
    filter(toupper(genes) %in% toupper(filtered_genes))
  #%>%
  #filter(group != "n.s") 
  
  return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
}


# Initialize an empty list to store results

results_up <- map(names(pathways), function(pathway) {
  get_top_genes_up(
    pathway_name = pathway,
    pathway_genes = pathways[[pathway]],
    limma_results = limmaRes_NTC
  )
})
results_down <- map(names(pathways), function(pathway) {
  get_top_genes_down(
    pathway_name = pathway,
    pathway_genes = pathways[[pathway]],
    limma_results = limmaRes_NTC
  )
})
# Optionally, name the elements of the list using pathway names
names(results_up) <- names(pathways)
names(results_down) <- names(pathways)

# results <- map(names(pathways), function(pathway) {
#   get_top_genes(
#     pathway_name = pathway,
#     pathway_genes = pathways[[pathway]],
#     limma_results = limmaRes_NTC
#   )
# })
# # Optionally, name the elements of the list using pathway names
# names(results) <- names(pathways)

# Access the results for each pathway
combined_genes_filtered <- bind_rows(
  results_up$mTORC1_or_Cholesterol$pathway_plot %>% mutate(gene_set = "mTORC1/Cholesterol"),
  #results$Cholesterol$pathway_plot%>%mutate(gene_set = "Cholesterol"),
  results_down$ISG_core$pathway_plot %>% mutate(gene_set = "ISG core")
  # results$ROC$pathway_plot%>%mutate(gene_set = "ROC")
  #results$Glycolysis$pathway_plot%>%mutate(gene_set = "Glycolysis")
)
# 
# genes_fig1 <- map_df(names(pathways), function(pathway) {
#   res <- get_top_genes(
#     pathway_name = pathway,
#     pathway_genes = pathways[[pathway]],
#     limma_results = limmaRes_NTC
#   )
#   
#   # Convert top genes to a dataframe and add pathway name
#   data.frame(
#     pathway = pathway,
#     genes = res$top_genes
#   )
# }) %>% filter(pathway %in% c("mTORC1_or_Cholesterol","ISG_core"))
genes_fig1 <- combined_genes_filtered %>%
  dplyr::select(genes, gene_set)%>%
  distinct()
colnames(genes_fig1) <- c("genes","pathway")
#add some genes of interest
mTORC1_or_Cholesterol <- c("Idi1", "Cyp51","Stard4", "Scd2","Mthfd2","Sqle",
                           "Fads2","Dhcr24",
                           "Hmgcs1","Ldlr","Plscr1","Acat1",
                           "Acat2","Msmo1")
cholesterol_df <- data.frame(
  # Cholesterol genes
  pathway = rep("mTORC1_or_Cholesterol", length(mTORC1_or_Cholesterol)),
  genes = mTORC1_or_Cholesterol # Cholesterol pathway
)
#%>%
genes_fig1 <- rbind(genes_fig1, cholesterol_df) %>% distinct()
genes_fig1 %>% write_rds(outdir("genes_fig1.rds"))

#fig 1.1 ---------------
limmaRes_NTC <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
dataVoom_NTC_in_ex <- read_rds(basedir("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(basedir("NTC_meta.rds"))


#fig-----

combined_genes_filtered$celltype <- factor(combined_genes_filtered$celltype,
                                   levels =c("HSC","MEP.early","MkP" ,
                                             "GMP", "Gran.P", 
                                             "Gran.", "Mono","Eo.Ba" ))
# Plotting the dot plot
Fig1.5 <- ggplot(combined_genes_filtered, aes(x = celltype, y = genes,
                                              color = pmin(3, pmax(-3, logFC)),
                                              size = pmin(5,-log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E")+#,
  #name=TeX("log_{2}(FC)"))+
  geom_point() +
  scale_size_continuous(
    range = c(0, 2),
    #breaks = c(0,2,5,10),
    #limits = c(0, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  labs(title = "Divergent Regulation of ISGs vs growth/metabolic genes",
       y = "Genes",
       x = "Cell Type",
       color = "logFC",
       size = "-log10(padj)") +
  facet_grid(rows = vars(gene_set), scales = "free", space = "free") +
  #coord_flip() +
  optimized_theme_fig()

Fig1.5
ggsave(outdir("Fig1.5.png"), plot = Fig1.5,dpi = 600,width = 1217, height = 5.5, units = "cm")
ggsave(outdir("Fig1.5_long.png"), plot = Fig1.5,dpi = 600,width = 126, height = 12, units = "cm")
###################################################
#Figure 1.6
example <- c("Idi1","Oas2","Tbp")
example <- intersect(example, rownames(dataVoom_NTC_in_ex$E))
dat.list<-list()
for(gg in unique(example)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta_in_ex %>%
    mutate(E = dataVoom_NTC_in_ex$E[gg,]) %>%
    rownames_to_column("samples") %>%
    remove_rownames()
  
  dat.list[[gg]] <- gene_data
}
dat.list <- bind_rows(dat.list,.id="gene")
head(dat.list)
######################################

create_gene_plots_NTC <- function(data, geneset, remove_guides = FALSE) {
  # Create the base plot with gene as a facet row and cell types as facet columns
  plot <- ggplot(data, aes(x = celltype, y = E, fill = tissue)) + 
    geom_boxplot(
      outlier.shape = NA,
      position = position_dodge(width = 0.8),  # Fixed dodge width
      color = "black",
      size = 0.2
    ) + 
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.3,     # Reduced jitter width
        dodge.width = 0.8       # Match dodge width with boxplot
      ), 
      alpha = 0.3,
      size = 0.8
    ) +
    facet_grid(rows = vars(gene), cols = vars(celltype), space = "free", scales = "free_x") +  
    labs(
      legend = "Experimental model",
      title = "Representative gene expression patterns"
    ) +
    xlab(NULL) +
    ylab("Scaled Gene Expression") +
    scale_fill_manual(
      values = c("ex.vivo" = "#6F987B", "in.vivo" = "#764BAB"),
      name = "Experimental model"
    ) + 
    theme_bw() +
    optimized_theme_fig() +  # Custom theme assumed to be defined elsewhere
    theme(
      axis.text.x = element_blank(),
      panel.spacing = unit(0.1, "lines")
    )
  
  # Optionally remove the guides/legends
  if (remove_guides) {
    plot <- plot + theme(legend.position = "none")
  }
  
  # Add Wilcoxon rank-sum test significance
  plot <- plot + stat_compare_means(
    aes(group = tissue), 
    method = "wilcox.test", 
    label = "p.signif",
    label.y = max(data$E, na.rm = TRUE) * 0.8
  )
  
  return(plot)
}

# Initialize empty list to store individual plots
all_gene_plots <- list()

# Filter data for selected genes and cell types
combined_data <- data.frame()

for (gg in example) {
  data <- dat.list %>% 
    filter(gene == gg) %>%
    filter(celltype %in% c("Eo.Ba", "HSC", "MkP", "Mono", "GMP", "Gran.P"))
  
  if (nrow(data) == 0) {
    warning(paste("No data available for gene:", gg))
    next
  }
  
  data$gene <- gg
  combined_data <- rbind(combined_data, data)
}

# Custom labels for genes
custom_labels <- c(
  "Oas2" = "Oas2 (ISG)",
  "Idi1" = "Idi1 (Cholesterol)",
  "Tbp" = "Tbp (Housekeeping)"
)

# Apply custom gene labels
combined_data <- combined_data %>%
  mutate(gene = dplyr::recode(gene, !!!custom_labels))

# Set the factor levels for cell type ordering
combined_data$celltype <- factor(combined_data$celltype,
                                 levels = c("HSC", "MEP.early", "MkP", 
                                            "GMP", "Gran.P", "Gran.", 
                                            "Mono", "Eo.Ba"))

# Generate the figure with guides
Fig1.6_with_guides <- create_gene_plots_NTC(combined_data, "example_NTC", remove_guides = FALSE)
# Create the plot without guides
Fig1.6_without_guides <- create_gene_plots_NTC(combined_data, "example_NTC", remove_guides = TRUE)

# Save both plots
ggsave(outdir("Fig1.6_with_guides.pdf"), Fig1.6_with_guides,dpi = 600,,width = 12, height = 8, units = "cm")
ggsave(outdir("Fig1.6_without_guides.pdf"), Fig1.6_without_guides, ,dpi = 600,width = 12, height = 11, units = "cm")

##########################
row3 <- Fig1.5 + Fig1.6_with_guides + plot_layout(widths = c(1, 3))

################

# Combine the three plots in a single column with no spacing
combined <- row1 / row2 / row3 + plot_layout(ncol = 1, heights = c(1, 0.5, 1.5)) & theme(plot.margin = margin(0, 0, 0, 0))
combined
ggsave(outdir("trial_combined.png"), w = 18.3, h = 24.7, units = "cm")

#####################
#extra
#Fig.1_supplementary 2
Indir7 <-  dirout("Ag_ScRNA_19_invivo_exvivo_izzo_zscore/")
longer_dataVoom_zscore <- read_rds(basedir("zscore_plot_izzo.rds"))
# View the resulting dataframe with z-scores

unique(longer_dataVoom_zscore$tissue)
# Create the plot


# Cap zscore at ±2
longer_dataVoom_zscore <- longer_dataVoom_zscore %>%
  mutate(zscore_capped = pmin(zscore, 2))
unique(longer_dataVoom_zscore$sample)
ggplot(longer_dataVoom_zscore, aes(x = gsub("_"," ",sample), y = genes , fill = zscore_capped)) +
  geom_tile(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) +  # Scatter plot with jitter for better visualization
  facet_grid(cols = vars(tissue_celltype), scales = "free", space = "free") +  # Separate by tissue, free_x ensures that each tissue has its own x-axis range
  labs(title = "Z-score of Gene Expression Across Tissues", 
       x = "Genes", 
       y = "Z-score") +
  scale_fill_gradient2(low = "#4C889C", high = "#D0154E", mid = "white", midpoint = 0) +
  optimized_theme_fig()+
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90,hjust = 0))+
  theme(panel.spacing = unit(0.2, "lines")) # Rotate x-axis labels for readability
# Remove minor gridlines

ggsave(outdir(paste0("Supplementary.Fig.2B.pdf")), w=18,
       h=8, units = "cm")


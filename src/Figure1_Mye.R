#figure 1
###############
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification_Mye.R")
source("src/Ag_enrichR_mouse_genes.R")
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

##################################################################################

basedir <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide_Mye")

out <- "Figure1_Mye"
outdir <- dirout("Figure1_Mye")

#load data

#Fig1A-----------
InDir1 <- dirout("FIG_02_scRNA_UMAPs_ar_Mye/")
pDT.labels <- read_rds(InDir1("pDT.labels.rds"))
# color coding
cluster_colors <- c(
  "Mono" = "#E69F00",      # Orange
  "Eo.Ba" = "#56B4E9",     # Sky Blue
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
cluster_labels <- c(
  "Mono" = "Mono: Monocytes",
  "Eo.Ba" = "Eo.Ba: Eosinophil/Basophil",
  "GMP" = "GMP: Granulocyte-Macrophage Progenitor",
  "MEP (early)" = "MEP (early): Early Mega-Erythroid Progenitor",
  "MkP" = "MkP: Megakaryocyte Progenitor",
  "Gran. P" = "Gran. P: Granulocyte Progenitor",
  "Gran." = "Gran.: Granulocyte",
  "HSC" = "HSC: Hematopoietic Stem Cell",
  "unclear" = "Unclear",
  "CLP" = "CLP: Common Lymphoid Progenitor",
  "MEP" = "MEP: Mega-Erythroid Progenitor",
  "Ery" = "Ery: Erythrocytes",
  "Imm.B.cell" = "Imm.B.cell: Immature B-Cell"
)

merged_data <- read_rds(InDir1("Cross_projected_on_in.vivo.rds"))
merged_data$functional.cluster <- gsub("Eo/Ba","Eo.Ba",merged_data$functional.cluster)
merged_data$functional.cluster <- factor(merged_data$functional.cluster, 
                                         levels = names(cluster_colors))

exclude <- merged_data[is.na(functional.cluster),]
merged_data <- merged_data %>%
filter(!(sample.x %in% exclude$sample.x))
Fig1A <- ggplot(merged_data[tissue != "leukemia"], aes(x = UMAP_1, y = UMAP_2)) + 
  
  geom_point(aes(color = functional.cluster), size = 0.00000001 ) + 
  
  geom_text_repel(data = pDT.labels %>%
                    filter(functional.cluster %in% c("Mono", "Eo.Ba", "GMP", "MEP (early)",
                                                     "MkP", "Gran. P", "Gran.", "HSC"
                    )),
                  aes(x = hex.x, y = hex.y, label = functional.cluster),
                  size = 2,                  # Adjust text size
                  box.padding = 0.21,         # Distance from points
                  point.padding = 0.21,       # Distance from label anchor
                  segment.color = "black",   # Line color
                  segment.size = 0.004,        # Line thickness
                  force = 20,                # Repelling force
                  max.overlaps = Inf) +
  facet_grid(cols = vars(tissue),
             labeller = labeller(tissue = c("ex.vivo" = "Ex vivo", "in.vivo" = "In vivo"))) + 
 
  scale_color_manual(name = "Celltype",
                     values = cluster_colors,
                     labels = cluster_labels,
                     guide = guide_legend(override.aes = list(size = 3))) + 
  labs(title = "Cell type composition") +
  # Adjusting the alpha scale for fraction
 
  labs(x =" UMAP 1", y = "UMAP 2")+
  # Adjusting the theme
  optimized_theme_fig() +
  theme(
    legend.position = "right",  # Color legend at the bottom
    legend.box = "horizontal",   # Horizontal alignment of legends
    legend.text = element_text(size = 5),     # Adjust legend text font size
    legend.title = element_blank(), # Remove title for color legend
    
  ) +
  optimized_theme_fig()+
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Fig1A
ggsave(outdir("Fig1A.png"),Fig1A,dpi=500, w=10.5, h=5, units = "cm")
ggsave(outdir("Fig1A.pdf"),Fig1A,dpi=300, w=10.5, h= 5, units = "cm")
#Fig1B ---------------
# Ensure Regulation is correctly factored and recoded
ex_in_NTC_per_ct <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
filtered_data <- ex_in_NTC_per_ct %>%
  mutate(Regulation = group) %>%
  filter(group != "n.s")%>%
  filter(abs(logFC)>1)

# Count the number of up and down genes for each cell type
gene_counts <- filtered_data %>%
  group_by(celltype, Regulation) %>%
  summarise(count = n()) %>%
  ungroup()

# Log10-transform the counts, adding 1 to avoid log(0)
gene_counts <- gene_counts %>%
  mutate(log10_count = log10(count + 1))
gene_counts <- gene_counts %>%
  mutate(Regulation = recode(Regulation,
                             "down" = "Higher expr in vivo",
                             "up" = "Higher expr ex vivo"))

# Convert to bidirectional values based on Regulation
gene_counts <- gene_counts %>%
  mutate(x_val = ifelse(Regulation == "Higher expr in vivo", -log10(count + 1), log10(count + 1)))

# Plot
Fig1B <- ggplot(gene_counts, aes(
  y = celltype,
  x = x_val,
  fill = Regulation
)) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "Higher expr in vivo" = "#4C889C",
      "Higher expr ex vivo" = "#D0154E"
    )
  ) +
  labs(
    y = "Cell type",
    x = expression(atop("Number of genes", 
                        paste(log[10](n + 1)))),
    title = "No. of DEGs"
  ) +
  coord_flip() +
  optimized_theme_fig() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

Fig1B

# Save
ggsave(outdir("Fig1B.pdf"), plot = Fig1B, width = 4, height = 5, units = "cm")
################################################################################
##Fig1C-------------
celltype_order <- c("HSC","MEP.early","MkP" ,"GMP", "Gran.P", "Gran.", "Mono","Eo.Ba" )
InDir1 <- dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide_Mye/")
gsea.res <- read_rds(InDir1("NTC_fgsea.rds"))
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge,
                                      function(vec) paste(vec[1:10], collapse = ","))
dbx<-"MSigDB_Hallmark_2020"
pDT <- gsea.res[db == dbx]
## Splitting the task to handle both ends of the NES spectrum-positive and negative
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=4),
                                                       by=c("celltype")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=4),
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
Fig1C <- ggplot(pDT, aes(y=celltype, x=pathway, color = pmin(pmax(NES, -2), 2), size=pmin(5, -log10(padj)))) +
         
         scale_color_gradient2(low = "#4C889C",
                               mid = "white",
                               high = "#D0154E",
                               name=TeX("NES"))+
         #name=TeX("log_{2}(FC)"))+
         geom_point() +
         scale_size_continuous(
           range = c(0, 1.8),
           #limits = c(0, 5),
           name=TeX("$-\\log_{10}(p_{adj})$"))+
         
         
         #xRot() +
         #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
         labs(y = "Cell type",
              x = "Pathways",
              title = "Enriched pathways")+
              
    #coord_flip()+
  optimized_theme_fig()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,
                                   ),
        legend.position = "right", legend.direction = "vertical",
        legend.justification = "bottom")
  
Fig1C
ggsave(outdir("Fig1C.pdf"),plot = Fig1C, w = 11,h = 5.5, units = "cm")

# Fig1C_vert <- ggplot(pDT, aes(x=celltype, y=pathway, color=pmin(2,NES), size=pmin(5, -log10(padj)))) +
#   
#   scale_color_gradient2(low = "#4C889C",
#                         mid = "white",
#                         high = "#D0154E",
#                         name=TeX("NES"))+
#   #name=TeX("log_{2}(FC)"))+
#   geom_point() +
#   scale_size_continuous(
#     range = c(0, 1.8),
#     limits = c(0, 5),
#     name=TeX("$-\\log_{10}(p_{adj})$"))+
#   
#   
#   #xRot() +
#   #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
#   labs(y = "Cell type",
#        x = "Pathways",
#        title = "Enriched pathways")+
#   
#   #coord_flip()+
#   optimized_theme_fig()+
#   theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
#         legend.position = "right", legend.direction = "vertical",
#         legend.justification = "bottom")
# 
# 
# ggsave(outdir("Fig1C_vert.pdf"),plot = Fig1C_vert, w = 7.5, h = 10, units = "cm")

#Fig1D-------------

#ISG and Cholesterol
# Define your pathways of interest
limmaRes_NTC <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
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
  
 
  pathway_plot <- limma_results %>%
    filter(toupper(genes) %in% toupper(filtered_genes))
  
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

# Access the results for each pathway
combined_genes_filtered <- bind_rows(
  results_up$mTORC1_or_Cholesterol$pathway_plot %>% mutate(gene_set = "mTORC1/Cholesterol"),
  results_down$ISG_core$pathway_plot %>% mutate(gene_set = "ISG core")

)

genes_fig1 <- combined_genes_filtered %>%
  dplyr::select(genes, gene_set)%>%
  distinct()
colnames(genes_fig1) <- c("genes","pathway")
#add some genes of interest
mTORC1_or_Cholesterol <- c("Idi1", "Cyp51","Stard4", "Scd2","Mthfd2","Sqle",
                           "Fads2","Dhcr24",
                           "Hmgcs1","Ldlr","Plscr1","Acat1",
                           "Acat2","Msmo1","Myc")
cholesterol_df <- data.frame(
  # Cholesterol genes
  pathway = rep("mTORC1_or_Cholesterol", length(mTORC1_or_Cholesterol)),
  genes = mTORC1_or_Cholesterol # Cholesterol pathway
)
#%>%
genes_fig1 <- rbind(genes_fig1, cholesterol_df) %>% distinct()
genes_fig1 %>% write_rds(outdir("genes_fig1.rds"))



combined_genes_filtered$celltype <- factor(combined_genes_filtered$celltype,
                                   levels =c("HSC","MEP.early","MkP" ,
                                             "GMP", "Gran.P", 
                                             "Gran.", "Mono","Eo.Ba" ))
# Plotting the dot plot
Fig1D <- ggplot(combined_genes_filtered, aes(x = celltype, y = genes,
                                              color = pmin(3, pmax(-3, logFC)),
                                              size = pmin(5,-log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
  name=TeX("log_{2}(FC)"))+
  geom_point() +
  scale_size_continuous(
    range = c(0, 1.8),
    #breaks = c(0,2,5,10),
    #limits = c(0, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  labs(title = "Downregulation of ISGs and upregulation
       of growth/metabolic genes in ex vivo ",
       y = "Genes",
       x = "Cell Type") +
  facet_grid(rows = vars(gene_set), scales = "free", space = "free") +
 # coord_flip() +
  optimized_theme_fig()

Fig1D
ggsave(outdir("Fig1D_horizontal.pdf"), plot = Fig1D, height = 4.8,width = 13.5, units = "cm")
ggsave(outdir("Fig1D.pdf"), plot = Fig1D, height = 12,width = 6, units = "cm")
#ggsave(outdir("Fig1D_long.png"), plot = Fig1D,dpi = 600,width = 126, height = 12, units = "cm")
###################################################
#Fig1E ---------------

dataVoom_NTC_in_ex <- read_rds(basedir("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(basedir("NTC_meta.rds"))

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
#function
create_gene_plots_NTC <- function(data, geneset, remove_guides = FALSE) {
  # Create the base plot with gene as a facet row and cell types as facet columns
  plot <- ggplot(data, aes(x = celltype, y = E, color = tissue, group = tissue)) + 
    geom_boxplot(
      fill = NA,  # Hollow boxplots
      outlier.shape = NA,
      position = position_dodge(width = 0.8),  # Match jitter
      size = 0.2
    ) + 
    # geom_jitter(
    #   position = position_jitterdodge(
    #     jitter.width = 0.3,
    #     dodge.width = 0.8
    #   ), 
    #   alpha = 0.3,
    #   size = 0.5,
    #   show.legend = FALSE
    # )+
    facet_grid(rows = vars(gene), cols = vars(celltype), space = "free_x",
               scales = "free",
               labeller = labeller(gene = label_wrap_gen(width = 18))) +  
    labs(
      legend = "Experimental model",
      title = "Representative gene expression patterns"
    ) +
    xlab(NULL) +
    ylab("Scaled Gene Expression") +
    scale_color_manual(
      values = c("ex.vivo" = "#6a3d9aff", "in.vivo" = "#d38d5fff"),
      name = "Experimental model",
      labels = c("ex.vivo" = "Ex vivo", "in.vivo" = "In vivo")
    ) +
    optimized_theme_fig() +
    theme(
      axis.text.x = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x = NULL
      
    )+
   theme( panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank())
  
  
  # Optionally remove the guides/legends
  if (remove_guides) {
    plot <- plot + theme(legend.position = "none")
  }
  
  # Add Wilcoxon rank-sum test significance
  plot <- plot + stat_compare_means(
    aes(group = tissue), 
    method = "wilcox.test", 
    label = "p.signif",
    label.y = max(data$E, na.rm = TRUE) * 0.5
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
Fig1E_with_wo_jitter_guides <- create_gene_plots_NTC(combined_data, "example_NTC", remove_guides = FALSE)
# Create the plot without guides
Fig1E_without_wo_jitter_guides <- create_gene_plots_NTC(combined_data, "example_NTC", remove_guides = TRUE)

# Save both plots
ggsave(outdir("Fig1E_with_wo_jitter_guides.pdf"), Fig1E_with_wo_jitter_guides,
       dpi = 600,width = 11, height = 6, units = "cm")



#combine----------------
first_row <- (Fig1A | plot_spacer() | Fig1B) +
  plot_layout(widths = c(1.8, 0.7, 1.5)) & 
  theme(plot.margin = margin(0, 0, 0, 0))

ggsave(outdir("first_row.pdf"),first_row, w= 18, h = 6.5, units = "cm")


# Left column: Fig1C on top, spacer below (50/50)
left <- Fig1C / plot_spacer() +
  plot_layout(heights = c(1, 2.5))  # Equal halves

# Right column: Fig1D spans full height
right <- Fig1D

second_row <- (left | right) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

ggsave(outdir("second_row.pdf"),second_row, w = 18, h = 11, units = "cm")


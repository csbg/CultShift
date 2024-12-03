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

##################################################################################
base<-"Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/"
basedir<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
out <- "Figure1"
outdir <- dirout("Figure1")
source("src/Ag_Optimized_theme.R")
########################
ENRICHR<-dirout(paste0("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/ENRICHR"))
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

#fig 1.1 ---------------
limmaRes_NTC <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
dataVoom_NTC_in_ex <- read_rds(basedir("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(basedir("NTC_meta.rds"))

top_genes <- limmaRes_NTC[limmaRes_NTC$genes %in% c("Idi1","Msmo1","Iigp1","Oas2"),] %>%
  unique()
#fig-----
Fig1.1 <- ggplot() +
  # Hexbin plot for the "others" group
  stat_bin_hex(data = filter(limmaRes_NTC, group == "n.s"), 
               aes(x = logFC, y = -log10(adj.P.Val), fill = ..count..), 
               bins = 20, color = NA, alpha = 0.7) +
  scale_fill_gradient(low = "lightgrey", high = "black",
                      limits = c(1, 5000), name = "Gene Count") +
  
  # Overlay points for the main groups
  geom_point(data = filter(limmaRes_NTC, group != "n.s"), 
             aes(x = logFC, y = -log10(adj.P.Val), color = group), 
             alpha = 0.5, size = 1.5) +
  
  # Add text labels for top genes with ggrepel
  geom_text_repel(
    data = top_genes,
    aes(x = logFC, y = -log10(adj.P.Val), label = genes),
    size = 3, 
    color = "black",
    max.overlaps = Inf,  # Ensure no labels are omitted
    box.padding = 0.4,
    point.padding = 0.4,
    segment.color = 'black',  # Color for the line pointing to the gene
    segment.size = 0.4,       # Thickness of the line
    force = 2,                # Increase the repelling force to reduce overlap
    force_pull = 0.4,         # Pull force for the label towards the point
    min.segment.length = 0,   # Ensure segments are always drawn, even for close points
    arrow = arrow(length = unit(0.02, "npc"), type = "closed", angle = 15)  # Add arrows/lines pointing to the gene
  )+
  
  # Manually setting colors for groups
  scale_color_manual(values = c(
    "up" = "#D0154E", 
    "down" = "#4C889C"
  ), name = "Group") +
  
  labs(title = "Ex-vivo vs in-vivo differentially expressed genes",
       x = "logFC",
       y = "-log10(adj.P)") +
  facet_wrap(vars(celltype), scales = "free") +
  theme_bw() +
  optimized_theme_fig()
Fig1.1
ggsave(outdir(paste0("Fig1.1.pdf")), plot=Fig1.1,w=12,h=9,units = "cm")
#fig1.2-----

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

# Fig
Fig1.2 <- ggplot(gene_counts,aes(y = celltype,
                                 x = ifelse(Regulation == "down",
                                            -log10(count), log10(count)), 
                                 fill = Regulation)) + geom_col() +
  scale_fill_manual(values = c("down" =  "#4C889C", "up" = "#D0154E")) +
  
  labs(y = "Interaction-KO",
       x = TeX("log_{10}(No. of genes+1)"),
       title = paste("DEGs")) +
  
  coord_flip()+
  optimized_theme_fig()+
  theme(legend.position = "right") 

ggsave(outdir(paste0("Fig1.2.pdf")), plot=Fig1.2,w=8,h=5, units = "cm")
#Figure 1.3
################################################################################
#FGSEA
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
pDT <- hierarch.ordering(pDT, "celltype", "pathway", "NES", TRUE)
Fig1.3 <- ggplot(pDT, aes(x=celltype, y=pathway, color=pmin(2,NES), size=pmin(5, -log10(padj)))) +
  
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E")+#,
  #name=TeX("log_{2}(FC)"))+
  geom_point() +
  scale_size_continuous(
    range = c(0, 3),
    limits = c(0, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  
  theme_bw()+
  xRot() +
  #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
  labs(x = "celltype",
       y = "Pathways",
       title = "Enriched pathways",
       size = "log10(padj)")+
  optimized_theme_fig()

Fig1.3
ggsave(outdir("Fig1.3fGSEA_plot_",dbx,".pdf"),plot = Fig1.3, w=8,h=8, units = "cm")

################################################################################

#######################################
#ISG and Cholesterol
# Define your pathways of interest
pathways <- list(
  ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
    filter(L1=="ISG_Core")%>%pull(value),
  mTORC1_or_Cholesterol = union(enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`,
                                enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`),
  ROC = enr.terms$MSigDB_Hallmark_2020$`Reactive Oxygen Species Pathway`,
  
  #Cholesterol = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  Hypoxia = enr.terms$MSigDB_Hallmark_2020$`Hypoxia`,
  Glycolysis = enr.terms$MSigDB_Hallmark_2020$`Glycolysis`
  # Add other pathways here as needed
)
# Define the function to get top genes for a given pathway
#here top 5 genes
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

results <- map(names(pathways), function(pathway) {
  get_top_genes(
    pathway_name = pathway,
    pathway_genes = pathways[[pathway]],
    limma_results = limmaRes_NTC
  )
})

# Optionally, name the elements of the list using pathway names
names(results) <- names(pathways)
# Access the results for each pathway
combined_genes_filtered <- bind_rows(
  results$mTORC1_or_Cholesterol$pathway_plot %>% mutate(gene_set = "mTORC1/Cholesterol"),
  #results$Cholesterol$pathway_plot%>%mutate(gene_set = "Cholesterol"),
  results$ISG_core$pathway_plot %>% mutate(gene_set = "ISG core")
  # results$ROC$pathway_plot%>%mutate(gene_set = "ROC")
  #results$Glycolysis$pathway_plot%>%mutate(gene_set = "Glycolysis")
)

genes_fig1 <- map_df(names(pathways), function(pathway) {
  res <- get_top_genes(
    pathway_name = pathway,
    pathway_genes = pathways[[pathway]],
    limma_results = limmaRes_NTC
  )
  
  # Convert top genes to a dataframe and add pathway name
  data.frame(
    pathway = pathway,
    genes = res$top_genes
  )
}) %>% filter(pathway %in% c("mTORC1_or_Cholesterol","ISG_core"))
#add some genes of interest
mTORC1_or_Cholesterol <- c("Idi1", "Cyp51","Stard4", "Scd2","Mthfd2","Sqle","Fads2","Dhcr24",
                           "Hmgcs1","Ldlr","Plscr1","Acat1",
                           "Acat2","Msmo1")
cholesterol_df <- data.frame(
  # Cholesterol genes
  pathway = rep("mTORC1_or_Cholesterol", length(mTORC1_or_Cholesterol)),
  genes = mTORC1_or_Cholesterol # Cholesterol pathway
)
#%>%
genes_fig1 <- rbind(genes_fig1, cholesterol_df) %>% distinct()
genes_fig1 %>%write_rds(outdir("genes_fig1.rds"))

# Plotting the dot plot
Fig1.4 <- ggplot(combined_genes_filtered, aes(x = celltype, y = genes,
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
  labs(title = "Differentially expressed genesets",
       x = "Cell Type",
       y = "Genes",
       color = "logFC",
       size = "-log10(padj)") +
  facet_grid(rows = vars(gene_set), scales = "free_y", space = "free") +
  theme_bw() +
  optimized_theme_fig()

Fig1.4
ggsave(outdir("Fig1.4.pdf"), plot = Fig1.4,width = 8, height = 12, units = "cm")
###################################################
#Figure 1.5
example <- c("Idi1","Oas2","Tbp")
example <- intersect(example, rownames(dataVoom_NTC_in_ex$E))
dat.list<-list()
for(gg in unique(example)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta_in_ex %>%
    mutate(E = scale(dataVoom_NTC_in_ex$E[gg,])) %>%
    rownames_to_column("samples") %>%
    remove_rownames()
  
  dat.list[[gg]] <- gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")
head(dat.list)
######################################
create_gene_plots_NTC <- function(data, geneset, remove_guides = FALSE) {
  # Create the base plot with gene as a facet row and cell types as facet columns
  plot <- ggplot(data, aes(x = celltype, y = E, fill = tissue)) + 
    geom_boxplot(
      outlier.shape = NA,
      position = position_dodge(width = 0.8),  # Adjust the width to create space between tissues
      color = "black",
      size = 0.2
    ) + 
    # Boxplot with tissue color fill
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.05, 
        dodge.width = 0.8  # Adjust the dodge width to match the boxplot
      ), 
      alpha = 0.4
    ) +  # Jittered points
    facet_grid(rows = vars(gene), cols = vars(celltype), space = "free", scales = "free_x") +  
    # Facet genes by row, cell types by column
    labs(
      legend = "Experimental model",
      title = "Gene expression ex-vivo vs in-vivo"
    ) +
    xlab(NULL) +  # Remove x-axis label
    ylab("Scaled Gene Expression") +
    scale_fill_manual(
      values = c("ex.vivo" = "#C1A0AC", "in.vivo" = "#87B1D6"),
      name = "Experimental model"
    ) + 
    theme_bw() +
    optimized_theme_fig() +  # Assuming `optimized_theme_fig()` is defined elsewhere
    theme(
      axis.text.x = element_blank(),
      panel.spacing = unit(0.1, "lines")  # Decrease space between facets
    ) + NULL
    #coord_fixed(ratio = 0.8)  # Adjust the aspect ratio to compress plot height or width
  
  # Optionally remove the guides/legends
  if (remove_guides) {
    plot <- plot + theme(legend.position = "none")  # Hide the legends
  }
  
  # Perform Wilcoxon rank-sum test and add significance stars
  plot <- plot + stat_compare_means(
    aes(group = tissue), 
    method = "wilcox.test", 
    label = "p.signif",  # Add significance stars
    label.y = max(data$E, na.rm = TRUE) * 0.8  # Adjust the position of the stars
  )
  
  return(plot)
}

# Initialize an empty list to store individual plots
all_gene_plots <- list()

# Loop over the example genes to combine the data
combined_data <- data.frame()  # Empty data frame to combine all data

for (gg in example) {
  data <- dat.list %>% 
    filter(gene == gg) %>%
    filter(celltype %in% c("Eo.Ba", "HSC", "MkP", "Mono", "GMP"))
  
  if (nrow(data) == 0) {
    warning(paste("No data available for gene:", gg))
    next  # Skip to the next gene if the data is empty
  }
  
  data$gene <- gg  # Add gene name as a column
  combined_data <- rbind(combined_data, data)  # Combine data for all genes
}
custom_labels <- c(
  "Oas2" = "Oas2 (ISG)",
  "Idi1" = "Idi1 (Cholesterol)",
  #"Brd9" = "Brd9",
  # Add more gene mappings as needed
  "Tbp" = "Tbp (Housekeeping)"  # Example for Tbp
)

# Step 2: Replace `gene` values with custom labels in `combined_data`
combined_data <- combined_data %>%
  mutate(gene = recode(gene, !!!custom_labels))

# Create the plot with guides
Fig1.5_with_guides <- create_gene_plots_NTC(combined_data, "example_NTC", remove_guides = FALSE)

# Create the plot without guides
Fig1.5_without_guides <- create_gene_plots_NTC(combined_data, "example_NTC", remove_guides = TRUE)

# Save both plots
ggsave(outdir("Fig1.5_with_guides.pdf"), Fig1.5_with_guides, width = 12, height = 10, units = "cm")
ggsave(outdir("Fig1.5_without_guides.pdf"), Fig1.5_without_guides, width = 11, height = 11, units = "cm")


combine_plots <- function(plot_list) {
  if (is.list(plot_list)) {
    return(wrap_plots(plot_list, ncol = 1))
  } else {
    return(plot_list)
  }
}

# Combine each plot into a single ggplot object
f1_1 <- combine_plots(Fig1.1)
f1_2 <- combine_plots(Fig1.2)
f1_3 <- combine_plots(Fig1.3)
f1_4 <- combine_plots(Fig1.4)
f1_5 <- combine_plots(Fig1.5_with_guides)

height_ratios <- c(1, 1.5)

# Define the width ratios for the first and second rows
width_ratios_row1 <- c(2.5, 0.8)   # Equal width for first row
width_ratios_row2 <- c(1,1, 1.8)   # Second plot in second row gets double width

# Create the first row layout
first_row <- wrap_plots(f1_1, f1_2, widths = width_ratios_row1)

# Create the second row layout
second_row <- wrap_plots(f1_3,f1_4, f1_5, widths = width_ratios_row2)

# Combine the two rows with the specified height ratios
combined_plot <- wrap_plots(first_row, second_row, ncol = 1) +
  plot_layout(heights = height_ratios)

# Save the combined plot to a PDF
ggsave(outdir("Figure1_combined1.pdf"), plot = combined_plot, width = 20, height = 20)
####################################
#
# Define the width ratios for the first row
width_ratios_row1 <- c(2)  # Adjusted for a small square Fig1.2
width_ratios_row2 <- c(0.8, 1,1.5)    # Adjusted for two plots in the second row

# Create the first row layout (3 plots)
first_row <- wrap_plots(f1_1,  widths = width_ratios_row1)

# Create the second row layout (2 plots)
second_row <- wrap_plots(f1_3,f1_4, f1_5, widths = width_ratios_row2)

# Combine both rows vertically using plot_layout
final_plot <- (first_row / second_row) + 
  plot_layout(heights = c(1, 1))  # Define height ratios

# Display the final plot
final_plot
ggsave(outdir("Figure1_combined1.pdf"), plot = final_plot, width = 20, height = 25)


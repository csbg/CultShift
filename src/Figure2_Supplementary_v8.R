source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
source("src/Ag_enrichR_mouse_genes.R")
library(tidyverse)
library(ggplot2)
library(dplyr)
library(latex2exp)
out <- "Figure2_Supplementary_v8"
basedir <- dirout("Figure2_Supplementary_v8")

Indir3 <- dirout("Ag_ScRNA_14_invivo_exvivo_external_zscore/")
Indir2 <- dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide")
Indir4 <- dirout("Ag_ScRNA_11_limma_all_ko_ex.vivo_vs_in.vivo_guide")
Indir5 <- dirout("Figure1")
#Supp_Fig1c----------------
##########
limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
dataVoom_NTC_in_ex <- read_rds(InDir_NTC("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(InDir_NTC("NTC_meta.rds"))
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

### ---------------------------------------------------------
### Load and clean FGSEA results
### ---------------------------------------------------------

gsea.res <- read_rds(Indir2("NTC_fgsea.rds"))
gsea.res[is.nan(NES), NES := 0]

gsea.res.export <- gsea.res[padj < 0.05][, -c("log2err", "size", "pval"), with = FALSE]
gsea.res.export$leadingEdge <- sapply(
  gsea.res.export$leadingEdge,
  function(vec) paste(vec[1:10], collapse = ",")
)

dbx <- "KEGG_2019_Mouse"

# Output directory
dat <- dirout(paste0(out, "FGSEA/", dbx))
output_file <- dat("GSEA_significant_", dbx, ".tsv")

if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

# Convert list columns to character before saving
df <- gsea.res[db == dbx]
df[] <- lapply(df, function(x) if (is.list(x)) sapply(x, toString) else x)

write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)


### ---------------------------------------------------------
### Select enriched pathways (top NES positive & negative)
### ---------------------------------------------------------

pDT <- df
pw.display.pos <- unique(
  pDT[padj < 0.05][order(-NES)][, head(.SD, 5), by = "celltype"]$pathway
)
pw.display.neg <- unique(
  pDT[padj < 0.05][order(NES)][, head(.SD, 5), by = "celltype"]$pathway
)

pw.display <- unique(c(pw.display.pos, pw.display.neg))
pDT <- pDT[pathway %in% pw.display]

if (nrow(pDT) == 0) stop("No pathways to display.")


### ---------------------------------------------------------
### Compute average NES per pathway for ordering
### ---------------------------------------------------------

pDT_agg <- pDT %>%
  group_by(pathway) %>%
  summarise(average_NES = mean(NES, na.rm = TRUE)) %>%
  arrange(desc(average_NES))

pDT <- pDT %>%
  mutate(pathway = factor(pathway, levels = pDT_agg$pathway)) %>%
  arrange(pathway)


### ---------------------------------------------------------
### Combine leading-edge genes per pathway
### ---------------------------------------------------------

combined_genes <- pDT %>%
  group_by(pathway) %>%
  summarise(Genes = list(unique(unlist(leadingEdge))), .groups = "drop")


### ---------------------------------------------------------
### Helper function: extract top DE genes per pathway
### ---------------------------------------------------------

get_top_genes <- function(pathway_name, limma_results,
                          logFC_threshold = 1,
                          pval_threshold = 0.05,
                          top_n = 5) {
  
  pathway_genes <- combined_genes %>%
    filter(pathway == pathway_name) %>%
    pull(Genes) %>%
    str_split(",") %>%                 # 🔥 split strings
    unlist() %>%
    str_trim() %>%                    # remove spaces
    unique()
  limma_results %>%
    filter(group != "n.s") %>%
    filter(adj.P.Val < pval_threshold) %>%
    filter(genes %in% pathway_genes) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    group_by(celltype) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(pathway = pathway_name)
}


### ---------------------------------------------------------
### Collect top genes across all pathways
### ---------------------------------------------------------

results <- map_dfr(unique(pDT$pathway), ~get_top_genes(.x, limmaRes_NTC)) %>%
  dplyr::select(genes, pathway) %>%
  distinct()


### ---------------------------------------------------------
### Select genes with interaction significance
### ---------------------------------------------------------

limmaRes <- read_rds(Indir4("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds")) %>%
  mutate(coef = gsub("interaction", "", coef))

genes_fig1 <- read_rds(Indir5("genes_fig1.rds"))

top_int_genes <- limmaRes %>%
  mutate(genes = ensg) %>%
  inner_join(results, by = "genes") %>%
  filter(group != "n.s") %>%
  filter(!genes %in% genes_fig1) %>%
  pull(genes) %>%
  unique()


### ---------------------------------------------------------
### Final gene selection for plotting
### ---------------------------------------------------------

top_genes <- limmaRes_NTC %>%
  inner_join(results, by = "genes") %>%
  filter(genes %in% top_int_genes) %>%
  filter(!genes %in% genes_fig1$genes) %>%
  group_by(genes) %>%
  filter(n_distinct(celltype[group != "n.s"]) >= 2) %>%
  ungroup() %>%
  group_by(pathway, genes) %>%
  slice_max(order_by = abs(logFC), n = 1) %>%
  ungroup() %>%
  group_by(pathway) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  dplyr::select(genes, pathway)


top_genes_NTC <- limmaRes_NTC %>%
  filter(genes %in% top_genes$genes) %>%
  inner_join(top_genes, by = "genes")


### ---------------------------------------------------------
### Pathway grouping and filtering
### ---------------------------------------------------------

top <- top_genes_NTC %>%
  mutate(geneset = case_when(
    pathway %in% c("Ribosome", "Ribosome biogenesis in eukaryotes") ~ "Ribosome machinery",
    pathway %in% c("DNA replication", "Cell cycle") ~ "Replication/ cell cycle",
    pathway == "Oxidative phosphorylation" ~ "Oxphos/ Electron transport",
    pathway == "Cell adhesion molecules (CAMs)" ~ "Cell adhesion molecules (CAMs)",
    TRUE ~ as.character(pathway)
  )) %>%
  group_by(genes) %>%
  filter(sum(group != "n.s") >= 2) %>%
  ungroup() %>%
  dplyr::select(-pathway)

# Only keep selected gene sets
desired_sets <- c(
  "Replication/ cell cycle",
  "Oxphos/ Electron transport",
  "Ribosome machinery",
  "Cell adhesion molecules (CAMs)"
)


top <- top %>%
  filter(geneset %in% desired_sets) %>%
  filter(!genes %in% c("Cdh1","Cd40","Mag","Pvr","Itgb7","Cox7a1","E2f2"))
#Adding ROS to the list-----------
ros_pw <- "Reactive Oxygen Species Pathway"
unique(grep("Reactive",gsea.res$pathway,value = T))
# 2. Pull all leading-edge genes (union)
ros_genes <- gsea.res %>%
  filter(pathway == "Reactive Oxygen Species Pathway") %>%
  pull(leadingEdge)%>%
  unlist() %>%
  unique()
### 1) Select ROS leading-edge genes
ros_pw <- "Reactive Oxygen Species Pathway"

ros_dt <- gsea.res %>%
  filter(pathway == ros_pw)

ros_genes <- ros_dt$leadingEdge %>%
  unlist() %>%
  unique()

### 2) Pick top 5 ROS genes per cell type
ros_top5_selection <- limmaRes_NTC %>%
  filter(
    genes %in% ros_genes,
    group != "n.s."
  ) %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  group_by(celltype) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(genes) %>%
  unique()

### 3) Retrieve ALL celltypes for only these selected genes
ros_top5_full <- limmaRes_NTC %>%
  filter(genes %in% ros_top5_selection) %>%
  mutate(
    pathway = ros_pw,
    geneset = "ROS pathway"
  )

ros_top5 <- ros_top5_full
################
# Append MYC manually-----------
top <- rbind(
  top,
  limmaRes_NTC %>%
    filter(genes == "Myc") %>%
    mutate(
      group = ifelse(logFC > 0, "up", "down"),
      geneset = "growth/cellcycle"
    )
)
#----------------
unique(top$geneset)
top <- top %>% 
  bind_rows(ros_top5)   # 
top$geneset <- factor(
  top$geneset,
  levels = c("Cell adhesion molecules (CAMs)",
             "Oxphos/ Electron transport",
             "Ribosome machinery",
             "Replication/ cell cycle",
             "ROS pathway",
             "growth/cellcycle"
             )
)
table <- top[,c("genes","geneset")]
write_rds(table , basedir("supp.Fig.2A_genes.rds"))
### ---------------------------------------------------------
### Final plot
### ---------------------------------------------------------

Supp_fig_2a <- ggplot(
  top,
  aes(
    y = celltype,
    x = genes,
    color = pmin(2, pmax(-2, logFC)),
    size = pmin(3, -log10(adj.P.Val))
  )
) +
  geom_point() +
  scale_color_gradient2(low = "#4C889C", mid = "white", high = "#D0154E") +
  scale_size_continuous(
    range = c(0, 1.5),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    title = "Differentially Expressed Genesets",
    y = "Cell Type",
    x = "Genes",
    color = "logFC",
    size = "-log10(padj)"
  ) +
  facet_grid(cols = vars(geneset), scales = "free", space = "free") +
  optimized_theme_fig() +
  theme(strip.text.x = element_text(angle = 90, hjust = 0),
        axis.text.x = element_text(angle = 90)
        )

Supp_fig_2a



#save
ggsave(basedir("Supp.Fig.2A.pdf"), plot = Supp_fig_2a, w = 18.2, h = 7, units = "cm")

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

# p_box <- ggplot(zscore_external, aes(x = sample_label, y = zscore)) +
#       geom_boxplot(aes(color = tissue), outlier.shape = NA, alpha = 0.8, width = 0.5) +
#       labs(x = NULL, y = "Z-score") +
#       optimized_theme_fig() +
#       scale_fill_manual(values = dataset_color_map[[ds]] ) +
#       theme(axis.text.x = element_blank(),
#             axis.ticks = element_blank(),
#             legend.position = "none")  # legend is unnecessary now
#     
#     if(!is.na(boundary_index)){
#       p_box <- p_box + geom_vline(xintercept = boundary_index + 0.2, linetype = "dashed", color = "grey50")
#     }
# 
# p_box <- p_box +
#       scale_color_manual(values = dataset_color_map[[ds]], name = "culture model") +
#       scale_x_discrete(drop = FALSE) +
#       labs(x = NULL, y = "Z-score") +
#       optimized_theme_fig() +
#       theme(axis.text.x = element_blank(),
#             axis.ticks = element_blank(),
#             legend.position = "bottom")
#     
    # -----------------------------
    # Density panel
    # -----------------------------
    # DENSITY
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




# Step 1: Calculate median zscore per sample per tissue_celltype across genes
median_per_sample <- zscore_external %>%
  group_by(tissue_celltype, sample, tissue) %>%
  summarize(median_zscore = median(zscore, na.rm = TRUE)) %>%
  ungroup()

# Define your manual palettes
my_color <- c(
  "#4A2C69",  # ex.vivo → violet
  "#F9C319",  # in.vivo → yellow
  "#D4E062",  # izzo → green
  "#E5971E"   # Anna_et_al → orange
)

ex_vivo_palette <- c("#4A2C69", "#613488", "#882A73", "#AE77AF")
# Step 2: Plot boxplot per tissue_celltype and jitter points per sample

median_per_sample <- median_per_sample %>%
  mutate(
    tissue_celltype = factor(tissue_celltype, levels = unique(tissue_celltype)),
   # tissue = gsub("Anna_et_al","Izzo et. al"),
    tissue = factor(tissue, levels = c("ex.vivo", "in.vivo", "izzo", "Anna_et_al"))  # optional, controls facet order
  )

plot3 <- ggplot(median_per_sample, aes(x = tissue_celltype, y = median_zscore)) +
  geom_boxplot(aes(color = tissue), alpha = 0.7, outlier.shape = NA) +
  #geom_jitter(color = "darkgrey", width = 0.2, size = 0.8, alpha = 0.7) +
  labs(
    title = "Median Z-score per Sample across ISG-core genes",
    x = "Celltype",
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
  facet_grid(~tissue, scales = "free", space = "free")+
  scale_x_discrete(labels = function(x) gsub("^(ex\\.vivo_|in\\.vivo_|izzo_|Anna_et_al_)", "", x)) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
plot3
ggsave(basedir(paste0("Sup_Fig2B", "_per_Tissue_line.pdf")), plot = plot3, w = 18, h = 5, units = "cm")


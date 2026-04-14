source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
library(tidyverse)
library(enrichR)
library(purrr)
library("scales")
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library(readr)

basedir <- dirout("Figure7_Supplementary_v8")

InDir7  <-  dirout("Ag_ScRNA_12_fgsea_overlap")

gsea.res_pert <- read_rds(InDir7("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction_invivo.rds"))
gsea.res_pert <- gsea.res_pert %>%
  filter(coef %in% grep("interaction",gsea.res_pert$coef, value = T))
gsea.res_pert$coef <- gsub("interaction","",gsea.res_pert$coef )

db = "MSigDB_Hallmark_2020"


###########
InDir8 <- dirout("Ag_ScRNA_17_JAKSTAT_Ar")

gsea.res_spleen <- read_rds(InDir8("fgsea_hom_vs_ex.vivo_per_CT.rds"))

###########
InDir9 <- dirout("Fig.Glio")
gsea.res_glio <- read_rds(InDir9("enrichment_interactionKO_at_noRT.rds"))



#select columns

gsea.res_glio <- gsea.res_glio %>%
  dplyr::select("db","NES","padj","celltype","coef","pathway","leadingEdge")

gsea.res_pert <- gsea.res_pert %>%
  dplyr::select("db","NES","padj","celltype","coef","pathway","leadingEdge")

gsea.res_spleen <- gsea.res_spleen %>%
  dplyr::select("db","NES","padj","celltype","coef","pathway","leadingEdge")%>%
  #filter(coef %in%)
  mutate(genotype = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(celltype = gsub("M", "Macrophage", celltype)) %>%
  mutate(celltype = gsub("T8", "T-cells", celltype))%>%
  filter(coef != "(Intercept)")%>%
  filter(celltype == "T-cells")
unique(gsea.res_spleen$coef)  
#gsea.res_spleen$genotype <- gsub("Interaction_","",gsea.res_spleen$genotype)
gsea.res_spleen$genotype <- gsub("treatmentex_vivo","WT",gsea.res_spleen$genotype )
# Step 2: Summarize to find KOs with at least one valid cell type
gsea.res_spleen <- gsea.res_spleen %>%
  filter(genotype %in% grep("Interaction",gsea.res_spleen$genotype, value = T))
gsea.res_spleen$coef <- tolower(gsub("genotype|:treatmentex_vivo","",gsea.res_spleen$coef))
gsea.res_spleen <- gsea.res_spleen %>%
  dplyr::select("db","NES","padj","celltype","coef","pathway","leadingEdge")
gsea.res_glio$dataset <- "Glioblastoma"
gsea.res_pert$dataset <- "Perturb-seq(Hematiopoietic)"
gsea.res_spleen$dataset <- "Spleen"

combined <- rbind(gsea.res_pert,gsea.res_spleen,gsea.res_glio)%>%
  filter(db == "MSigDB_Hallmark_2020")


write_rds(combined,basedir("GSEA_all_datasets_KO_interaction.rds"))

pathways_int <- c("Interferon Gamma Response", "Interferon Alpha Response",
                  "Cholesterol Homeostasis",
                  "mTORC1 Signaling","Myc Targets V2", "Myc Targets V1")

combined_selected <- combined %>%
  filter(pathway %in% pathways_int)


df <- combined_selected %>%
  mutate(logpadj = -log10(padj))

ggplot(df, aes(y = str_wrap(pathway,50), x = coef)) +
  geom_point(aes(size = pmin(logpadj, 3), color = pmax(-2,pmin(NES,2)))) +
  facet_grid(cols = vars(celltype), rows = vars(dataset), scales = "free", space = "free") +
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        midpoint = 0, name = "NES") +
  scale_size_continuous(
    range = c(0,1),
    name = "-log10(padj)")+
  labs(y =NULL,
       x = "Pathway",
       title = "Enrichment Dot Plot"
  )+
  
  optimized_theme_fig()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(angle = 0, hjust = 0),
    legend.position = "bottom")
ggsave(basedir("Enrichment across_dataset and celltype.pdf"), w =42, h = 7.5, units = "cm")
#ISG_core_Sup.Fig.7B---------------------
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))

#gene_based analysis

# Step 1: Filter for significant genes (adj.P.Val < 0.05 and abs(logFC) > 1)
limmaRes_significant <- limmaRes %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)  # Only significantly altered genes

data <- summary_df %>%
  filter(coef %in% koi) %>%
  filter(Count >= 10) %>%
  inner_join(ko_flags, by = c("celltype", "coef")) %>%
  filter(valid_ko)
data %>% write_rds(basedir("kos_fig3.rds"))

# grouping gene sets specifically for ISG and cholesterol
genes_fig_1 <- read_rds(InDir_genes("genes_fig1.rds"))
colnames(genes_fig_1) <- c("ensg","pathway")
genes_fig_1$pathway <- gsub("mTORC1_or_Cholesterol","mTORC1/Cholesterol", genes_fig_1$pathway)

filtered_genes <- limmaRes %>%
  inner_join(summary_df,by =c("coef","celltype"))%>%
  filter(Count >10)%>%
  #filter(group != "n.s")%>%
  filter(ensg %in% genes_fig_1$ensg, coef %in% koi) %>%
  group_by(celltype, coef) %>%
  merge(genes_fig_1, by = "ensg") %>%
  left_join(ko_flags, by = c("coef" = "genotype", "celltype")) %>%  # Merge with KO flags per cell type
  filter(valid_ko == TRUE)  # Keep only valid KOs for the specific cell type
#
isg_count_summary <- filtered_genes %>%
  filter(pathway == "ISG core") %>%
  group_by(coef) %>%
  summarize(
    n_sig_genes = sum(adj.P.Val < 0.05 & abs(logFC) > 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_genes))

ggplot(isg_count_summary, aes(x = reorder(coef, n_sig_genes), y = n_sig_genes)) +
  geom_segment(aes(xend = coef, y = 0, yend = n_sig_genes), color = "#D0154E") +
  geom_point(color = "#D0154E", size = 4) +
  coord_flip() +
  labs(
    title = "Number of ISGs with significant interaction effect per KO",
    x = "KO", y = "Number of significant ISGs"
  ) +
  optimized_theme_fig()


# fig---------
Sup.Fig.7B <- ggplot(filtered_genes %>%
                  filter(pathway == "ISG core"), aes(x = coef, y = ensg,
                                                     color = pmin(2, pmax(-2, logFC)) ,
                                                     size = pmin(3, -log10(adj.P.Val))
                  )) +  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0,1.5),
    #limits = c(0,5),
    #breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(title = "Interaction effect of ISG core genes",
       x = "KOs",
       y = "Genes")+
  facet_grid(cols = vars(celltype), rows = vars(pathway), scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()+theme(
    legend.position = "bottom",
    strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
    axis.text.x = element_text(angle = 45))


Sup.Fig.7B
#paper--------------
ggsave(
  filename = basedir("Sup.Fig.7B.pdf"),
  plot = Sup.Fig.7B,
  width = 16,
  height = 8.2,
  units = "cm"
)
#UMAP---Sup.Fig.7C------------------

inDir1 <- dirout_load("Ag_SCRNA_05_01_UMAPs_and_celltypes")

xu <- xlab("UMAP 1")
yu <- ylab("UMAP 2")

# SETTINGS
annList <- readRDS(inDir1("ProjVivo_celltypes.RDS"))

# Other projections
umap.proj <- list(
  original = readRDS(inDir1("ProjMonocle.RDS")),
  #izzo = readRDS(dirout_load("Ag_SCRNA_UMAPs")("ProjIzzo.RDS")),
  in.vivo = readRDS(inDir1("ProjVivo.RDS")),
  in.vivo.X = readRDS(inDir1("ProjVivoX.RDS"))
)

# SIMPLE-

mobjs <- list()

tissue<-c("ex.vivo","in.vivo")
for(tissuex in tissue){
  (base::load(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis//Ag_SCRNA_02_01_Integration/",tissuex,"/soupx/MonocleObject.RData")))
  mobjs[[tissuex]] <- monocle.obj
}

# Remove duplets
for(tissuex in names(mobjs)){
  mobjs[[tissuex]] <- mobjs[[tissuex]][,colnames(mobjs[[tissuex]]) %in% annList$rn]
}

# Initialize the guide and tissue columns in annList
annList$guide <- NA
annList$tissue <- NA

# Match and extract guide and tissue information from `ex.vivo`
match_ex_vivo <- match(annList$rn, colnames(mobjs$ex.vivo))
idx <- which(is.na(annList$guide) & !is.na(match_ex_vivo))
annList$guide[idx] <- mobjs$ex.vivo@colData$guide[match_ex_vivo[idx]]
#annList$guide[is.na(annList$guide) &!is.na(match_ex_vivo)] <- mobjs$ex.vivo@colData$guide[match_ex_vivo[(match_ex_vivo)]]
annList$tissue[!is.na(match_ex_vivo)] <- "ex.vivo"  # Assign "ex.vivo" for matched rows

# Match and extract guide and tissue information from `in.vivo`
match_in_vivo <- match(annList$rn, colnames(mobjs$in.vivo))
annList$guide[is.na(annList$guide) & !is.na(match_in_vivo)] <- mobjs$in.vivo@colData$guide[match_in_vivo[!is.na(match_in_vivo)]]
annList$tissue[is.na(annList$tissue) & !is.na(match_in_vivo)] <- "in.vivo"  # Assign "in.vivo" for matched rows

# Filter annList to keep only rows where guide is "NTC"
annList <- annList[annList$guide %in% c("NTC","Brd9_AS_70306", "Brd9_BR_48004", "Brd9_BR_48005"), ]
unique(annList$guide)

annList <- annList %>%
  mutate(
    functional.cluster = case_when(
      functional.cluster %in% c("MEP (G1)", "MEP (pert.)", "MEP (S)") ~ "MEP",  # Grouping all MEP categories to "MEP"
      functional.cluster == "Imm. B-cell" ~ "Imm.B.cell",  # Correct "Imm. B-cell" to "Imm.B.cell"
      is.na(functional.cluster) ~ "Unknown",  # Handle NA values if necessary
      TRUE ~ functional.cluster  # Leave other values unchanged
    )
  )%>%
  mutate(genotype = gsub("_.*","",guide))
unique(annList$genotype)
#annList$tissue <- gsub("ex.vivo", "ex.vivo", annList$tissue)
in.vivo <- umap.proj$in.vivo
# Filter in.vivo based on annList rn to keep only NTC samples

in.vivo <- inner_join(in.vivo, annList, by = c("rn", "tissue"))

# Generate hexbin object based on filtered in.vivo (only NTC samples)
hex.obj <- hexbin::hexbin(x = in.vivo$UMAP_1, y = in.vivo$UMAP_2, xbins = 100, IDs = TRUE)
in.vivo <- cbind(in.vivo, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
pDT <- in.vivo
pDT <- pDT[, .N, by = c("hex.x", "hex.y", "functional.cluster")]
pDT[, sum := sum(N), by = c("hex.x", "hex.y")]
pDT[, frac := N / sum]

# Filter clusters with significant fractions
pDT <- pDT[frac > 0.25]

# Merge summary data back with the original dataset
merged_data <- inner_join(in.vivo, pDT, by = c("hex.x", "hex.y", "functional.cluster"), all.x = TRUE)


# Check the unique values in functional.cluster to ensure it's working

# Generate cluster labels for significant clusters
pDT.labels <- pDT[frac > 0.25, .(hex.x = median(hex.x), hex.y = median(hex.y)), by = c("functional.cluster")]
pDT.labels %>% write_rds(basedir("pDT.labels.rds"))
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
  "MEP" = "pink",        # Lighter gray for MEP
  "Ery" = "red",        # Slightly darker gray for Ery
  "Imm.B.cell" = "gray"     # Other
)

unique(merged_data$guide)
merged_data %>% write_rds(basedir("Cross_projected_on_in.vivo.rds"))
# Ensure factor ordering for correct label display
merged_data$functional.cluster <- factor(merged_data$functional.cluster, 
                                         levels = names(cluster_colors))

ggplot(merged_data[tissue != "leukemia"], aes(x = UMAP_1, y = UMAP_2)) + 
  
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
  facet_grid(cols = vars(tissue), rows = vars(genotype)) + 
  # Defining color manual scale for clusters
  # Defining color manual scale for clusters
  scale_color_manual(name = "functional.cluster", values = cluster_colors) + 
  
  
  optimized_theme_fig() +
  
  # Positioning legends separately
  # Positioning legends separately
  theme(
    axis.text.x = element_text(angle = 0),
    #,
    #legend.position = "none")
    legend.position = "bottom",  # Color legend at the bottom
    legend.box = "horizontal",   # Horizontal alignment of legends
    legend.text = element_text(size = 5),     # Adjust legend text font size
    legend.title = element_blank(), # Remove title for color legend
    legend.spacing = unit(0.5, "cm"),  # Adjust spacing between legends
    legend.key.size = unit(0.5, "lines")  # Adjust size of the legend keys
  )
# Correct axis labels (you can define xu and yu separately if needed)
#ggsave(basedir("UMAP_InvivoX_NTC_Brd9.pdf"), w=3,h=2)
ggsave(basedir("UMAP_InvivoX_NTC_Brd9.pdf"), w = 5,h = 4)
ggsave(basedir("UMAP_InvivoX_NTC_Brd9.png"), w = 3,h = 2)

ggplot(merged_data[tissue != "leukemia"], aes(x=UMAP_1, UMAP_2)) + 
  themeNF() +
  stat_binhex(aes(fill=log10(..count..)), bins=100) + 
  scale_fill_gradientn(colors = c("lightgrey", "#1f78b4", "#e31a1c", "#ff7f00")) +
  facet_grid(cols = vars(tissue), rows = vars(genotype)) + 
  xu + yu+
  
  optimized_theme_fig()+
  theme(legend.position = "none")
ggsave(basedir("UMAP_InvivoX_NTC_Brd9_distribution.png"), w=3,h=2)
#













# ###############
# InDir4 <- dirout("Figure1")
# genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
# colnames(genes_fig_1) <- c("gene","pathways")
# goi <- genes_fig_1[genes_fig_1$pathways == "ISG core",]$gene
# 
# #Int pert
# Int_pert <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
#   mutate(KO = gsub("interaction","",coef))%>%
#   mutate(genes = ensg)%>%
#   filter(genes %in% goi)%>%
#   dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes"))
# Int_pert$dataset <- "Perturb-seq(Hematopoiesis)"
# 
# 
# #Int_spleen
# InDir1 <- dirout("JAK_STAT_splenic_M_T")
# Int_spleen <- read_rds(InDir1("combined_jakstat_diff_exp.rds"))
# Int_spleen <- Int_spleen %>% filter(grepl("treatmentex_vivo",Int_spleen$coef))
# 
# Int_spleen$group <- ifelse(Int_spleen$logFC >= 1 & 
#                              Int_spleen$adj.P.Val <= 0.05, "up", 
#                            ifelse(Int_spleen$logFC <= -1 & 
#                                     Int_spleen$adj.P.Val <= 0.05, "down", "n.s"))
# unique(genes_fig_1[genes_fig_1$pathways == "ISG core",]$pathways)
# 
# # Modify the 'coef' column for any KO
# Int_spleen <- Int_spleen %>%
#   mutate(celltype = gsub("M", "Macrophage", cell_type)) %>%
#   mutate(celltype = gsub("T8", "T-cells", cell_type)) %>%
#   filter(celltype == "T-cells")%>%
#   filter(coef != "treatmentex_vivo")%>%
#   mutate(KO = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
#   mutate(KO = gsub("Interaction_","",KO))%>%
#   mutate(genes = gene)%>%
#   #filter(genes %in% goi)%>%
#   dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
# Int_spleen$dataset <-"Spleen"
# 
# summary_spleen_int <- Int_spleen %>%
#   group_by(KO) %>%
#   summarise(signif_count = sum(group %in% c("up", "down")),
#             .groups = "drop")
# selected_spleen_KOs <- summary_spleen_int %>%
#   filter(signif_count > 10) %>%
#   pull(KO) %>%
#   unique()
# Fig_interaction <-  ggplot(summary_spleen_int %>% filter(KO %in% selected_spleen_KOs),aes(
#   x = reorder(KO,signif_count),
#   y = log10(signif_count)
# )) +
#   geom_col(color = "darkgrey", fill = NA, width = 0.6) +
#   
#   labs(
#     title = "No. of genes with interaction effects per KO",
#     x = NULL,
#     y = expression(atop("Number of genes with", 
#                         paste("interaction effects ", log[10](n))))
#   ) +
#   # Custom theme with no legend if not needed
#   optimized_theme_fig() + 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#     panel.spacing = unit(0.02, "cm"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#     
#   )
# # Display the plot
# Fig_interaction
# 
# # ggsave(basedir("interaction_KO_number_spleen.pdf"),
# #        w= 2.8, h = 3, units = "cm")
# #########
# InDir <- dirout("Glioblastoma_limmaRes.3way.mod")
# 
# limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))
# 
# limmaRes_int <-  limmaRes %>% filter(coef %in% grep("tissueex.vivo.genotype[A-Za-z0-9]*$",
#                                                     limmaRes$coef, value = T))
# limmaRes_int$celltype <- "glio"
# 
# Int_glio <- limmaRes_int %>%
#   mutate(KO = gsub("tissueex.vivo.genotype","", coef))%>%
#   mutate(genes = ensg)%>%
#  
#   dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
# Int_glio$dataset <- "Glioblastoma"
# 
# summary_interaction_KO <- Int_glio%>%
#   # group by interaction term
#   group_by(KO) %>%
#   summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")
# selected_KOs <-summary_interaction_KO %>%
#   filter(signif_count> 10)%>%
#   pull(KO)%>%
#   unique()
# Fig_interaction <-  ggplot(summary_interaction_KO %>% filter(KO %in% selected_KOs),aes(
#   x = reorder(KO,signif_count),
#   y = log10(signif_count)
# )) +
#   geom_col(color = "darkgrey", fill = NA, width = 0.6) +
#   
#   labs(
#     title = "No. of genes with interaction effects per KO",
#     x = NULL,
#     y = expression(atop("Number of genes with", 
#                         paste("interaction effects ", log[10](n))))
#   ) +
#   # Custom theme with no legend if not needed
#   optimized_theme_fig() + 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#     panel.spacing = unit(0.02, "cm"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#     
#   )
# # Display the plot
# Fig_interaction
# # ggsave(basedir("interaction_KO_number.pdf"),
# #        w= 8, h = 3, units = "cm")
# 
# 
# 
# 
# limma_subset <- rbind(Int_glio,Int_pert,Int_spleen) 
# 
# # limma_subset$cell_type <- recode(limma_subset$cell_type,"Macrophage" = "M",
# #                                  "T-cells"  = "T")
# Fig <- limma_subset %>%
#   ggplot(aes(
#     x = KO, 
#     y = genes, 
#     color = pmin(2, pmax(-2, logFC)),
#     size = pmin(5, -log10(adj.P.Val))
#   )) +
#   geom_point() +
#   
#   scale_color_gradient2(
#     low = "#4C889C",
#     mid = "white",
#     high = "#D0154E",
#     name = TeX("log_{2}(FC)")
#   ) +
#   scale_size_continuous(
#     range = c(0, 1.8),
#     name = TeX("$-\\log_{10}(p_{adj})$")
#   ) +
#   labs(
#     title = "Differentially expressed ISGs and\ngrowth/metabolic genes",
#     x = "Cell type",
#     y = "Genes"
#     #caption = "M: Macrophages    T: T-cells"
#   ) +
#   facet_wrap(
#     ~dataset*celltype,      # << put facets side-by-side
#     scales = "free_y", 
#     nrow = 3
#   ) +
#   optimized_theme_fig()+
#   theme(
#     legend.position = "bottom",
#     axis.text.x = element_text(angle = 90),
#     plot.caption = element_text(
#       hjust = 0,          # Left-align
#       size = 5, color = "black"    # Or "Times", "Courier", etc. (must be installed)
#     ))
# Fig
# ggsave(basedir("Combined_dataset_ISG_int.pdf"), w = 27,
#        h = 8, units = "cm")
# 
# library(dplyr)
# library(ggplot2)
# library(stringr)
# 
# # List of dataset names
# datasets <- unique(limma_subset$dataset)
# 
# for (d in datasets) {
#   
#   df <- limma_subset %>% filter(dataset == d)
#   
#   Fig_d <- ggplot(df, aes(
#     x = KO, 
#     y = genes, 
#     color = pmin(2, pmax(-2, logFC)),
#     size = pmin(5, -log10(adj.P.Val))
#   )) +
#     geom_point() +
#     scale_color_gradient2(
#       low = "#4C889C",
#       mid = "white",
#       high = "#D0154E",
#       name = TeX("log_{2}(FC)")
#     ) +
#     scale_size_continuous(
#       range = c(0, 1.8),
#       name = TeX("$-\\log_{10}(p_{adj})$")
#     ) +
#     labs(
#       title = paste0("Differentially expressed ISGs and growth/metabolic genes Dataset: ", d),
#       x = "Cell type",
#       y = "Genes"
#     ) +
#     facet_grid(
#       . ~ celltype,                # only facet by celltype now
#       scales = "free",
#       space = "free"
#     ) +
#     optimized_theme_fig() +
#     theme(
#       legend.position = "bottom",
#       axis.text.x = element_text(angle = 90),
#       plot.caption = element_text(hjust = 0, size = 5, color = "black")
#     )
#   
#   # Print to screen
#   print(Fig_d)
#   
#   # Save file with dataset name
#   outfile <- basedir(paste0("ISG_plot_", d, ".pdf"))
#   ggsave(outfile, Fig_d, width = 27, height = 8, units = "cm")
# }
# 
# 
# summary_invivo_KO <- Int_glio%>%
#   # group by interaction term
#   mutate(KO = gsub("_in.vivo_noRT","", coef))%>%
#   group_by(KO) %>%
#   summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")
# 
# Fig_interaction <-  ggplot(summary_interaction_KO,aes(
#   x = reorder(KO,signif_count),
#   y = log10(signif_count)
# )) +
#   geom_col(color = "darkgrey", fill = NA, width = 0.6) +
#   
#   labs(
#     title = "No. of genes with interaction effects per KO",
#     x = NULL,
#     y = expression(atop("Number of genes with", 
#                         paste("interaction effects ", log[10](n))))
#   ) +
#   # Custom theme with no legend if not needed
#   optimized_theme_fig() + 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#     panel.spacing = unit(0.02, "cm"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#     
#   )
# # Display the plot
# Fig_interaction
# ggsave(out("interaction_KO_number.pdf"))
#        
# ###############
# 
# 
# #############
# #
# InDir4 <- dirout("Figure1")
# genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
# colnames(genes_fig_1) <- c("gene","pathways")
# 
# #Int pert
# Int_pert <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
#   mutate(KO = gsub("interaction","",coef))%>%
#   mutate(genes = ensg)%>%
#   #filter(genes %in% genes_fig_1$gene)%>%
#   dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes"))
# Int_pert$dataset <- "Perturb-seq(Hematopoiesis)"
# 
# 
# 
# 
# 
# 
# 
# #Int_spleen
# InDir1 <- dirout("Figure4")
# Int_spleen <- read_rds(InDir1("combined_jakstat_diff_exp.rds"))
# Int_spleen <- Int_spleen %>% filter(grepl("treatmentex_vivo",Int_spleen$coef))
# 
# Int_spleen$group <- ifelse(Int_spleen$logFC >= 1 & 
#                              Int_spleen$adj.P.Val <= 0.05, "up", 
#                            ifelse(Int_spleen$logFC <= -1 & 
#                                     Int_spleen$adj.P.Val <= 0.05, "down", "n.s"))
# 
# # Modify the 'coef' column for any KO
# Int_spleen <- Int_spleen %>%
#   mutate(celltype = gsub("M", "Macrophage", cell_type)) %>%
#   mutate(celltype = gsub("T8", "T-cells", cell_type)) %>%
#   filter(celltype == "T-cells")%>%
#   filter(coef != "treatmentex_vivo")%>%
#   mutate(KO = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
#   mutate(KO = gsub("Interaction_","",KO))%>%
#   mutate(genes = gene)%>%
#   #filter(genes %in% genes_fig_1$gene)%>%
#   dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
# Int_spleen$dataset <-"Spleen"
# 
# summary_spleen_int <- Int_spleen %>%
#   group_by(KO) %>%
#   summarise(signif_count = sum(group %in% c("up", "down")),
#             .groups = "drop")
# selected_spleen_KOs <- summary_spleen_int%>%
#   filter(signif_count > 10) %>%
#   pull(KO) %>%
#   unique()
# Fig_interaction <-  ggplot(summary_spleen_int %>% filter(KO %in% selected_spleen_KOs),aes(
#   x = reorder(KO,signif_count),
#   y = log10(signif_count)
# )) +
#   geom_col(color = "darkgrey", fill = NA, width = 0.6) +
#   
#   labs(
#     title = "No. of genes with interaction effects per KO",
#     x = NULL,
#     y = expression(atop("Number of genes with", 
#                         paste("interaction effects ", log[10](n))))
#   ) +
#   # Custom theme with no legend if not needed
#   optimized_theme_fig() + 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#     panel.spacing = unit(0.02, "cm"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#     
#   )
# # Display the plot
# Fig_interaction
# 
# ggsave(basedir("interaction_KO_number_spleen.pdf"),
#        w= 2.8, h = 3, units = "cm")
# #########
# InDir <- dirout("Glioblastoma_limmaRes.3way.mod")
# 
# limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))
# 
# limmaRes_int <-  limmaRes %>% filter(coef %in% grep("tissueex.vivo.genotype[A-Za-z0-9]*$",
#                                                     limmaRes$coef, value = T))
# limmaRes_int$celltype <- "glio"
# 
# Int_glio <- limmaRes_int %>%
#   mutate(KO = gsub("tissueex.vivo.genotype","", coef))%>%
#   mutate(genes = ensg)%>%
#   #filter(genes %in% genes_fig_1$gene)%>%
#   dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
# Int_glio$dataset <- "Glioblastoma"
# 
# summary_invivo_KO <- Int_glio%>%
#   # group by interaction term
#   group_by(KO) %>%
#   summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")
# selected_KOs <-summary_interaction_KO%>%
#   filter(signif_count> 10)%>%
#   pull(KO)%>%
#   unique()
# Fig_interaction <-  ggplot(summary_interaction_KO %>% filter(KO %in% selected_KOs),aes(
#   x = reorder(KO,signif_count),
#   y = log10(signif_count)
# )) +
#   geom_col(color = "darkgrey", fill = NA, width = 0.6) +
#   
#   labs(
#     title = "No. of genes with interaction effects per KO",
#     x = NULL,
#     y = expression(atop("Number of genes with", 
#                         paste("interaction effects ", log[10](n))))
#   ) +
#   # Custom theme with no legend if not needed
#   optimized_theme_fig() + 
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
#     panel.spacing = unit(0.02, "cm"),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#     
#   )
# # Display the plot
# Fig_interaction
# ggsave(basedir("interaction_KO_number.pdf"),
#        w= 8, h = 3, units = "cm")
# 
# 
# colnames(Int_glio)
# 
# colnames(Int_spleen)
# 
# 
# colnames(Int_pert)
# Int_glio$group <- NULL
# Int_spleen$group <- NULL
# limma_subset <- rbind(Int_glio,Int_pert,Int_spleen) 
# 
# # limma_subset$cell_type <- recode(limma_subset$cell_type,"Macrophage" = "M",
# #                                  "T-cells"  = "T")
# Fig <- limma_subset %>%
#   ggplot(aes(
#     x = KO, 
#     y = genes, 
#     color = pmin(2, pmax(-2, logFC)),
#     size = pmin(5, -log10(adj.P.Val))
#   )) +
#   geom_point() +
#   
#   scale_color_gradient2(
#     low = "#4C889C",
#     mid = "white",
#     high = "#D0154E",
#     name = TeX("log_{2}(FC)")
#   ) +
#   scale_size_continuous(
#     range = c(0, 1.8),
#     name = TeX("$-\\log_{10}(p_{adj})$")
#   ) +
#   labs(
#     title = "Differentially expressed ISGs and\ngrowth/metabolic genes",
#     x = "Cell type",
#     y = "Genes"
#     #caption = "M: Macrophages    T: T-cells"
#   ) +
#   facet_grid(
#     dataset~celltype,      # << put facets side-by-side
#     scales = "free", 
#     space = "free"
#   ) +
#   optimized_theme_fig()+
#   theme(
#     legend.position = "bottom",
#     axis.text.x = element_text(angle = 90),
#     plot.caption = element_text(
#       hjust = 0,          # Left-align
#       size = 5, color = "black"    # Or "Times", "Courier", etc. (must be installed)
#     ))
# Fig
# de
# ggsave(basedir("Combined_dataset_ISG_int.pdf"), w = 27,
#        h = 8, units = "cm")
# 
# library(dplyr)
# library(ggplot2)
# library(stringr)
# 
# # List of dataset names
# datasets <- unique(limma_subset$dataset)
# 
# for (d in datasets) {
#   
#   df <- limma_subset %>% filter(dataset == d)
#   
#   Fig_d <- ggplot(df, aes(
#     x = KO, 
#     y = genes, 
#     color = pmin(2, pmax(-2, logFC)),
#     size = pmin(5, -log10(adj.P.Val))
#   )) +
#     geom_point() +
#     scale_color_gradient2(
#       low = "#4C889C",
#       mid = "white",
#       high = "#D0154E",
#       name = TeX("log_{2}(FC)")
#     ) +
#     scale_size_continuous(
#       range = c(0, 1.8),
#       name = TeX("$-\\log_{10}(p_{adj})$")
#     ) +
#     labs(
#       title = paste0("Differentially expressed ISGs and growth/metabolic genes Dataset: ", d),
#       x = "Cell type",
#       y = "Genes"
#     ) +
#     facet_grid(
#       . ~ celltype,                # only facet by celltype now
#       scales = "free",
#       space = "free"
#     ) +
#     optimized_theme_fig() +
#     theme(
#       legend.position = "bottom",
#       axis.text.x = element_text(angle = 90),
#       plot.caption = element_text(hjust = 0, size = 5, color = "black")
#     )
#   
#   # Print to screen
#   print(Fig_d)
#   
#   # Save file with dataset name
#   outfile <- basedir(paste0("ISG_plot_", d, ".pdf"))
#   ggsave(outfile, Fig_d, width = 27, height = 8, units = "cm")
# }
# 
# ################
# Fig <- limma_subset %>%
#   filter(genes %in% goi)%>%
#   ggplot(aes(
#     x = KO,
#     y = genes,
#     color = pmin(2, pmax(-2, logFC)),
#     size = pmin(5, -log10(adj.P.Val))
#   )) +
#   geom_point() +
#   
#   scale_color_gradient2(
#     low = "#4C889C", mid = "white", high = "#D0154E",
#     name = TeX("log_{2}(FC)")
#   ) +
#   scale_size_continuous(
#     range = c(0, 1.5),
#     name = TeX("$-\\log_{10}(p_{adj})$")
#   ) +
#   
#   facet_grid(
#     ~dataset*celltype,
#     scales = "free_x",
#     space  = "free_x"     # <<< equal spacing, but dropping unused KOs
#   ) +
#   
#   labs(
#     title = "Differentially expressed ISGs and\ngrowth/metabolic genes",
#     x = "Cell type",
#     y = "Genes"
#   ) +
#   
#   optimized_theme_fig() +
#   theme(
#     legend.position = "bottom",
#     axis.text.x = element_text(angle = 90),
#     plot.caption = element_text(hjust = 0, size = 5, color = "black")
#   )
# Fig
# ggsave(basedir("ISGs_across_datasets_2.pdf"), w = 35.6, h = 8.2, units = "cm")  
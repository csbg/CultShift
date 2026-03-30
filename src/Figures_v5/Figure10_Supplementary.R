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

basedir <- dirout("Figure10_Supplementary")
###############
InDir4 <- dirout("Figure1")
genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
colnames(genes_fig_1) <- c("gene","pathways")
goi <- genes_fig_1[genes_fig_1$pathways == "ISG core",]$gene

#Int pert
Int_pert <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(KO = gsub("interaction","",coef))%>%
  mutate(genes = ensg)%>%
  filter(genes %in% goi)%>%
  dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes"))
Int_pert$dataset <- "Perturb-seq(Hematopoiesis)"


#Int_spleen
InDir1 <- dirout("Figure4")
Int_spleen <- read_rds(InDir1("combined_jakstat_diff_exp.rds"))
Int_spleen <- Int_spleen %>% filter(grepl("treatmentex_vivo",Int_spleen$coef))

Int_spleen$group <- ifelse(Int_spleen$logFC >= 1 & 
                             Int_spleen$adj.P.Val <= 0.05, "up", 
                           ifelse(Int_spleen$logFC <= -1 & 
                                    Int_spleen$adj.P.Val <= 0.05, "down", "n.s"))
unique(genes_fig_1[genes_fig_1$pathways == "ISG core",]$pathways)

# Modify the 'coef' column for any KO
Int_spleen <- Int_spleen %>%
  mutate(celltype = gsub("M", "Macrophage", cell_type)) %>%
  mutate(celltype = gsub("T8", "T-cells", cell_type)) %>%
  filter(celltype == "T-cells")%>%
  filter(coef != "treatmentex_vivo")%>%
  mutate(KO = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(KO = gsub("Interaction_","",KO))%>%
  mutate(genes = gene)%>%
  filter(genes %in% goi)%>%
  dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
Int_spleen$dataset <-"Spleen"

summary_spleen_int <- Int_spleen %>%
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")),
            .groups = "drop")
selected_spleen_KOs <- summary_spleen_int%>%
  filter(signif_count > 10) %>%
  pull(KO) %>%
  unique()
Fig_interaction <-  ggplot(summary_spleen_int %>% filter(KO %in% selected_spleen_KOs),aes(
  x = reorder(KO,signif_count),
  y = log10(signif_count)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
# Display the plot
Fig_interaction

# ggsave(basedir("interaction_KO_number_spleen.pdf"),
#        w= 2.8, h = 3, units = "cm")
#########
InDir <- dirout("Glioblastoma_limmaRes.3way.mod")

limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))

limmaRes_int <-  limmaRes %>% filter(coef %in% grep("tissueex.vivo.genotype[A-Za-z0-9]*$",
                                                    limmaRes$coef, value = T))
limmaRes_int$celltype <- "glio"

Int_glio <- limmaRes_int %>%
  mutate(KO = gsub("tissueex.vivo.genotype","", coef))%>%
  mutate(genes = ensg)%>%
  filter(genes %in% goi)%>%
  dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
Int_glio$dataset <- "Glioblastoma"

summary_invivo_KO <- Int_glio%>%
  # group by interaction term
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")
selected_KOs <-summary_interaction_KO%>%
  filter(signif_count> 10)%>%
  pull(KO)%>%
  unique()
Fig_interaction <-  ggplot(summary_interaction_KO %>% filter(KO %in% selected_KOs),aes(
  x = reorder(KO,signif_count),
  y = log10(signif_count)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
# Display the plot
Fig_interaction
# ggsave(basedir("interaction_KO_number.pdf"),
#        w= 8, h = 3, units = "cm")




limma_subset <- rbind(Int_glio,Int_pert,Int_spleen) 

# limma_subset$cell_type <- recode(limma_subset$cell_type,"Macrophage" = "M",
#                                  "T-cells"  = "T")
Fig <- limma_subset %>%
  ggplot(aes(
    x = KO, 
    y = genes, 
    color = pmin(2, pmax(-2, logFC)),
    size = pmin(5, -log10(adj.P.Val))
  )) +
  geom_point() +
  
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = TeX("log_{2}(FC)")
  ) +
  scale_size_continuous(
    range = c(0, 1.8),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    title = "Differentially expressed ISGs and\ngrowth/metabolic genes",
    x = "Cell type",
    y = "Genes"
    #caption = "M: Macrophages    T: T-cells"
  ) +
  facet_wrap(
    ~dataset*celltype,      # << put facets side-by-side
    scales = "free_y", 
    nrow = 3
  ) +
  optimized_theme_fig()+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90),
    plot.caption = element_text(
      hjust = 0,          # Left-align
      size = 5, color = "black"    # Or "Times", "Courier", etc. (must be installed)
    ))
Fig
ggsave(basedir("Combined_dataset_ISG_int.pdf"), w = 27,
       h = 8, units = "cm")

library(dplyr)
library(ggplot2)
library(stringr)

# List of dataset names
datasets <- unique(limma_subset$dataset)

for (d in datasets) {
  
  df <- limma_subset %>% filter(dataset == d)
  
  Fig_d <- ggplot(df, aes(
    x = KO, 
    y = genes, 
    color = pmin(2, pmax(-2, logFC)),
    size = pmin(5, -log10(adj.P.Val))
  )) +
    geom_point() +
    scale_color_gradient2(
      low = "#4C889C",
      mid = "white",
      high = "#D0154E",
      name = TeX("log_{2}(FC)")
    ) +
    scale_size_continuous(
      range = c(0, 1.8),
      name = TeX("$-\\log_{10}(p_{adj})$")
    ) +
    labs(
      title = paste0("Differentially expressed ISGs and growth/metabolic genes Dataset: ", d),
      x = "Cell type",
      y = "Genes"
    ) +
    facet_grid(
      . ~ celltype,                # only facet by celltype now
      scales = "free",
      space = "free"
    ) +
    optimized_theme_fig() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90),
      plot.caption = element_text(hjust = 0, size = 5, color = "black")
    )
  
  # Print to screen
  print(Fig_d)
  
  # Save file with dataset name
  outfile <- basedir(paste0("ISG_plot_", d, ".pdf"))
  ggsave(outfile, Fig_d, width = 27, height = 8, units = "cm")
}


###############
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

summary_invivo_KO <- Int_glio%>%
  # group by interaction term
  mutate(KO = gsub("_in.vivo_noRT","", coef))%>%
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")

Fig_interaction <-  ggplot(summary_interaction_KO,aes(
  x = reorder(KO,signif_count),
  y = log10(signif_count)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
# Display the plot
Fig_interaction
ggsave(out("interaction_KO_number.pdf"),
       w= 12, h = 4, units = "cm")


###########

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
head(combined)
unique(combined$dataset)
write_rds(combined,basedir("GSEA_all_datasets_KO_interaction.rds"))

pathways_int <- c("Interferon Gamma Response", "Interferon Alpha Response",
                  "Cholesterol Homeostasis",
                  "mTORC1 Signaling","Myc Targets V2", "Myc Targets V1")

combined_selected <- combined %>%
  filter(pathway %in% pathways_int)

head(combined_selected)


library(ggplot2)
library(dplyr)

df <- combined_selected %>%
  mutate(logpadj = -log10(padj))

ggplot(df, aes(y = str_wrap(pathway,50), x = coef)) +
  geom_point(aes(size = pmin(logpadj, 3), color = pmax(-2,pmin(NES,2)))) +
  facet_grid(~celltype*dataset, scales = "free", space = "free") +
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
ggsave(basedir("Enrichment across_dataset and celltype.pdf"), w =38.5, h = 5.5, units = "cm")

#############
#
InDir4 <- dirout("Figure1")
genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
colnames(genes_fig_1) <- c("gene","pathways")

#Int pert
Int_pert <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(KO = gsub("interaction","",coef))%>%
 mutate(genes = ensg)%>%
  #filter(genes %in% genes_fig_1$gene)%>%
  dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes"))
Int_pert$dataset <- "Perturb-seq(Hematopoiesis)"







#Int_spleen
InDir1 <- dirout("Figure4")
Int_spleen <- read_rds(InDir1("combined_jakstat_diff_exp.rds"))
Int_spleen <- Int_spleen %>% filter(grepl("treatmentex_vivo",Int_spleen$coef))

Int_spleen$group <- ifelse(Int_spleen$logFC >= 1 & 
                           Int_spleen$adj.P.Val <= 0.05, "up", 
                         ifelse(Int_spleen$logFC <= -1 & 
                                  Int_spleen$adj.P.Val <= 0.05, "down", "n.s"))

# Modify the 'coef' column for any KO
Int_spleen <- Int_spleen %>%
  mutate(celltype = gsub("M", "Macrophage", cell_type)) %>%
  mutate(celltype = gsub("T8", "T-cells", cell_type)) %>%
  filter(celltype == "T-cells")%>%
  filter(coef != "treatmentex_vivo")%>%
  mutate(KO = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(KO = gsub("Interaction_","",KO))%>%
  mutate(genes = gene)%>%
  #filter(genes %in% genes_fig_1$gene)%>%
  dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
Int_spleen$dataset <-"Spleen"

summary_spleen_int <- Int_spleen %>%
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")),
            .groups = "drop")
selected_spleen_KOs <- summary_spleen_int%>%
  filter(signif_count > 10) %>%
  pull(KO) %>%
  unique()
Fig_interaction <-  ggplot(summary_spleen_int %>% filter(KO %in% selected_spleen_KOs),aes(
  x = reorder(KO,signif_count),
  y = log10(signif_count)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
# Display the plot
Fig_interaction

ggsave(basedir("interaction_KO_number_spleen.pdf"),
       w= 2.8, h = 3, units = "cm")
#########
InDir <- dirout("Glioblastoma_limmaRes.3way.mod")

limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))

limmaRes_int <-  limmaRes %>% filter(coef %in% grep("tissueex.vivo.genotype[A-Za-z0-9]*$",
                                                    limmaRes$coef, value = T))
limmaRes_int$celltype <- "glio"

Int_glio <- limmaRes_int %>%
  mutate(KO = gsub("tissueex.vivo.genotype","", coef))%>%
  mutate(genes = ensg)%>%
  #filter(genes %in% genes_fig_1$gene)%>%
  dplyr::select(c("KO","celltype","logFC","adj.P.Val","genes","group"))
Int_glio$dataset <- "Glioblastoma"

summary_invivo_KO <- Int_glio%>%
  # group by interaction term
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")
selected_KOs <-summary_interaction_KO%>%
  filter(signif_count> 10)%>%
  pull(KO)%>%
  unique()
Fig_interaction <-  ggplot(summary_interaction_KO %>% filter(KO %in% selected_KOs),aes(
  x = reorder(KO,signif_count),
  y = log10(signif_count)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
# Display the plot
Fig_interaction
ggsave(basedir("interaction_KO_number.pdf"),
       w= 8, h = 3, units = "cm")




limma_subset <- rbind(Int_glio,Int_pert,Int_spleen) 

# limma_subset$cell_type <- recode(limma_subset$cell_type,"Macrophage" = "M",
#                                  "T-cells"  = "T")
Fig <- limma_subset %>%
  ggplot(aes(
    x = KO, 
    y = genes, 
    color = pmin(2, pmax(-2, logFC)),
    size = pmin(5, -log10(adj.P.Val))
  )) +
  geom_point() +
  
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = TeX("log_{2}(FC)")
  ) +
  scale_size_continuous(
    range = c(0, 1.8),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    title = "Differentially expressed ISGs and\ngrowth/metabolic genes",
    x = "Cell type",
    y = "Genes"
    #caption = "M: Macrophages    T: T-cells"
  ) +
  facet_grid(
    dataset~celltype,      # << put facets side-by-side
    scales = "free", 
    space = "free"
  ) +
  optimized_theme_fig()+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90),
    plot.caption = element_text(
      hjust = 0,          # Left-align
      size = 5, color = "black"    # Or "Times", "Courier", etc. (must be installed)
    ))
Fig
ggsave(basedir("Combined_dataset_ISG_int.pdf"), w = 27,
       h = 8, units = "cm")

library(dplyr)
library(ggplot2)
library(stringr)

# List of dataset names
datasets <- unique(limma_subset$dataset)

for (d in datasets) {
  
  df <- limma_subset %>% filter(dataset == d)
  
  Fig_d <- ggplot(df, aes(
    x = KO, 
    y = genes, 
    color = pmin(2, pmax(-2, logFC)),
    size = pmin(5, -log10(adj.P.Val))
  )) +
    geom_point() +
    scale_color_gradient2(
      low = "#4C889C",
      mid = "white",
      high = "#D0154E",
      name = TeX("log_{2}(FC)")
    ) +
    scale_size_continuous(
      range = c(0, 1.8),
      name = TeX("$-\\log_{10}(p_{adj})$")
    ) +
    labs(
      title = paste0("Differentially expressed ISGs and growth/metabolic genes Dataset: ", d),
      x = "Cell type",
      y = "Genes"
    ) +
    facet_grid(
      . ~ celltype,                # only facet by celltype now
      scales = "free",
      space = "free"
    ) +
    optimized_theme_fig() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90),
      plot.caption = element_text(hjust = 0, size = 5, color = "black")
    )
  
  # Print to screen
  print(Fig_d)
  
  # Save file with dataset name
  outfile <- basedir(paste0("ISG_plot_", d, ".pdf"))
  ggsave(outfile, Fig_d, width = 27, height = 8, units = "cm")
}

################
Fig <- limma_subset %>%
  filter(genes %in% goi)%>%
  ggplot(aes(
    x = KO,
    y = genes,
    color = pmin(2, pmax(-2, logFC)),
    size = pmin(5, -log10(adj.P.Val))
  )) +
  geom_point() +
  
  scale_color_gradient2(
    low = "#4C889C", mid = "white", high = "#D0154E",
    name = TeX("log_{2}(FC)")
  ) +
  scale_size_continuous(
    range = c(0, 1.5),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  
  facet_grid(
    ~dataset*celltype,
    scales = "free_x",
    space  = "free_x"     # <<< equal spacing, but dropping unused KOs
  ) +
  
  labs(
    title = "Differentially expressed ISGs and\ngrowth/metabolic genes",
    x = "Cell type",
    y = "Genes"
  ) +
  
  optimized_theme_fig() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90),
    plot.caption = element_text(hjust = 0, size = 5, color = "black")
  )
Fig
ggsave(basedir("ISGs_across_datasets_2.pdf"), w = 35.6, h = 8.2, units = "cm")  
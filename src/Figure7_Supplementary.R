source("src/00_init.R")

# renv::snapshot(lockfile = "renv_NF.lock")
source("~/code/resources/RFunctions/Basics.R")
source("src/Ag_Optimized_theme_fig.R")
require(tidyverse)
require(data.table)
require(enrichR)
library(dplyr)
require(latex2exp)
library(patchwork)


out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"
InDir <- dirout("Figure4")
InDir2 <- dirout("../Ag_ScRNA_22_JAKSTAT_Ar/")
basedir <- dirout("Figure7_Supplementary")


selected_coef <- c("Interaction_STAT1KO","Interaction_STAT2KO",
                   "Interaction_TYK2CMV", "Interaction_IRF9KO")
selected_KO <- gsub("Interaction_","",selected_coef)


#############################################
gsea.res <- read_rds(InDir2("fgsea_hom_vs_ex.vivo_per_CT.rds"))
unique(gsea.res$coef)

FGSEA <- gsea.res
pDT <- FGSEA %>%
  #filter(coef %in%)
  mutate(genotype = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(celltype = gsub("M", "Macrophage", celltype)) %>%
  mutate(celltype = gsub("T8", "T-cells", celltype))%>%
  filter(coef != "(Intercept)")
unique(pDT$genotype)  
#pDT$genotype <- gsub("Interaction_","",pDT$genotype)
pDT$genotype <- gsub("treatmentex_vivo","WT",pDT$genotype )
# Step 2: Summarize to find KOs with at least one valid cell type
pDT <- pDT %>%
  filter(genotype %in% c(grep("Interaction",pDT$genotype, value = T),"WT"))
db = "MSigDB_Hallmark_2020"
pDT <- pDT %>%
  filter(db == "MSigDB_Hallmark_2020") 

# Step 3 continued: Keep only valid KOs for the specific cell type
pDT <- pDT %>% filter( padj < 0.05)

# mutate(alpha_value = if_else(valid_ko, 1, 0))  # Set alpha based on validity
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5),by=c("coef", "celltype","pathway")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("coef", "celltype","pathway")]$pathway)
# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(union(pw.display.pos, pw.display.neg))
pDT <- pDT[pathway %in% pw.display]
# Remove duplicate rows
pDT <- pDT %>% distinct()
pDT_agg <- pDT %>%
  group_by(pathway) %>%
  summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
  arrange(desc(average_NES))  # Ordering pathways by the average NES, highest first

# Step 2: Create a factor for pathway that reflects the aggregated NES order
pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)


#pDT$coef <- gsub("Interaction_","",pDT$genotype)
pDT <- pDT %>%
  filter(genotype %in% c(selected_coef, "WT"))
fig7A_supplementary <- ggplot(pDT, aes(x = gsub("Interaction_","",genotype), y = pathway,
                                       color = NES,
                                       size = pmin(5, -log10(padj)))) +
  geom_point() + 
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name=TeX("log_{2}(FC)")) +
  geom_point(data = pDT[padj < 0.05], shape = 1) +
  scale_size_continuous(
    range = c(0, 2),
    limits = c(0, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  theme_bw() +
  xRot() +
  labs(x="KOs",
       title = "Enriched pathways in culture effect and interaction effect")+
  facet_grid(cols= vars(celltype),scales = "free",space = "free") +  # Create facets for each cell type
  theme(strip.text.y = element_text(angle = 0)) +
  optimized_theme_fig() +
  theme(
    strip.text.x = element_text(angle = 90)
  )
fig7A_supplementary
ggsave(basedir("Fig7A_supplementary.pdf"), 
       plot = fig6A_supplementary, 
       width = 8, 
       height = 11.5, 
       units = "cm")

#Sup.Fig.6A-----------
filtered_data <- read_rds(InDir("filtered_data.rds"))
goi <- "Oas1a"
boxplot_jitter <- ggplot(filtered_data,
                         aes(x = genotype, y = E, color = tissue)) + 
  # Boxplot with tissue color fill
  geom_boxplot(
    outlier.colour = NA,
    position = position_dodge(width = 0.8) # Adjust width to create space between tissues
    
  ) + 
  # geom_jitter(data = filtered_data%>%filter(genotype != "WT"),
  #             
  #             position = position_jitterdodge(
  #               jitter.width = 0.2,
  #               dodge.width = 0.8  # Adjust dodge width to match the boxplot
  #             ),
  #             alpha = 0.5)+
  
  facet_grid(rows = vars(tissue),cols = vars(cell_type), scales = "free",
             labeller = labeller(
               tissue = c("ex_vivo" = "Ex vivo",
                          "in_vivo" = "In vivo")
             ))+
  
  scale_color_manual(
    values = c("ex_vivo" = "#6a3d9aff", "in_vivo" = "#d38d5fff"),
    name = "Experimental model",
    labels = c("ex_vivo" = "Ex vivo", "in_vivo" = "In vivo")
  ) +
  optimized_theme_fig() +
  labs(
       x = "KOs",
       y = "Expression",
       title = "Difference in KO-effects in ex vivo cultured cells\n(Gene: Oas1a)") +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#supplementaryfig4B----
ggsave(basedir("SupplementaryFig7B_no_jitter.pdf"), 
       plot = boxplot_jitter, 
       width = 12, 
       height = 6, 
       units = "cm")

###################################################################



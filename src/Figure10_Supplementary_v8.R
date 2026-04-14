###############
source("src/00_init.R")
#source("src/Ag_enrichR_mouse_genes.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(biomaRt)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
library(GEOquery)
library(R.utils)
library(msigdbr)
library(fgsea)
library(readxl)
library(latex2exp)
library(Biobase)
library(Matrix)
library(circlize)
library(grid)

InDir1 <- dirout("Ag_top_pathway_genes")
InDir <- dirout("Glioblastoma_limmaRes.3way.mod")
basedir <- dirout("Figure10_Supplementary")
limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))


limmaRes_NTC <- limmaRes %>% filter(coef == "tissueex.vivo")
limmaRes_int <-  limmaRes %>% filter(coef %in% grep("tissueex.vivo.genotype[A-Za-z0-9]*$",
                                                    limmaRes$coef, value = T))

limmaRes_NTC$KO <- "NTC"
limmaRes_NTC$coef <- "NTC"

head(limmaRes_int)

top_NTC <- limmaRes_NTC %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(desc(abs(logFC)))%>%
  slice_head(n = 50)

top_genes <-  limmaRes_NTC %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(desc(abs(logFC)))%>%
  slice_head(n = 50)%>%
  pull(ensg)%>%
  unique() 

colnames(limmaRes_NTC) %in% colnames(limmaRes_int)

#   

limmaRes_RT_int <- limmaRes %>%
  filter(coef == "tissueex.vivo.RT_statusRT")
top_NTC$celltype <-"glio"
limmaRes_RT_int$KO <- "Radiotherapy"
limmaRes_RT_int$celltype <- "glio"
limmaRes_int$celltype <- "glio"
limmaRes_int$KO <- gsub("tissueex.vivo.genotype","",limmaRes_int$coef)

top_NTC_int <- top_NTC %>%
  rbind(limmaRes_int %>% filter(ensg %in% top_NTC$ensg))%>%
  rbind(limmaRes_RT_int %>% filter(ensg %in% top_NTC$ensg))

top_NTC_int <- annotate_pathway(top_NTC_int, pathways)
top_NTC_int$KO <- factor(top_NTC_int$KO, 
                         levels = c("NTC",
                                    setdiff(top_NTC_int$KO, "NTC")))

# fig---------
plot <- ggplot(top_NTC_int , aes(x = KO, y = ensg,
                                 color = pmin(2, pmax(-2, logFC)),
                                 size = pmin(3, -log10(adj.P.Val))))+  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0,2),
    #limits = c(0,5),
    #breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(
    title = str_wrap("Interaction effect (KO/RT) - culture model) of top culture effect genes", width = 80),
    x = "NTC/Interaction",
    y = "Genes")+
  #facet_grid(cols = vars(pathway), scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()+theme(
    legend.position = "right",
    strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
    axis.text.x = element_text(angle = 45)
  )

plot
ggsave(basedir("top_NTCs_interction.pdf"), w = 13, h = 12.5, units = "cm")
######################

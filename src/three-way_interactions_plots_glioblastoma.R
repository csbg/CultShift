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
basedir <- dirout("three-way_interactions_plots_glioblastoma")
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


######################
#plot
Fig_interaction_threeway <-  ggplot(summary_interaction_threeway,aes(
  x = KO,
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
Fig_interaction_threeway
ggsave(out("No.of_genes_3way_interaction.pdf"))
#############
#correlation 3way------------
unique(limmaRes_3way$coef)
merged_data <- limmaRes_3way %>%
  inner_join(limmaRes_NTC, by = c("ensg"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

# Step 1: Calculate correlation with p-value for each KO and celltype
correlation_results <- merged_data %>%
  mutate(KO = gsub("tissueex.vivo.genotype|.RT_statusRT","",coef.x))%>%
  inner_join(summary_interaction_KO, by = c("KO")) %>%
  #filter(signif_count > 10) %>%
  group_by(KO)%>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results <- correlation_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Step 3: Add significance labels based on p-value thresholds
correlation_results <- correlation_results %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )
#order
correlation_results <- correlation_results %>%
  mutate(KO = factor(KO, levels = KO_order))
Correlation_3way_to_NTC <- ggplot(correlation_results, aes(x = KO, y = cor_abs)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  geom_text(aes(label = significance), 
            y =  0.25, 
            color = "black",
            size = 1.5) + 
  labs(
    title = "Correlation of 3-way interaction effect to culture effect",
    x = NULL)+
  # fill = "Cell Type")+
  ylab(expression(atop("Correlation of interaction effect", "to culture effects")))+
  
  
  optimized_theme_fig()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
Correlation_3way_to_NTC
ggsave(out("correlation_3way_to_ntc.pdf"),
       w= 8, h =4)

#########
# ISGs with interaction effects-----------
InDir <- dirout("Ag_top_pathway_genes")

genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))
genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

####
source("src/Ag_enrichR_mouse_genes.R")
pathways <- list(
  ISG_core_all = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
    filter(L1=="ISG_Core")%>%pull(value),
  Cholesterol_all = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  mTORC1_all = enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`,
  ISGs_all = union(enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`,
                   enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`),
  ROC_all = enr.terms$MSigDB_Hallmark_2020$`Reactive Oxygen Species Pathway`,
  
  #Cholesterol = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  Hypoxia_all = enr.terms$MSigDB_Hallmark_2020$`Hypoxia`,
  Glycolysis_all = enr.terms$MSigDB_Hallmark_2020$`Glycolysis`,
  Protein_loc_all = Reduce(union, list(
    enr.terms$MSigDB_Hallmark_2020$`Myc Targets V1`,
    enr.terms$GO_Biological_Process_2023$`Protein Import (GO:0017038)`,
    enr.terms$GO_Biological_Process_2023$`Protein Insertion Into ER Membrane (GO:0045048)`,
    enr.terms$GO_Biological_Process_2023$`Protein Insertion Into Membrane (GO:0051205)`,
    enr.terms$GO_Biological_Process_2021$`RNA biosynthetic process (GO:0032774)`
  ))
)
annotate_pathway <- function(df, pathways_list) {
  df %>%
    mutate(
      pathway = case_when(
        ensg %in% pathways_list$ISGs_all ~ "ISGs_all",
        ensg %in% pathways_list$ISG_core_all ~ "ISG_core_all",
        ensg %in% pathways_list$mTORC1_all ~ "mTORC1_all",
        ensg %in% pathways_list$Cholesterol_all ~ "Cholesterol_all",
        ensg %in% pathways_list$ROC_all ~ "ROS_all",
        ensg %in% pathways_list$Hypoxia_all ~ "Hypoxia_all",
        ensg %in% pathways_list$Glycolysis_all ~ "Glycolysis_all",
        ensg %in% pathways_list$Protein_loc_all ~ "Protein_localization_all",
        TRUE ~ "Other"
      )
    )
}
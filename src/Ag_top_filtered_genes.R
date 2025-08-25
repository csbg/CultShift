source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification_Mye.R")
source("src/Ag_enrichR_mouse_genes.R")
library(tidyverse)
library(enrichR)
library(purrr)
basedir <- dirout("Ag_top_filtered_genes")
limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
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
genes_fig1 %>% write_rds(basedir("genes_fig1.rds"))


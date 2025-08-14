###############
source("src/00_init.R")
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
#library(fgsea)
#library(msigdbr)
#library(dplyr)
library(ggplot2)

outdir <- dirout("Figure1_2/")
enrich <- dirout("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/Enrichr")
#load data
inDir1<-dirout_load("Ag_ScRNA_09_Normalized_counts_ex_in_izzo_NTC")
inDir2<-dirout_load("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/ENRICHR")
inDir3<-dirout_load("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
base<-"Ag_ScRNA_16_per_celltype_gene_plots_guide/"
basedir<-dirout("Ag_ScRNA_16_per_celltype_gene_plots_guide")

################################################################################
#load limma results from NTC

limmaRes_NTC<-read_rds(inDir3("limma_perCTex.vivovsin.vivo.rds"))
#load enrichment_results from NTC
enr.res.all_up<-read_rds(inDir2("enr.res.all_NTC_up.rds")) 
enr.res.all_down<-read_rds(inDir2("enr.res.all_NTC_down.rds"))
#####################
#NTC
################################################################################
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

pathway_names <- list(
  "Cholesterol Homeostasis",
  "mTORC1 Signaling",
  "Interferon Alpha Response",
  "Interferon Gamma Response",
  "Inflammatory Response",
  "p53 Pathway",
  "Myc Targets V1",
  "IL−6/JAK/STAT3 Signaling",
  "Oxidative Phosphorylation",
  "TGF−beta Signaling",
  "Unfolded Protein Response",
  "E2F Targets",
  "Protein Secretion",
  "ISG_core"  # Include ISG_core as a pathway
)

#pathway_names <- names(enr.terms$MSigDB_Hallmark_2020)
pathway_genes <- list()

# Loop through each pathway name and extract genes
for (pathway_name in pathway_names) {
  # Check if pathway_name exists in enr.terms$MSigDB_Hallmark_2020
  if (pathway_name %in% names(enr.terms$MSigDB_Hallmark_2020)) {
    # Extract genes for the pathway
    pathway_genes[[pathway_name]] <- enr.terms$MSigDB_Hallmark_2020[[pathway_name]]
  } else {
    # Handle case where pathway_name does not exist
    warning(paste("Pathway", pathway_name, "not found in enrichment results. Skipping."))
  }
}
ISG_core<-read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)
pathway_genes[["ISG_Core"]] <- ISG_core

# Function to calculate overlap percentage
calculate_overlap_percentage <- function(pathway_name, pathway_genes, limmaRes) {
  # Get the genes for the current pathway
  genes_in_pathway <- pathway_genes[[pathway_name]]
  
  # Initialize a list to store results for each group
  overlap_counts_list <- list()
  
  # Calculate overlap for each group ('up', 'down')
  for (group_name in c("up", "down")) {
    # Filter genes for the specific group and pathway
    overlap_genes <- limmaRes %>%
      filter(group == group_name) %>%
      filter(toupper(genes) %in% toupper(genes_in_pathway))
    
    # Calculate the percentage of overlap per cell type
    overlap_percentage <- overlap_genes %>%
      group_by(celltype) %>%
      summarise(count = n()) %>%
      mutate(
        total_genes = length(genes_in_pathway),
        percentage_overlap = (count / total_genes) * 100,
        pathway = pathway_name,
        group = group_name
      )
    
    # Store the overlap percentage in the list
    overlap_counts_list[[group_name]] <- overlap_percentage
  }
  
  # Combine the results for all groups into one data frame
  overlap_counts_df <- bind_rows(overlap_counts_list)
  return(overlap_counts_df)
}
# Calculate overlaps for all pathways
overlap_results_list <- map(names(pathway_genes), 
                            calculate_overlap_percentage, 
                            pathway_genes = pathway_genes, 
                            limmaRes = limmaRes_NTC)
overlap_all <- bind_rows(overlap_results_list)

# Plot the overlap percentage
ggplot(overlap_all, aes(x = celltype, y = percentage_overlap, fill= group)) +
  geom_bar(stat = "identity", 
           width = 0.8) +
  scale_fill_manual(values = c("down" =  "#4C889C", "up" = "#D0154E")) +
  labs(#title = "Overlap Percentage for MSigDB Hallmark 2020 Pathways",
    x = "Cell Type",
    y = "Percentage of Overlapping Genes",
    fill = "Group") +
  facet_grid(rows = vars(group), cols = vars(pathway)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 90),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text.x = element_text(size = 14, angle = 90, hjust = 0),
        strip.text.y = element_text(size = 14, angle = 360, hjust = 1))
ggsave(outdir("percentage_genes_all.pdf"),width = 30,h = 10)

################################################################################
#actual counts
################################################################################
# Function to calculate overlap percentage
calculate_overlap_count <- function(pathway_name, pathway_genes, limmaRes) {
  # Get the genes for the current pathway
  genes_in_pathway <- pathway_genes[[pathway_name]]
  
  # Initialize a list to store results for each group
  overlap_counts_list <- list()
  
  # Calculate overlap for each group ('up', 'down')
  for (group_name in c("up", "down")) {
    # Filter genes for the specific group and pathway
    overlap_genes <- limmaRes %>%
      filter(group == group_name) %>%
      filter(toupper(genes) %in% toupper(genes_in_pathway))
    
    # Calculate the percentage of overlap per cell type
    overlap_count <- overlap_genes %>%
      group_by(celltype) %>%
      summarise(count = n()) %>%
      mutate(
        total_genes = length(genes_in_pathway),
        pathway = pathway_name,
        group = group_name
      )
    
    # Store the overlap percentage in the list
    overlap_counts_list[[group_name]] <- overlap_count
  }
  
  # Combine the results for all groups into one data frame
  overlap_counts_df <- bind_rows(overlap_counts_list)
  return(overlap_counts_df)
}

# Calculate overlaps for all pathways
overlap_results_list <- map(names(pathway_genes), 
                            calculate_overlap_count, 
                            pathway_genes = pathway_genes, 
                            limmaRes = limmaRes_NTC)
overlap_all <- bind_rows(overlap_results_list)

# Plot the overlap percentage
ggplot(overlap_all, aes(x = celltype, y = count, fill= group)) +
  geom_bar(stat = "identity", 
           width = 0.8) +
  scale_fill_manual(values = c("down" =  "#4C889C", "up" = "#D0154E")) +
  labs(#title = "Overlap Percentage for MSigDB Hallmark 2020 Pathways",
    x = "Cell Type",
    y = "Percentage of Overlapping Genes",
    fill = "Group") +
  facet_grid(rows = vars(group), cols = vars(pathway)) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text.x = element_text(size = 15, angle = 90, hjust = 1),
        strip.text.y = element_text(size = 15, angle = 360, hjust = 1))
ggsave(outdir("count_genes.pdf"),width = 20,h = 10)

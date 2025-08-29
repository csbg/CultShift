
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
require(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
enrich<-dirout("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/Enrichr")
################################################################################
#function1
#Function to extract genes from terms
extract_unique_genes <- function(enrichr_data,
                                 terms,
                                 db) {
  unique_genes <- enrichr_data %>%
    filter(db == db & Term %in% terms) %>%
    pull(Genes) %>%
    unique() %>%
    strsplit(";") %>%
    unlist() %>%
    toupper() %>%
    unique()
  return(unique_genes)
}
################################################################################
#function2

# General plotting function for specific terms and datasets
plot_genes_by_term <- function(limma_data, 
                               enrichr_data_up,
                               enrichr_data_down,
                               term,
                               termname,
                               db,
                               output_file_prefix,
                               list_coef,dir_name) {
  # Combine upregulated and downregulated enrichment results
  combined_enrichr_data <- bind_rows(enrichr_data_up, enrichr_data_down)
  
  # Extract unique genes for the given term
  unique_genes <- extract_unique_genes(combined_enrichr_data, term, db)
  
  # Filter limma data for the extracted genes and the given coefficients
  genes_KO <- limma_data %>%
    filter(coef %in% list_coef) %>%
    filter(toupper(ensg) %in% unique_genes)
  
  # Check if genes_KO has more than one gene
  if (nrow(genes_KO) == 0) {
    cat("No data to plot for term:", term, "\n")
    return(NULL)
  }
  
  # Create the plot
  plot <- ggplot(genes_KO, aes(y = gsub("ex.vivo:", "", coef), x = ensg, size = -log10(adj.P.Val), color = logFC)) +
    geom_point(alpha = 0.6) +
    scale_color_gradient2(high = "red", low = "blue") +
    theme(axis.text = element_text(size = 15)) +
    facet_wrap(vars(celltype), scales = "free_x") +
    theme_bw(12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 8),
          axis.text.y = element_text(size = 8)) +
    scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
    ggtitle(paste(term, "- Interaction(ex.vivo-in.vivo)")) +
    ylab(paste("Knockouts"))+
    xlab("genes")
  plot
  # Create directory if it doesn't exist
  
  dir<-dirout(paste0(base,dir_name,termname))

  # Save the plot
  ggsave(dir(paste0(output_file_prefix,"_",termname,"_summary.pdf")),
         plot = plot, width = length(unique_genes) * 0.2 + 5,
         height = length(unique_genes) * 0.18 + 1,
         limitsize = FALSE)
}


#s
################################################################################

#load data
inDir<-dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
inDir1<-dirout_load("Ag_ScRNA_09_Normalized_counts_ex_in_izzo_NTC")
inDir2<-dirout_load("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/ENRICHR")
inDir3<-dirout_load("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
base<-"Ag_ScRNA_16_per_celltype_gene_plots_guide/"
basedir<-dirout("Ag_ScRNA_16_per_celltype_gene_plots_guide")

IFN_inter<-dirout(paste0(base,"IFN_inter"))
IFN_base<-"Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/"
IFN_NTC<-dirout(paste0(base,"IFN_NTC"))
IFN_NTC_base<-"Ag_ScRNA_16_per_celltype_gene_plots_guide/IFN_NTC"
################################################################################
#load limma results from NTC
dataVoom_NTC<-read_rds(inDir1("dataVoom_perCT_NTC_izzo_ex_in.rds"))
dataVoom_NTC_in_ex<-read_rds(inDir3("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaRes_NTC<-read_rds(inDir3("limma_perCTex.vivovsin.vivo.rds"))
#load enrichment_results from NTC
enr.res.all_up<-read_rds(inDir2("enr.res.all_NTC_up.rds"))%>% 
enr.res.all_down<-read_rds(inDir2("enr.res.all_NTC_down.rds"))
#load enrichment_results from interaction
down_enrichr_filtered_KO<-read_rds(enrich("count_above10_down_logFC_1_enrichr.rds"))
up_enrichr_filtered_KO<-read_rds(enrich("count_above10_up_logFC_1_enrichr.rds"))
#load_limma across KO
limmaRes<-read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
limmaRes<-limmaRes%>%filter(celltype != "MEP")
dataVoom_Eo.Ba<-read_rds(inDir("Eo.Ba_dataVoom.rds"))
dataVoom_Mono<-read_rds(inDir("Mono_dataVoom.rds"))
dataVoom_MkP<-read_rds(inDir("MkP_dataVoom.rds"))
dataVoom_GMP<-read_rds(inDir("GMP_dataVoom.rds"))
dataVoom_HSC<-read_rds(inDir("HSC_dataVoom.rds"))
dataVoom_MEP<-read_rds(inDir("MEP_dataVoom.rds"))
NTC_meta<-read_rds(inDir1("meta_all_NTC.rds"))
NTC_meta_in_ex<-read_rds(inDir3("NTC_meta.rds"))
meta<-read_rds(inDir("meta.rds"))
################################################################################
#NTC
################################################################################


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
  "Cholesterol"  # Include Cholesterol as a pathway
)

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
Cholesterol<-read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))
Cholesterol<-Cholesterol%>%filter(L1=="Cholesterol")%>%pull(value)
pathway_genes[["Cholesterol"]] <- Cholesterol
pathway_names
###########################
#from NTC
##########################

################################################################################
#NTC including izzo
############################
#set directory
#fig_1.3
NTC_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots_guide/Pathways_genes_NTC/"))
#
# Assuming you have the genes for Interferon Gamma and Interferon Beta responses
interferon_gamma_genes <- enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`
interferon_beta_genes <- enr.terms$MSigDB_Hallmark_2020$`Interferon Beta Response`

# Combine both gene sets (if you want to consider overlaps between both responses)
interferon_genes <- union(interferon_gamma_genes, interferon_beta_genes)
length(interferon_genes)
length(Cholesterol)
#overlap_ifn <- toupper(interferon_genes)[toupper(interferon_genes) %in% toupper(Cholesterol)]
# Filter limmaRes for logFC < 1 and overlap with interferon genes

Cholesterol<-read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))
Cholesterol<-Cholesterol%>%filter(L1=="Cholesterol")%>%pull(value)
#####################
# Assuming you have already loaded Cholesterol genes and Interferon genes

# Filter limmaRes_NTC for Cholesterol genes
Cholesterol_genes <- limmaRes_NTC %>% 
  filter(group != "n.s")%>%
  group_by(celltype) %>%
  filter(logFC < 1, adj.P.Val < 0.05) %>%
  filter(toupper(genes) %in% toupper(Cholesterol))%>%
  arrange(logFC)%>%
  slice_head(n = 5) %>%
  pull(genes)
ISG_plot <-  limmaRes_NTC %>% 
  #filter(group != "n.s")%>%
  #filter(logFC < 1, adj.P.Val < 0.05) %>%
  filter(toupper(genes) %in% toupper(Cholesterol_genes))
# Extract Cholesterol Homeostasis genes
cholesterol_homeostasis_genes <- enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`
cholesterol_genes <- limmaRes_NTC %>% 
  filter(group != "n.s")%>%
  group_by(celltype) %>%
  filter(logFC > 1, adj.P.Val < 0.05) %>%
  filter(toupper(genes) %in% toupper(cholesterol_homeostasis_genes))%>%
  arrange(desc(logFC))%>%
  slice_head(n = 5) %>%
  pull(genes)
cholesterol_plot <-  limmaRes_NTC %>% 
  # filter(group != "n.s")%>%
  # filter(logFC > 1, adj.P.Val < 0.05) %>%
  filter(toupper(genes) %in% toupper(cholesterol_genes))
combined_genes_filtered <- bind_rows(
  cholesterol_plot %>% mutate(gene_set = "Cholesterol Homeostasis"),
  ISG_plot %>% mutate(gene_set = "ISG Core")
)

# Filter limmaRes_NTC for Interferon genes

# Plotting the dot plot
ggplot(combined_genes_filtered, aes(x = celltype, y = genes,
                                    color = logFC,
                                    size = pmin(30,-log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(
    range = c(1, 5),  # Set the actual size range from 2 to 30
    breaks = c(1, 5, 10, 20, 30)  # Set specific breaks to create distinct point sizes
  ) +
  labs(title = "ISG Core and Cholesterol Homeostasis Genes",
       x = "Cell Type",
       y = "Genes",
       color = "logFC",
       size = "-log10(adj.P.Val)") +
  facet_grid(rows = vars(gene_set), scales = "free_y", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10),
    strip.text.y.right = element_text(angle = 90, size = 12, face = "bold"),  # Rotate gene set labels to 90 degrees
    strip.placement = "outside",  # Move the facet labels outside the plot
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.background = element_blank()  # Remove the background from the strip
  ) +
  theme(strip.switch.pad.grid = unit(0.1, "cm"))  # Adjust padding between strip and plot

# Save the plot
ggsave(NTC_dir("Cholesterol_and_Cholesterol_homeostasis_logFC.pdf"), width = 12, height = 8)

#figure 1.4
#example genes
example <- c("Idi1","Irf7","Ubc","Gbp2","Lmna","Tbp","Tubb")
example <- intersect(example, rownames(dataVoom_NTC_in_ex$E))
dat.list<-list()
for(gg in unique(example)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta_in_ex %>%
    mutate(E = scale(dataVoom_NTC_in_ex$E[gg,])) %>%
    rownames_to_column("samples") %>%
    remove_rownames()
  
  # Group the data by tissue and celltype, and calculate the average E for each group
  #avg_gene_data <- gene_data %>%
  #group_by(tissue, celltype) %>%
  #summarise(avg_E = mean(E))
  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")
######################################
geneset<-"example_NTC"
create_gene_plots_NTC <- function(data, gene, geneset) {
  # Define the directory to save the plots
  dir <- dirout(paste0("Ag_ScRNA_16_per_celltype_gene_plots_guide/NTC/"))
  
  # Create the plot
  ggplot(data[data$gene == gene, ], aes(x = celltype, y = E, fill = tissue)) + 
    geom_boxplot(
      outlier.shape = NA, 
      position = position_dodge(width = 0.8),  # Adjust the width to create space between tissues
      color = "black", 
      size = 0.8
    ) +  # Boxplot with tissue color fill
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.2, 
        dodge.width = 0.8  # Adjust the dodge width to match the boxplot
      ), 
      alpha = 0.6
    ) +  # Jittered points
    facet_grid(cols = vars(celltype), scales = "free") +
    labs(title = gene) +
    xlab(NULL) +  # Remove x-axis label
    ylab("Scaled Gene Expression") +
    scale_fill_manual(values = c("#214C73", "#9E7F9C")) +  # Custom colors for tissues
    theme_minimal() +  # A cleaner theme
    theme(
      axis.text.x = element_blank(),  # Hide x-axis text
      axis.ticks.x = element_blank(),  # Hide x-axis ticks
      legend.title = element_text(size = 10),  # Legend title size
      legend.text = element_text(size = 8),  # Legend text size
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.y = element_text(size = 12)
    )
  
  # Save the plot
  ggsave(dir(paste0(geneset, "_", gene, ".pdf")), width = 10, height = 6)
}


for (gg in example){
  data<-dat.list%>%filter(gene==gg)%>%
    filter(celltype %in% c("Eo.Ba","HSC","MkP","Mono","GMP"))
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene,"example_NTC")
    
  })
}
################################################################################
# Define the function to get top genes for a given pathway
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


# Define your pathways of interest
pathways <- list(
  mTORC1 = enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`,
  #Cholesterol = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  Hypoxia = enr.terms$MSigDB_Hallmark_2020$`Hypoxia`,
  Glycolysis = enr.terms$MSigDB_Hallmark_2020$`Glycolysis`
  # Add other pathways here as needed
)

# Initialize an empty list to store results
results <- list()

# Loop through each pathway and apply the function
for (pathway in names(pathways)) {
  results[[pathway]] <- get_top_genes(
    pathway_name = pathway,
    pathway_genes = pathways[[pathway]],
    limma_results = limmaRes_NTC
  )
}

# Access the results for each pathway
combined_genes_filtered <- bind_rows(
  results$mTORC1$pathway_plot %>% mutate(gene_set = "mTORC1"),
  #results$Cholesterol$pathway_plot%>%mutate(gene_set = "Cholesterol"),
  results$Hypoxia$pathway_plot%>%mutate(gene_set = "Hypoxia"),
  results$Glycolysis$pathway_plot%>%mutate(gene_set = "Glycolysis")
)
# Plotting the dot plot
ggplot(combined_genes_filtered, aes(x = celltype, y = genes,
                                    color = logFC,
                                    size = pmin(30,-log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(
    range = c(1, 5),  # Set the actual size range from 2 to 30
    breaks = c(1, 5, 10, 20, 30)  # Set specific breaks to create distinct point sizes
  ) +
  labs(title = "ISG Core and Cholesterol Homeostasis Genes",
       x = "Cell Type",
       y = "Genes",
       color = "logFC",
       size = "-log10(adj.P.Val)") +
  facet_grid(rows = vars(gene_set), scales = "free_y", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10),
    strip.text.y.right = element_text(angle = 90, size = 12, face = "bold"),  # Rotate gene set labels to 90 degrees
    strip.placement = "outside",  # Move the facet labels outside the plot
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.background = element_blank()  # Remove the background from the strip
  ) +
  theme(strip.switch.pad.grid = unit(0.1, "cm"))  # Adjust padding between strip and plot

################################################################################

# Call the function with your data
plot_genes_by_term(limmaRes, down_enrichr_filtered_KO,
                   term, db, output_file_prefix, 
                   list_coef)



######################################
# Iterate over each unique gene
Cholesterol <- intersect(Cholesterol, rownames(dataVoom_NTC$E))
dat.list<-list()
for(gg in unique(Cholesterol)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta %>%
    mutate(E = scale(dataVoom_NTC[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
  
  # Group the data by tissue and celltype, and calculate the average E for each group
  #avg_gene_data <- gene_data %>%
    #group_by(tissue, celltype) %>%
    #summarise(avg_E = mean(E))
  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")
head(dat.list)

##########################################
create_gene_plots_NTC <- function(data, gene) {
  ggplot(data[data$gene == gene ,], aes(x = tissue, y = E)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(rows = vars(celltype), scales = "free") +
    labs(title = gene)+
    xlab(paste0(gene))+
    ylab(paste0("scaled-Normalized gene exp"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                                      hjust = 1))
  
  
    
  ggsave(ko_dir(paste0(gene,".pdf")))
}
head(dat.list)
#For  
for (gg in IFN_genes){
  data<-dat.list%>%filter(gene==gg)%>%
  filter(celltype %in% c("Eo.Ba","HSC","MkP","Mono","GMP"))
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene)
    
  })
}
################################################################################

# Call the function with your data
plot_genes_by_term(limmaRes, down_enrichr_filtered_KO,
                   term, db, output_file_prefix, 
                   list_coef)

################################################################################
#Interaction
ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots_guide/IFN_Interaction_int_genes/"))
IFN_genes <- extract_unique_genes(down_enrichr_filtered_KO, c("Interferon Alpha Response",
                                                      "Interferon Gamma Response"),
                                  db= "MSigDB_Hallmark_2020")
IFN_genes<-limmaRes%>%filter(toupper(ensg)%in% IFN_genes)%>%
  pull(ensg)%>%
  unique()
# Calculate the number of up and downregulated genes for each coefficient and cell type
# Define your cutoffs and filter the dataframe
adj_p_cutoff <- 0.05
logfc_cutoff <- 1

summary_df <- limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")
######################
count_threshold =10
coefficients <- summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)%>%
  pull(coef)%>%
  unique()
#
list_ko <- limmaRes%>%
  filter(coef %in% coefficients)%>%
  pull(coef)%>%
  unique()

for (KO in list_ko){
  IFN_genes_KO<-limmaRes %>%
    filter(coef == KO) %>%
    #filter(adj.P.Val < 0.05) %>%
    filter(logFC > 1) %>%
    filter(toupper(ensg) %in% IFN_genes)
  
  
  ggplot(IFN_genes_KO,aes(x=celltype,y=ensg,size=pmin(-log10(adj.P.Val),5),
                          color=logFC))+
    
    geom_point()+scale_color_gradient2(high="red", low="blue")+
    theme(axis.text = element_text(size = 15)) +
    #facet_wrap(vars(celltype), scales = "free_x")+
    theme_bw(12)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
    scale_size_continuous(range=c(0,5), limits = c(0,5))+
    ggtitle("IFN genes-Interaction(ex.vivo-in.vivo)")+
    ylab("IFN response genes")
  ko_dir<-dirout(paste0(base,"/IFN_interaction_down/",KO))
  ggsave(ko_dir(paste0(KO,"-interaction_upregulated_IFN_genes.pdf")))
  
}

################################################################################
#interaction_pathways and genes
list_coef <- coefficients  # Replace with your actual list of coefficients

#term <- "Cholesterol Homeostasis"  # Replace with your desired term
#termname<- "Cholesterol Homeostasis"
db <- "MSigDB_Hallmark_2020"  # Replace with your desired database
limma_data <- limmaRes
enrichr_data_up <- up_enrichr_filtered_KO
enrichr_data_down <- down_enrichr_filtered_KO
#output_file_prefix <- "Cholesterol"  # Replace with your desired output file prefix
###############
# List of terms and corresponding conditions
pathways <- list(
  "Cholesterol Homeostasis",
  "mTORC1 Signaling",
  c("Interferon Alpha Response",
    "Interferon Gamma Response"),
  "Inflammatory Response",
  "p53 Pathway",
  "Myc Targets V1",
  "IL−6/JAK/STAT3 Signaling",
  "Oxidative Phosphorylation",
  "TGF−beta Signaling",
  "Unfolded Protein Response",
  "E2F Targets",
  "Protein Secretion"
)
termnames<- c(
  "Cholesterol Homeostasis",
  "mTORC1 Signaling",
  "IFN_response",
  "Inflammatory Response",
  "p53 Pathway",
  "Myc Targets V1",
  "IL−6/JAK/STAT3 Signaling",
  "Oxidative Phosphorylation",
  "TGF−beta Signaling",
  "Unfolded Protein Response",
  "E2F Targets",
  "Protein Secretion"
)
termnames <- gsub(" ", "_", termnames)
map2(pathways,
     termnames, ~ plot_genes_by_term(limmaRes, 
                                     up_enrichr_filtered_KO, 
                                     down_enrichr_filtered_KO, 
                                     .x, 
                                     .y, 
                                     db, 
                                     .y,
                                     list_coef,
                                     dir_name = "Pathways_genes_in_interaction_effect/"))
###################


################################################################################
#excluding izzo-----------------------
################################################################################
# Iterate over each unique gene
Cholesterol <- intersect(Cholesterol, rownames(dataVoom_NTC_in_ex$E))
dat.list<-list()
for(gg in unique(Cholesterol)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta_in_ex %>%
    mutate(E = scale(dataVoom_NTC_in_ex$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
  
  # Group the data by tissue and celltype, and calculate the average E for each group
  #avg_gene_data <- gene_data %>%
  #group_by(tissue, celltype) %>%
  #summarise(avg_E = mean(E))
  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")

##########################################
create_gene_plots_NTC <- function(data, gene) {
  ggplot(data[data$gene == gene ,], aes(x = tissue, y = E)) + 
    geom_boxplot() +
    geom_jitter(aes(colour=tissue)) +
    facet_grid(rows = vars(celltype), scales = "free") +
    labs(title = gene)+
    xlab(paste0(gene))+
    ylab(paste0("scaled-Normalized gene exp"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots_guide/IFN_NTC_ex_in/"))
  
  ggsave(ko_dir(paste0(gene,".pdf")))
}
#For each KO where there were significantly up genes, make the plot
for (gg in IFN_genes){
  data<-dat.list%>%filter(gene==gg)%>%
    filter(celltype %in% c("Eo.Ba","HSC","MkP","Mono","GMP"))
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene)
    
  })
}
################################################################################
calculate_median <- function(x) {
  return(data.frame(y = median(x)))
}

mean_data <- dat.list %>%
    group_by(celltype, tissue) %>%
    summarise(mean_scaled_E = mean(E))
  
# Calculate median scaled_E for each celltype, tissue, and genotype
median_data <- dat.list %>%
    group_by(celltype, tissue) %>%
    summarise(median_scaled_E = median(E))
# Create the violin plot
ggplot(dat.list, aes(x = tissue, y = E, fill = tissue)) +
    #geom_violin(trim = FALSE) +
    geom_boxplot(aes(middle = median(E)))+
    #geom_jitter()+
    facet_wrap(~ celltype, scales = "free") +
    labs(x = "Tissue", y = "Scaled E", fill = "Genotype", title = paste0("IFN genes in upregulated_",comp,"_interaction")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/UP/",comp,"Summary_IFN/"))
ggsave(ko_dir("Boxplot_IFN_genes_in_Brd9_up.pdf"),h=min(20,nrow(data)*0.04+1))
  
}  
################################################################################
#Interaction--------------------------------------------------------------------
################################################################################
filter_limma_results <- function(limmaRes,
                                 direction,
                                 logFC_threshold,
                                 adj_p_val_threshold,
                                 gene_list) {
  filter_direction <-
    if (direction == "up")
      logFC > logFC_threshold
  else
    logFC < -logFC_threshold
  limmaRes %>%
    filter(filter_direction & adj.P.Val < adj_p_val_threshold) %>%
    filter(toupper(ensg) %in% gene_list)
}
#########################
#Interaction plot
#########################
create_interaction_plot <-
  function(data, direction, base_dir, gene_set_name) {
    ggplot(data, aes(
      x = celltype,
      y = ensg,
      size = pmin(-log10(adj.P.Val), 5),
      color = logFC
    )) +
      geom_point() +
      scale_color_gradient2(high = "red", low = "blue") +
      theme(axis.text = element_text(size = 15)) +
      facet_wrap(vars(celltype), scales = "free_x") +
      theme_bw(12) +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
      scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
      ggtitle(paste(
        gene_set_name,
        "Interaction (ex.vivo-in.vivo) -",
        direction,
        "regulated"
      )) +
      ylab("Response genes")
    
    ggsave(base_dir(paste0("Summary_of_genes_", direction, ".png")), width = 27, height = 12)
  }
##########################
#
##########################
create_gene_plots <- function(data, gene, KO, base_dir) {
  ggplot(data[data$gene == gene, ], aes(x = genotype, y = scaled_E)) +
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype ~ tissue, scales = "free") +
    labs(title = gene) +
    xlab(paste0(KO, " KO"))
  
  ko_dir <- file.path(base_dir, KO)
  if (!dir.exists(ko_dir))
    dir.create(ko_dir, recursive = TRUE)
  ggsave(file.path(ko_dir, paste0(gene, ".pdf")))
}

process_genes <-
  function(limmaRes,
           direction,
           logFC_threshold,
           adj_p_val_threshold,
           gene_list,
           base_dir) {
    filtered_data <-
      filter_limma_results(limmaRes,
                           direction,
                           logFC_threshold,
                           adj_p_val_threshold,
                           gene_list)
    
    # Summarized Plot
    create_interaction_plot(filtered_data, direction, base_dir)
    
    # Detailed Plots for Each KO
    list_ko <- filtered_data %>%
      pull(coef) %>%
      str_replace("ex.vivo:", "") %>%
      unique()
    
    dat.list <- list()
    
    for (KO in list_ko) {
      for (ct in unique(meta$celltype)) {
        dataVoom_ct <- get(paste0("dataVoom_", ct))
        
        list_of_genes <- filtered_data %>%
          filter(coef == paste0("ex.vivo:", KO)) %>%
          pull(ensg)
        
        if (any(rownames(dataVoom_ct$E) %in% list_of_genes)) {
          for (goi in list_of_genes) {
            if (goi %in% rownames(dataVoom_ct$E)) {
              gene_data <- meta[meta$celltype == ct, ] %>%
                mutate(E = dataVoom_ct$E[goi, ]) %>%
                rownames_to_column("sample1") %>%
                filter(genotype %in% c(KO, "NTC")) %>%
                mutate(scaled_E = scale(E)) %>%
                mutate(gene = goi) %>%
                mutate(celltype = ct) %>%
                mutate(comparison = KO)
              
              dat.list[[paste0(ct, "_", goi, KO)]] <- gene_data
            }
          }
        }
      }
    }
    
    limma_KO_genes <- bind_rows(dat.list, .id = "celltype_gene_KO")
    write_rds(limma_KO_genes, file.path(
      base_dir,
      paste0(
        "genes_in_each_ko_",
        direction,
        "regulated_for_all_celltypes.rds"
      )
    ))
    
    for (comp in unique(limma_KO_genes$comparison)) {
      data <- limma_KO_genes %>%
        filter(comparison == comp)
      lapply(unique(data$gene), function(gene) {
        create_gene_plots(data, gene, comp, base_dir)
      })
    }
  }


#############******************################################################
#############******************################################################
#*#Interaction upregulated------------------------------------------------------
 ################################################################################
 # logFC-----------
 ################################################################################
IFN_genes <- extract_unique_genes(up_enrichr_filtered_KO, c("Interferon Alpha Response",
                                                              "Interferon Gamma Response"),
                                  db= "MSigDB_Hallmark_2020")
list_ko <- limmaRes%>%
  filter(logFC>1 & adj.P.Val<0.05)%>%
  pull(coef)%>%
  str_replace("ex.vivo:","")%>%
  unique()

for (KO in list_ko){
  IFN_genes_KO<-limmaRes %>%
    filter(coef == paste0("ex.vivo:",KO)) %>%
    #filter(adj.P.Val < 0.05) %>%
    #filter(logFC > 1) %>%
    filter(toupper(ensg) %in% IFN_genes)
  
  
  ggplot(IFN_genes_KO,aes(x=celltype,y=ensg,size=-log10(adj.P.Val),
                            color=logFC))+
    
    geom_point()+scale_color_gradient2(high="red", low="blue")+
    theme(axis.text = element_text(size = 15)) +
    #facet_wrap(vars(celltype), scales = "free_x")+
    theme_bw(12)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
    scale_size_continuous(range=c(0,5), limits = c(0,5))+
    ggtitle("IFN genes-Interaction(ex.vivo-in.vivo)")+
    ylab("IFN response genes")
  ko_dir<-dirout(paste0(base,"/IFN_inter/",KO))
  ggsave(ko_dir(paste0(KO,"-interaction_upregulated_IFN_genes.pdf")))
  
}
###################
#all of them
list_coef <- limmaRes%>%
  filter(logFC>1 & adj.P.Val<0.05)%>%
  pull(coef)%>%
  unique()

IFN_genes_KO<-limmaRes %>%
  filter(coef %in% list_coef) %>%
  #Note:not filtering for log adj p value again
  filter(logFC > 1) %>%
  filter(toupper(ensg) %in% IFN_genes)
ggplot(IFN_genes_KO,aes(y=gsub("ex.vivo:","",coef),x=ensg,size=pmin(-log10(adj.P.Val),5),
                        color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5))+
  ggtitle("IFN genes-Interaction(ex.vivo-in.vivo)")+
  ylab("IFN response genes")
ggsave(IFN_inter("Summary_of_IFN_genes.png"),w=27,h=12)
###################
#Expression Upregulated
##############
#NOTE:CHECK IF CAPITALIZED
gene_list <-c("Idi1", "Cyp51","Stard4", "Scd2","Mthfd2","Sqle","Fads2","Dhcr24",
                           "Hmgcs1","Ldlr","Plscr1","Acat1",
                           "Acat2")

#Theses are the KO with significant upregulated interaction terms
list_ko <- limmaRes%>%
  filter(logFC>1 & adj.P.Val<0.05,coef %in% coefficients)%>%
  pull(coef)%>%
  str_replace("interaction","")%>%
  unique()

dat.list <- list()
for (KO in list_ko){
  list_of_genes<-limmaRes%>%
    #FILTERING ONLY SIGNIFICANT-UP IFN genes for that particular KO
    filter(coef==paste0("ex.vivo:",KO))%>%
             # ONLY UPREGULATED GENES
             filter(logFC>1 & adj.P.Val<0.05)%>%
    filter(toupper(ensg) %in% gene_list)%>%
    pull(ensg)
  for (ct in unique(meta$celltype)) {
    # Get the dataVoom object corresponding to the current cell type
    dataVoom_ct <- get(paste0("dataVoom_", ct))
    #CAPITALIZE
  
    # Check if goi exists in the row names of dataVoom_ct$E
    if (any(rownames(dataVoom_ct$E) %in% unique(list_of_genes))){
      for (goi in unique(list_of_genes)) {
        # Proceed only if goi exists in the row names of dataVoom_ct$E
        if (goi %in% rownames(dataVoom_ct$E)) {
          # Subset the metadata and E values for the current gene and cell type
          gene_data <- meta[meta$celltype == ct,] %>%
            mutate(E = dataVoom_ct$E[goi,]) %>%
            rownames_to_column("sample1") %>%
            filter(genotype %in% c(KO, "NTC")) %>%
            mutate(scaled_E = scale(E)) %>%
            mutate(gene = goi)%>%
            mutate(celltype=ct)%>%
            mutate(comparison=KO)
          
          # Store the gene data in the list
          dat.list[[paste0(ct, "_", goi,KO)]] <- gene_data
        }
      }
    }
  }
}
#Here we have sign.upregulated genes for each KO in all celltypes(not only the one with sig.up) of that KO.
limma_KO_up_IFN_genes<-bind_rows(dat.list,.id = "celltype_gene_KO")
limma_KO_up_IFN_genes%>%write_rds(basedir("IFN_genes_in_each_ko_upregulated_for all_celltypes"))
head(limma_KO_up_IFN_genes)
######################################################
# IFN genes-------------------------------------------
######################################################
basedir

IFN_UP<-dirout(paste0(base,"IFN_inter/UP/"))
limma_KO_up_IFN_genes$tissue<-factor(limma_KO_up_IFN_genes$tissue,levels=c("in.vivo","ex.vivo"))
unique(limma_KO_up_IFN_genes$genotype)
limma_KO_up_IFN_genes$genotype<-factor(limma_KO_up_IFN_genes$genotype,
                                       levels=c("NTC",setdiff(limma_KO_up_IFN_genes$genotype,"NTC")))

colnames(limma_KO_up_IFN_genes)
# Function to create plots for each gene
create_gene_plots <- function(data, gene,KO) {
  ggplot(data[data$gene == gene,], aes(x = genotype, y = scaled_E,color=)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype~tissue, scales = "free") +
    labs(title = gene)+
    xlab(paste0(KO,"KO"))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots_guide/IFN_inter/UP/",KO))
  ggsave(ko_dir(paste0(gene,".pdf")))
}
#For each KO where there were significantly up genes, make the plot
for (comp in unique(limma_KO_up_IFN_genes$comparison)){
  data<-limma_KO_up_IFN_genes%>%
    filter(comparison==comp)
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots(data, gene,comp)
    
  })
}
# Generate plots for each gene


#####################################################################
#Barplopt
# Group data by celltype, tissue, and genotype, calculate the mean scaled_E
# Define a custom function to calculate median and return a data frame
calculate_median <- function(x) {
  return(data.frame(y = median(x)))
}
for (comp in unique(limma_KO_up_IFN_genes$comparison)){
  data<-limma_KO_up_IFN_genes%>%
    filter(comparison==comp)
  #calculate mean  
  mean_data <- data %>%
    group_by(celltype, tissue, genotype) %>%
    summarise(mean_scaled_E = mean(scaled_E))
    
  # Calculate median scaled_E for each celltype, tissue, and genotype
  median_data <- data %>%
      group_by(celltype, tissue, genotype) %>%
      summarise(median_scaled_E = median(scaled_E))
  # Create the violin plot
  ggplot(data, aes(x = tissue, y = scaled_E, fill = genotype)) +
    #geom_violin(trim = FALSE) +
    geom_boxplot(aes(middle = mean(scaled_E)))+
    #geom_jitter()+
    facet_wrap(~ celltype, scales = "free") +
    labs(x = "Tissue", y = "Scaled E", fill = "Genotype", title = paste0("IFN genes in upregulated_",comp,"_interaction")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/UP/",comp,"Summary_IFN/"))
  ggsave(ko_dir("Boxplot_IFN_genes_in_Brd9_up.pdf"),h=min(20,nrow(data)*0.04+1))
    
}   

# Plot mean scaled_E for each celltype in each tissue
ggplot(mean_data, aes(x = celltype, y = mean_scaled_E, fill = genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ tissue, scales = "free") +
  labs(x = "Cell Type", y = "Mean Scaled E", fill = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(IFN_inter(paste0("mean_of_scaled_E_across_IFN_gene",".pdf")))
###############################################################
#####################################################################
#
#
#####################################################################
#library(ggplot2)

# Define a custom function to calculate median and return a data frame
calculate_median <- function(x) {
  return(data.frame(y = median(x)))
}

# Calculate median scaled_E for each celltype, tissue, and genotype
median_data <- data_Brd9 %>%
  group_by(celltype, tissue, genotype) %>%
  summarise(median_scaled_E = median(scaled_E))


# Create the violin plot
ggplot(data_Brd9, aes(x = tissue, y = scaled_E, fill = genotype)) +
  geom_violin(trim = FALSE) +
  #geom_hline(data = median_data, aes(yintercept = median_scaled_E),
  #color = "black", linetype = "solid", size = 1) +
  #geom_point(data = median_data, aes(y = median_scaled_E), 
  # color = "black", size = 3, shape = 95) +
  #geom_point(data = median_data, aes(x = celltype, y = median_scaled_E),
  #color = "black", size = 3, shape = 95) +
  facet_wrap(~ celltype, scales = "free") +
  labs(x = "Cell Type", y = "Scaled E", fill = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(IFN_inter(""))

#####
# Boxplot
# Calculate median scaled_E for each genotype, tissue, and celltype
median_data <- data_Brd9 %>%
  group_by(genotype, tissue, celltype) %>%
  summarise(median_scaled_E = median(scaled_E))

# Create the violin plot
ggplot(data_Brd9, aes(x = tissue, y = scaled_E, fill = genotype)) +
  #geom_violin(trim = FALSE) +
  geom_boxplot(aes(middle = mean(scaled_E)))+
  #geom_jitter()+
  facet_wrap(~ celltype, scales = "free") +
  labs(x = "Tissue", y = "Scaled E", fill = "Genotype", title = "IFN genes in upregulated Brd9KO interaction.pdf") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(IFN_inter("Boxplot_IFN_genes_in_Brd9_up.pdf"))
################################################################################
#
Cholesterol_inter<-dirout(paste0(base,"Cholesterol_inter"))
#Cholesterol homeostasis
Cholesterol_genes <- extract_unique_genes(, c("Cholesterol Homeostasis"))
################################################################################
#Interaction upregulated------------------------------------------------------
################################################################################
# logFC-----------

###################
#Expression Upregulated
##############
#NOTE:CHECK IF CAPITALIZED
gene_list<-Cholesterol_genes
#Theses are the KO with significant upregulated interaction terms
list_ko <- limmaRes%>%
  filter(logFC< -1 & adj.P.Val<0.05,coef %in% coefficients)%>%
  pull(coef)%>%
  str_replace("ex.vivo:","")%>%
  unique()

dat.list <- list()
for (KO in list_ko){
  list_of_genes<-limmaRes%>%
    #FILTERING ONLY SIGNIFICANT-UP IFN genes for that particular KO
    filter(coef==paste0("ex.vivo:",KO))%>%
    # ONLY UPREGULATED GENES
    filter(logFC< -1 & adj.P.Val<0.05)%>%
    filter(toupper(ensg) %in% gene_list)%>%
    pull(ensg)
  for (ct in unique(meta$celltype)) {
    # Get the dataVoom object corresponding to the current cell type
    dataVoom_ct <- get(paste0("dataVoom_", ct))
    #CAPITALIZE
    
    # Check if goi exists in the row names of dataVoom_ct$E
    if (any(rownames(dataVoom_ct$E) %in% unique(list_of_genes))){
      for (goi in unique(list_of_genes)) {
        # Proceed only if goi exists in the row names of dataVoom_ct$E
        if (goi %in% rownames(dataVoom_ct$E)) {
          # Subset the metadata and E values for the current gene and cell type
          gene_data <- meta[meta$celltype == ct,] %>%
            mutate(E = dataVoom_ct$E[goi,]) %>%
            rownames_to_column("sample1") %>%
            filter(genotype %in% c(KO, "NTC")) %>%
            mutate(scaled_E = scale(E)) %>%
            mutate(gene = goi)%>%
            mutate(celltype=ct)%>%
            mutate(comparison=KO)
          
          # Store the gene data in the list
          dat.list[[paste0(ct, "_", goi,KO)]] <- gene_data
        }
      }
    }
  }
}
#Here we have sign.upregulated genes for each KO in all celltypes(not only the one with sig.up) of that KO.
limma_KO_down_Cholesterol_genes<-bind_rows(dat.list,.id = "celltype_gene_KO")
limma_KO_down_Cholesterol_genes%>%write_rds(basedir("Cholesterol_genes_in_each_ko_downregulated_for all_celltypes"))

######################################################
# Cholesterol_genes-------------------------------------------
######################################################
Cholesterol_down<-dirout(paste0(base,"Cholesterol_Inter/DOWN/"))

# Function to create plots for each gene
create_gene_plots <- function(data, gene,KO) {
  ggplot(data[data$gene == gene,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype~tissue, scales = "free") +
    labs(title = gene)+
    xlab(paste0(KO,"KO"))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/Cholesterol_Inter/DOWN/",KO))
  ggsave(ko_dir(paste0(gene,".pdf")))
}
#For each KO where there were significantly up genes, make the plot
for (comp in unique(limma_KO_down_Cholesterol_genes$comparison)){
  data<-limma_KO_down_Cholesterol_genes%>%
    filter(comparison==comp)
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots(data, gene,comp)
    
  })
}
# Generate plots for each gene


#####################################################################
#Barplopt
# Group data by celltype, tissue, and genotype, calculate the mean scaled_E
# Define a custom function to calculate median and return a data frame
calculate_median <- function(x) {
  return(data.frame(y = median(x)))
}
for (comp in unique(limma_KO_down_Cholesterol_genes$comparison)){
  data<-limma_KO_down_Cholesterol_genes%>%
    filter(comparison==comp)
  #calculate mean  
  mean_data <- data %>%
    group_by(celltype, tissue, genotype) %>%
    summarise(mean_scaled_E = mean(scaled_E))
  
  # Calculate median scaled_E for each celltype, tissue, and genotype
  median_data <- data %>%
    group_by(celltype, tissue, genotype) %>%
    summarise(median_scaled_E = median(scaled_E))
  # Create the violin plot
  ggplot(data, aes(x = tissue, y = scaled_E, fill = genotype)) +
    #geom_violin(trim = FALSE) +
    geom_boxplot(aes(middle = mean(scaled_E)))+
    #geom_jitter()+
    facet_wrap(~ celltype, scales = "free") +
    labs(x = "Tissue", y = "Scaled E", fill = "Genotype", title = paste0("Cholesterol genes in downregulated_",comp,"_interaction")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/UP/",comp,"Summary_IFN/"))
  ggsave(ko_dir("Boxplot_Cholesterol_in",comp,"down.pdf"),h=min(20,nrow(data)*0.04+1))
  
}   
##################################################################
##################################################################
basedir<-dirout_load("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype/FGSEA")
gsea.res<-read_rds(basedir("FGSEA_Interaction_across_KO.RDS"))
head(gsea.res)
gsea.res%>%filter(pathway %in% grep("Cholesterol Homeostasis",gsea.res$pathway,value = T))%>%
  pull(leadingEdge)%>% unlist()%>%unique()

gsea.res%>%filter(pathway %in% grep("Cholesterol Homeostasis",gsea.res$pathway,value = T))%>%
  pull(leadingEdge)%>% unlist()%>%unique()
##################################################################
limmaRes%>%filter(coef=="ex.vivo:Brd9")%>%
  filter(ensg %in% c("Actb", "Gapdh", "Rpl13a", "Hprt", "Tbp"))
list_of_genes<-c("Actb", "Gapdh", "Rpl13a", "Hprt", "Tbp")
dat.list <- list()
for (ct in unique(meta$celltype)) {
  # Get the dataVoom object corresponding to the current cell type
  dataVoom_ct <- get(paste0("dataVoom_", ct))
  #CAPITALIZE
  
  # Check if goi exists in the row names of dataVoom_ct$E
  if (any(rownames(dataVoom_ct$E) %in% unique(list_of_genes))){
    for (goi in unique(list_of_genes)) {
      # Proceed only if goi exists in the row names of dataVoom_ct$E
      if (goi %in% rownames(dataVoom_ct$E)) {
        # Subset the metadata and E values for the current gene and cell type
        gene_data <- meta[meta$celltype == ct,] %>%
          mutate(E = dataVoom_ct$E[goi,]) %>%
          rownames_to_column("sample1") %>%
          filter(genotype %in% c("Brd9", "NTC")) %>%
          mutate(scaled_E = scale(E)) %>%
          mutate(gene = goi)%>%
          mutate(celltype=ct)%>%
          mutate(comparison="Brd9")
        
        # Store the gene data in the list
        dat.list[[paste0(ct, "_", goi,"Brd9")]] <- gene_data
      }
    }
  }}
control_genes<-bind_rows(dat.list,.id = "celltype_gene_KO")
create_gene_plots <- function(data, gene,KO) {
  ggplot(data[data$gene == gene,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype~tissue, scales = "free") +
    labs(title = gene)+
    xlab(paste0(KO,"KO"))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/Brd9_Control_genes",KO))
  ggsave(ko_dir(paste0(gene,".pdf")))
}
head(control_genes)
#For each KO where there were significantly up genes, make the plot
for (comp in unique(control_genes$comparison)){
  data<-control_genes%>%
    filter(comparison==comp)
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots(data, gene,comp)
    
  })
}
# Generate plots for each gene

IFN_genes_logFC<-limmaRes%>%
  filter(toupper(ensg)%in% IFN_genes)%>%
  filter(abs(logFC)>1)

ggplot(IFN_genes_logFC,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                        color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5))+
  ggtitle("IFN genes-Interaction(ex.vivo-in.vivo)")+
  ylab("IFN response genes")
ggsave(basedir("IFN_genes_in_ko_interaction.pdf"),w=25,h=10)

#############################################################################################

Cholesterol_Hom_logFC<-limmaRes%>%filter(toupper(ensg)%in% Cholesterol_Hom_genes )%>%
  filter(abs(logFC)>1)

ggplot(Cholesterol_Hom_logFC,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                           color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5))+
  ggtitle("Cholesterol_Hom-Interaction(ex.vivo-in.vivo)")+
  ylab("Cholesterol_Hom")
ggsave(basedir("Cholesterol_Hom_in_ko_interaction.pdf"),w=25,h=10)

############################################################################################
############################################
#Cholesterol_NTC
############################################
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, 
                                          c("Cholesterol Homeostasis"))
limmaRes_NTC$genes
Cholesterol_genes<-limmaRes_NTC%>%filter(toupper(genes)%in% Cholesterol_genes)%>%
  pull(genes)%>%
  unique()
Cholesterol_genes <- intersect(Cholesterol_genes, rownames(dataVoom_NTC$E))
dat.list<-list()
for(gg in unique(Cholesterol_genes)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta_in_ex %>%
    mutate(E = scale(dataVoom_NTC_in_ex$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
  
  # Group the data by tissue and celltype, and calculate the average E for each group
  #avg_gene_data <- gene_data %>%
  #group_by(tissue, celltype) %>%
  #summarise(avg_E = mean(E))
  
  # Store the average gene data in the list
  dat.list[[gg]] <- gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")

##########################################
create_gene_plots_NTC <- function(data, gene) {
  ggplot(data[data$gene == gene ,], aes(x = tissue, y = E)) + 
    geom_boxplot() +
    geom_jitter(aes(colour=tissue)) +
    facet_grid(cols = vars(celltype), scales = "free") +
    labs(title = gene)+
    xlab(paste0(gene))+
    ylab(paste0("scaled-Normalized gene exp"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots_guide/Cholesterol_NTC_ex_in/"))
  
  ggsave(ko_dir(paste0(gene,".pdf")),h=8)
}
#For each KO where there were significantly up genes, make the plot
for (gg in Cholesterol_genes){
  data<-dat.list%>%filter(gene==gg)%>%
    filter(celltype %in% c("Eo.Ba","HSC","MkP","Mono","GMP"))
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene)
    
  })
}
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, 
                                          c("Cholesterol Homeostasis"))
data<-limmaRes_NTC
#############################
create_gene_plots_NTC <- function(data, gene) {
  # Filter the data for the specific gene
  data_gene <- data[data$gene == gene,]
  
  # Ensure the 'tissue' column is a factor
  data_gene$tissue <- as.factor(data_gene$tissue)
  
  # Create the plot
  plot <- ggplot(data_gene, aes(x = tissue, y = E)) + 
    geom_boxplot() +
    geom_jitter(aes(colour = tissue), width = 0.2, height = 0) +
    facet_grid(cols = vars(celltype), scales = "free") +
    labs(title = gene) +
    xlab(paste0(gene)) +
    ylab("scaled-Normalized gene exp") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Create directory if it doesn't exist
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots_guide/Cholesterol_NTC_ex_in/"))
  # Save the plot
  ggsave(ko_dir(paste0(gene,".pdf")),h=8)
}

# For each KO where there were significantly up genes, make the plot
for (gg in Cholesterol_genes) {
  data <- dat.list %>%
    filter(gene == gg) %>%
    filter(celltype %in% c("Eo.Ba", "HSC", "MkP", "Mono", "GMP"))
  
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots_NTC(data, gene)
  })
}
head(enr.res.all_up)

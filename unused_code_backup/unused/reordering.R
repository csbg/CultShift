
#load data
inDir<-dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
inDir1<-dirout_load("Ag_ScRNA_09_Normalized_counts_ex_in_izzo_NTC")
inDir2<-dirout_load("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment/ENRICHR")
inDir3<-dirout_load("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC")
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
enr.res.all_up<-read_rds(inDir2("enr.res.all_NTC_up.rds"))
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
#meta<-read_rds(inDir("meta.rds"))
################################################################################
#function2: Plot overlapping genes from enriched pathways for relevant KOs
################################################################################
limma_data 
enrichr_data_up
enrichr_data_down
term
termname
db
output_file_prefix
list_coef
dir_name

# General plotting function for specific terms and datasets
plot_genes_by_term <- function(limma_data, 
                               enrichr_data_up,
                               enrichr_data_down,
                               term,
                               termname,
                               db,
                               output_file_prefix,
                               list_coef,
                               dir_name) {
  # Combine upregulated and downregulated enrichment results
  combined_enrichr_data <- bind_rows(enrichr_data_up, enrichr_data_down)
  
  # Extract unique genes for the given term
  unique_genes <- extract_unique_genes(combined_enrichr_data, term, db)
  
  # Filter limma data for the extracted genes and the given coefficients
  genes_KO <- limma_data %>%
    filter(celltype !="MEP")%>%
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
    ylab(paste(term, "genes"))
  plot
  # Create directory if it doesn't exist
  
  dir<-dirout(paste0(base,dir_name,termname))
  
  # Save the plot
  ggsave(dir(paste0(output_file_prefix,"_",termname,"_summary.pdf")),
         plot = plot, width = length(unique_genes) * 0.2 + 5,
         height = length(unique_genes) * 0.18 + 1,
         limitsize = FALSE)
}


###############################################################################
#function2_usage1
###############################################################################
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

#interaction_pathways and genes
list_coef <- coefficients  # Replace with your actual list of coefficients

#term <- "Cholesterol Homeostasis"  # Replace with your desired term
#termname<- "Cholesterol Homeostasis"
db <- "MSigDB_Hallmark_2020"  # Replace with your desired database
limma_data <- limmaRes
enrichr_data_up <- up_enrichr_filtered_KO
enrichr_data_down <- down_enrichr_filtered_KO
#output_file_prefix <- "Cholesterol"  # Replace with your desired output file prefix

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
                                     list_coef ="",
                                     dir_name = "Pathways_genes_in_interaction_effect/"))
###################
#function2_version2:for NTC
###################
# General plotting function for specific terms and datasets
limmaRes_NTC$genes<-gsub("\\.\\.\\.\\d+", "", rownames(limmaRes_NTC))
limma_data<-limmaRes_NTC
plot_genes_by_term <- function(limma_data, 
                               enrichr_data_up,
                               enrichr_data_down,
                               term,
                               termname,
                               db,
                               output_file_prefix,
                               dir_name) {
  # Combine upregulated and downregulated enrichment results
  combined_enrichr_data <- bind_rows(enrichr_data_up, enrichr_data_down)
  
  # Extract unique genes for the given term
  unique_genes <- extract_unique_genes(combined_enrichr_data, term, db)
  
  term="Cholesterol Homeostasis"
  db= "MSigDB_Hallmark_2020"
  genes_KO <- limma_data %>%
      filter(toupper(genes) %in% unique_genes)

  
  # Check if genes_KO has more than one gene
  if (nrow(genes_KO) == 0) {
    cat("No data to plot for term:", term, "\n")
    return(NULL)
  }
  # Check if genes_KO has more than one gene
  if (nrow(genes_KO) == 0) {
    cat("No data to plot for term:", term, "\n")
    return(NULL)
  }
  head(limmaRes)
  # Create the plot
  plot <- ggplot(genes_KO, aes(y = celltype, x = genes, size = -log10(adj.P.Val), color = logFC)) +
    geom_point(alpha = 0.6) +
    scale_color_gradient2(high = "red", low = "blue") +
    theme(axis.text = element_text(size = 15)) +
    facet_wrap(vars(celltype), scales = "free_x") +
    theme_bw(12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 8),
          axis.text.y = element_text(size = 8)) +
    scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
    ggtitle(paste(term, "- Interaction(ex.vivo-in.vivo)")) +
    ylab(paste(term, "genes"))
  plot
  # Create directory if it doesn't exist
  
  dir<-dirout(paste0(base,dir_name,termname))
  
  # Save the plot
  ggsave(dir(paste0(output_file_prefix,termname,"_summary.pdf")),
         plot = plot, width = length(unique_genes) * 0.2 + 5,
         height = length(unique_genes) * 0.18 + 1,
         limitsize = FALSE)
}



enr.res.all_up<-read_rds(inDir2("enr.res.all_NTC_up.rds"))
enr.res.all_down<-read_rds(inDir2("enr.res.all_NTC_down.rds"))

enrichr_data_up <- enr.res.all_up
enrichr_data_down <- enr.res.all_up
limmaRes<-limmaRes_NTC
limmaRes
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
     termnames, ~ plot_genes_by_term(
       limma_data =limmaRes_NTC, 
       enrichr_data_up,
       enrichr_data_down,
       term=.x,
       termname=.y,
       db,
       output_file_prefix = "NTC_",
       dir_name = "Pathways_genes_in_NTC/"))

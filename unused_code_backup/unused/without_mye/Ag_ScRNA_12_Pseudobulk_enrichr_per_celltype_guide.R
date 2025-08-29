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
#####################################################################
inDir  <-  dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
base  <-  "Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/"
basedir  <-  dirout("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/")
outdir <- dirout("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/")
source("src/Ag_Optimized_theme.R")
########################################################################
limmaRes  <-  read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
limmaRes  <-  limmaRes%>% filter(celltype != "MEP")
################################
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
#############################################################################
#functions-------------------------------------------------------------------
#############################################################################
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

# Filter rows where Count > 10
summary_df %>%
  filter(Count > 10) %>%
  # Plot with log10 scale on y-axis
  ggplot(aes(x = gsub("interaction", "", coef), y = Count, fill = Regulation)) +
  geom_col(position = position_dodge(preserve = "single")) +    # Separate bars for Upregulated and Downregulated
  scale_fill_manual(values = c("Downregulated" =  "#4C889C", "Upregulated" = "#D0154E")) +
  facet_grid(rows = vars(celltype)) +
  labs(title="No: of genes with interaction efect of culture condition in KOs",
         x="KOs")+
  scale_y_log10() +  # Apply log10 scale to the y-axis
  optimized_theme()
ggsave(inDir("Differential_genes_count.pdf"))

###############################################################################
#function for enrichment analysis for
#limmaRes:input limma results table
#databases: databases in Enrichr
#coefficients: coefficients in limmaRes, that has to be used for enrichment,
#if you have a 
#list of KOs for example, only do the enrichment analysis for those
#logFC_threshold
#directory for saving results
#file_prefix


#function
perform_enrichment_analysis  <-  function(limmaRes,
                                        databases,
                                        coefficients,
                                        logFC_threshold = 1,
                                        directory,
                                        file_prefix="") {
  enrichment_results_up  <-  list()
  enrichment_results_down  <-  list() 
  
  # Iterate over each cell type and coefficient
  for (ct in unique(limmaRes$celltype)) {
    cat("Processing:", ct, "\n")
    for (coefx in unique(coefficients)) {
      cat("Processing:", ct, coefx, "\n")
      
      # Filter for upregulated and downregulated ensg
      genes_up  <-  limmaRes %>%
        filter(celltype == ct & coef == coefx & group == "up") %>%
        pull(ensg)
      
      genes_down  <-  limmaRes %>%
        filter(celltype == ct & coef == coefx & group == "down") %>%
        pull(ensg)
      
      if (length(genes_up) > 0) {
        enr_res_up  <-  enrichr(genes_up, databases = databases)
        if (!is.null(enr_res_up)) {
          # Filter out list elements with zero rows
          enr_res_up_filtered  <-  Filter(function(x) nrow(x) > 0, enr_res_up)
          if (length(enr_res_up_filtered) > 0) {  # Check if there are any remaining list elements
            # Bind rows of the filtered list
            enr_res_up  <-  bind_rows(enr_res_up_filtered, .id = "db") %>%
              mutate(celltype = ct, coef = coefx)  # Add ct and coef as columns
            enrichment_results_up[[paste(ct, coefx, "up", sep = "_")]]  <-  enr_res_up
          }
        }
      }
      
      if (length(genes_down) > 0) {
        enr_res_down  <-  enrichr(genes_down, databases = databases)
        if (!is.null(enr_res_down)) {
          # Filter out list elements with zero rows
          enr_res_down_filtered  <-  Filter(function(x) nrow(x) > 0, enr_res_down)
          if (length(enr_res_down_filtered) > 0) {  # Check if there are any remaining list elements
            # Bind rows of the filtered list
            enr_res_down  <-  bind_rows(enr_res_down_filtered, .id = "db") %>%
              mutate(celltype = ct, coef = coefx)  # Add ct and coef as columns
            enrichment_results_down[[paste(ct, coefx, "down", sep = "_")]]  <-  enr_res_down
          }
        }
      }
    }
  }
  
  # Combine results 
  
  down_enrichr  <-  bind_rows(enrichment_results_down, .id = "ct_coef")
  up_enrichr  <-  bind_rows(enrichment_results_up, .id = "ct_coef")
  
  #create directory for enrichr and save rds files
  directory <-  dirout(paste0("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/","Enrichr_results_updated/"))
  #enrich <- dirout(paste0(base,"Enrichr"))
  write_rds(down_enrichr, directory(paste0(file_prefix,"down_logFC_1_enrichr.rds")))
  write_rds(up_enrichr, directory(paste0(file_prefix,"up_logFC_1_enrichr.rds")))
}
################################################################################
#make directory for output


# For all KOs


################################################################################
# Only for the knockouts with atleast 10 genes----------------------------------
################################################################################
#for selected KOs
#filtered based on KOs with atleast 10 differentially regulated genes between tissue types
count_threshold = 10
coefficients  <-  summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)%>%
  pull(coef)%>%
  unique()

perform_enrichment_analysis(limmaRes,
                            databases,
                            coefficients,
                            directory = directory,
                            file_prefix = "count_above10_")
################################################################################
#plots
################################################################################
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
    filter(coef %in% coefficients) %>%
    filter(toupper(ensg) %in% unique_genes)
  #%>%
    #filter(adj.P.Val <0.05)
  
  # Check if genes_KO has more than one gene
  if (nrow(genes_KO) == 0) {
    cat("No data to plot for term:", term, "\n")
    return(NULL)
  }
  
 
  
  # Create directory if it doesn't exist
  #termname<-"IFN_response"
  #dir_name<-"testdir"
 dir <- dirout(paste0(base,dir_name))
 # output_file_prefix<-"IFN_response"
  
  # Save the plot
 # Get the number of unique cell types to determine the number of rows
 # Set the number of genes
 num_genes <- length(unique_genes)
 
 # Set a maximum allowed width for the plot (for example, 20 units)
 max_width <- 10
 
 # Initial calculation for the number of rows (default is 2)
 #facet_rows <- ifelse(num_genes > 50, 4, 2)
 
 # Calculate the initial width based on the number of genes
 #initial_width <- num_genes * 0.45 + 5
 
 # If the initial width exceeds the maximum allowed width, increase the number of rows
 # if (initial_width > max_width) {
 #   additional_rows <- ceiling(initial_width / max_width)  # Calculate how many additional rows are needed
 #   facet_rows <- facet_rows * additional_rows  # Increase the number of rows accordingly
 #   initial_width <- max_width  # Set the width to the maximum allowed width
 # }
 # Create the plot
 
 # Adjust the plot by setting the number of rows in facet_wrap
 plot <- ggplot(genes_KO,
                aes(x = gsub("interaction", "", coef), y = ensg,
                    size = -log10(adj.P.Val), color = logFC)) +
   geom_point(alpha = 0.6) +
   scale_color_gradient2(high = "red", low = "blue") +
   theme(axis.text = element_text(size = 15)) +
   facet_wrap(vars(celltype), #nrow = facet_rows,
              scales = "free_x",nrow =1) +  # Dynamic number of rows
   theme_bw(12) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
         axis.text.y = element_text(size = 8)) +
   scale_size_continuous(range = c(2, 6)) +
   ggtitle(paste(term, "- Interaction(ex.vivo-in.vivo)")) +
   ylab(paste("Knockouts")) +
   xlab("genes")+optimized_theme()+
   theme(axis.text.x = element_text(angle = 90))
 
 
  ggsave(dir(paste0(output_file_prefix,termname,"_summary.pdf")),
         plot = plot, 
         width = length(unique_genes) * 0.2 + 1,
         height = length(unique_genes) * 0.2 + 5,
         limitsize = FALSE)
}


down_enrichr_filtered_KO <- read_rds(directory(paste0("count_above10_down_logFC_1_enrichr.rds")))
up_enrichr_filtered_KO <- read_rds(directory(paste0("count_above10_up_logFC_1_enrichr.rds")))


db <- "MSigDB_Hallmark_2020"  # Replace with your desired database
limma_data <- limmaRes
head(limmaRes)
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
                                     coefficients,
                                     dir_name ="Pathways_genes_in_interaction_effect/"))

coefficients <-c("interactionBrd9","interactionPhf10",
                 "interactionBcl11a","interactionRcor1","interactionWdr82","interactionSmc3")
map2(pathways,
     termnames, ~ plot_genes_by_term(limmaRes, 
                                     up_enrichr_filtered_KO, 
                                     down_enrichr_filtered_KO, 
                                     .x, 
                                     .y, 
                                     db, 
                                     .y,
                                     coefficients,
                                     dir_name ="Pathways_genes_in_interaction_effect_selected/"))


################################################################################

dfs<-list(up_enrichr_filtered_KO,down_enrichr_filtered_KO)
names(dfs)<-c("up","down")
#modify for KOs.
pmap(list(names(dfs), dfs), function(name, df) {
  walk(databases, ~{
    db <- .x
    df_filtered <- df[df$db == db & df$celltype != "MEP", ]
    
    # Check if subset is not empty
    if (nrow(df_filtered) > 0) {
      # Remove the prefix from column names
      colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
      
      # Perform filtering and mutation
      plotting_enr <- df_filtered %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
      
      ggplot(plotting_enr, aes(x = gsub("interaction","",coef), 
                               y = Term, color = log2(Odds.Ratio),
                               size = pmin(10,neg.log10.Adjusted.P.value))) +
        
        geom_point() +
        
        scale_size_continuous(
          range = c(2, 6)  # Set the sizes in the legend
        )+
        scale_color_gradientn(
          colors = c("pink", "red"),
          #breaks = c(0, 5),  # Set custom breaks for the color scale
          #labels = c("0", "2", "5"),  # Set custom labels for the breaks
          #limits = c(0, 10)  # Set the limits for the color scale
        ) +
        facet_wrap(vars(celltype))+
        
        ggtitle(paste0(db))+
        optimized_theme()
        # theme(axis.text = element_text(size = 15),
        #       axis.title = element_text(size=15),
        #       legend.text = element_text(size=12),
        #       legend.title = element_text(size=12),
        #       title = element_text(size=18))theme(axis.text = element_text(size = 15),
              # axis.title = element_text(size=15),
              # legend.text = element_text(size=12),
              # legend.title = element_text(size=12),
              # title = element_text(size=18))
      
      
      ggsave(outdir(paste0(name, "_per.celltype_", db, ".pdf")), 
             w = 10, h = length(unique(plotting_enr$Term)) * 0.2 + 3, limitsize = FALSE)
      
      # Do something with 'result'
    } else {
      warning("Subset is empty for iteration ")
    }
  })
})
#######################################

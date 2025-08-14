source("src/00_init.R")
require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)
require(enrichR)
# renv::snapshot(lockfile = "renv_NF.lock")

source("~/code/resources/RFunctions/Basics.R")

out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/"
base <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Enrichr"

base_fgsea<-"/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/FGSEA"


################################################################################
dirout_jak <- function(out, ext="", init=TRUE){
  out.dir <- paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/", "/JAKSTAT/", out, "/")
  if(init){
    dir.create(out.dir,showWarnings=FALSE); 
    message("Setting output directory: ", out.dir)
  }
  function(...){
    paste0(out.dir, paste0(...), ext)
  }
}

#function
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
res<-read.delim(paste0(out,"/DEG.tsv"))
res$probe<-res$rn
res$rn<-NULL
res<-as.data.frame(res)
gmap<-as.data.frame(read_rds(file = file.path(out, "DEG_GMP.RDS")))
res<-merge(res,gmap[,c("probe","gene")],by="probe")

# Visualize results ---------------------------------------------------------
limmaRes <- res %>% filter(grepl("treatmentex_vivo",res$genotype))
limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))


################################################################################
#function
perform_enrichment_analysis <- function(limmaRes,
                                        databases,
                                        coefficients,
                                        logFC_threshold = 1,
                                        directory,
                                        file_prefix="") {
  enrichment_results_up <- list()
  enrichment_results_down <- list() 
  
  # Iterate over each cell type and coefficient
  for (ct in unique(limmaRes$cell_type)) {
    cat("Processing:", ct, "\n")
    for (coefx in unique(coefficients)) {
      cat("Processing:", ct, coefx, "\n")
      
      # Filter for upregulated and downregulated genes
      genes_up <- limmaRes %>%
        filter(cell_type == ct & genotype == coefx & group == "up" & logFC > logFC_threshold) %>%
        pull(gene)
      
      genes_down <- limmaRes %>%
        filter(cell_type == ct & genotype == coefx & group == "down" & logFC < -logFC_threshold) %>%
        pull(gene)
      
      if (length(genes_up) > 0) {
        enr_res_up <- enrichr(genes_up, databases = databases)
        if (!is.null(enr_res_up)) {
          # Filter out list elements with zero rows
          enr_res_up_filtered <- Filter(function(x) nrow(x) > 0, enr_res_up)
          if (length(enr_res_up_filtered) > 0) {  # Check if there are any remaining list elements
            # Bind rows of the filtered list
            enr_res_up <- bind_rows(enr_res_up_filtered, .id = "db") %>%
              mutate(cell_type = ct, coef = coefx)  # Add ct and coef as columns
            enrichment_results_up[[paste(ct, coefx, "up", sep = "_")]] <- enr_res_up
          }
        }
      }
      
      if (length(genes_down) > 0) {
        enr_res_down <- enrichr(genes_down, databases = databases)
        if (!is.null(enr_res_down)) {
          # Filter out list elements with zero rows
          enr_res_down_filtered <- Filter(function(x) nrow(x) > 0, enr_res_down)
          if (length(enr_res_down_filtered) > 0) {  # Check if there are any remaining list elements
            # Bind rows of the filtered list
            enr_res_down <- bind_rows(enr_res_down_filtered, .id = "db") %>%
              mutate(cell_type = ct, coef = coefx)  # Add ct and coef as columns
            enrichment_results_down[[paste(ct, coefx, "down", sep = "_")]] <- enr_res_down
          }
        }
      }
    }
  }
  
  # Combine results 
  
  down_enrichr <- bind_rows(enrichment_results_down, .id = "ct_coef")
  up_enrichr <- bind_rows(enrichment_results_up, .id = "ct_coef")
  
  #create directory for enrichr and save rds files
  
  #
  
  write_rds(down_enrichr, paste0(directory,"/",file_prefix,"down_logFC_",logFC_threshold,"_enrichr.rds"))
  write_rds(up_enrichr, paste0(directory,"/",file_prefix,"up_logFC_",logFC_threshold,"_enrichr.rds"))
}
#function usage
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019")
directory<-"/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Enrichr"

coefficients<-unique(limmaRes$genotype)
perform_enrichment_analysis(limmaRes, 
                            databases, coefficients,
                            directory = directory,logFC_threshold = 3)
################################################################################

################################################################################
IFN_limmRes <- limmaRes %>%
  filter(gene %in% IFN_genes & genotype == "treatmentex_vivo")

# Calculate percentage of downregulated genes (logFC < -1)
downregulated <- IFN_limmRes %>%
  filter(logFC < -1 & adj.P.Val < 0.05)
nrow(downregulated)
length(IFN_genes)
unique(IFN_genes)

percent_downregulated <- nrow(downregulated) / length(IFN_genes) * 100

# Create a bar plot for percentage of downregulated genes
percent_plot <- ggplot() +
  geom_bar(aes(x = "Downregulated", y = percent_downregulated), stat = "identity", fill = "blue") +
  geom_text(aes(x = "Downregulated", y = percent_downregulated + 2, label = paste0(round(percent_downregulated, 1), "%")), vjust = -0.5) +
  geom_bar(aes(x = "Not Downregulated", y = 100 - percent_downregulated), stat = "identity", fill = "gray") +
  geom_text(aes(x = "Not Downregulated", y = 100 - percent_downregulated + 2, label = paste0(round(100 - percent_downregulated, 1), "%")), vjust = -0.5) +
  coord_flip() +
  labs(x = "", y = "Percentage") +
  ggtitle("Percentage of Downregulated IFN Genes") +
  theme_minimal()

percent_plot
# Subset limmaRes for specific genes
limma_subset <- limmaRes %>%
  filter(gene %in% IFN_genes & genotype == "treatmentex_vivo")

# Filter limma results for IFN genes and treatmentex_vivo genotype
IFN_limmRes <- limmaRes %>%
  filter(gene %in% IFN_genes & genotype == "treatmentex_vivo")

# Count total IFN genes
total_IFN_genes <- length(IFN_genes)

# Count upregulated IFN genes (assuming logFC > 1 as upregulated)
upregulated_IFN <- IFN_limmRes %>%
  filter(logFC < -1) # Adjust threshold as needed

total_upregulated_IFN <- nrow(upregulated_IFN)

# Create a summary dataframe for plotting
summary_df <- data.frame(
  Category = c("Total IFN Genes", "Upregulated IFN Genes"),
  Count = c(total_IFN_genes, total_upregulated_IFN)
)
##########################
# Create the bar plot
 ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.3) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +
  labs(
    title = "Summary of IFN Genes",
    x = "Category",
    y = "Count"
  ) +
  theme(legend.position = "none")
################################################################################
# Only for the knockouts with atleast 10 genes----------------------------------
################################################################################
#for selected KOs
# plotting function for number of up/down genes
plot_genes <- function(data, 
                       count_threshold = 0,
                       width = 0.2,
                       title_suffix = "",
                       file_suffix = "") {
  
  # Filter out rows with count equal to 0 and based on count threshold
  filtered_data <- data %>% 
    filter(Count != 0) %>% 
    filter(Count >= count_threshold)
  
  # Plot
  p <- ggplot(filtered_data, 
              aes(x = genotype,
                  y = ifelse(Regulation == "Downregulated",
                             -log10(Count), log10(Count)), 
                  fill = Regulation)) +
    geom_col(width = width) +
    scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
    labs(x = "Interaction-KO",
         y = "Log10(Number of Genes)",
         title = paste("Upregulated and downregulated genes", title_suffix)) +
    theme_bw(12) +
    theme(axis.text = element_text(size = 15)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    facet_grid(cols =vars(cell_type)) +
    coord_flip() +
    theme(strip.text = element_text(size = 15))
  
  # Save plot
  ggsave(enrich(paste0("Number_of_genes", file_suffix, ".png")), plot = p,
         height = length(coefficients)*0.25+3)
}


#for all KOs
# Calculate the number of up and downregulated genes for each coefficient and cell type
summary_df <- limmaRes %>%
  group_by(cell_type, genotype) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")

head(summary_df)
plot_genes(summary_df, count_threshold = 0, title_suffix = "", file_suffix = "")
#for Kos with above 10 genes differentially regulated between tissue type
plot_genes(summary_df, count_threshold = 10, title_suffix = " (Count >= 10)", 
           file_suffix = "_above10")

# Define your cutoffs and filter the dataframe
adj_p_cutoff <- 0.05
logfc_cutoff <- 1


#filtered based on KOs with atleast 10 differentially regulated genes between tissue types
count_threshold = 10
coefficients <- summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)%>%
  pull(genotype)%>%
  unique()
perform_enrichment_analysis(limmaRes,
                            databases,
                            coefficients,
                            directory = directory,
                            file_prefix = "count_above10_")

################################################################################
#plot and save function for enriched results

plot_and_save<- function(enr.res, db, output_file,file_prefix = "",base) {
  # Subset the enrichment results for the current database
  plotting_enr <- enr.res[enr.res$db == db,] %>%
    filter(Odds.Ratio > 5 & Adjusted.P.value < 0.05) %>%
    mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
  # Check the number of overlapping genes
  plotting_enr <- plotting_enr %>%
    mutate(overlap_count = sapply(strsplit(Genes, ";"), length)) %>%
    filter(overlap_count > 1)
  
  # Create the plot
  ggplot(plotting_enr, aes(x = gsub("ex.vivo:","",coef), 
                           y = Term, color = log2(Odds.Ratio),
                           size = neg.log10.Adjusted.P.value)) +
    labs(x = "KOs",
         y = "Terms-Odds.Ratio > 5 & Adjusted.P.value < 0.05", 
         title = paste0(output_file))+
    geom_point() +
    
    # coord_flip()+
    facet_wrap(vars(cell_type))+
    theme_bw(12)+
    scale_size(range = c(1,3))+
    #scale_size_continuous(1,4)+
    ggtitle(paste0(output_file,"_",db))+
    theme(axis.title = element_text(size=7)) +
    theme(strip.text = element_text(size = 7))+
    theme(legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))+
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size = 7),
          axis.text.y = element_text(size = 7))
  
  
  # Save the plot
  enrichr<-dirout_jak(paste0("enrichr_plot/",db))
  
  ggsave(filename = enrichr(paste0(file_prefix,output_file,"_", db, ".pdf")),
         width = length(unique(coefficients))*0.35+4,
         #h=18,
         height = (length(unique(plotting_enr$Term)) * 0.155 + 3),
         limitsize = FALSE)
}



#set directory/filepath

enrich<-dirout_jak("Enrichr")
# Apply the function to each database using purrr::map()
down_enrichr<-read_rds(enrich("/down_logFC_3_enrichr.rds"))
up_enrichr<-read_rds(enrich("/up_logFC_3_enrichr.rds"))

map(databases, ~ plot_and_save(enr.res = up_enrichr, db = .x, output_file = "Upregulated_enrichr"))
map(databases, ~ plot_and_save(enr.res = down_enrichr, db = .x, output_file = "Downregulated_enrichr"))

################################################################################
# Only for the knockouts with atleast 10 genes----------------------------------
################################################################################
#filteredKO
down_enrichr_filtered_KO<-read_rds(enrich("count_above10_down_logFC_1_enrichr.rds"))
up_enrichr_filtered_KO<-read_rds(enrich("count_above10_up_logFC_1_enrichr.rds"))
map(databases, ~ plot_and_save(enr.res = up_enrichr_filtered_KO, db = .x,
                               output_file = "Upregulated_enrichr",
                               file_prefix = "count_above10"))
map(databases, ~ plot_and_save(enr.res = down_enrichr_filtered_KO, db = .x,
                               output_file = "Downregulated_enrichr",
                               file_prefix = "count_above10"))

################################################################################
# Define your cutoffs and filter the dataframe
adj_p_cutoff <- 0.05
logfc_cutoff <- 1

# Calculate the number of up and downregulated genes for each coefficient and cell type
summary_df <- limmaRes %>%
  group_by(cell_type, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")
###############################################

# plotting function for number of up/down genes
plot_genes <- function(data, 
                       count_threshold = 0,
                       width = 0.2,
                       title_suffix = "",
                       file_suffix = "") {
  
  # Filter out rows with count equal to 0 and based on count threshold
  filtered_data <- data %>% 
    filter(Count != 0) %>% 
    filter(Count >= count_threshold)
  
  # Plot
  p <- ggplot(filtered_data, 
              aes(x = gsub("ex.vivo:", "", coef),
                  y = ifelse(Regulation == "Downregulated",
                             -log10(Count), log10(Count)), 
                  fill = Regulation)) +
    geom_col(width = width) +
    scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
    labs(x = "Interaction-KO",
         y = "Log10(Number of Genes)",
         title = paste("Upregulated and downregulated genes", title_suffix)) +
    theme_bw(12) +
    theme(axis.text = element_text(size = 15)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    facet_grid(cols =vars(cell_type)) +
    coord_flip() +
    theme(strip.text = element_text(size = 15))
  
  # Save plot
  ggsave(basedir(paste0("Number_of_genes", file_suffix, ".png")), plot = p,
         height = length(coefficients)*0.25+3)
}

#for all KOs
plot_genes(summary_df, count_threshold = 0, title_suffix = "", file_suffix = "")
#for Kos with above 10 genes differentially regulated between tissue type
plot_genes(summary_df, count_threshold = 10, title_suffix = "counts_above_10", 
           file_suffix = "_above10")
##########################
#do it by myself
# Filter out rows with count equal to 0 and based on count threshold
filtered_data <- summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)
head(summary_df)
# Plot
p <- ggplot(filtered_data, 
            aes(x = genotype,
                y = ifelse(Regulation == "Downregulated",
                           -log10(Count), log10(Count)), 
                fill = Regulation)) +
  geom_col(width = width) +
  scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
  labs(x = "Interaction-KO",
       y = "Log10(Number of Genes)",
       title = paste("Upregulated and downregulated genes", "counts_above_10")) +
  theme_bw(12) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid(cols =vars(cell_type)) +
  coord_flip() +
  theme(strip.text = element_text(size = 15))

# Save plot
basedir
ggsave(paste0(base,"Number_of_genes", "counts_above_10", ".pdf"), plot = p,
       )
################################################################################
#fgsea -------------------------------------------------------------------------
################################################################################
out <- dirout("EXT_02_EnrichR_Genesets/")
# # # Download gene sets ------------------------------------------------------
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)

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
save(enr.terms, file=out("Genesets_Mouse.RData"))
################################################################################
################################################################################
#IFN-genes
################################################################################
up_enrichr<-read_rds(paste0(directory,"/","up_logFC_1_enrichr.rds"))
down_enrichr<-read_rds(paste0(directory,"/","down_logFC_1_enrichr.rds"))

#from enrichment analysis-

IFN_down<-extract_unique_genes(enrichr_data = down_enrichr, term=c("Interferon Alpha Response",
                                                                   "Interferon Gamma Response"),
                               db="MSigDB_Hallmark_2020")
##################################################
#all_downregulated_IFN genes across condtions
limmaRes%>%
  filter(toupper(gene) %in% IFN_down)%>%
  filter(adj.P.Val < 0.05)%>%
  filter(logFC < -1)%>%
  ggplot(aes(x = genotype, y = gene, 
             size = pmax(4,-log10(adj.P.Val)),
             color = pmin(10,logFC))) +
  geom_point() +
  scale_color_gradient2(high = "red", low = "blue") +
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(cell_type), scales = "free_x") +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
  ggtitle(paste("ex.vivo-in.vivo_JAK-STAT")) 

# Create directory if it doesn't exist

dir<-dirout_jak(paste0("/IFN_genes"))

# Save the plot
ggsave(dir(paste0("IFN_genes_downregulated","_summary.pdf")),
       width = length(IFN_down) * 0.05 ,
       height = length(IFN_down) * 0.1 + 2)
#limitsize = FALSE)
##################################################
# Define a function to summarize and plot IFN gene data
plot_Gene_set_summary <- function(gene_set, limmaRes, output_dir,gene_set_name=names(gene_set)) {
  
  # Filter and get unique genes from limmaRes based on the provided gene_set
  Gene_set <- intersect(gene_set, unique(limmaRes$gene))
  
  # Summarize data
  plot_data <- limmaRes %>%
    filter(gene %in% Gene_set, adj.P.Val < 0.05) %>%
    mutate(Regulation = case_when(
      logFC > 0 ~ "Upregulated",
      logFC < -1 ~ "Downregulated",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Regulation)) %>%
    group_by(genotype, cell_type, Regulation) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = genotype, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_grid(. ~ cell_type, scales = "free_x", space = "free") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
    theme_minimal() +
    labs(
      title = paste("Count of Upregulated and Downregulated Genes by Genotype and Cell Type"),
      subtitle = gene_set_name,
      x = "Genotype",
      y = "Count"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")
  
  # Save the plot with a unique filename based on gene_set_name
  output_filename <- file.path(output_dir, paste0("Gene_set_summary_", gene_set_name, ".pdf"))
  ggsave(output_filename, plot = p, width = 10, height = 6)
  
  return(output_filename)
}

# Directory to save plots
output_dir <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Enrichr/"

# Iterate through each gene set name in pathway_names and create plots
plot_files <- lapply(pathway_genes, function(pathway_genes) {
  plot_Gene_set_summary(pathway_genes, limmaRes, output_dir)
})

# Print the list of saved plot files
print(plot_files)
# Example pathway names
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
  "Protein Secretion"
)


# Initialize a list to store genes for each pathway



# Define function to summarize and plot gene set data
plot_Gene_set_summary <- function(gene_set_name, genes, limmaRes, output_dir) {
  
  # Extract genes from gene set
  Gene_set <- genes[[gene_set_name]]
  
  # Filter and get unique genes from limmaRes based on the provided gene_set
  Gene_set <- intersect(Gene_set, unique(limmaRes$gene))
  
  # Check if Gene_set is empty
  if (length(Gene_set) == 0) {
    message(paste("No genes found for", gene_set_name))
    return(NULL)
  }
  
  # Summarize data for upregulated and downregulated genes
  plot_data <- limmaRes %>%
    filter(gene %in% Gene_set, adj.P.Val < 0.05) %>%
    mutate(Regulation = case_when(
      logFC > 0 ~ "Upregulated",
      logFC < -1 ~ "Downregulated",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Regulation)) %>%
    group_by(genotype, cell_type, Regulation) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Check if plot_data is empty
  if (nrow(plot_data) == 0) {
    message(paste("No data found for", gene_set_name))
    return(NULL)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = genotype, y = Count, fill = Regulation)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_grid(. ~ cell_type, scales = "free_x", space = "free") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
    theme_minimal() +
    labs(
      title = paste("Count of Upregulated and Downregulated Genes by Genotype and Cell Type"),
      subtitle = gene_set_name,
      x = "Genotype",
      y = "Count"
    ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom")
  
  # Save the plot with a unique filename based on gene_set_name
  output_filename <- file.path(output_dir, paste0("Gene_set_summary_", gene_set_name, ".pdf"))
  ggsave(output_filename, plot = p, width = 10, height = 6)
  
  return(output_filename)
}

# Iterate through each gene set in pathway_genes and create plots
plot_files <- lapply(names(pathway_genes), function(gene_set_name) {
  plot_Gene_set_summary(gene_set_name, pathway_genes, limmaRes, output_dir)
})
length(ISG_core)
# Print the list of saved plot files
print(plot_files)

plot_data <- limmaRes %>%
  filter(gene %in% Gene_set, adj.P.Val < 0.05) %>%
  mutate(Regulation = case_when(
    logFC > 0 ~ "Upregulated",
    logFC < -1 ~ "Downregulated",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Regulation)) %>%
  group_by(genotype, cell_type, Regulation) %>%
  summarise(Count = n(), .groups = 'drop')




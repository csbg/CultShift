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
########################################################################
limmaRes  <-  read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
limmaRes  <-  limmaRes%>% filter(celltype != "MEP")
# Define your cutoffs and filter the dataframe
adj_p_cutoff  <-  0.05
logfc_cutoff  <-  1

# Calculate the number of up and downregulated genes for each coefficient and cell type
summary_df  <-  limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")
count_threshold = 10
coefficients  <-  summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)%>%
  pull(coef)%>%
  unique()
###############################################
# plotting function for number of up/down genes
# Improved plotting function
plot_genes  <-  function(data, 
                         count_threshold = 0,
                         width = 0.2,
                         title_suffix = "",
                         file_suffix = "") {
  
  # Filter out rows with count equal to 0 and based on count threshold
  filtered_data  <-  data %>% 
    filter(Count != 0) %>% 
    filter(Count >= count_threshold) %>%
    mutate(LogCount = ifelse(Regulation == "Downregulated", -log10(abs(Count)), log10(Count)),
           CountLabel = abs(Count))
  
  # Determine y-axis limits for padding
  max_log_count  <-  max(abs(filtered_data$LogCount), na.rm = TRUE)
  y_limit  <-  c(-max_log_count * 1.5, max_log_count * 1.5)
  
  # Plot
  p  <-  ggplot(filtered_data, 
                aes(x = gsub("ex.vivo:", "", coef),
                    y = LogCount,
                    fill = Regulation))  +
    geom_col(width = width) +
    geom_text(aes(label = CountLabel), 
              position = position_dodge(width = width), 
              vjust = ifelse(filtered_data$LogCount > 0, -0.5, 1.5)) +
    scale_fill_manual(values = c("Upregulated" = "#8A264A", "Downregulated" = "#5782A7")) +
    labs(x = "Interaction-KO",
         y = "Log10(Number of Genes)",
         title = paste("Upregulated and Downregulated Genes", title_suffix)) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text = element_text(size = 15)) +
    facet_grid(rows = vars(celltype)) +
    #coord_flip() +
    scale_y_continuous(labels = abs, limits = y_limit) # Ensure y-axis labels are positive and add padding
  
  # Save plot
  ggsave(filename = basedir(paste0("Number_of_genes", file_suffix, ".png")), plot = p,
         height = length(unique(filtered_data$coef)) * 0.25 + 4, width = 10)
}


plot_genes(summary_df, count_threshold = 10, title_suffix = " (Count >= 10)", 
           file_suffix = "_above10")
#########################################
#plot_enrichment analysis----------------
#########################################
#plot and save function for enriched results
plot_and_save <-  function(enr.res, db, output_file,file_prefix = "") {
  # Subset the enrichment results for the current database
  plotting_enr  <-  enr.res[enr.res$db == db,] %>%
    filter(Odds.Ratio > 5 & Adjusted.P.value < 0.05) %>%
    mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
  # Check the number of overlapping genes
  plotting_enr  <-  plotting_enr %>%
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
    facet_wrap(vars(celltype))+
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
  dir <- dirout(paste0(base,"Enrichr/",db))
  ggsave(filename = dir(paste0(file_prefix,output_file,"_", db, ".pdf")),
         width = length(unique(coefficients))*0.35+4,
         #h=18,
         height = (length(unique(plotting_enr$Term)) * 0.155 + 3),
         limitsize = FALSE)
}

#set directory/filepath
directory <- dirout(paste0(base,"Enrichr"))
enrich <- directory
##############################
#All KOs
##############################
# # Apply the function to each database using purrr::map()
# down_enrichr <- read_rds(enrich("down_logFC_1_enrichr.rds"))
# up_enrichr <- read_rds(enrich("up_logFC_1_enrichr.rds"))
# 
# map(databases, ~ plot_and_save(enr.res = up_enrichr, db = .x, output_file = "Upregulated_enrichr"))
# map(databases, ~ plot_and_save(enr.res = down_enrichr, db = .x, output_file = "Downregulated_enrichr"))

################################################################################
# Only for the knockouts with atleast 10 genes----------------------------------
################################################################################

#filteredKO
down_enrichr_filtered_KO <- read_rds(enrich("count_above10_down_logFC_1_enrichr.rds"))
up_enrichr_filtered_KO <- read_rds(enrich("count_above10_up_logFC_1_enrichr.rds"))
map(databases, ~ plot_and_save(enr.res = up_enrichr_filtered_KO, db = .x,
                               output_file = "Upregulated_enrichr",
                               file_prefix = "count_above10"))
map(databases, ~ plot_and_save(enr.res = down_enrichr_filtered_KO, db = .x,
                               output_file = "Downregulated_enrichr",
                               file_prefix = "count_above10"))
#################################################################################
#Plot_genes_from specific pathways from enrichr
#################################################################################
# Define the function
#filtered:refers whether data being plotted is filtered for adj.pvalue and log2 odds ratio

plot_specific_pathway <- function(data, term, db, direction, filtered = TRUE) {
  # Filter the data for the specified term and conditions
  plotting_data <- data %>%
    filter(db == db) %>%
    filter(Term == term)
  
  # Check the number of overlapping genes
  plotting_data <- plotting_data %>%
    mutate(overlap_count = sapply(strsplit(Genes, ";"), length)) %>%
    filter(overlap_count > 1)
  
  # Check if there is data to plot
  if (nrow(plotting_data) == 0) {
    stop("No data to plot. Please check your filters and input dataset.")
  }
  
  if (filtered) {
    plotting_data <- plotting_data %>%
      filter(Odds.Ratio > 5 & Adjusted.P.value < 0.05) %>%
      mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
  } else {
    plotting_data <- plotting_data %>%
      mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
  }
  
  # Generate the plot
  p <- ggplot(plotting_data, aes(x = gsub("ex.vivo:", "", coef), 
                                 y = celltype, 
                                 color = log2(Odds.Ratio), 
                                 size = neg.log10.Adjusted.P.value)) +
    labs(x = "KOs", title = paste0(direction, "_", term)) +
    geom_point() +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(size = 7),
          strip.text = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 7))
  
  filt <- ifelse(filtered, "filtered", "unfiltered")
  
  dir <- dirout(paste0(base, "Enrichr/Specific_pathways/", term))
  ggsave(dir(paste0(term, "_", direction, "_", db, "_", filt, ".png")),
         plot = p, width = 8, height = 8)
}

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

# Conditions with corresponding enrichment data
conditions <- list(
  "Up" = up_enrichr_filtered_KO,
  "Down" = down_enrichr_filtered_KO
)

count_threshold = 10
coefficients <- summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold) %>%
  pull(coef) %>%
  unique()

# Database name
db <- "MSigDB_Hallmark_2020"

# Function to apply plot_specific_pathway
apply_plot_specific_pathway <- function(data, term, db, direction, filtered) {
  filtered_data <- data %>%
    filter(db == db) %>%
    filter(Term == term)
  
  if (filtered) {
    filtered_data <- filtered_data %>%
      filter(Odds.Ratio > 5 & Adjusted.P.value < 0.05)
  }
  
  if (nrow(filtered_data) == 0) {
    cat("No data to plot for term:", term, "in direction:", direction,
        "with filtering set to", filtered, "\n")
  } else {
    plot_specific_pathway(data = data, term = term, db = db,
                          direction = direction, filtered = filtered)
  }
}

# Example of applying the function to each condition
for (direction in names(conditions)) {
  data <- conditions[[direction]]
  for (term in unlist(pathways)) {
    apply_plot_specific_pathway(data, term, db, direction, filtered = TRUE)
  }
}
################################################################################
#fgsea plots -------------------------------------------------------------------
################################################################################
gsea.res <- read_rds(out("FGSEA_Interaction_across_KO.RDS"))
# cleanup / export results
gsea.res[is.nan(NES), NES := 0]
gsea.res.export  <-  gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge  <-  sapply(gsea.res.export$leadingEdge,
                                        function(vec) paste(vec[1:10], collapse = ","))

for(dbx in unique(gsea.res$db)){
  dat <- dirout(paste0(base,"FGSEA/",dbx))
  write.tsv(gsea.res.export[db == dbx], dat("GSEA_significant_",dbx,".tsv"))
}

dbx <- "MSigDB_Hallmark_2020"
# Prepare for plotting
process_db <- function(dbx, gsea_res, base,size_n = 0,coefs) {
  pDT <- gsea_res[gsea_res$db == dbx & size > size_n & padj < 0.05 &
                    coef %in% coefs,]
  head(gsea.res)
  ## Splitting the task to handle both ends of the NES spectrum - positive and negative
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5), by=c("celltype", "coef")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("celltype", "coef")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  pDT <- pDT[pathway %in% pw.display]
  
  # Hierarchical ordering
  pDT <- hierarch.ordering(pDT, "celltype", "pathway", "NES", TRUE)
  
  # Plotting
  ggplot(pDT, aes(x=gsub("ex.vivo:", "", coef), y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
    scale_color_gradient2(low="blue", mid="white", high="red") +
    geom_point(data=pDT[padj < 0.05]) +
    scale_size_continuous(range=c(0, 5), limits = c(0, 5)) +
    theme_bw(12) +
    xRot() +
    facet_wrap(vars(celltype)) +
    labs(x = "interaction-KO") +
    theme(axis.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
  
  # Save the plot
  dat <- dirout(paste0(base, "FGSEA/", dbx))
  ggsave(dat("GSEA_plot_", dbx,"_",size_n, ".pdf"), width = 30,
         height = length(unique(pDT$pathway)) * 0.3 + 3, limitsize = FALSE)
}

# Apply the function to each unique `db` using purrr::walk
unique_dbs <- unique(gsea.res$db)
walk(unique_dbs, process_db, gsea_res = gsea.res, base = base,size_n =5,coefs=coefficients)
################################################################################
#specific_pathways_from fgsea
################################################################################
dbx <- "MSigDB_Hallmark_2020"
#MSigDB
pDT  <-  gsea.res[db == dbx]
## Splitting the task to handle both ends of the NES spectrum-positive and negative
pw.display.pos  <-  unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5), by=c("celltype", "coef")]$pathway)
pw.display.neg  <-  unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("celltype", "coef")]$pathway)

# Combine and remove duplicates across both positive and negative selections
pw.display  <-  unique(c(pw.display.pos, pw.display.neg))
pDT  <-  pDT[pathway %in% pw.display]
#pDT  <-  hierarch.ordering(pDT, "pathway", "celltype", "NES", TRUE)
pDT  <-  hierarch.ordering(pDT, "celltype", "pathway", "NES", TRUE)

Interested_pathways <- c("IL-6/JAK/STAT3 Signaling","Inflammatory Response",
                         "Interferon Alpha Response","Interferon Gamma Response",
                         "Cholesterol Homeostasis","mTORC1 Signaling",
                         "TNF−alpha Signaling via NF−kB")
pDT <- pDT[pathway %in% Interested_pathways & coef %in% coefficients ]
ggplot(pDT, aes(x=gsub("ex.vivo:","",coef), 
                y=celltype, color=NES, 
                size=pmin(5, -log10(padj)))) +
  
  scale_color_gradient2(low="blue", mid="white", high="red") +
  geom_point(data=pDT[padj < 0.05]) +
  scale_size_continuous(range=c(1.3,5), limits = c(1.3,5)) +
  theme_bw(12) +
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))+
  facet_wrap(vars(pathway),)+
  labs(x = "interaction-KO")+#space="free", scales="free") +
  theme(strip.text = element_text(size = 15))
#theme(strip.text=element_text(angle=0),strip.placement = "outside")
dat <- dirout(paste0(base,"FGSEA/",dbx))
ggsave(dat("GSEA_plot_",dbx,"selected_pathways.png"),w=15,h=18)
################################################################################
#function
extract_unique_genes  <-  function(dataframe, term,dbx) {
  # Filter the dataframe for rows where the Term is in the specified list of terms
  pDT  <-  subset(gsea.res, pathway ==term & db==dbx)
  
  # Extract the genes column from the filtered dataframe
  genes  <-  filtered_table$Genes
  
  # Split the genes column by ";" and flatten the resulting list
  all_genes  <-  unlist(strsplit(genes, ";"))
  
  # Select only unique genes
  unique_genes  <-  unique(all_genes)
  
  return(unique_genes)
}

Interested_pathways <- c("IL-6/JAK/STAT3 Signaling","Inflammatory Response",
                         "Interferon Alpha Response","Interferon Gamma Response",
                         "Cholesterol Homeostasis","mTORC1 Signaling",
                         "TNF−alpha Signaling via NF−kB")

#IFN alpha response

################################################################################
list_ko  <-  limmaRes%>%
  filter(coef %in% coefficients)%>%
  pull(coef)%>%
  unique()

for (KO in list_ko){
  IFN_genes_KO <- limmaRes %>%
    filter(coef == KO) %>%
    #filter(adj.P.Val < 0.05) %>%
    filter(logFC > 1) %>%
    filter(toupper(ensg) %in% IFN_genes)
  
  
  ggplot(IFN_genes_KO,aes(x=celltype,
                          y=ensg,
                          size=pmin(-log10(adj.P.Val),5),
                          color=logFC))+
    
    geom_point()+
    scale_color_gradient2(high="red", low="blue")+
    theme(axis.text = element_text(size = 15)) +
    #facet_wrap(vars(celltype), scales = "free_x")+
    theme_bw(12)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
    scale_size_continuous(range=c(0,5), limits = c(0,5))+
    ggtitle("IFN genes-Interaction(KO.ex.vivo-KO.in.vivo)")+
    ylab("IFN response genes")
  ko_dir <- dirout(paste0(base,"/IFN_interaction_down/",KO))
  ggsave(ko_dir(paste0(KO,"-interaction_downregulated_IFN_genes.pdf")))
  
}
enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`

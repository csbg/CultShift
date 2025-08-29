
###############
source("src/00_init.R")
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
require(fgsea)
library(latex2exp)
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification_Mye.R")
#####################################################################
inDir  <-  dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide_Mye/")

InDir5 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide_Mye/")
base  <-  "Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_Mye/"
basedir  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_Mye")
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
################################################################################
adj_p_cutoff <- 0.05
logfc_cutoff <- 1
# data
meta <- fread(InDir5("meta_cleaned.tsv")) # Read data
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   

meta <- meta[, -1, drop = FALSE] 
colnames(meta) <- gsub("rowname","sample1", colnames(meta))
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop")%>%
  mutate(coef = genotype)


selected_KOs <- meta %>%
  group_by(genotype, tissue, celltype) %>%                  # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample1), .groups = 'drop') %>% # Count distinct samples for each group
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>% # Spread tissue to separate columns (in.vivo and ex.vivo)
  group_by(genotype) %>%                                    # Regroup by genotype
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) %>%              # Keep genotypes that have at least one celltype with 3+ samples in both tissues
  pull(genotype) %>% unique()


summary_df <- limmaRes %>%
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
correlation_deg <- read_rds(basedir("correlation_deg.rds"))
KO_list <- correlation_deg %>% filter(correlation < 0.5, num_degs >= 10) %>%
  pull(genotype)%>%
  unique()
koi <- Reduce(intersect, list(selected_KOs,  KO_list))#, coefficients)) #only valid_k
################################################################################
#fgsea--------------------------------------------------------------------------
################################################################################

enr.terms <- enrichrGetGenesets(databases)
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

################################################################################
# fgsea -------------------------------------------------------------------
# Initialize the result table
################################################################################

# Initialize the result table
gsea.res <- data.table() 

run_gsea <- function(limmaRes, enr.terms, celltypes = NULL, coefs = NULL) {
  
  # Initialize the result table
  gsea_res <- data.table() 
  
  # Determine cell types to process
  if (is.null(celltypes)) {
    celltypes <- unique(limmaRes$celltype)
  }
  
  # Determine coefficients to process
  if (is.null(coefs)) {
    coefs <- unique(limmaRes$coef)
  }
  
  # Loop through each cell type
  for (ct in celltypes) {
    
    # Loop through each coefficient
    for (de_grp in coefs) {
      
      # Loop through each database in the enrichment terms
      for (dbx in names(enr.terms)) {
        
        # Subset the limma results based on the current cell type and coefficient
        subset_limmaRes <- limmaRes[limmaRes$celltype == ct & limmaRes$coef == de_grp, ]
        
        # Extract statistics (logFC) and assign gene names as names
        stats <- with(subset_limmaRes, setNames(logFC, nm = ensg))
        
        # Skip this iteration if there are missing values in stats
        if (any(is.na(stats))) {
          next
        }
        
        # Perform fgsea analysis
        fgsea_output <- fgsea(
          pathways = enr.terms[[dbx]],
          stats = stats
          #minSize = 15,   # Example additional arguments, adjust as necessary
          #maxSize = 500,  # Example additional arguments, adjust as necessary
          #nperm = 1000    # Example additional arguments, adjust as necessary
        )
        
        # Check if fgsea output is not empty and append the results to gsea_res
        if (length(fgsea_output) > 0) {
          gsea_res <- rbind(gsea_res, data.table(fgsea_output,
                                                 coef = de_grp,
                                                 celltype = ct,
                                                 db = dbx))
        }
      }
    }
  }
  
  # Return the combined GSEA results
  return(gsea_res)
}
gsea.res <- run_gsea(limmaRes, enr.terms, celltypes = unique(limmaRes$celltype),
                     coefs =unique(limmaRes$coef))
gsea.res %>% write_rds(basedir("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))

#############################################
#for all databases
# Assuming 'databaselist' is a vector of databases you want to plot
  # replace with actual database names


#}

#selected pathways


######################################################
# gsea.res <- read_rds(InDir7("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
gsea.res$coef <- gsub("interaction","",gsea.res$coef )
# Summarize to find KOs with at least one valid cell type
valid_ko_summary <- ko_flags %>%
  group_by(genotype) %>%
  summarize(has_valid_celltype = any(valid_ko), .groups = 'drop')




create_gsea_plot <- function(db) {
  
  # Filter the data for the given database
  pDT <- gsea.res %>%
    filter(coef %in% koi) %>%
    filter(db == !!db) %>%  # Correctly reference the current database
    left_join(ko_flags, by = c("coef", "celltype")) %>%
    left_join(summary_df, by = c("coef", "celltype")) %>%
    filter(Count > 5) %>%
    filter(valid_ko)%>%
    dplyr::select(-c("Count", "Regulation"))
  # Merge with KO flags
  
  # Step 3 continued: Keep only valid KOs for the specific cell type
  pDT <- pDT %>% filter(valid_ko == TRUE, padj < 0.05)
  
  # Select the pathways for plotting (both positive and negative)
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n = 7), by = c("coef", "celltype", "pathway")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n = 7), by = c("coef", "celltype", "pathway")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(union(pw.display.pos, pw.display.neg))
  # Remove duplicate rows
  pDT <- pDT %>% distinct()
  if (nrow(pDT) > 0) {
    # Step 1: Aggregate NES values across all celltypes (mean of NES per pathway)
    pDT_agg <- pDT %>%
      group_by(pathway) %>%
      summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
      arrange(desc(abs(average_NES)))%>%  # Ordering pathways by the average NES, highest first
      slice_head( n=30)
    # Step 2: Create a factor for pathway that reflects the aggregated NES order
    pDT <- pDT %>%
      filter(pathway %in% pDT_agg$pathway)
    pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)}
  
  # Filter pDT to include only selected pathways
  # **Convert list-columns to character format**
  dat <- pDT %>%
    mutate(across(where(is.list), ~sapply(., toString)))  # Converts lists to comma-separated strings
  
  # Save the filtered data table
  write.table(dat, file = basedir(paste0("fgsea_", db, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # Create the plot for the current database
  fig <- ggplot(pDT, aes(x = coef, y = pathway, color = NES, size = pmin(5, -log10(padj)))) +
    geom_point() + 
    scale_color_gradient2(low = "#4C889C", mid = "white", high = "#D0154E", name = TeX("log_{2}(FC)")) +
    geom_point(data = pDT[padj < 0.05], shape = 1) +
    scale_size_continuous(range = c(0, 2), limits = c(0, 5), name = TeX("$-\\log_{10}(p_{adj})$")) +
    theme_bw() +
    xRot() +
    labs(x = "KOs") +
    facet_grid(cols = vars(celltype), scales = "free", space = "free") +  # Create facets for each cell type
    optimized_theme_fig()+
    theme(strip.text.x = element_text(angle = 90))
  #strip.text.y = element_text(angle = 0)) +
  
  # Save the plot for the current database
  ggsave(basedir("fig2.2_supplementary_fgsea", db, "_3rep.pdf"), fig,
         w =30, h =30, units = "cm")
  
  
  
}

# Example usage: 
# Create the plot for "MSigDB_Hallmark_2020"
create_gsea_plot("MSigDB_Hallmark_2020")
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
# You can also loop through the databases to generate plots for all of them
lapply(databases, create_gsea_plot)
#############################################################################
############################################
enr.terms$MSigDB_Hallmark_2020$
  pathways <- c(
    "Cholesterol Homeostasis",
    "mTORC1 Signaling",
    "Interferon Alpha Response",
    "Interferon Gamma Response",
    #"Inflammatory Response",
    #"p53 Pathway",
    "Myc Targets V1",
    "KRAS Signaling Dn",
    #"IL−6/JAK/STAT3 Signaling",
    "Oxidative Phosphorylation",
    "TGF−beta Signaling",
    #"Unfolded Protein Response",
    "E2F Targets")
#"Protein Secretion"

#pathways <- names(enr.terms$MSigDB_Hallmark_2020)
dbx <- "MSigDB_Hallmark_2020"
pDT  <-  gsea.res[db == dbx & pathway %in% pathways]
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5),by=c("coef", "celltype","pathway")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("coef", "celltype","pathway")]$pathway)
# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(c(pw.display.pos, pw.display.neg))

pDT <- pDT[pathway %in% pw.display] %>%
  filter(coef %in% paste0("interaction",koi))

# Ensure that 'celltype' and 'NES' are not logical and have correct data types
if (!is.character(pDT$celltype) && !is.factor(pDT$celltype)) {
  stop("Column 'celltype' must be character or factor.")
}
if (!is.numeric(pDT$NES)) {
  stop("Column 'NES' must be numeric.")
}



if (nrow(pDT) > 0) {
  # Aggregate NES values across all cell types
  pDT_agg <- pDT %>%
    group_by(pathway) %>%
    summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
    arrange(desc(average_NES)) 
  
  #pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)
  # Ensure pDT follows the same order as pDT_agg
  pDT <- pDT %>%
    mutate(pathway = factor(pathway, levels = pDT_agg$pathway)) %>%
    arrange(factor(pathway, levels = pDT_agg$pathway))
  
  ggplot(pDT, aes(y=gsub("interaction","",coef), x=celltype, color=NES, size=pmin(5, -log10(padj)))) +
    geom_point() + scale_color_gradient2(low="blue", mid="white", high="red") +
    geom_point(data=pDT[padj < 0.05], shape=1) +
    scale_size_continuous(range=c(0,3),
                          limits = c(0,5)) +
    theme_bw(12) +
    xRot() +
    facet_grid(cols=vars(pathway))+#,space="free", scales="free") +
    theme(strip.text.y=element_text(angle=0))+
    optimized_theme_fig()
  ggsave(basedir("GSEA_plot_","selected_pathways",dbx,".pdf"), w=30,
         h=length(unique(pDT$pathway)) * 0.2 + 3, limitsize = FALSE)}


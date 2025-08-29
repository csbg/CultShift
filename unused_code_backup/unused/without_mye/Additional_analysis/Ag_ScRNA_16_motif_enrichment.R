###############
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(msigdbr)
library(httr)
library(jsonlite)
require(latex2exp)
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10) 
library(biomaRt)
# Load the library
library(rtracklayer)# Replace with appropriate genome version
#####################################################################

base  <-  "Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/"
basedir  <-  dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
outdir <- dirout("Ag_ScRNA_16_motif_enrichment/")

meta <- read_rds(basedir("NTC_meta.rds"))
meta$sample1 <- rownames(meta)

limmaRes <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))#%>%
 # mutate(coef = gsub("interaction","",coef))

########################################################################

ENRICHR<- dirout(paste0("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/ENRICHR"))
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


#############################################################################
#functions-------------------------------------------------------------------
#############################################################################


adj_p_cutoff <- 0.05
logfc_cutoff <- 1
summary_df <- limmaRes %>%
  group_by(celltype) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")


# Filter rows where Count > 10
summary_df %>%
  filter(Count > 10) %>%
  ggplot(aes(x=celltype,y = Count, fill = Regulation)) +
  #geom_col()+
  geom_col(position = position_dodge(preserve = "single")) +    # Separate bars for Upregulated and Downregulated
  scale_fill_manual(values = c("Downregulated" =  "#4C889C", "Upregulated" = "#D0154E")) +
  #facet_grid(rows = vars(celltype)) +
  labs(title="No: of genes with interaction efect of culture condition in KOs",
       x="KOs")+
  scale_y_log10() +  # Apply log10 scale to the y-axis
  optimized_theme_fig()
ggsave(outdir("Differential_genes_count.pdf"))
#function


# Separate upregulated (logFC > 1) and downregulated (logFC < -1) genes by celltype
top_genes_list_with_celltype <- limmaRes %>%
  filter(adj.P.Val < 0.05) %>%                # Filter by adjusted P-value < 0.05 (significant genes)
  mutate(Direction = case_when(               # Create a Direction column to classify genes
    logFC > 1  ~ "Upregulated",               # Genes with logFC > 1 are upregulated
    logFC < -1 ~ "Downregulated",             # Genes with logFC < -1 are downregulated
    TRUE      ~ NA_character_)) %>%           # Ignore genes with 1 <= logFC <= -1 (if needed)
  filter(!is.na(Direction)) %>%               # Remove NA values (genes with logFC between -1 and 1)
  group_by(celltype, Direction) %>%           # Group by both celltype and Direction
  arrange(desc(abs(logFC))) %>%               # Sort by absolute value of logFC in descending order
  summarise(genes = list(unique(genes)), .groups = "drop") %>%  # Create a list of unique genes
  mutate(genotype = "NTC") %>%                # Add genotype info
  split(.$Direction) %>%                      # Split the data frame by Direction (Upregulated and Downregulated)
  lapply(function(df) split(df, df$celltype)) # Split further by celltype within each Direction


enrichr_dbs = c(
              "TRRUST_Mouse_2020"
              )
# Perform enrichment for each gene list (upregulated and downregulated) by celltype
enrichment_results <- lapply(names(top_genes_list_with_celltype$Upregulated), function(name) {
  
  # Extract the upregulated and downregulated gene lists for the current celltype
  up_genes <- top_genes_list_with_celltype$Upregulated[[name]]$genes[[1]]
  down_genes <- top_genes_list_with_celltype$Downregulated[[name]]$genes[[1]]
  
  # Log the length of the gene lists to see if they are empty
  message(paste("Number of Upregulated genes for", name, ":", length(up_genes)))
  message(paste("Number of Downregulated genes for", name, ":", length(down_genes)))
  
  enrichment_list <- list()
  
  # Perform enrichment analysis for upregulated genes if they exist
  if (length(up_genes) > 0) {
    enrich_result_up <- enrichr(up_genes, enrichr_dbs)
    if (!is.null(enrich_result_up[[enrichr_dbs]]) && nrow(enrich_result_up[[enrichr_dbs]]) > 0) {
      result_df_up <- enrich_result_up[[enrichr_dbs]]
      result_df_up$celltype <- name
      result_df_up$direction <- "Upregulated"
      enrichment_list$upregulated <- result_df_up
    } else {
      message(paste("No enrichment results for Upregulated genes in:", name))
    }
  }
  
  # Perform enrichment analysis for downregulated genes if they exist
  if (length(down_genes) > 0) {
    enrich_result_down <- enrichr(down_genes, enrichr_dbs)
    if (!is.null(enrich_result_down[[enrichr_dbs]]) && nrow(enrich_result_down[[enrichr_dbs]]) > 0) {
      result_df_down <- enrich_result_down[[enrichr_dbs]]
      result_df_down$celltype <- name
      result_df_down$direction <- "Downregulated"
      enrichment_list$downregulated <- result_df_down
    } else {
      message(paste("No enrichment results for Downregulated genes in:", name))
    }
  }
  
  # Return combined enrichment results for the current celltype (up and downregulated)
  if (length(enrichment_list) > 0) {
    return(enrichment_list)
  } else {
    return(NULL)
  }
})

# Combine all non-NULL enrichment results into a single data frame
enrichment_results_df <- do.call(rbind, unlist(enrichment_results, recursive = FALSE))

# View the combined results
if (!is.null(enrichment_results_df)) {
  head(enrichment_results_df)
} else {
  message("No significant enrichment results found.")
}

# Save results to a CSV file if data is available
if (!is.null(enrichment_results_df)) {
  write.csv(enrichment_results_df, "TRUST_enrichment_results.csv", row.names = FALSE)
} else {
  message("No enrichment results to save.")
}

# Visualizing the top enriched transcription factors using dot plot
if (!is.null(enrichment_results_df)) {
  top_terms <- enrichment_results_df %>%
    filter(Odds.Ratio > 10) %>%
    filter(Adjusted.P.value < 0.05) %>%
    filter(Combined.Score > 5) %>%
    #slice_head(n = 30) %>%
    pull(Term)
  
  top_results <- enrichment_results_df %>%
    filter(Term %in% top_terms)
  
  # Dot plot to visualize the top enriched transcription factors
  ggplot(top_results, aes(y = Term, 
                          x = celltype, 
                          size = -log10(Adjusted.P.value), 
                          color = log2(Odds.Ratio))) +
    geom_point() +
    scale_color_gradient2(low = "#4C889C",
                          mid = "white",
                          high = "#D0154E") +
    labs(title = "Top 10 Enriched Transcription Factors (Dot Plot)",
         x = "Celltype",
         y = "Transcription Factor",
         size = "-log10(Adjusted P-value)",
         color = "log2(Odds.Ratio)") +
    theme_minimal() +
    facet_grid(rows = vars(direction), space = "free", scales = "free") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
} else {
  message("No results to visualize.")
}
##########################
#promoter motif
# Install and load the necessary libraries

# Load required libraries
library(dplyr)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(org.Mm.eg.db)

# Function to extract promoters
extract_promoters <- function(gene_symbols, genome, upstream = 2000, downstream = 200) {
  # Map gene symbols to Entrez IDs
  gene_mapping <- select(org.Mm.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns = "ENTREZID")
  
  # Remove NA mappings
  gene_mapping <- gene_mapping[!is.na(gene_mapping$ENTREZID), ]
  
  # Load TxDb for mouse genome
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  
  # Extract promoter regions for each gene
  promoter_coords <- promoters(genes(txdb, filter = list(gene_id = gene_mapping$ENTREZID)),
                               upstream = upstream, downstream = downstream)
  
  # Extract promoter sequences
  promoter_seqs <- getSeq(genome, promoter_coords)
  
  # Add gene symbols to promoter sequences
  names(promoter_seqs) <- gene_mapping$SYMBOL[match(names(promoter_seqs), gene_mapping$ENTREZID)]
  
  return(promoter_seqs)
}

# Function to process celltype promoters
process_celltype_promoters <- function(top_genes_list_with_celltype, genome) {
  promoter_list <- lapply(names(top_genes_list_with_celltype$Upregulated), function(celltype) {
    message(paste("Processing celltype:", celltype))
    
    # Extract upregulated and downregulated gene lists
    up_genes <- if (!is.null(top_genes_list_with_celltype$Upregulated[[celltype]])) {
      top_genes_list_with_celltype$Upregulated[[celltype]]$genes[[1]]
    } else {
      character(0)
    }
    
    down_genes <- if (!is.null(top_genes_list_with_celltype$Downregulated[[celltype]])) {
      top_genes_list_with_celltype$Downregulated[[celltype]]$genes[[1]]
    } else {
      character(0)
    }
    
    message(paste("Upregulated genes:", length(up_genes)))
    message(paste("Downregulated genes:", length(down_genes)))
    
    # Extract promoters if there are genes
    up_promoters <- if (length(up_genes) > 0) extract_promoters(up_genes, genome) else NULL
    down_promoters <- if (length(down_genes) > 0) extract_promoters(down_genes, genome) else NULL
    
    list(Upregulated = up_promoters, Downregulated = down_promoters)
  })
  
  names(promoter_list) <- names(top_genes_list_with_celltype$Upregulated)
  return(promoter_list)
}

# Processed input data: top_genes_list_with_celltype
# Ensure this object is structured properly before running the function

# Example: Generate promoter list for mouse genome (MM10)
promoter_list <- process_celltype_promoters(top_genes_list_with_celltype, genome = BSgenome.Mmusculus.UCSC.mm10)

# Save the extracted promoters to FASTA files
save_promoters_to_fasta <- function(promoter_list) {
  
  
  lapply(names(promoter_list), function(celltype) {
    up_fasta <- outdir(paste0(celltype, "_Upregulated_promoters.fasta"))
    down_fasta <- paste0(celltype, "_Downregulated_promoters.fasta")
    
    if (!is.null(promoter_list[[celltype]]$Upregulated)) {
      writeXStringSet(promoter_list[[celltype]]$Upregulated, outdir(paste0(celltype, "_Upregulated_promoters.fasta")))
      message(paste("Saved Upregulated promoters for", celltype, "to", up_fasta))
    } else {
      message(paste("No Upregulated promoters for", celltype))
    }
    
    if (!is.null(promoter_list[[celltype]]$Downregulated)) {
      writeXStringSet(promoter_list[[celltype]]$Downregulated, outdir(paste0(celltype, "_Downregulated_promoters.fasta")))
      message(paste("Saved Downregulated promoters for", celltype, "to", down_fasta))
    } else {
      message(paste("No Downregulated promoters for", celltype))
    }
  })
}

# Save promoter sequences to FASTA files
save_promoters_to_fasta(promoter_list)

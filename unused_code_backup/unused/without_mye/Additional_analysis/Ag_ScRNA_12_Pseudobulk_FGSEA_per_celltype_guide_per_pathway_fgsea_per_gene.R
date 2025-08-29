###############
source("src/00_init.R")
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
require(fgsea)
library(latex2exp)

source("src/Ag_Optimized_theme_fig.R")
#####################################################################
inDir  <-  dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
InDir1 <- dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide")
InDir5 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
base  <-  "Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_per_pathway_fgsea_per_gene"
basedir  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_per_pathway_fgsea_per_gene")
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


#data
gsea.res <- read_rds(InDir1("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))

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
correlation_deg <- read_rds(InDir1("correlation_deg.rds"))
KO_list <- correlation_deg %>% filter(correlation < 0.5, num_degs >= 10) %>%
  pull(genotype)%>%
  unique()
koi <- Reduce(intersect, list(selected_KOs,  KO_list))#, coefficients)) #only valid_k

# Loop over each database and generate the plot
for (dbx in databases) {
  # Filter for the current database and the pathways of interest
  pDT <- gsea.res[db == dbx]
  
  # Select top pathways based on NES
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5), by = c("coef", "celltype", "pathway")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by = c("coef", "celltype", "pathway")]$pathway)
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  
  # Filter for displayed pathways and KOI
  pDT <- pDT[pathway %in% pw.display] %>%
    filter(coef %in% paste0("interaction", koi))
  
  # Convert logical columns if necessary
  logical_columns <- sapply(pDT, is.logical)
  if (any(logical_columns)) {
    pDT[logical_columns] <- lapply(pDT[logical_columns], as.character)
  }
  
  # Ensure 'celltype' and 'NES' columns are appropriate types
  if (!is.character(pDT$celltype) && !is.factor(pDT$celltype)) {
    stop("Column 'celltype' must be character or factor.")
  }
  if (!is.numeric(pDT$NES)) {
    stop("Column 'NES' must be numeric.")
  }
  for (col_name in names(pDT)) {
    if (is.list(pDT[[col_name]])) {
      pDT[[col_name]] <- sapply(pDT[[col_name]], toString)  # Convert list to string
    }
  }
  # Set up the output directory and write table
  dir <- dirout(paste0(base, "/", dbx, "/"))
  
  write.csv(pDT, file = dir(paste0(dbx, ".csv")))
  
  #Plot GSEA results for the current database
  plot <- ggplot(pDT, aes(x = gsub("interaction", "", coef), y = celltype, color = NES, size = pmin(5, -log10(padj)))) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    geom_point(data = pDT[padj < 0.05], shape = 1) +
    scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
    theme_bw(12) +
    xRot() +
    facet_wrap(vars(pathway)) +
    theme(strip.text.y = element_text(angle = 0))
  
  
  
  # Save the plot with a unique filename for each database
  ggsave(dir(paste0( "GSEA_plot_", dbx, "_selected_pathways.pdf")),
         plot = plot, width = 30,
         height = length(unique(pDT$pathway)) * 0.2 + 3,
         limitsize = FALSE)
}

dbx <-"Reactome_2022"

for (dbx in unique(gsea.res$db)) {
  #for (genotype in koi) {
  dat <- dirout(paste0(base, "FGSEA/", dbx))
  # Define file path
  output_file <- dat("GSEA_significant_", dbx, ".tsv")
  
  # Ensure the directory exists before writing
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)  # Create directory if it doesn't exist
  }
  # Convert list columns to character before writing
  df <- gsea.res[db == dbx]
  df[] <- lapply(df, function(x) if (is.list(x)) sapply(x, toString) else x)
  
  # Now, write the table
  write.table(df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  
  
  
  # Get the pathways for this database
  pDT <- gsea.res[db == dbx ]
  
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n = 5), by = c("celltype")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n = 5), by = c("celltype")]$pathway)
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  
  pDT <- pDT[pathway %in% pw.display]
  
  if (nrow(pDT) > 0) {
    # Aggregate ABSOLUTE NES values across all cell types- change can be in any direction
    pDT_agg <- pDT %>%
      group_by(pathway) %>%
      summarize(average_NES = mean(abs(NES), na.rm = TRUE)) %>%
      arrange(desc(average_NES)) 
    
    #pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)
    # Ensure pDT follows the same order as pDT_agg
    pDT <- pDT %>%
      mutate(pathway = factor(pathway, levels = pDT_agg$pathway)) %>%
      arrange(factor(pathway, levels = pDT_agg$pathway))  # Explicitly reorder pDT
    # Step 3: Plot with the new pathway order (highest NES first)
    pDT$genotype <- gsub("interaction","",pDT$coef)
    pDT <- inner_join(pDT,ko_flags[,c("celltype","genotype","valid_ko")],by=c("genotype", "celltype"))%>%
      filter(valid_ko)
    
    #pDT <- pDT %>%
    #  filter(valid_ko)
    ggplot(pDT, aes(y=genotype, x=pathway, color=pmin(2,NES), size=pmin(3, -log10(padj)))) +
      
      scale_color_gradient2(low = "#4C889C",
                            mid = "white",
                            high = "#D0154E"
      )+#,
      #name=TeX("log_{2}(FC)"))+
      geom_point() +
      scale_size_continuous(
        range = c(0, 2.5),
        limits = c(0, 5),
        name=TeX("$-\\log_{10}(p_{adj})$"))+
      
      
      #xRot() +
      facet_grid(cols=vars(celltype),space="free", scales="free")+
      labs(y = "celltype",
           x = "Pathways",
           title = "Enriched pathways",
           size = "log10(padj)")+
      
      optimized_theme_fig()+
      coord_flip()
    
    ggsave(dat(paste0(dbx, "_top_pathways.pdf")))
    # Extract pathways that were plotted
    pathways <- unique(pDT$pathway)
    # Split the leadingEdge column into separate genes and create a combined list for each pathway
    combined_genes <- pDT %>%
      group_by(pathway) %>%
      summarise(
        Genes = list(unique(unlist(leadingEdge))),  # Combine and take the union of gene lists
        .groups = "drop"
      )
    
    # View the cleaned-up result
    print(combined_genes)
    
    
    
    # Function to extract top genes for a pathway (based on the db)
    get_top_genes <- function(pathway_name, 
                              dbx,
                              limma_results,
                              logFC_threshold = 1,
                              pval_threshold = 0.05,
                              top_n = 5) {
      # Dynamically get pathway genes based on the database
      pathway_genes <- combined_genes[combined_genes$pathway == pathway_name,]$Genes
      pathway_genes <- pathway_genes[[1]]
      #pathway_genes <- pathway_genes[pathway_genes != ""]
      # Filter limma results for the pathway genes
      filtered_genes <- limma_results %>%
        filter(group != "n.s") %>%
        group_by(celltype,coef) %>%
        filter(adj.P.Val < pval_threshold) %>%
        filter(ensg %in% pathway_genes) %>%
        arrange(adj.P.Val) %>%
        arrange(desc(abs(logFC))) %>%
        slice_head(n = top_n) %>%
        pull(ensg) %>%
        unique()
      
      
      # Extract and return data for the filtered genes
      pathway_plot <- limma_results %>%
        filter(toupper(ensg) %in% toupper(filtered_genes))
      
      #%>%
      #filter(group != "n.s") 
      
      return(pathway_plot)
    }
    
    pathways <- unique(combined_genes$pathway)
    # Apply function to all pathways in this db
    results <- map(pathways, function(pathway) {
      get_top_genes(
        pathway_name = pathway,
        dbx = dbx,
        limma_results = limmaRes
      )
    })
    
    # Optionally, name the elements of the list using pathway names
    names(results) <- pathways
    # Combine results into one dataframe
    combined_genes_filtered <- bind_rows(results,.id = "gene_set")
    
    # Convert cell type to a factor for correct ordering
    combined_genes_filtered$celltype <- factor(combined_genes_filtered$celltype,
                                               levels = c("HSC", "MEP.early", "MkP", 
                                                          "GMP", "Gran.P", "Gran.", 
                                                          "Mono", "Eo.Ba"))
    
    # Create a directory for the database
    dat1 <- dirout(paste0(base, "FGSEA/", dbx, "/GenePlots/"))
    #dir.create(file.path(base, "GenePlots", dbx), showWarnings = FALSE, recursive = TRUE)
    gene_sets <- unique(combined_genes_filtered$gene_set)
    walk(unique(gene_sets), function(pathway) {
      combined_genes_filtered$genotype <- gsub("interaction","",combined_genes_filtered$coef)
      pathway_data <- combined_genes_filtered %>% filter(gene_set == pathway)%>%
        inner_join(ko_flags, by =c("genotype","celltype"))%>%
        filter(valid_ko)
      
      plot <- ggplot(pathway_data, aes(x = genotype, y = ensg,
                                       color = pmin(3, pmax(-3, logFC)),
                                       size = pmin(5, -log10(adj.P.Val)))) +
        geom_point() +  
        facet_grid(cols = vars(celltype),space = "free",scales = "free")+
        scale_color_gradient2(low = "#4C889C", mid = "white", high = "#D0154E") +
        scale_size_continuous(range = c(0, 3), name = TeX("$-\\log_{10}(p_{adj})$")) +
        labs(title = paste("Differentially Expressed Genes -", pathway),
             y = "Gene",
             x = "Cell Type",
             color = "logFC",
             size = "-log10(padj)") +
        optimized_theme_fig()+
        theme(strip.text.x = element_text(angle = 90))
      
      # 🔹 Sanitize pathway name (replace special characters)
      safe_pathway <- gsub("[^A-Za-z0-9_-]", "_", pathway)  # Replaces problematic characters
      
      # 🔹 Construct path correctly
      # path1 <- file.path(base, "FGSEA", dbx, safe_pathway)
      # dir.create(path1, showWarnings = FALSE, recursive = TRUE)  # Ensure directory exists
      
      # 🔹 Save plot
      ggsave(dat1(paste0(safe_pathway, "_top_pathways.pdf")), 
             plot = plot)
    })
    
  }
}

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

################################################################################
# Set up directories and file paths
inDir <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
base <- "Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/"
basedir <- dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/")

FGSEA <- dirout(paste0(base, "FGSEA/"))
ENRICHR <- dirout(paste0(base, "ENRICHR/"))

# Load limma results
dataVoom_NTC <- read_rds(inDir("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaRes_NTC <- read_rds(inDir("limma_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(inDir("NTC_meta.rds"))

################################################################################
# Download gene sets from EnrichR
ENRICHR.DBS <-union(ENRICHR.DBS,
                    c("GO_Biological_Process_2021",
                      "TRRUST_Transcription_Factors_2019",
                      "Reactome_2022",
                      "GO_Molecular_Function_2023",
                      "GO_Biological_Process_2023",
                      "CellMarker_2024"))
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
# Convert to mouse
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "", c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl) {
  dbl <- lapply(dbl, function(gs) {
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})

################################################################################
# EnrichR Enrichment Analysis
perform_enrichment_analysis <- function(limmaRes, direction) {
  enr.res.list <- list()
  
  for (ct in unique(limmaRes$celltype)) {
    # Extract genes of interest (GOI)
    goi <- limmaRes %>%
      filter(group == direction & celltype == ct) %>%
      pull(genes)
    
    # Perform enrichment analysis
    enr.res <- enrichr(goi, databases = c("KEGG_2019_Mouse",
                                          "MSigDB_Hallmark_2020",
                                          "WikiPathways_2019_Mouse",
                                          "GO_Biological_Process_2021",
                                          "TRRUST_Transcription_Factors_2019",
                                          "Reactome_2022",
                                          "GO_Molecular_Function_2023",
                                          "GO_Biological_Process_2023",
                                          "CellMarker_2024"))
    
    # Combine results
    enr.res <- bind_rows(enr.res, .id = "db")
    enr.res.list[[ct]] <- enr.res
  }
  
  # Combine results for all cell types
  bind_rows(enr.res.list, .id = "celltype")
}

# Perform enrichment analysis for up and down-regulated genes
limmaRes_NTC <- limmaRes_NTC[limmaRes_NTC$celltype != "MEP",]
enr.res.all_up <- perform_enrichment_analysis(limmaRes_NTC, "up")
write_rds(enr.res.all_up, ENRICHR("enr.res.all_NTC_up.rds"))

enr.res.all_down <- perform_enrichment_analysis(limmaRes_NTC, "down")
write_rds(enr.res.all_down, ENRICHR("enr.res.all_NTC_down.rds"))

################################################################################
# Plot Enrichment Results
enr.res.all_up <- read_rds(ENRICHR("enr.res.all_NTC_up.rds"))
enr.res.all_down <- read_rds(ENRICHR("enr.res.all_NTC_down.rds"))
dfs <- list(enr.res.all_up, enr.res.all_down)
names(dfs) <- c("up", "down")
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
pmap(list(names(dfs), dfs), function(name, df) {
  walk(databases, ~{
    db <- .x
    df_filtered <- df[df$db == db & df$celltype != "MEP", ]
    
    if (nrow(df_filtered) > 0) {
      # Remove the prefix from column names
      colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
      
      # Perform filtering and mutation
      plotting_enr <- df_filtered %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
      
      ggplot(plotting_enr, aes(x = celltype, y = Term, color = log2(Odds.Ratio), size = pmin(10, neg.log10.Adjusted.P.value))) +
        geom_point() +
        scale_size_continuous(range = c(2, 6)) +
        scale_color_gradientn(colors = c("pink", "red")) +
        ggtitle(paste0(db)) +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              title = element_text(size = 18))
      
      ggsave(ENRICHR(paste0(name, "_per.celltype_", db, ".pdf")), w = 10, h = length(unique(plotting_enr$Term)) * 0.2 + 3, limitsize = FALSE)
    }
  })
})

##################
# Without Filter Plot
##################
pmap(list(names(dfs), dfs), function(name, df) {
  walk(databases, ~{
    db <- .x
    df_filtered <- df[df$db == db, ]
    df_filtered <- df_filtered %>% mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
    
    if (nrow(df_filtered) > 0) {
      # Remove the prefix from column names
      colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
      
      plotting_enr <- df_filtered %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        pull(Term)
      enr <- df_filtered %>% filter(Term %in% plotting_enr)
      
      ggplot(enr, aes(x = celltype, y = Term, color = log2(Odds.Ratio), size = pmin(10, neg.log10.Adjusted.P.value))) +
        geom_point() +
        scale_size_continuous(range = c(2, 6)) +
        scale_color_gradientn(colors = c("pink", "red"))
      
      ggsave(ENRICHR(paste0(name, "_per.celltype_", db, "without_filter.pdf")), w = 10, h = length(unique(enr$Term)) * 0.2 + 3, limitsize = FALSE)
    }
  })
})

################################################################################
# FGSEA
################################################################################
gsea.res <- data.table()
for (ct in unique(limmaRes_NTC$celltype)) {
  for (dbx in names(enr.terms)) {
    subset_limmaRes_NTC <- limmaRes_NTC[limmaRes_NTC$celltype == ct, ]
    stats <- with(subset_limmaRes_NTC, setNames(logFC, nm = genes))
    
    if (any(is.na(stats))) {
      next
    }
    
    fgsea_output <- fgsea(pathways = enr.terms[[dbx]], stats = stats)
    
    if (length(fgsea_output) > 0) {
      gsea.res <- rbind(gsea.res, data.table(fgsea_output, celltype = ct, db = dbx))
    }
  }
}
gsea.res %>% write_rds(basedir("NTC_fgsea.rds"))

################################################################################
# Plot FGSEA Results
################################################################################
gsea.res <- read_rds(basedir("NTC_fgsea.rds"))
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][, -c("log2err", "NES", "size", "pval"), with = F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, function(vec) paste(vec[1:10], collapse = ","))

for (dbx in unique(gsea.res$db)) {
  dat <- dirout(paste0(base, "FGSEA/", dbx))
  write.tsv(gsea.res.export[db == dbx], dat("GSEA_significant_", dbx, ".tsv"))
  pDT <- gsea.res[db == "MSigDB_Hallmark_2020"]
  pDT <- gsea.res[db == dbx]
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n = 5), by = c("celltype")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n = 5), by = c("celltype")]$pathway)
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  
  pDT <- pDT[pathway %in% pw.display]
  
  
  
  if (nrow(pDT) > 0) {
    # Step 1: Aggregate NES values across all celltypes (mean of NES per pathway)
    pDT_agg <- pDT %>%
      group_by(pathway) %>%
      summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
      arrange(desc(average_NES))  # Ordering pathways by the average NES, highest first
    
    # Step 2: Create a factor for pathway that reflects the aggregated NES order
    pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)
    
    # Step 3: Plot with the new pathway order (highest NES first)
    ggplot(pDT, aes(x = celltype, y = pathway, size = pmin(10, -log10(padj)), color = NES)) +
      geom_point() +
      scale_size_continuous(range = c(3, 6)) +
      scale_color_gradientn(colors = c("blue", "white", "red")) +
      theme(axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            title = element_text(size = 18))
    
    
  }
}



# Iterate over each database
for (dbx in unique(gsea.res$db)) {
  dat <- dirout(paste0(base, "FGSEA/", dbx))
  write.tsv(gsea.res.export[db == dbx], dat("GSEA_significant_", dbx, ".tsv"))
  
  # Get the pathways for this database
  pDT <- gsea.res[db == dbx]
  
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n = 5), by = c("celltype")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n = 5), by = c("celltype")]$pathway)
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  
  pDT <- pDT[pathway %in% pw.display]
  
  if (nrow(pDT) > 0) {
    # Aggregate NES values across all cell types
    pDT_agg <- pDT %>%
      group_by(pathway) %>%
      summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
      arrange(desc(average_NES)) 
    
    pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)
    
    # Step 3: Plot with the new pathway order (highest NES first)
    ggplot(pDT, aes(y=celltype, x=pathway, color=pmin(2,NES), size=pmin(5, -log10(padj)))) +
      
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
      #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
      labs(y = "celltype",
           x = "Pathways",
           title = "Enriched pathways",
           size = "log10(padj)")+
      
      optimized_theme_fig()+
      coord_flip()
    ggsave(outdir("Fig1.4fGSEA_plot_",dbx,"._long.pdf"),plot = Fig1.4, w=12,h=7, units = "cm")
    # Extract pathways that were plotted
    pathways <- unique(pDT$pathway)
    
    # Function to extract top genes for a pathway (based on the db)
    get_top_genes <- function(pathway_name, dbx, limma_results, logFC_threshold = 1, pval_threshold = 0.05, top_n = 5) {
      # Dynamically get pathway genes based on the database
      pathway_genes <- unique(unlist(enr.terms[[dbx]][[pathway_name]]))
      
      # Filter limma results for the pathway genes
      filtered_genes <- limma_results %>%
        filter(group != "n.s") %>%
        group_by(celltype) %>%
        filter(adj.P.Val < pval_threshold) %>%
        filter(toupper(genes) %in% toupper(pathway_genes)) %>%
        arrange(adj.P.Val) %>%
        slice_head(n = top_n) %>%
        pull(genes) %>%
        unique()
      
      # Extract and return data for the filtered genes
      pathway_plot <- limma_results %>%
        filter(toupper(genes) %in% toupper(filtered_genes)) %>%
        mutate(gene_set = pathway_name)
      
      return(pathway_plot)
    }
    
    # Apply function to all pathways in this db
    results <- map(pathways, function(pathway) {
      get_top_genes(
        pathway_name = pathway,
        dbx = dbx,
        limma_results = limmaRes_NTC
      )
    })
    
    # Combine results into one dataframe
    combined_genes_filtered <- bind_rows(results)
    
    # Convert cell type to a factor for correct ordering
    combined_genes_filtered$celltype <- factor(combined_genes_filtered$celltype,
                                               levels = c("HSC", "MEP.early", "MkP", 
                                                          "GMP", "Gran.P", "Gran.", 
                                                          "Mono", "Eo.Ba"))
    
    # Create a directory for the database
    dir.create(file.path(base, "GenePlots", dbx), showWarnings = FALSE, recursive = TRUE)
    
    # Generate separate plots for each pathway
    walk(unique(combined_genes_filtered$gene_set), function(pathway) {
      pathway_data <- combined_genes_filtered %>% filter(gene_set == pathway)
      
      pathway_plot <- ggplot(pathway_data, aes(x = celltype, y = genes,
                                               color = pmin(3, pmax(-3, logFC)),
                                               size = pmin(5, -log10(adj.P.Val)))) +
        geom_point() +  
        scale_color_gradient2(low = "#4C889C", mid = "white", high = "#D0154E") +
        scale_size_continuous(range = c(0, 2), name = TeX("$-\\log_{10}(p_{adj})$")) +
        labs(title = paste("Differentially Expressed Genes -", pathway),
             y = "Gene",
             x = "Cell Type",
             color = "logFC",
             size = "-log10(padj)") +
        optimized_theme_fig()
      
      # Save plot to respective folder
      ggsave(filename = file.path(base, "GenePlots", dbx, paste0(pathway, ".png")), 
             plot = pathway_plot, width = 8, height = 6)
    })
  }
}

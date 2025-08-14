##################################
#enrichment function using enrichr
#################################

perform_enrichment_analysis <- function(data,
                                        logFC_cutoff_up,
                                        logFC_cutoff_down,
                                        databases,
                                        output_file_prefix) {
  # Filter data based on coefficient and group and extract genes
  genes_up <- data %>%
    filter(group == "up" & logFC > logFC_cutoff_up) %>%
    pull(ensg)
  
  genes_down <- data %>%
    filter(group == "down" & logFC < logFC_cutoff_down) %>%
    pull(ensg)
  
  # Perform enrichment analysis for upregulated genes
  perform_enrichment <- function(genes, output_file) {
    enr.res <- enrichr(genes, databases = databases)
    enr.res <- bind_rows(enr.res, .id = "db")
    
    map(databases, ~{
      db <- .x
      plotting_enr <- enr.res[enr.res$db == db,] %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
      
      ggplot(plotting_enr, aes(x = db, y = Term, size = log2(Odds.Ratio),
                               color = neg.log10.Adjusted.P.value)) +
        geom_point() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
        theme_bw(12)
      
      ggsave(out(paste0(output_file, "_", db, ".pdf")), w = 20, h = length(unique(plotting_enr$Term)) * 0.2 + 3, limitsize = FALSE)
    })
  }
  
  perform_enrichment(genes_up, paste0(output_file_prefix, "_upregulated.pdf"))
  perform_enrichment(genes_down, paste0(output_file_prefix, "_downregulated.pdf"))
}

####################################################

# Example usage:

# For all cell types together
perform_enrichment_analysis(data = data,
                            logFC_cutoff_up = 1, 
                            logFC_cutoff_down = -1,
                            databases = c("KEGG_2019_Mouse", "MSigDB_Hallmark_2020", "WikiPathways_2019_Mouse", "GO_Biological_Process_2021"),
                            output_file_prefix = "all_celltypes")

# For each cell type separately
for (celltype in unique(meta$celltype)) {
  celltype_data <- limmaRes[limmaRes$celltype == celltype, ]
  perform_enrichment_analysis(data = celltype_data,
                              logFC_cutoff_up = 1, 
                              logFC_cutoff_down = -1,
                              databases = c("KEGG_2019_Mouse", "MSigDB_Hallmark_2020", "WikiPathways_2019_Mouse", "GO_Biological_Process_2021"),
                              output_file_prefix = paste0(celltype, "_celltype"))
}


# Example usage:

# For all cell types together
perform_enrichment_analysis(data = limmaRes,
                            coef_col = logFC,
                            group_col = group,
                            logFC_cutoff_up = 1, 
                            logFC_cutoff_down = -1,
                            databases = c("KEGG_2019_Mouse", "MSigDB_Hallmark_2020", "WikiPathways_2019_Mouse", "GO_Biological_Process_2021"),
                            output_file_prefix = "all_celltypes")

# For each cell type separately
limmaRes<-top
for (celltype in unique(meta$celltype)) {
  celltype_data <- limmaRes[limmaRes$celltype == celltype, ]
  perform_enrichment_analysis(data = celltype_data,
                              coef_col = logFC,
                              group_col = group,
                              logFC_cutoff_up = 1, 
                              logFC_cutoff_down = -1,
                              databases = c("KEGG_2019_Mouse", "MSigDB_Hallmark_2020", "WikiPathways_2019_Mouse", "GO_Biological_Process_2021"),
                              output_file_prefix = paste0(celltype, "_celltype"))
}

# Example usage:

# For all cell types together
perform_enrichment_analysis(data = limmaRes,
                            coef = logFC,
                            group_col = group,
                            logFC_cutoff_up = 1, 
                            logFC_cutoff_down = -1,
                            databases = c("KEGG_2019_Mouse", "MSigDB_Hallmark_2020", "WikiPathways_2019_Mouse", "GO_Biological_Process_2021"),
                            output_file_prefix = "all_celltypes")

# For each cell type separately
for (celltype in unique(meta$celltype)) {
  celltype_data <- limmaRes[limmaRes$celltype == celltype, ]
  perform_enrichment_analysis(data = celltype_data,
                              coef = logFC,
                              group_col = group,
                              logFC_cutoff_up = 1, 
                              logFC_cutoff_down = -1,
                              databases = c("KEGG_2019_Mouse", "MSigDB_Hallmark_2020", "WikiPathways_2019_Mouse", "GO_Biological_Process_2021"),
                              output_file_prefix = paste0(celltype, "_celltype"))
}

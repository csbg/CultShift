# Obtains top n (based on adj. P.Val from limma Results) for each celltype 
#genes from a given pathway

get_top_genes <- function(pathway_name,
                          pathway_genes,
                          limma_results, logFC_threshold = 1,
                          pval_threshold = 0.05, top_n = 5) {
  # Filter limma results for the pathway genes
  filtered_genes <- limma_results %>%
    filter(group != "n.s") %>%  # Skip non-significant genes
    filter(adj.P.Val < pval_threshold, abs(logFC) > logFC_threshold) %>%
    filter(toupper(ensg) %in% toupper(pathway_genes)) %>%  # Match pathway genes
    group_by(celltype, coef) %>%
    arrange(adj.P.Val) %>%  # Sort by adjusted p-value
    slice_head(n = top_n) %>%  # Get the top 5 genes per group (celltype)
    pull(ensg) %>% 
    unique()
  
  # Filter the main table to return data for plotting
  pathway_plot <- limma_results %>%
    filter(toupper(ensg) %in% toupper(filtered_genes), coef %in% koi)%>%
    mutate(pathway = pathway_name)
  
  return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
}
get_top_genes_pathway <- function(pathway_name,
                          pathway_genes,
                          limma_results, 
                          logFC_threshold = 1,
                          pval_threshold = 0.05, 
                          top_n = 10) {
  # Filter limma results for pathway genes and significance thresholds
  filtered_genes <- limma_results %>%
    filter(group != "n.s") %>%  # Exclude non-significant genes
    filter(adj.P.Val < pval_threshold, abs(logFC) > logFC_threshold) %>%
    filter(toupper(ensg) %in% toupper(pathway_genes))  # Match pathway genes
  
  # Select top 10 genes across all cell types and coefficients by adjusted p-value
  top_genes <- filtered_genes %>%
    arrange(adj.P.Val) %>%  # Sort by adjusted p-value
    slice_head(n = top_n) %>%  # Take the top 10 genes
    pull(ensg) %>%
    unique()
  
  # Filter the main table to get data for plotting
  pathway_plot <- limma_results %>%
    filter(toupper(ensg) %in% toupper(top_genes)) %>%
    mutate(pathway = pathway_name)
  
  return(list(top_genes = top_genes, pathway_plot = pathway_plot))
}


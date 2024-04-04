#########################################################
#**Ag_ScRNA_09_pseudobulk_per_celltype_limma_function**#
#########################################################
#Aarathy
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)

########################
performDE <- function(meta, counts, tissue_type1, tissue_type2) {
  # Filtering
  d0 <- DGEList(counts)
  d0 <- calcNormFactors(d0)
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  
  group <- interaction(meta$celltype, meta$tissue)
  unique(group)
  
  mm <- model.matrix(~0 + group)
  rownames(mm) <- rownames(meta)
  
  # Normalization
  dataVoom <- voom(d, mm)
  limmaFit <- lmFit(dataVoom, mm)
  limmaFit <- eBayes(limmaFit)
  
  # Generate comparisons based on tissue types
  
  comparisons <- lapply(unique(meta$celltype), function(celltype) {
    comp1 <- paste0("group", celltype, ".", tissue_type1)
    comp2 <- paste0("group", celltype, ".", tissue_type2)
    list(comp1, comp2, celltype)
  })
  
  # Initialize list to store top table results
  top_table <- list()
  
  # Perform differential expression analysis for each comparison
  for (comp_list in comparisons) {
    group1 <- comp_list[[1]]
    group2 <- comp_list[[2]]
    id <- comp_list[[3]]
    
    # Create contrast matrix
    contrast <- paste0(group1, "-", group2)
    
    # Debugging: Print contrast matrix
    print(contrast)
    
    # Ensure that contrast inputs are accepted
    if (!(all(c(group1, group2) %in% colnames(limmaFit)))) {
      stop("Invalid column names for comparison:", group1, " and ", group2)
    }
    
    # Fit model and perform differential expression analysis
    tmp <- contrasts.fit(limmaFit, makeContrasts(contrast, levels = colnames(coef(limmaFit))))
    tmp <- eBayes(tmp)
    top_table[[id]] <- topTable(tmp, sort.by = "P", n = Inf)
  }
  
  # Bind results and add metadata
  top_table_res <- bind_rows(top_table, .id = 'celltype') %>%
    mutate(genes = gsub("\\...\\d+", "", rownames(.)))
  
  # Create a column 'group' to group the genes as up, down, or ns
  top_table_res$group <- ifelse(top_table_res$logFC >= 1 & 
                                  top_table_res$adj.P.Val <= 0.05, "up", 
                                ifelse(top_table_res$logFC <= -1 & 
                                         top_table_res$adj.P.Val <= 0.05, "down", "n.s"))
  
 
}

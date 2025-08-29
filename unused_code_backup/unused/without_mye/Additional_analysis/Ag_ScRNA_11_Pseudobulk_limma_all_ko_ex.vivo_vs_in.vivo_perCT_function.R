source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
####################################################################
#function
####################################################################

performDE <- function(meta, counts, model_formula) {
  # Check if model_formula is a string, convert to formula object
  if (!is.character(model_formula)) {
    stop("Model formula must be provided as a string.")
  }
  model_formula <- as.formula(model_formula)
  
  # Extract variable names from the formula
  formula_terms <- all.vars(model_formula)
  
  # Check if variables in the formula are present in the metadata
  if (!all(formula_terms %in% names(meta))) {
    missing_vars <- setdiff(formula_terms, names(meta))
    stop(paste("Variables in the model formula are not present in the metadata:", paste(missing_vars, collapse = ", ")))
  }
  
  # edgeR based normalization
  d0 <- counts[, rownames(meta)]
  d0 <- DGEList(d0)
  d0 <- calcNormFactors(d0)
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  
  # Initialize an empty list to store results per cell type
  de_results_list <- list()
  
  # Iterate over each unique cell type
  unique_celltypes <- unique(meta$celltype)
  for (celltype in unique_celltypes) {
    cat("Performing DE analysis for cell type:", celltype, "\n")
    
    # Subset metadata and counts data for the current cell type
    meta_subset <- meta[meta$celltype == celltype, ]
    counts_subset <- counts[, rownames(meta_subset)]
    
    # Define the model matrix internally based on the provided formula and metadata subset
    modelMatrix <- model.matrix(model_formula, data = meta_subset)
    
    # voom
    dataVoom <- voom(counts_subset, modelMatrix, plot = T)  # Disable plot for each cell type
    limmaFit <- lmFit(dataVoom, modelMatrix)
    limmaFit <- eBayes(limmaFit)
    
    # Get unique coefficients for "ex.vivo" and "in.vivo"
    unique_coefs_exvivo <- grep("^ex.vivo:", colnames(limmaFit$coef), value = TRUE)
    unique_coefs_invivo <- setdiff(colnames(limmaFit$coef), unique_coefs_exvivo)
    
    # Check if there are any coefficients present in both "ex.vivo" and "in.vivo"
    common_coefs <- intersect(unique_coefs_exvivo, unique_coefs_invivo)
    if (length(common_coefs) == 0) {
      cat("No common coefficients found for cell type:", celltype, "\n")
      next  # Skip to the next iteration if no common coefficients found
    }
    
    # Perform DE analysis for the common coefficients
    limmaRes <- map_dfr(common_coefs, ~{
      coefx <- .x
      topTable(limmaFit, coef = coefx, number = Inf) %>%
        rownames_to_column("ensg") %>%
        mutate(coef = coefx)
    }) %>%
      filter(coef != "(Intercept)")
    limmaRes <- limmaRes %>%
      mutate(coef = str_replace(coef, "genotype", "")) %>%
      mutate(coef = str_replace(coef, "tissue", ""))
    
    # Up and downregulated genes for the common coefficients
    limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                               limmaRes$adj.P.Val <= 0.05, "up", 
                             ifelse(limmaRes$logFC <= -1 & 
                                      limmaRes$adj.P.Val <= 0.05, "down", "n.s"))
    
    # Store results for the current cell type
    de_results_list[[celltype]] <- limmaRes
  }
  
  return(de_results_list)
}


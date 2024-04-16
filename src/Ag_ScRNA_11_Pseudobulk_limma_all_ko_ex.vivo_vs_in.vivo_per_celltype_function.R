source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
####################################################################
#function

performDE <- function(meta, counts, model_formula, output_dir = "dataVoom_plots/") {
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
  
  # Initialize an empty list to store results per cell type
  de_results_list <- list()
  
  # Iterate over each unique cell type
  unique_celltypes <- unique(meta$celltype)
  for (ct in unique_celltypes) {
    cat("Performing DE analysis for cell type:", ct, "\n")
    
    # Subset metadata and counts data for the current cell type
    meta_subset <- meta %>% filter(celltype == ct)
    counts_subset <- counts[, rownames(meta_subset)]
    
    # Prepare the data for differential expression analysis
    d <- DGEList(counts_subset)
    d <- calcNormFactors(d)
    cutoff <- 1
    drop <- which(apply(cpm(d), 1, max) < cutoff)
    d <- d[-drop, ]
    
    # Define the model matrix internally based on the provided formula and metadata subset
    modelMatrix <- model.matrix(model_formula, data = meta_subset)
    
    # Perform voom transformation and fit the model using limma
    
    
    # Save the dataVoom plot
    pdf(basedir(paste0(ct,"_dataVoom.pdf")))
    dataVoom <- voom(d, modelMatrix,plot = T)
    dev.off()
    
    dataVoom <- voom(d, modelMatrix)
    
    limmaFit <- lmFit(dataVoom, modelMatrix)
    limmaFit <- eBayes(limmaFit)
    
    # Extract and format results for each coefficient
    limmaRes <- map_dfr(colnames(limmaFit$coef), function(coef_name) {
      topTable(limmaFit, coef = coef_name, number = Inf) %>%
        rownames_to_column("ensg") %>%
        mutate(coef = coef_name, 
               group = case_when(
                 logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
                 logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
                 TRUE ~ "n.s"),
               celltype = ct) %>% # Include cell type in the results
        select(-contains("Intercept")) # Optionally remove intercept results if unnecessary
    }, .id = "coefficient")
    
    # Store results for the current cell type
    de_results_list[[ct]] <- limmaRes
  }
  
  # Combine all results into a single data frame
  final_results_table <- bind_rows(de_results_list)
  return(final_results_table)
}


# Usage example


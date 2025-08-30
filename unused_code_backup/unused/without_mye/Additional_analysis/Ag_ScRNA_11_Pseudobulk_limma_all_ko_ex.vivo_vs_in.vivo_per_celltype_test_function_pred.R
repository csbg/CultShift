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
  ct <- "Mono"
  for (ct in unique_celltypes) {
    cat("Performing DE analysis for cell type:", ct, "\n")
    
    # Subset metadata and counts data for the current cell type
    meta_subset <- meta %>% filter(celltype == ct & !is.na(genotype))
    
    counts_subset <- counts[, rownames(meta_subset)]
    
    # Prepare the data for differential expression analysis
    d0 <- DGEList(counts_subset)
    d0 <- calcNormFactors(d0,method = "TMM")
    threshold <- 30
    drop <- which(apply(cpm(d0), 1, max) < threshold)
    d <- d0[-drop, ]
    
    #########################
    cpm_values_before <- cpm(d0)
    
    # Calculate row sums of CPM values before filtering
    row_sums_cpm_before <- rowSums(cpm_values_before)
    
    # Convert row sums to a data frame for ggplot2
    row_sums_cpm_before_df <- data.frame(RowSumsCPM = row_sums_cpm_before)
    
    # Plot the density of row sums of CPM values before filtering
    p_before <- ggplot(row_sums_cpm_before_df, aes(x = log10(RowSumsCPM + 1))) +
      geom_density(fill = "blue", alpha = 0.4) +
      labs(
        title = "Density Plot of Row Sums of CPM Values (Log Scale) Before Filtering",
        x = "Row Sums of CPM (Log Scale)",
        y = "Density"
      ) +
      theme_minimal()
    
    # Save the plot before filtering with default filename
    ggsave(basedir(paste0(ct,"counts_density_before_filtering.pdf")), plot = p_before)
    
    # Report the number of genes remaining after filtering
    cat("Number of genes remaining after threshold filtering:",nrow(d$counts), "\n")
    
    # Calculate CPM for the filtered DGEList
    cpm_values_after <- cpm(d)
    
    # Calculate row sums of CPM values after filtering
    row_sums_cpm_after <- rowSums(cpm_values_after)
    
    # Convert row sums to a data frame for ggplot2
    row_sums_cpm_after_df <- data.frame(RowSumsCPM = row_sums_cpm_after)
    
    # Plot the density of row sums of CPM values after filtering
    p_after <- ggplot(row_sums_cpm_after_df, aes(x = log10(RowSumsCPM + 1))) +
      geom_density(fill = "blue", alpha = 0.4) +
      labs(
        title = "Density Plot of Row Sums of CPM Values (Log Scale) After Filtering",
        x = "Row Sums of CPM (Log Scale)",
        y = "Density"
      ) +
      theme_minimal()
    
    # Save the plot after filtering with threshold in filename
    ggsave(filename = basedir(paste0(ct,"counts_density_after_filtering_threshold_", threshold, ".pdf")), plot = p_after)
    ###############################################################################
    #model matrix
    modelMatrix <- model.matrix(model_formula, data = meta_subset)
    dim(d)           # Should return (genes x samples)
    dim(modelMatrix) #
    
    
    # Save the dataVoom plot
    pdf(basedir(paste0(ct,"_dataVoom.pdf")))
    dataVoom <- voom(d, modelMatrix,plot = T)
    dev.off()
    
    
    # Optionally save the dataVoom object as an RDS file
    saveRDS(dataVoom, basedir(paste0(ct, "_dataVoom.rds")))
    ############
    # Ensure the modelMatrix column names are syntactically valid
    colnames(modelMatrix) <- make.names(colnames(modelMatrix))
    limmaFit <- lmFit(dataVoom, modelMatrix)
    
    # Check which coefficients are not estimable
    non_estimable <- is.na(limmaFit$coefficients)
    non_estimable_coefs <- colnames(modelMatrix)[colSums(non_estimable) > 0]
    
    # Remove non-estimable coefficients from the model matrix
    estimable_coefs <- setdiff(colnames(modelMatrix), non_estimable_coefs)
    modelMatrix_estimable <- modelMatrix[, estimable_coefs, drop = FALSE]
    
    # Re-fit the model without the non-estimable coefficients
    limmaFit_estimable <- lmFit(dataVoom, modelMatrix_estimable)
    limmaFit_estimable <- eBayes(limmaFit_estimable)
    
    # Generate valid contrasts
    genotypes_estimable <- grep("^genotype", colnames(modelMatrix_estimable), value = TRUE)
    
    contrasts_list <- lapply(genotypes_estimable, function(genotype) {
      interaction_term <- paste0("tissueex.vivo.", genotype)
      if (interaction_term %in% colnames(modelMatrix_estimable)) {
        return(paste0(interaction_term, " + ", genotype))
      } else {
        warning(paste("Skipping non-existent interaction term:", interaction_term))
        return(NULL)
      }
    })
    
    contrasts_list <- Filter(Negate(is.null), contrasts_list)
    if (length(contrasts_list) > 0) {
      contrast_matrix <- makeContrasts(contrasts = contrasts_list, levels = modelMatrix_estimable)
    
    # Fit the contrast with the updated model
      limmaFit.contrast <- contrasts.fit(limmaFit_estimable, contrast_matrix)
      limmaFit.contrast <- eBayes(limmaFit.contrast)
    # Fit the contrast with the updated model
        
    
    # Extract and format results for each contrast coefficient
      limmaRes.contrast <- map_dfr(colnames(contrast_matrix), function(contrast_name) {
      topTable(limmaFit.contrast, coef = contrast_name, number = Inf) %>%
        rownames_to_column("ensg") %>%
        mutate(coef = contrast_name,
               group = case_when(
                 logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
                 logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
                 TRUE ~ "n.s"),
               celltype = ct) %>% # Include cell type in the results
        select(-contains("Intercept")) # Optionally remove intercept results if unnecessary
      }, .id = "coefficient")
        
    
    # Extract and format results for each coefficient
      limmaRes <- map_dfr(colnames(limmaFit_estimable$coef), function(coef_name) {
      topTable(limmaFit_estimable, coef = coef_name, number = Inf) %>%
        rownames_to_column("ensg") %>%
        mutate(coef = coef_name, 
               group = case_when(
                 logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
                 logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
                 TRUE ~ "n.s"),
               celltype = ct) %>% # Include cell type in the results
        select(-contains("Intercept")) # Optionally remove intercept results if unnecessary
     }, .id = "coefficient")
    
     
      # Add contrast results to the main results
      limmaRes <- rbind(limmaRes.contrast, limmaRes) 
    } else {
      # If no contrasts are available, just continue with the main results without contrast
      cat("No valid contrasts found, proceeding with the main results only.\n")
    }
    
    # Store results for the current cell type
    de_results_list[[ct]] <- limmaRes
      # Store results for the current cell type
      de_results_list[[ct]] <- limmaRes
      
    }
    
  
  # Combine all results into a single data frame
  final_results_table <- bind_rows(de_results_list)
  return(final_results_table)
}

# Usage example
############################################################

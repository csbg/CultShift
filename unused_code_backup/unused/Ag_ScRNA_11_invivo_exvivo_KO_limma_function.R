##################################################
#**Ag_ScRNA_11_invivo_exvivo_KO_limma_function**#
##################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#06-12-23
#Aarathy
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)


##########################################################
# Create a function for differential expression analysis
##########################################################
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
  
  # Define the model matrix internally based on the provided formula and metadata
  modelMatrix <- model.matrix(model_formula, data = meta)
  
  #A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
  #A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
  #The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs
  # voom
  dataVoom <- voom(d, modelMatrix, plot = T)
  dataVoom%>%write_rds(basedir("dataVoom.rds"))
  limmaFit <- lmFit(dataVoom, modelMatrix)
  limmaFit <- eBayes(limmaFit)
  
  # limma Results
   limmaRes <- list() # start an empty list
  
  limmaRes <- map_dfr(unique(colnames(limmaFit$coef)), ~{
    coefx <- .x
    topTable(limmaFit, coef = coefx, number = Inf) %>%
      rownames_to_column("ensg") %>%
      mutate(coef = coefx)
  }) %>%
    filter(coef != "(Intercept)")
  limmaRes <- limmaRes %>%
    mutate(coef = str_replace(coef, "genotype", "")) %>%
    #mutate(coef = str_replace(coef, "tissueex.vivo:", "interaction_")) %>%
    mutate(coef = str_replace(coef, "tissue", ""))
  #mutate(coef = str_replace(coef, "exvivo:", "interaction_")) 
  #up and downregulated
  limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                             limmaRes$adj.P.Val <= 0.05, "up", 
                           ifelse(limmaRes$logFC <= -1 & 
                                    limmaRes$adj.P.Val <= 0.05, "down", "n.s"))
    
  
  return(limmaRes)
}

# 


####################################################################
#***Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01***#
####################################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#14-03-24
#Aarathy
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
#####################################################################
inDir<-dirout("/Ag_ScRNA_08_Pseudobulk")
base<-"Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01/"
basedir<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01/")
#####################################################################
#source DE function
#####################################################################
source("src/Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_function.R")
################################################
#load data and clean metadata
################################################
#metadata
meta<-read.csv(inDir("metadata.csv"),row.names=1)

#exclude given celltypes
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","Gran.")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]

#only select the genotypes present in both tissue conditions
meta<-meta[meta$genotype %in% meta[meta$tissue=="ex.vivo",]$genotype,]

# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)-]", ".", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta<-meta%>%filter(!grepl("NA",rownames(meta)))
#counts
counts <- read.csv(inDir("ex_counts.csv"), row.names = 1)
counts<-cbind(counts,read.csv(inDir("in_counts.csv"), row.names = 1))
counts<-counts[,rownames(meta)]
stopifnot(all(colnames(counts)==rownames(meta)))
stopifnot(all(unique(meta[meta$tissue=="ex.vivo",]$celltype)==unique(meta[meta$tissue=="in.vivo",]$celltype)))

################################################
#factors and levels
################################################
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))
################################################
#perform DE independent of celltype
model_formula <- "~tissue*genotype"
result <- performDE(meta, counts,model_formula)
limmaRes<-result$limmaRes
dataVoom<-result$dataVoom
limmaRes%>%write_rds(basedir("limma_ex.vivo_vs_in.vivo_all_CT.rds"))
dataVoom%>%write_rds(basedir("dataVoom_ex.vivo_vs_in.vivo_all_CT.rds"))


  
  
  
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
    dataVoom <- voom(counts_subset, modelMatrix, plot = FALSE)  # Disable plot for each cell type
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
      mutate(coef = str_replace(coef, "tissue", ""))
    
    # Up and downregulated genes
    limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                               limmaRes$adj.P.Val <= 0.05, "up", 
                             ifelse(limmaRes$logFC <= -1 & 
                                      limmaRes$adj.P.Val <= 0.05, "down", "n.s"))
    
    # Store results for the current cell type
    de_results_list[[celltype]] <- limmaRes
  }
  
  return(de_results_list)
}


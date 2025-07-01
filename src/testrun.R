
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(purrr)
library(gridExtra)
require(fgsea)
library(msigdbr)
#####################################################################
inDir  <-  dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
base  <-  "Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/"
basedir  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide.R")
########################################################################
limmaRes  <-  read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
limmaRes  <-  limmaRes%>% filter(celltype != "MEP")
################################
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019")################################################################################
#fgsea--------------------------------------------------------------------------
################################################################################


# # # Download gene sets ------------------------------------------------------
enr.terms  <-  enrichrGetGenesets(ENRICHR.DBS)
# # save(enr.terms, file=out("Genesets_Human.RData"))

# Convert to mouse --------------------------------------------------------
hm.map  <-  fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm  <-  unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm)  <-  c("Mouse", "Human")
enr.terms  <-  lapply(enr.terms, function(dbl){
  dbl  <-  lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})

################################################################################
# fgsea -------------------------------------------------------------------
# Initialize the result table
################################################################################

# Initialize the result table
gsea.res <- data.table() 

run_gsea <- function(limmaRes, enr.terms, celltypes = NULL, coefs = NULL) {
  
  # Initialize the result table
  gsea_res <- data.table() 
  
  # Determine cell types to process
  if (is.null(celltypes)) {
    celltypes <- unique(limmaRes$celltype)
  }
  
  # Determine coefficients to process
  if (is.null(coefs)) {
    coefs <- unique(limmaRes$coef)
  }
  
  # Loop through each cell type
  for (ct in celltypes) {
    
    # Loop through each coefficient
    for (de_grp in coefs) {
      
      # Loop through each database in the enrichment terms
      for (dbx in names(enr.terms)) {
        
        # Subset the limma results based on the current cell type and coefficient
        subset_limmaRes <- limmaRes[limmaRes$celltype == ct & limmaRes$coef == de_grp, ]
        
        # Extract statistics (logFC) and assign gene names as names
        stats <- with(subset_limmaRes, setNames(logFC, nm = ensg))
        
        # Skip this iteration if there are missing values in stats
        if (any(is.na(stats))) {
          next
        }
        
        # Perform fgsea analysis
        fgsea_output <- fgsea(
          pathways = enr.terms[[dbx]],
          stats = stats
          #minSize = 15,   # Example additional arguments, adjust as necessary
          #maxSize = 500,  # Example additional arguments, adjust as necessary
          #nperm = 1000    # Example additional arguments, adjust as necessary
        )
        
        # Check if fgsea output is not empty and append the results to gsea_res
        if (length(fgsea_output) > 0) {
          gsea_res <- rbind(gsea_res, data.table(fgsea_output,
                                                 coef = de_grp,
                                                 celltype = ct,
                                                 db = dbx))
        }
      }
    }
  }
  
  # Return the combined GSEA results
  return(gsea_res)
}
gsea.res <- run_gsea(limmaRes, enr.terms, celltypes = unique(limmaRes$celltype), coefs =unique(limmaRes$coef))
gsea.res%>%write_rds(basedir("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))

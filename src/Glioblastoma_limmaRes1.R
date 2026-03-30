##############################
# Modular GEO RNA-seq pipeline
##############################
#Aarathy
###############
source("src/00_init.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(biomaRt)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
library(GEOquery)
library(R.utils)
library(msigdbr)
library(fgsea)
library(readxl)
library(latex2exp)
library(Biobase)
library(Matrix)
source("src/Ag_ScRNA_11_invivo_exvivo_KO_limma_function.R")
InDir1 <- dirout("Ag_top_pathway_genes")

InDir <- dirout("Glioblastoma_NTC")

out <- dirout("Glioblastoma_limmaRes")
meta <- read_rds(InDir("metadata.rds"))
counts <- read_rds(InDir("counts.rds"))
meta$genotype <- gsub("non-targeting","NTC",meta$genotype)
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
unique(meta$condition)
#subsetting
meta$tissue <- ifelse(grepl("invitro", meta$condition), "ex.vivo",
                      ifelse(grepl("preinf", meta$condition), "in.vivo",
                             "CED"))
meta <- meta %>%
  filter(tissue != "CED")
meta$genotype <-   paste(meta$genotype, meta$RT_status, sep = "_")

counts <- counts[, rownames(meta)]
#########
meta$tissue <- factor(meta$tissue)
meta$genotype <- factor(meta$genotype)

meta$genotype <- relevel(meta$genotype, ref = "NTC_noRT")
meta$tissue <- relevel(meta$tissue, ref = "in.vivo")
####################
model_formula <- "~ tissue * genotype"
design <- model.matrix(~ tissue * genotype, data = meta)
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")
keep <- filterByExpr(d0,design)
d0 <- d0[keep,]

# voom transformation
dataVoom <- voom(d0, design, plot=TRUE)

colnames(design) <- make.names(colnames(design))
limmaFit <- lmFit(dataVoom, design)
colnames(design)
non_estimable_coefs <-  nonEstimable(design)

# Remove non-estimable coefficients from the model matrix

estimable_coefs <- setdiff(colnames(design), non_estimable_coefs)
design_estimable <- design[, estimable_coefs, drop = FALSE]
# Re-fit model using estimable coefficients only
limmaFit_estimable <- lmFit(dataVoom, design_estimable)
limmaFit_estimable <- eBayes(limmaFit_estimable)


# Generate valid contrasts
genotypes_estimable <- grep("^genotype", colnames(design_estimable), value = TRUE)

contrasts_list <- lapply(genotypes_estimable, function(genotype) {
  interaction_term <- paste0("tissueex.vivo.", genotype)
  if (interaction_term %in% colnames(design_estimable)) {
    return(paste0(interaction_term, " + ", genotype))
  } else {
    warning(paste("Skipping non-existent interaction term:", interaction_term))
    return(NULL)
  }
})

contrasts_list <- Filter(Negate(is.null), contrasts_list)
if (length(contrasts_list) > 0) {
  contrast_matrix <- makeContrasts(contrasts = contrasts_list, levels = design_estimable)
  
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
               TRUE ~ "n.s")) %>% # Include cell type in the results
      dplyr::select(-contains("Intercept")) # Optionally remove intercept results if unnecessary
  }, .id = "coefficient")
  
  
  # Extract and format results for each coefficient
  limmaRes <- map_dfr(colnames(limmaFit_estimable$coef), function(coef_name) {
    topTable(limmaFit_estimable, coef = coef_name, number = Inf) %>%
      rownames_to_column("ensg") %>%
      mutate(coef = coef_name, 
             group = case_when(
               logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
               logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
               TRUE ~ "n.s")
      ) %>% 
      dplyr::select(-contains("Intercept")) # Optionally remove intercept results if unnecessary
  }, .id = "coefficient")
  
  
  # Add contrast results to the main results
  limmaRes <- rbind(limmaRes.contrast, limmaRes) 
} else {
  # If no contrasts are available, just continue with the main results without contrast
  cat("No valid contrasts found, proceeding with the main results only.\n")
} #


limmaRes <-limmaRes%>%
  filter(coef != "(Intercept)")
unique(limmaRes$coef)
limmaRes <- limmaRes %>%
  mutate(coef = str_replace(coef, "^genotype", "in.vivo")) %>%
  mutate(coef = str_replace(coef, "tissue", ""))

limmaRes <- limmaRes %>%
  mutate(coef = str_replace(coef, "ex.vivo.genotype", "interaction"))
# Modify the coefficients
limmaRes <- limmaRes %>%
  mutate(coef = case_when(
    # Replace "interactionX + genotypeX" with "ex.vivoX"
    str_detect(coef, "interaction") & str_detect(coef, "\\+ genotype") ~ 
      str_replace(coef, "interaction", "ex.vivo") %>% 
      str_replace(" \\+ genotype.*", ""),  # Escape + to treat it as a literal character
    # Keep other coefficients unchanged
    TRUE ~ coef
  ))
limmaRes$RT_status <- ifelse(grepl("noRT", limmaRes$coef), "noRT", "RT")

limmaRes %>%
  write_rds(out("Glioblastoma_limmaRes.rds"))

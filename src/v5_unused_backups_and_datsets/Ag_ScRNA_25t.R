source("src/00_init.R")
source("src/Ag_ScRNA_11_invivo_exvivo_KO_limma_function.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(stringr)
library(data.table)
basedir <- dirout("Ag_ScRNA_25_PerturbNet")
inDir <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")
inDir1 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")



pred <- read.csv(basedir("predicted_expression_SUM_direct.csv"))
pred <- as.data.frame(pred)

# Make rownames
rownames(pred) <- paste(
  pred$genotype,
  pred$cell_type,
  pred$tissue,
  sep = "_"
)

pred <- pred[, -(1:3)]

pred_counts <- t(pred)

#build metadata for pred------
meta_pred <- data.frame(
  sample = colnames(pred_counts)
)

meta_pred$genotype <- sub("_.*", "", meta_pred$sample)
meta_pred$celltype <- sub("^[^_]+_([^_]+)_.*", "\\1", meta_pred$sample)
meta_pred$tissue   <- sub(".*_([^_]+)$", "\\1", meta_pred$sample)

rownames(meta_pred) <- meta_pred$sample
meta_pred$rowname <- rownames(meta_pred)

meta_pred <- meta_pred %>%
  mutate(
    # Check for discrepancies based on rowname and correct celltype
    celltype = case_when(
      grepl("GMP \\(early\\)", rowname) & celltype != "GMP.early" ~ "GMP.early", 
      grepl("GMP \\(late\\)", rowname) & celltype != "GMP.late" ~ "GMP.late",
      grepl("Gran\\. P", rowname) & celltype != "Gran.P" ~ "Gran.P",
      grepl("MEP \\(G1\\)" , rowname) & celltype != "MEP.G1"  ~ "MEP.G1" ,
      grepl("MEP \\(pert\\.\\)" , rowname) & celltype != "MEP.pert."  ~ "MEP.pert." ,
      grepl("MEP \\(S\\)"  , rowname) & celltype != "MEP.S"   ~ "MEP.S" ,
      grepl("MEP \\(early\\)"  , rowname) & celltype != "MEP.early" ~ "MEP.early" ,
      grepl("Imm. B-cell"  , rowname) & celltype == "Imm\\. B-cell"  ~ "Imm.B.cell", 
      TRUE ~ celltype
    )
  )

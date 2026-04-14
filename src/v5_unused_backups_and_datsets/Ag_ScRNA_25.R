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
      grepl("Eo/Ba"  , rowname) & celltype != "Eo.Ba"  ~ "Eo.Ba",
      TRUE ~ celltype
    )
  )


colnames(pred_counts) <- meta_pred$rowname[match(colnames(pred_counts),meta_pred$sample)]
rownames(meta_pred) <- meta_pred$rowname
all(rownames(meta_pred) == colnames(pred_counts))
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell",
                          "Imm.B.cell","MEP","MEP.pert.","MEP.S", "MEP.G1","GMP.late",
                          "GMP.early","Imm. B-cell","nan")

meta_pred <- meta_pred%>%
  filter(!(celltype %in% celltypes_to_exclude))
pred_counts <- pred_counts[,rownames(meta_pred)]
meta_pred$rowname <- NULL
meta_pred$sample <- NULL

unique(meta_pred$celltype)
#meta data
counts <- read.delim(inDir("combined_in_ex.vivo_with_Mye_counts_guide.tsv"), row.names = 1)

common_genes <- intersect(rownames(counts), rownames(pred_counts))


counts_all <- cbind(
  counts[common_genes, ],
  pred_counts[common_genes, ]
)
#observed metadata
meta_obs  <- fread(inDir1("meta_cleaned.tsv"))
meta_obs <- as.data.frame(meta_obs)               # Convert to dataframe (optional)
rownames(meta_obs) <- meta_obs[[1]]   
meta <- meta_obs
meta <- meta[, -1, drop = FALSE]
#ex.vivo from opriginal meta
meta_ex <- meta %>%
  dplyr::select(genotype,tissue,celltype) %>%
  filter(tissue == "ex.vivo")
#in vivo NTC from original meta
meta_in_NTC <- meta %>%
  dplyr::select(genotype,tissue,celltype) %>%
  filter(tissue == "in.vivo", genotype == "NTC")
meta <- rbind(meta_ex,meta_in_NTC)
t <- meta[meta$tissue == "ex.vivo",]
meta_pred <- meta_pred%>%
  filter(genotype %in% unique(t$genotype))
#combine
meta_all <- rbind(meta,meta_pred)

meta_all <- meta_all %>%
  filter(celltype %in% unique(t$celltype))
#filter data
counts_all <- counts_all[rownames(meta_all)]
all(rownames(meta_all) == colnames(counts_all))
###############
#baselines
meta_all$genotype <- factor(meta_all$genotype, levels=c("NTC", 
                                                        unique(setdiff(meta_all$genotype,"NTC"))))
meta_all$tissue <- factor(meta_all$tissue, levels=c("in.vivo", "ex.vivo"))

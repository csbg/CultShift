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
InDir <- dirout("Ag_top_pathway_genes")

out <- dirout("Glioblastoma")
# -----------------------------
# USER PARAMETERS
# -----------------------------

data <-  read_rds("/media/AGFORTELNY/PROJECTS/TfCf_AG/GL261_integrated_20230619.Rds")
# Filter out rows with NA in the sgRNA column
# 1. Identify cells to keep (non-NA sgRNA)
cells_keep <- rownames(data@meta.data)[
  !is.na(data@meta.data$sgRNA) &  # Keep only non-NA sgRNA
    (
      (data@meta.data$source %in% c("CED", "preinf") & data@meta.data$sorted == "MACSFACS") |
        !(data@meta.data$source %in% c("CED", "preinf"))  # keep all other sources
    )
]


# Subset metadata
data@meta.data <- data@meta.data[cells_keep, ]

# Subset RNA assay counts and data
data@assays$RNA@counts <- data@assays$RNA@counts[, cells_keep]
data@assays$RNA@data <- data@assays$RNA@data[, cells_keep]

# If using SCT assay
if ("SCT" %in% names(data@assays)) {
  data@assays$SCT@counts <- data@assays$SCT@counts[, cells_keep]
  data@assays$SCT@data <- data@assays$SCT@data[, cells_keep]
}
table(data@meta.data$orig.ident,data@meta.data$source)

###############
##################################################
#***ScRNA_08_Aggregate_pseudobulk_Ar.R***********#
##################################################

# 1. Create a "condition" column by pasting sgRNA and orig.ident
data@meta.data$condition <- paste(data@meta.data$sgRNA, data@meta.data$orig.ident, sep = "_")
unique(data@meta.data$condition)
# 2. Split cells by condition
cells_by_condition <- split(rownames(data@meta.data), data@meta.data$condition)

# 3. Sum counts per condition
# We'll sum the RNA assay counts (can change to SCT if needed)
count_matrix <- data@assays$RNA@counts

summed_counts <- sapply(cells_by_condition, function(cells) {
  Matrix::rowSums(count_matrix[, cells, drop = FALSE])
})

# Convert to sparse matrix if needed
summed_counts <- as(summed_counts, "dgCMatrix")

summed_counts %>%
  write_rds(out("counts.rds"))
# 4. Create metadata for the summed counts
meta <- data.frame(
  sgRNA = sapply(strsplit(names(cells_by_condition), "_"), `[`, 1),
  orig.ident = sapply(strsplit(names(cells_by_condition), "_"), `[`, 2),
  RT_status = sapply(strsplit(names(cells_by_condition), "_"), function(x) if(any(grepl("RT|noRT", x))) grep("RT|noRT", x, value = TRUE) else NA),
  source = sapply(strsplit(names(cells_by_condition), "_"), function(x) {
    if(any(grepl("48hit", x))) {
      "invitro"
    } else {
      match <- intersect(x, c("CED", "preinf"))
      if(length(match) == 1) match else NA
    }
  }),
  row.names = names(cells_by_condition)
)

# Now you can safely create the condition column
meta$condition <- paste(meta$source, meta$RT_status, sep = "_")
meta$genotype <- meta$sgRNA
meta %>%
  write_rds(out("metadata.rds"))

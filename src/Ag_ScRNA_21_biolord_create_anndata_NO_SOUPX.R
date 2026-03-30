
###############
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
# Load your Monocle object
library(zellkonverter)
library(SingleCellExperiment)
library(monocle3)
library(Seurat)
library(SeuratDisk)
basedir = dirout("Ag_ScRNA_17_proj_celltype_in_monocle_obj_NO_SOUPX/")
# outdir = dirout("/vscratch/aarathy/perturb-net/data/NO_SOUPX")
# Load your Monocle object


# Load Monocle objects
monocle_obj_ex <- readRDS(basedir("ex.vivo_monocle_proj.rds"))
monocle_obj_in <- readRDS(basedir("in.vivo_monocle_proj.rds"))
# Load the Monocle objects
# Ensure the genes match between the two objects
common_genes <- intersect(rownames(monocle_obj_ex), rownames(monocle_obj_in))

# Subset the expression data to match the common genes
monocle_obj_ex <- monocle_obj_ex[common_genes, ]
monocle_obj_in <- monocle_obj_in[common_genes, ]

# Concatenate the expression matrices from both Monocle objects

# Safe extraction
ex_counts <- as(counts(monocle_obj_ex), "dgCMatrix")
in_counts <- as(counts(monocle_obj_in), "dgCMatrix")
common_genes <- intersect(rownames(ex_counts), rownames(in_counts))
ex_counts <- ex_counts[common_genes, ]
in_counts <- in_counts[common_genes, ]

combined_expression_data <- cbind(ex_counts, in_counts)

# Concatenate the cell metadata and add suffixes to differentiate them
colData(monocle_obj_ex)$dataset <- "exvivo"  # Add a dataset column for exvivo
colData(monocle_obj_in)$dataset <- "invivo"  # Add a dataset column for invivo
library(zellkonverter)  # for conversion to/from AnnData
library(SingleCellExperiment)
library(Matrix)

# Extract sparse counts
ex_counts <- as(counts(monocle_obj_ex), "dgCMatrix")
in_counts <- as(counts(monocle_obj_in), "dgCMatrix")

# Keep only common genes
common_genes <- intersect(rownames(ex_counts), rownames(in_counts))
ex_counts <- ex_counts[common_genes, , drop = FALSE]
in_counts <- in_counts[common_genes, , drop = FALSE]

# Combine counts
combined_counts <- cbind(ex_counts, in_counts)

ex_meta <- as.data.frame(colData(monocle_obj_ex))
in_meta <- as.data.frame(colData(monocle_obj_in))

# Add dataset label
ex_meta$dataset <- "exvivo"
in_meta$dataset <- "invivo"

# Combine
combined_cell_data <- rbind(ex_meta, in_meta)
rownames(combined_cell_data) <- colnames(combined_counts)


sce <- SingleCellExperiment(
  assays = list(counts = combined_counts),
  colData = combined_cell_data
)



writeH5AD(sce, "/vscratch/aarathy/perturb-net/data/combined_monocle_NO_SOUPX.h5ad")

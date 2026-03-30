
###############
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
# Load your Monocle object
library(zellkonverter)
library(SingleCellExperiment)
library(monocle3)
library(Seurat)
library(SeuratDisk)
basedir = dirout("Ag_ScRNA_20_proj_celltype_in_monocle_obj/")
outdir = dirout("Ag_ScRNA_21_biolord_create_h5ad")
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
combined_expression_data <- cbind(
  assay(monocle_obj_ex),  # Expression data from exvivo
  assay(monocle_obj_in)   # Expression data from invivo
)

# Concatenate the cell metadata and add suffixes to differentiate them
colData(monocle_obj_ex)$dataset <- "exvivo"  # Add a dataset column for exvivo
colData(monocle_obj_in)$dataset <- "invivo"  # Add a dataset column for invivo

all_cols <- union(colnames(colData(monocle_obj_ex)), colnames(colData(monocle_obj_in)))

ex_meta <- as.data.frame(colData(monocle_obj_ex))
in_meta <- as.data.frame(colData(monocle_obj_in))

# Add missing columns with NA
ex_meta[, setdiff(all_cols, names(ex_meta))] <- NA
in_meta[, setdiff(all_cols, names(in_meta))] <- NA

# Reorder columns to match
ex_meta <- ex_meta[, all_cols]
in_meta <- in_meta[, all_cols]

combined_cell_data <- rbind(ex_meta, in_meta)



# # Combine cell metadata
# combined_cell_data <- rbind(colData(monocle_obj_ex), colData(monocle_obj_in))

# Combine the two objects into a single Monocle object
combined_monocle_obj <- new_cell_data_set(
  expression_data = combined_expression_data,
  cell_metadata = combined_cell_data,
  gene_metadata = rowData(monocle_obj_ex)  # Assuming gene metadata is the same
)


# Convert the combined Monocle object to a Seurat object
seurat_obj <- CreateSeuratObject(
  counts = combined_expression_data, 
  meta.data = combined_cell_data
)
combined_cell_data <- as.data.frame(combined_cell_data)
# Create the Seurat object from the combined data
seurat_obj <- CreateSeuratObject(
  counts = combined_expression_data, 
  meta.data = combined_cell_data
)



# Save Seurat object as .h5Seurat file
SaveH5Seurat(seurat_obj, filename = outdir("combined_monocle_seurat.h5Seurat"))

Convert(outdir("combined_monocle_seurat.h5Seurat"), dest = "h5ad")

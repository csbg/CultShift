library(monocle3)
library(data.table)
library(Matrix)


source("src/00_init.R")
library(Matrix)
library(monocle3)
require("sceasy")
basedir <- "SCRNA_02_01_Integration/"
out <- dirout(paste0(basedir, "/Bet_et_al/"))
out.base <- dirout(basedir)
inDir <- dirout_load("Ag_ScRNA_23_external_dataset3")
dir <- dirout("Ag_ScRNA_23_external_dataset3/")
# Load raw and normalized data
raw_dt <- fread(dir("GSM2388072/GSM2388072_basal_bone_marrow.raw_umifm_counts.csv.gz"))
norm_dt <- fread(dir("GSM2388072/GSM2388072_basal_bone_marrow.filtered_normalized_counts.csv.gz"))

# Extract cell metadata (first 5 columns)
cell_metadata <- raw_dt[, 1:5]
setnames(cell_metadata, c("cell_id", "barcode", "library_id", "run_id", "passed_filter"))

# Set rownames as cell ID
cell_metadata <- as.data.frame(cell_metadata)
rownames(cell_metadata) <- cell_metadata$cell_id

# Extract gene expression matrix (columns 6+)
expr_mat_raw <- as.matrix(raw_dt[, 6:ncol(raw_dt)])
rownames(expr_mat_raw) <- raw_dt[[1]]  # cell_id in column 1
expr_mat_raw <- t(expr_mat_raw)        # gene x cell for monocle3

# Create gene metadata (fake, required)
gene_metadata <- data.frame(
  gene_short_name = rownames(expr_mat_raw),
  row.names = rownames(expr_mat_raw)
)

# Create Monocle3 CellDataSet (CDS)
cds <- new_cell_data_set(
  expr_mat_raw,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)
#pass filter
cds <- cds[, cds@colData$passed_filter == 1]
#Remove cells with zero total counts
cds <- cds[, Matrix::colSums(SingleCellExperiment::counts(cds)) != 0]
# Estimate size factors
cds <- estimate_size_factors(cds)  # 
# Proceed to preprocessing
cds <- preprocess_cds(cds)

#UMAP
cds <- reduce_dimension(cds,preprocess_method = "PCA", verbose = TRUE)
set.seed(42)

cds <- reduce_dimension(cds,
                        reduction_method = "UMAP",
                        
                        verbose = TRUE)


# Clustering
set.seed(12121)
cds = cluster_cells(cds)

monocle.file = out("MonocleObject.RData")
# Ensure the target directory exists
dir.create(dirname(monocle.file), recursive = TRUE, showWarnings = FALSE)

# Save the object
save(cds, file = monocle.file)





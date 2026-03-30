#!/usr/bin/env Rscript

source("src/00_init.R")

library(Seurat)
library(Matrix)
library(data.table)
library(sceasy)

out <- dirout("Ag_SCRNA_01_01_Seurat")

# ------------------------------------------------------------------
# Load sample annotation
# ------------------------------------------------------------------
SANN <- fread(out("SampleAnnotation.tsv"), sep = "\t", header = TRUE)

# Keep relevant tissues only
SANN <- SANN[tissue %in% c("in.vivo", "ex.vivo")]

# ------------------------------------------------------------------
# Helper: read CellRanger filtered matrix
# ------------------------------------------------------------------
read_cellranger <- function(sample) {
  path <- file.path(
    "/media/AGFORTELNY/PROJECTS/TfCf/Data",
    sample,
    "outs",
    "filtered_feature_bc_matrix.h5"
  )
  
  if (!file.exists(path)) {
    stop("Missing filtered matrix for ", sample)
  }
  
  mat <- Read10X_h5(path)
  if (is.list(mat)) mat <- mat[["Gene Expression"]]
  stopifnot(class(mat) == "dgCMatrix")
  mat
}

# ------------------------------------------------------------------
# Create Seurat objects per sample WITHOUT normalization
# ------------------------------------------------------------------
seurat.list <- list()

for (dsx in unique(SANN$sample)) {
  message("Processing: ", dsx)
  
  counts <- read_cellranger(dsx)
  
  sobj <- CreateSeuratObject(
    counts = counts,          # keep raw counts
    project = dsx,
    min.cells = 3,
    min.features = 200
  )
  
  # Add sample-level metadata
  meta <- SANN[sample == dsx]
  for (cn in setdiff(colnames(meta), "sample")) {
    sobj[[cn]] <- meta[[cn]][1]
  }
  
  # QC metrics (percent.mt calculation is fine)
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")
  
  # Conservative filtering
  sobj <- subset(
    sobj,
    subset =
      nFeature_RNA > 300 &
      nCount_RNA > 500 &
      percent.mt < 15
  )
  
  # **No normalization, no scaling, no variable features, no PCA/UMAP**
  # Only keep raw counts as integers
  
  seurat.list[[dsx]] <- sobj
}

stopifnot(length(seurat.list) > 0)

# ------------------------------------------------------------------
# Merge ALL samples into ONE object
# ------------------------------------------------------------------
message("Merging all samples into one Seurat object")

seurat.all <- merge(
  seurat.list[[1]],
  y = seurat.list[-1],
  add.cell.ids = names(seurat.list),
  project = "TfCf_All"
)

# Reduce size for AnnData compatibility, keep only raw counts
seurat.all <- DietSeurat(
  seurat.all,
  assays = "RNA",
  counts = TRUE,
  data = FALSE,         # remove normalized data slot
  scale.data = FALSE,   # remove scaled data
  dimreducs = NULL      # no PCA/UMAP
)

# ------------------------------------------------------------------
# Export to single AnnData file
# ------------------------------------------------------------------
out.file <- "/vscratch/aarathy/perturb-net/data/TfCf_All_rawcounts.h5ad"

sceasy::convertFormat(
  seurat.all,
  from = "seurat",
  to = "anndata",
  outFile = out.file
)

message("Single AnnData file with raw counts written: ", out.file)

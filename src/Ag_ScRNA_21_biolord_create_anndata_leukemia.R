
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
outdir = dirout("Ag_ScRNA_21_biolord_create_h5ad_leukenmia")
# Load your Monocle object

mobjs <- list()
tissue <- "leukemia"
base::load("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_02_01_Integration/leukemia/MonocleObject.RData")
mobjs[["leukemia"]] <- monocle.obj
mobjs$leukemia@colData
ann <- read.delim("/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis/GSE213511_CellAnnotation_leukemia.tsv")
rownames(ann) <- paste(ann$Barcode, ann$Sample,sep = "_")
all(rownames(ann) %in% rownames(mobjs$leukemia@colData))

# Make sure both have unique rownames
stopifnot(!any(duplicated(rownames(ann))))
stopifnot(!any(duplicated(rownames(mobjs$leukemia@colData))))

# Match annotation rows to colData rows
matched_idx <- match(rownames(mobjs$leukemia@colData), rownames(ann))

# Add new column 'celltype' — assign where match is found
mobjs$leukemia@colData$celltype <- NA_character_
mobjs$leukemia@colData$celltype[!is.na(matched_idx)] <- ann$Clusters[matched_idx[!is.na(matched_idx)]]


# Convert Monocle colData to a regular data.frame
meta_df <- as.data.frame(mobjs$leukemia@colData)
expr_matrix <- mobjs$leukemia@assays
################
library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)

# 0) Coerce to SCE
sce <- as(mobjs$leukemia, "SingleCellExperiment")

# 1) Make names unique/non-empty (AnnData requirement)
rownames(sce) <- make.unique(rownames(sce))
colnames(sce) <- make.unique(colnames(sce))

# 2) Ensure sparse matrices (smaller .h5ad)
to_sparse <- function(m) {
  if (is.null(m)) return(NULL)
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  m
}
if ("counts" %in% assayNames(sce))    assay(sce, "counts")    <- to_sparse(assay(sce, "counts"))
if ("logcounts" %in% assayNames(sce)) assay(sce, "logcounts") <- to_sparse(assay(sce, "logcounts"))

# 3) Choose what becomes AnnData's X (change to "counts" if you prefer raw in X)
X_name <- if ("logcounts" %in% assayNames(sce)) "logcounts" else "counts"

# 4) Make sure output dir exists
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 5) Write .h5ad (other assays become layers automatically)
outfile <- file.path(outdir("leukemia.h5ad"))
zellkonverter::writeH5AD(
  sce,
  file = outfile,
  X_name = X_name
)

message("Wrote: ", outfile, "  (X = ", X_name, ")")



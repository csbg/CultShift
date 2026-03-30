# Modified from David Lara et al, Nature Genetics 2023
# Purpose: Projection of ex vivo scRNA-seq onto in vivo reference
# ------------------------------------------------------------------------------

# --- Initialization ------------------------------------------------------------
#package_version("Seurat")
# CODEBASE is where the code is (this code and also other functions)
if(Sys.getenv("CODEBASE") == ""){
  print("Setting CODEBASE")
  Sys.setenv(CODEBASE=paste0(Sys.getenv("HOME"), "/code/"))
}

# GFS points to where the data is stored and results will be written to
if(dir.exists("/media/AGFORTELNY")){
  print("Setting GFS within singularity to /media/AGFORTELNY")
  Sys.setenv(GFS="/media/AGFORTELNY")
}
#source("../../tfcf_AG/src/FUNC_ProjecTILs_PLUS.R")
#source("src/00_init.R")
require(ProjecTILs)
require(umap)
library(readr)
library(Seurat)
library(data.table)
library(readr)
#require(biomaRt)

base.dir <- "Ag_SCRNA_04_01_proj_ex.vivo_D2/"
out <- dirout(base.dir)
InDir <- dirout("Invivo_Invitro_DAVID_new")



# --- Load Monocle objects ------------------------------------------------------
sobjs <- list(
  ex.vivo = NULL,
  in.vivo = NULL
)

paths <- c(
  in.vivo = InDir("Invivo4_5/rna_integrated.rds"),
  ex.vivo = InDir("aarathy_Invitro1_20260310/rna_integrated.rds")
)

for (nm in names(paths)) {
  message("Loading: ", paths[nm])

  sobjs[[nm]] <- read_rds(paths[nm])
}

# singleR cell types ------------------------------------------------------
singleR.cell.types <- readRDS(InDir("Invivo4_5/cell_types_singler.rds"))


#reference
sobjs$in.vivo@meta.data$celltype <- singleR.cell.types$tfcf$labels

# Make sure the assay exists
Assays(sobjs$in.vivo)       # check what assays are present
DefaultAssay(sobjs$in.vivo) <- "RNA"  # or whichever assay exists


# Normalize, find variable features, scale
sobjs$in.vivo <- NormalizeData(sobjs$in.vivo)
sobjs$in.vivo <- FindVariableFeatures(sobjs$in.vivo)
sobjs$in.vivo <- ScaleData(sobjs$in.vivo)

# Run PCA
sobjs$in.vivo <- RunPCA(sobjs$in.vivo, features = VariableFeatures(sobjs$in.vivo))

# Run UMAP
sobjs$in.vivo <- RunUMAP(sobjs$in.vivo, dims = 1:30)

# UMAP coordinates
umap_coords <- sobjs$in.vivo@reductions$umap@cell.embeddings
head(umap_coords)

# Optional: visualize
DimPlot(sobjs$in.vivo, reduction = "umap", group.by = "celltype")

ref <- sobjs$in.vivo
query <- sobjs$ex.vivo
DefaultAssay(ref) <- "RNA"
DefaultAssay(query) <- "RNA"


fix_assay_v5_to_v4 <- function(seurat_obj, assay = "RNA") {
  
  counts <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")
  data   <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "data")
  
  # Step 1: create assay with counts only
  a4 <- SeuratObject::CreateAssayObject(counts = counts)
  
  # Step 2: manually add normalized data
  a4@data <- data
  
  seurat_obj[[assay]] <- a4
  return(seurat_obj)
}

ref   <- fix_assay_v5_to_v4(ref)
query <- fix_assay_v5_to_v4(query)
# Ensure normalized (safe to rerun)
ref   <- NormalizeData(ref)
query <- NormalizeData(query)
query.projected <- Run.ProjecTILs(query, ref = ref)

saveRDS(query.projected, file = out("query_projected.rds"))

query.projected <- Run.ProjecTILs(query, ref = ref)
# Ensure matching cells
matching_cells <- intersect(colnames(ref), singleR.cell.types$cellname)
ref <- subset(ref, cells = matching_cells)

# Extract UMAP embeddings from monocle and align with final ref
ref.umap.original <- reducedDims(ref.monocle)$UMAP
common_cells <- intersect(colnames(ref), rownames(ref.umap.original))
ref <- subset(ref, cells = common_cells)
ref.umap.original <- ref.umap.original[common_cells, , drop = FALSE]
stopifnot(all(colnames(ref) == rownames(ref.umap.original)))

# Find variable features
ref <- FindVariableFeatures(ref)

# PCA
set.seed(1234)
which.assay <- "integrated"
varfeat <- ref@assays[[which.assay]]@var.features
refdata <- data.frame(t(ref@assays[[which.assay]]@data[varfeat,]))
refdata <- refdata[, sort(colnames(refdata))]
ref.pca <- prcomp(refdata, rank. = 50, scale. = TRUE, center = TRUE, retx = TRUE)

# UMAP
seed <- 1234
n.neighbors <- 30
min.dist <- 0.3
metric <- "cosine"
ndim <- 10
umap.config <- umap.defaults
umap.config$n_neighbors <- n.neighbors
umap.config$min_dist <- min.dist
umap.config$metric <- metric
umap.config$n_components <- 2
umap.config$random_state <- seed
umap.config$transform_state <- seed
ref.umap <- umap(ref.pca$x[, 1:ndim], config = umap.config)
colnames(ref.umap$layout) <- c("UMAP_1", "UMAP_2")

# Add PCA and UMAP to object
ref@reductions$UMAP@cell.embeddings <- ref.umap$layout
ref@reductions$PCA@cell.embeddings <- ref.pca$x
ref@reductions$PCA@feature.loadings <- ref.pca$rotation
colnames(ref@reductions$PCA@cell.embeddings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl = TRUE)
colnames(ref@reductions$PCA@feature.loadings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl = TRUE)

# Store PCA and UMAP objects
ref@misc$pca_object <- ref.pca
ref@misc$umap_object <- ref.umap
ref@misc$projecTILs <- "in.vivo"

# Add metadata labels
stopifnot(all(colnames(ref) %in% singleR.cell.types$cellname))
ref <- AddMetaData(ref, as.factor(singleR.cell.types[match(colnames(ref), singleR.cell.types$cellname), ]$labels), col.name = "functional.cluster")

write.tsv(
  merge(
    data.table(ref@meta.data, keep.rownames = "cell_id"),
    data.table(ref@reductions$UMAP@cell.embeddings, keep.rownames = "cell_id"),
    by = "cell_id"
  ),
  out("Output_in.vivo", ".tsv")
)

# Prepare for projection
ref.use <- ref
ref.use@reductions$umap <- ref.use@reductions$UMAP
ref.use@reductions$pca <- ref.use@reductions$PCA

# Projection and prediction loop
for (tx in c("ex.vivo")) {
  query <- as.Seurat.NF(mobjs[[tx]])
  for (sx in unique(query$sample)) {
    query.use <- subset(query, cells = colnames(query)[query$sample == sx])
    query.use@reductions$umap <- query.use@reductions$UMAP
    query.use@reductions$pca <- query.use@reductions$PCA
    
    proj <- make.projection(
      query = query.use,
      ref = ref.use,
      filter.cells = FALSE,
      fast.mode = TRUE,
      seurat.k.filter = 100
    )
    
    pred <- cellstate.predict(ref = ref.use, query = proj)
    
    # Cross-map UMAP with aligned cells
    stopifnot(identical(colnames(ref.use), rownames(ref.umap.original)))
    proj.umap.original <- ref.umap.predict(ref = ref.use, query = proj, ref.umap = ref.umap.original)
    
    # Export results
    meta_dt <- data.table(pred@meta.data, keep.rownames = "cell_id")
    umap_dt <- data.table(proj@reductions$umap@cell.embeddings, keep.rownames = "cell_id")
    
    write.tsv(
      merge(meta_dt, umap_dt, by = "cell_id"),
      out("Output_", tx, "_", sx, ".tsv")
    )
    
    
    meta_dt <- data.table(pred@meta.data, keep.rownames = "cell_id")
    umap_orig_dt <- data.table(proj.umap.original, keep.rownames = "cell_id")
    
    write.tsv(
      merge(meta_dt, umap_orig_dt, by = "cell_id"),
      out("OutputCrossprojection_", tx, "_", sx, ".tsv")
    )
    
  }
}


# Summarize results --------------------------------
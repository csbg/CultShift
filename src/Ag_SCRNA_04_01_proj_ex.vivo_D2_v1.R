# Modified from David Lara et al, Nature Genetics 2023
# Purpose: Projection of ex vivo scRNA-seq onto in vivo reference
# ------------------------------------------------------------------------------

# --- Initialization ------------------------------------------------------------
#package_version("Seurat")
# CODEBASE is where the code is (this code and also other functions)

#source("../../tfcf_AG/src/FUNC_ProjecTILs_PLUS.R")
source("src/00_init.R")
require(ProjecTILs)
require(umap)

library(readr)
#require(biomaRt)
library(Matrix)
library(SeuratObject)
library(Seurat)
options(Seurat.object.assay.version = "v3")
base.dir <- "Ag_SCRNA_04_01_proj_ex.vivo_D2/"
out <- dirout(base.dir)
InDir <- dirout("Invivo_Invitro_DAVID_new")


ex.vivo   <- readRDS(out("counts_ex_vivo.rds"))
ex_meta   <- readRDS(out("ex.vivo_meta.rds"))

in.vivo   <- readRDS(out("counts_in_vivo.rds"))
in_meta   <- readRDS(out("in.vivo_meta.rds"))

singleR.cell.types <- readRDS(InDir("Invivo4_5/cell_types_singler.rds"))


ex.vivo <- as(ex.vivo, "dgCMatrix")
in.vivo <- as(in.vivo, "dgCMatrix")
str(ex.vivo)
# Ensure metadata alignment
ex_meta <- ex_meta[colnames(ex.vivo), ]
in_meta <- in_meta[colnames(in.vivo), ]
library(SeuratObject)

# Create assay
assay <- CreateAssayObject(counts = ex.vivo)

# Create identities
idents <- factor(rep("unknown", ncol(ex.vivo)))
names(idents) <- colnames(ex.vivo)

# Build Seurat object properly
query <- new(
  Class = "Seurat",
  assays = list(RNA = assay),
  meta.data = ex_meta,
  active.assay = "RNA",
  active.ident = idents,
  graphs = list(),
  reductions = list(),
  images = list(),
  project.name = "query"
)

# Set identities explicitly
Idents(query) <- idents

# --- Load seurat objects ------------------------------------------------------
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

ex.vivo <- readRDS(out("counts_ex_vivo.rds"))
ex.vivo <- as(ex.vivo, "dgCMatrix")
ex_meta <- readRDS(out("ex.vivo_meta.rds"))
dim(ex.vivo)
head(colnames(ex.vivo))
head(rownames(ex.vivo))

dim(ex_meta)
head(rownames(ex_meta))

query <- CreateSeuratObject(
  counts = ex.vivo,
  meta.data = ex_meta
)

#referenceout
sobjs$in.vivo@meta.data$celltype <- singleR.cell.types$tfcf$labels


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
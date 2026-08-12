# ============================================================
# Project ex vivo cells onto in vivo reference cell types
# ============================================================

# --- Initialization ------------------------------------------------------------

source("src/00_init.R")

library(Seurat)
library(SeuratObject)
library(ProjecTILs)
library(umap)
library(readr)
library(data.table)
library(Matrix)
library(ggplot2)
library(dplyr)

base.dir <- "Ag_SCRNA_04_01_proj_ex.vivo_D2/"
out <- dirout(base.dir)

InDir <- dirout("Invivo_Invitro_DAVID_new")


# ============================================================================
# 1. Load Seurat objects
# ============================================================================

paths <- c(
  in.vivo  = InDir("Invivo4_5/rna_integrated.rds"),
  ex.vivo  = InDir("aarathy_Invitro1_20260310/rna_integrated.rds")
)

sobjs <- lapply(paths, function(path) {
  message("Loading: ", path)
  readRDS(path)
})


# ============================================================================
# 2. Load in vivo SingleR annotations
# ============================================================================

singleR.cell.types <- readRDS(
  InDir("Invivo4_5/cell_types_singler.rds")
)

# Add cell types to in vivo object
sobjs$in.vivo$celltype <- singleR.cell.types$tfcf$labels


# ============================================================================
# 3. Prepare in vivo reference
# ============================================================================

ref <- sobjs$in.vivo

# Extract RNA counts
ref_counts <- ref[["RNA"]]@layers$counts
ref_counts <- Matrix(ref_counts, sparse = TRUE)

rownames(ref_counts) <- rownames(ref)
colnames(ref_counts) <- colnames(ref)

# Rebuild clean Seurat object
ref.clean <- CreateSeuratObject(
  counts = ref_counts,
  meta.data = ref@meta.data
)

# Standard RNA processing
ref.clean <- NormalizeData(ref.clean)
ref.clean <- FindVariableFeatures(ref.clean)
ref.clean <- ScaleData(ref.clean)
ref.clean <- RunPCA(ref.clean)
ref.clean <- RunUMAP(
  ref.clean,
  dims = 1:30
)


# ============================================================================
# 4. Create ProjecTILs reference
# ============================================================================

ref.use <- make.reference(
  ref = ref.clean,
  annotation.column = "celltype",
  recalculate.umap = TRUE
)


# ============================================================================
# 5. Prepare ex vivo query
# ============================================================================

query <- sobjs$ex.vivo

# Extract RNA counts
query_counts <- query[["RNA"]]@layers$counts
query_counts <- Matrix(query_counts, sparse = TRUE)

rownames(query_counts) <- rownames(query)
colnames(query_counts) <- colnames(query)

# Rebuild clean query object
query.clean <- CreateSeuratObject(
  counts = query_counts,
  meta.data = query@meta.data
)

# Normalize query
query.clean <- NormalizeData(query.clean)

DefaultAssay(query.clean) <- "RNA"


# ============================================================================
# 6. Project ex vivo cells onto in vivo reference
# ============================================================================

query.projected <- Run.ProjecTILs(
  query.clean,
  ref = ref.use,
  query.assay = "RNA",
  filter.cell = FALSE,
  skip.normalize = TRUE
)


# ============================================================================
# 7. Add predicted cell types back to ex vivo object
# ============================================================================

# Check that cell order is identical
stopifnot(
  all(colnames(sobjs$ex.vivo) == colnames(query.projected))
)

# Add projected annotation
sobjs$ex.vivo$celltype <-
  query.projected$functional.cluster

sobjs$ex.vivo$functional.cluster.conf <-
  query.projected$functional.cluster.conf


# ============================================================================
# 8. Visualize projection
# ============================================================================

# Query projected UMAP
DimPlot(
  query.projected,
  group.by = "functional.cluster",
  reduction = "umap"
) +
  ggtitle("Ex vivo cells projected onto in vivo reference")


# Reference UMAP
DimPlot(
  ref.use,
  group.by = "celltype",
  reduction = "umap"
) +
  ggtitle("In vivo reference")


# ============================================================================
# 9. Combined reference + query UMAP
# ============================================================================

# Reference
ref_umap <- Embeddings(
  ref.use,
  "umap"
)

ref_df <- as.data.frame(ref_umap)

ref_df$dataset <- "Reference"
ref_df$celltype <- ref.use$celltype


# Query
query_umap <- Embeddings(
  query.projected,
  "umap"
)

query_df <- as.data.frame(query_umap)

query_df$dataset <- "Query"
query_df$celltype <-
  query.projected$functional.cluster


# Plot
ggplot() +
  
  geom_point(
    data = ref_df,
    aes(
      x = UMAP_1,
      y = UMAP_2
    ),
    color = "lightgrey",
    alpha = 0.5,
    size = 0.3
  ) +
  
  geom_point(
    data = query_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = celltype
    ),
    size = 0.5,
    alpha = 0.8
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Ex vivo cells projected onto in vivo reference",
    color = "Predicted cell type"
  )


# ============================================================================
# 10. Save annotated objects
# ============================================================================

saveRDS(
  sobjs$ex.vivo,
  InDir("aarathy_Invitro1_20260310/rna_integrated_annotated_proj.rds")
)

saveRDS(
  sobjs$in.vivo,
  InDir("Invivo4_5/rna_integrated_annotated_invivo.rds")
)


# ============================================================================
# 11. Final checks
# ============================================================================

message(
  "Ex vivo cells: ",
  ncol(sobjs$ex.vivo)
)

message(
  "In vivo cells: ",
  ncol(sobjs$in.vivo)
)

message(
  "Ex vivo cells with predicted cell type: ",
  sum(!is.na(sobjs$ex.vivo$celltype))
)

message(
  "Ex vivo cells without predicted cell type: ",
  sum(is.na(sobjs$ex.vivo$celltype))
)

table(
  sobjs$ex.vivo$celltype,
  useNA = "ifany"
)
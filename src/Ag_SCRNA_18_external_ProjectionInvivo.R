#modifie_from_David Lara et al 2023 nature genetics
source("src/00_init.R")

require(ProjecTILs)
require(umap)
require(biomaRt)
basedir <- "Ag_SCRNA_18_external_ProjectionInvivo/"
out <- dirout(basedir)
source("src/FUNC_ProjecTILs_PLUS.R")
InDir <- dirout("Ag_SCRNA_02_01_Integration/")

# Annotation --------------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)

# Read in vivo data and perform differnetial expression -------------------
mobjs <- list()

# List of file paths and names
paths <- c("ex.vivo_with_Mye" = InDir("ex.vivo_with_Mye/soupx/MonocleObject.RData"),
           
           "in.vivo"    = InDir("in.vivo/soupx/MonocleObject.RData"))

base::load(paths["ex.vivo_with_Mye"])
mobjs[["ex.vivo_with_Mye"]] <- monocle.obj
unique(mobjs$ex.vivo_with_Mye@colData$functional.cluster)
base::load(paths["in.vivo"])
mobjs[["in.vivo"]] <- monocle.obj

# singleR cell types ------------------------------------------------------
singleR.cell.types <- readRDS(dirout_load("SCRNA_06_02_MergeMarkers")("CellTypes_in.vivo.RDS"))


# Function to transform monocle3 to Seurat objects ----------------------------------
x <- mobjs$in.vivo
as.Seurat.NF <- function(x){
  logcounts(x) <- counts(x)
  x <- as.Seurat(x)
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e6)
  x <- RenameAssays(x, RNA = "integrated")
  x
}

# Prepare reference -------------------------------------------------------

# the following code and code at the end of this section enables loading of previously stored data, this is dangerous if data gets updated
# ref.file <- out("reference.rds")
# if(file.exists(ref.file)){
#   ref <- readRDS(ref.file)
# } else {

ref.monocle <- mobjs$in.vivo
ref <- as.Seurat.NF(ref.monocle)
ref.umap.original <- reducedDims(ref.monocle)$UMAP
ref <- FindVariableFeatures(ref)

# PCA
set.seed(1234)
which.assay="integrated"
varfeat <- ref@assays[[which.assay]]@var.features
refdata <- data.frame(t(ref@assays[[which.assay]]@data[varfeat,]))
refdata <- refdata[, sort(colnames(refdata))]
ref.pca <- prcomp(refdata, rank. = 50, scale. = TRUE, center = TRUE, retx=TRUE)

# UMAP
seed=1234
n.neighbors=30
min.dist=0.3
metric="cosine"
ndim=10
umap.config <- umap.defaults
umap.config$n_neighbors = n.neighbors
umap.config$min_dist = min.dist
umap.config$metric = metric
umap.config$n_components = 2
umap.config$random_state = seed
umap.config$transform_state = seed
ref.umap <- umap(ref.pca$x[,1:ndim], config=umap.config)
colnames(ref.umap$layout) <- c("UMAP_1","UMAP_2")

# add to object
ref@reductions$UMAP@cell.embeddings <- ref.umap$layout
ref@reductions$PCA@cell.embeddings <- ref.pca$x
ref@reductions$PCA@feature.loadings <- ref.pca$rotation
colnames(ref@reductions$PCA@cell.embeddings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl=TRUE)
colnames(ref@reductions$PCA@feature.loadings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl=TRUE)
#Store the complete PCA and UMAP object in @misc
ref@misc$pca_object <- ref.pca
ref@misc$umap_object <- ref.umap
ref@misc$projecTILs="in.vivo"

# Add labels
#stopifnot(all(colnames(ref) %in% singleR.cell.types$cellname))
colnames(ref)[!colnames(ref) %in% singleR.cell.types$cellname]
ref <- AddMetaData(ref, as.factor(singleR.cell.types[match(colnames(ref), cellname),]$labels), col.name = "functional.cluster")

# Export table
write.tsv(
  merge(
    data.table(ref@meta.data, keep.rownames = TRUE),
    data.table(ref@reductions$UMAP@cell.embeddings, keep.rownames = TRUE),
    by="rn"
  ), out("Output_in.vivo",".tsv"))
    # Save
 # saveRDS(ref, ref.file)
# }


# Run projection and cell type prediction ----------------------------------------------------------
ref.use <- ref
ref.use@reductions$umap <- ref.use@reductions$UMAP
ref.use@reductions$pca <- ref.use@reductions$PCA

tx <- "Anna_et_al"

mobjs[[tx]]
for (tx in c( "Anna_et_al")){
  query <- as.Seurat.NF(mobjs[[tx]])
  query@meta.data
  query@reductions$umap <- query@reductions$UMAP
  query@reductions$pca <- query@reductions$PCA
#rownames to gene names

  # Check if gene names are Ensembl IDs in query and gene symbols in reference
  is_ensembl <- grepl("^ENSMUSG", rownames(query)[1])
  is_symbol  <- !grepl("^ENSMUSG", rownames(ref.use)[1])
  
  if (is_ensembl && is_symbol) {
    message("Mapping Ensembl IDs to gene symbols for ", tx)
    
    # Remove Ensembl version numbers (e.g., ENSMUSG00000000001.1 → ENSMUSG00000000001)
    ensembl_ids <- gsub("\\..*", "", rownames(query))
    
    # Query biomaRt
    mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    gene.map <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "mgi_symbol"),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = mart
    )
    
    # Filter gene map: remove empty symbols and duplicates
    gene.map <- gene.map[gene.map$mgi_symbol != "", ]
    gene.map <- gene.map[!duplicated(gene.map$ensembl_gene_id), ]
    
    # Map Ensembl to gene symbols
    symbol_map <- setNames(gene.map$mgi_symbol, gene.map$ensembl_gene_id)
    ens_ids <- gsub("\\..*", "", rownames(query[["integrated"]]))
    new_symbols <- symbol_map[ens_ids]
    keep <- !is.na(new_symbols)
    
    # Subset query to mapped genes only
    query <- subset(query, features = rownames(query[["integrated"]])[keep])
    new_symbols <- unname(new_symbols[keep])
    
    # Rename assay matrices
    DefaultAssay(query) <- "integrated"
    assay <- query[["integrated"]]
    
    assay@counts <- assay@counts[keep, , drop = FALSE]
    assay@data   <- assay@data[keep, , drop = FALSE]
    if (!is.null(assay@scale.data) && nrow(assay@scale.data) > 0) {
      assay@scale.data <- assay@scale.data[keep, , drop = FALSE]
    }
    
    # Rename rows to gene symbols
    rownames(assay@counts)     <- new_symbols
    rownames(assay@data)       <- new_symbols
    if (!is.null(assay@scale.data) && nrow(assay@scale.data) > 0) {
      rownames(assay@scale.data) <- new_symbols
    }
    
    query[["integrated"]] <- assay
    message("Renaming completed for ", tx)
    
  } else {
    message("Gene format already matched — no renaming applied for ", tx)
  }
  
    
  
  
  # Harmonize features
  common.genes <- intersect(rownames(query), rownames(ref.use))
  query.sub <- subset(query, features = common.genes)
  ref.sub <- subset(ref.use, features = common.genes)
  
  # Make projection
  
  proj <- make.projection(
    query = query.sub,
    ref = ref.sub,
    filter.cells = FALSE,
    fast.mode = FALSE,
    seurat.k.filter = 100
  )
  
  # Prediction
  pred <- cellstate.predict(ref=ref.use, query=proj)
  
  # Cross-map UMAP
  proj.umap.original <- ref.umap.predict(ref=ref.use, 
                                         query=proj, ref.umap = ref.umap.original)
  message("Projection complete for ", tx)
  # Export results
  # Normal
  write.tsv(
    merge(
      data.table(pred@meta.data, keep.rownames = TRUE),
      data.table(proj@reductions$umap@cell.embeddings, keep.rownames = TRUE),
      by="rn"
    ), out("Output_", tx, ".tsv"))
  
  # Cross-projected to original UMAP
  write.tsv(
    merge(
      data.table(pred@meta.data, keep.rownames = TRUE),
      data.table(proj.umap.original, keep.rownames = TRUE),
      by="rn"
    ), out("OutputCrossprojection_", tx, ".tsv"))
}

#Update annotation


# Human/Mouse gene mapping ------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
SANN <- fread(PATHS$SCRNA$ANN)


######################################################
# Load Anna dataset -------------------------------------------------------
dataset <-"Anna_et_al"
for (dataset in c("Anna_et_al", "Bet_et_al")){
  annot <- fread(out("OutputCrossprojection_", dataset, ".tsv"))
  head(annot)
  path <- dirout(paste0("Ag_SCRNA_02_01_Integration/",dataset))
  (base::load(path("MonocleObject.RData")))
  # Step 1: Add rownames of colData to a new column
  cds_coldata <- as.data.frame(colData(cds))
  cds_coldata$rn <- rownames(cds_coldata)
  
  # Step 2: Merge with annot using "rn" as the key
  merged_coldata <- merge(cds_coldata, annot[, .(rn, functional.cluster)], by = "rn", all.x = TRUE)
  
  # Step 3: Set rownames again (important!)
  rownames(merged_coldata) <- merged_coldata$rn
  
  # Step 4: Remove the temporary "rn" column if you want (optional)
    # Step 5: Assign back to cds
  colData(cds) <- S4Vectors::DataFrame(merged_coldata)
  monocle.file <- path("MonocleObject.RData")
  # Store full dataset
  save(cds, file=monocle.file)
  
}
cds@colData
# # Summarize results -------------------------------------------------------
# ff <- list.files(out(""), pattern="Output_")
# ff <- lapply(ff, function(fx) fread(out(fx)))
# 
# ff2 <- list.files(out(""), pattern="OutputCrossprojection_")
# ff2 <- lapply(ff2, function(fx) fread(out(fx)))
# 
# pDT.list <- list(
#   original=rbindlist(ff, fill=TRUE),
#   crossproject=rbindlist(ff2, fill=TRUE)
# )
# pred@meta.data
# 
# for(xx in names(pDT.list)){
#   pDT <- pDT.list[[xx]]
#   
#   # Hex plot
#   ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) + 
#     theme_bw(12) +
#     geom_hex(data=pDT[tissue == tx], bins=100) + 
#     geom_density_2d(data=pDT[tissue != "in.vivo"], aes(color=tissue))
#   ggsave(out(xx, "_Hexplot.pdf"), w=6,h=5)
#   
#   # Hexplot by tissue
#   ggplot(pDT, aes(x=UMAP_1, y=UMAP_2)) + 
#     theme_bw(12) +
#     stat_binhex(aes(fill=log10(..count..)), bins=100) + 
#     facet_grid(. ~ tissue)
#   ggsave(out(xx, "_Hexplot.byTissue.pdf"), w=16,h=5)
#   
#   # Celltypes
#   ggplot(pDT, aes(x=UMAP_1, y=UMAP_2, color=functional.cluster)) + 
#     theme_bw(12) +
#     scale_color_brewer(palette = "Paired") +
#     geom_point(shape=1, alpha=0.3)
#   ggsave(out(xx, "_Celltypes.png"), w=7,h=5)
#   
#   # Numbers of cell types
#   ggplot(pDT, aes(x=functional.cluster)) + 
#     theme_bw(12) +
#     geom_bar() + 
#     facet_grid(. ~ tissue) +
#     xRot() +
#     scale_y_log10()
#   ggsave(out(xx, "_Celltypes.numbers.pdf"), w=9,h=4)
#   
#   
#   # Compare to SingleR ------------------------------------------------------
#   xDT <- merge(pDT, singleR.cell.types, by.x="rn", by.y="cellname")
#   
#   if(nrow(xDT) == 0) next
#   
#   jDT <- data.table()
#   for(tx in unique(xDT$tissue)){
#     jMT <- jaccard.twolists(
#       l1=with(xDT[tissue == tx], split(rn, functional.cluster)),
#       l2=with(xDT[tissue == tx], split(rn, labels))
#     )
#     
#     jDT <- rbind(jDT, data.table(melt(data.table(jMT, keep.rownames = TRUE), id.vars = "rn"),
#                                  tissue=tx))
#   }
#   ggplot(jDT, aes(x=rn, y=variable, fill=value)) +
#     theme_bw(12) +
#     geom_tile() +
#     facet_grid(. ~ tissue) +
#     xRot() +
#     scale_fill_gradient(low="white", high="blue")
#   ggsave(out(xx, "_ComparisonToSingleR.pdf"), w=13,h=5)
# }
# 

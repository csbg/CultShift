
#############
source("src/00_init.R")
library("dplyr")
library("stringr")
library("tidyverse")
library("monocle3")
require(ProjecTILs)
require(umap)
source("src/00_init.R")
base.dir <- "SCRNA_08_01_ProjectionInvivo/"
out <- dirout(base.dir)

source("src/FUNC_ProjecTILs_PLUS.R")

#############
#paths------
#############
out <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")
inDir1<- dirout_load("/SCRNA_10_collect_UMAPs")
inDir <- dirout_load("Ag_SCRNA_02_01_Integration/Anna_et_al")

# Human/Mouse gene mapping ------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
SANN <- fread(PATHS$SCRNA$ANN)


######################################################
# Load Anna dataset -------------------------------------------------------
(base::load(inDir("MonocleObject.RData")))
cds <- NULL
monocle.obj <- cds
######################################################


# Annotation --------------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)

# Read in vivo data and perform differnetial expression -------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


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
stopifnot(all(colnames(ref) %in% singleR.cell.types$cellname))
ref <- AddMetaData(ref, as.factor(singleR.cell.types[match(colnames(ref), cellname),]$labels), col.name = "functional.cluster")

# Export table
write.tsv(
  merge(
    data.table(ref@meta.data, keep.rownames = TRUE),
    data.table(ref@reductions$UMAP@cell.embeddings, keep.rownames = TRUE),
    by="rn"
  ), out("Output_in.vivo",".tsv"))

#   # Save
#   saveRDS(ref, ref.file)
# }


monocle.obj <-
  preprocess_cds(cds, verbose = TRUE) %>%
  reduce_dimension(preprocess_method = "PCA", verbose = TRUE)
set.seed(42)

#why? what alignment group for Anna? 
# monocle.obj <- align_cds(monocle.obj,
# alignment_group = "sample_broad",
# residual_model_formula_str = "~Phase",
# verbose = TRUE)

# The alignment_group parameter is used to specify the grouping variable based on 
# which the alignment will be performed. This variable typically represents the 
# batches or groups that you want to align.
# 
# For example, if you have scRNA-seq data from different experimental conditions 
# (e.g., control vs. treatment), you might want to align cells based on these 
# conditions. In this case, you would specify "condition" as the alignment_group.
monocle.obj <- reduce_dimension(monocle.obj,
                                reduction_method = "UMAP",
                                #what for Anna?
                                #preprocess_method = "Aligned",
                                verbose = TRUE)
# Clustering
set.seed(12121)
monocle.obj = cluster_cells(monocle.obj)

####################
#process_monocle_data <- function(monocle.obj, annotations, tissue){
  # Extract annotations
annot <- Anna.ann[rn %in% colnames(monocle.obj)][, c("rn", "celltype","sample")]


rownames(monocle.obj@colData)
unique(annot$celltype)
#Replace Eo. and Ba. with Eo.Ba
annot$celltype<-gsub("Eo|Ba", "Eo.Ba", annot$celltype)
#annot$celltype<-gsub("IMP[12]","GMP",annot$celltype)

# Add cell type to colData
monocle.obj@colData["rn"] <- rownames(monocle.obj@colData)
monocle.obj <- monocle.obj[, annot$rn]
stopifnot(all(annot$rn == colnames(monocle.obj)))
monocle.obj@colData["celltype"] <- annot$celltype
monocle.obj@colData["sample"]<-annot$sample
#metadata 
data<-as.data.frame(monocle.obj@colData)
head(data)
celltype_samp <- data %>%
    group_by(sample,celltype,) %>%
    summarize(num_celltypes_per_sample = n_distinct(rn))
data<-celltype_samp
# Plotting using ggplot
ggplot(data, aes(x = celltype, y = num_celltypes_per_sample, fill = sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Cell Type", y = "Number of Cell Types per Sample", fill = "Sample") +
    ggtitle("Number of Cell Types per Sample for Each Cell Type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(out("celltype_abundance_Anna.png"))
# Create genotype column
# monocle.obj@colData["genotype"] <- gsub("_.+$", "", monocle.obj@colData$guide)
  
  # Split metadata based on celltype, and sample
  # also split based on sample because when doing DE analysis, otherwise you only get one sample per condition.
mat_sam_ct <- with(monocle.obj@colData, split(rn, paste(celltype, sample)))
sort(sapply(mat_sam_ct, length))
# sapply here applies the function rowsums to the given list of vectors-mat_gt_ct_ex. 
#Each of the elements in mat_gt_ct_ex is a vector specifying the sample names in that group. 
# Using sapply, the function rowsums is applied to the counts matrix of 
#the exvivo with those selected columns.
  
# Apply function rowSums to the given list of vectors
result <- sapply(mat_sam_ct, function(bcs) Matrix::rowSums(counts(monocle.obj)[, bcs, drop = FALSE]))
  
#########
#metadata
#########
meta <- data.frame(
    cell = gsub("\\s+$", "", names(mat_sam_ct)),
    #genotype = "NA",
    sample = "NA",
    celltype = "NA",
    tissue = "Anna"
  )
rownames(meta)<-meta$cell
# This loop iterates over each unique cell type in monocle.obj@colData$celltype, 
# then for each one, it assigns that cell type to meta$celltype wherever 
#meta$cell has in the name, the current cell type ct

for (sam in unique(monocle.obj@colData$sample)) {
  meta$sample[grepl(sam, meta$cell)] <- sam
  }
  
for (ct in unique(monocle.obj@colData$celltype)) {
    meta$celltype[grepl(ct, meta$cell)] <- ct
  }
  
df_counts_Anna<-as.data.frame(result)
meta_Anna<- meta
colnames(meta_Anna)
# Merge and write combined counts

write.table(df_counts_Anna,out("Anna_counts.tsv"), row.names = TRUE,col.names = T)


# Combine metadata

write.table(meta_Anna, file = out("Anna_metadata.tsv"), row.names = F,col.names = T)


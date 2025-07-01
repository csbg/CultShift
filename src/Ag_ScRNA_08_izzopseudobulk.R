
#############
source("src/00_init.R")
library("dplyr")
library("stringr")
library("tidyverse")
#############
#paths------
#############
out <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")
inDir1<- dirout_load("/SCRNA_10_collect_UMAPs")
inDir <- dirout_load("SCRNA_02_01_Integration")

# Human/Mouse gene mapping ------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
SANN <- fread(PATHS$SCRNA$ANN)


######################################################
# Load izzo dataset -------------------------------------------------------
(base::load(dirout_load("SCRNA_05_01_SingleR")("izzo.RData")))

izzoCDS <- new_cell_data_set(expression_data = izzoMT, 
                             cell_metadata = data.frame(row.names=colnames(izzoMT),
                                                        tissue=rep("Izzo", ncol(izzoMT))))

head(izzo.ann)
#izzo.ann <- fread("metadata/IzzoEtAl.metadata.csv")
izzo.ann <- izzo.ann[grepl("^WT\\d$", orig.ident),]
izzo.ann[, V1 := gsub("WT6", "WT4", V1)]
izzo.ann[, rn := V1]
izzo.ann[, orig.ident := gsub("WT6", "WT4", orig.ident)]
izzo.ann[,clusterName := gsub("-[0-9]+", "", clusterName)]
stopifnot(length(setdiff(colnames(izzoMT), izzo.ann$V1))==0)
izzoMT <- izzoMT[,izzo.ann$V1]
head(izzoMT)
head(izzo.ann)
stopifnot(all(colnames(izzoMT) == izzo.ann$V1))
izzo.ann$sample <- izzo.ann$orig.ident
izzo.ann$celltype <- izzo.ann$clusterName
izzo.ann$rn <- izzo.ann$V1
izzo.ann <- izzo.ann[,c("rn","sample","celltype")]
head(izzo.ann)




monocle.obj <-
  preprocess_cds(izzoCDS, verbose = TRUE) %>%
  reduce_dimension(preprocess_method = "PCA", verbose = TRUE)
set.seed(42)

#why? what alignment group for izzo? 
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
                                #what for izzo?
                                #preprocess_method = "Aligned",
                                verbose = TRUE)
# Clustering
set.seed(12121)
monocle.obj = cluster_cells(monocle.obj)

####################
#process_monocle_data <- function(monocle.obj, annotations, tissue){
  # Extract annotations
annot <- izzo.ann[rn %in% colnames(monocle.obj)][, c("rn", "celltype","sample")]


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
ggsave(out("celltype_abundance_izzo.png"))
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
    tissue = "izzo"
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
  
df_counts_izzo<-as.data.frame(result)
meta_izzo<- meta
colnames(meta_izzo)
# Merge and write combined counts

write.table(df_counts_izzo,out("izzo_counts.tsv"), row.names = TRUE,col.names = T)


# Combine metadata

write.table(meta_izzo, file = out("izzo_metadata.tsv"), row.names = F,col.names = T)


##################################################
#***ScRNA_08_Aggregate_pseudobulk_Ar.R***********#
##################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#06-12-23
#Aarathy
#############
#setup------
#############
source("src/00_init.R")
library("dplyr")
library("stringr")
library("tidyverse")
#############
#paths------
#############
out <- dirout("/AG_01_test")
inDir1<- dirout_load("/SCRNA_10_collect_UMAPs")
inDir <- dirout_load("SCRNA_02_01_Integration")

#############
#Functions----
#############
process_monocle_data <- function(monocle_obj, annotations, tissue){
  # Extract annotations
  annot <- annotations[rn %in% colnames(monocle_obj)][, c("rn", "functional.cluster")]
  colnames(annot) <- c("rn", "celltype")
  
  # Add cell type to colData
  monocle_obj@colData["rn"] <- rownames(monocle_obj@colData)
  monocle_obj <- monocle_obj[, annot$rn]
  stopifnot(all(annot$rn == colnames(monocle_obj)))
  monocle_obj@colData["celltype"] <- annot$celltype
  
  # Create genotype column
  monocle_obj@colData["genotype"] <- gsub("_.+$", "", monocle_obj@colData$guide)
  
  # Split metadata based on celltype, genotype, and sample
  # also split based on sample because when doing DE analysis, otherwise you only get one sample per condition.
  mat_gt_ct <- with(monocle_obj@colData, split(rn, paste(genotype, celltype, orig.ident)))
  sort(sapply(mat_gt_ct, length))
  # sapply here applies the function rowsums to the given list of vectors-mat_gt_ct_ex. 
  #Each of the elements in mat_gt_ct_ex is a vector specifying the sample names in that group. 
  # Using sapply, the function rowsums is applied to the counts matrix of 
  #the exvivo with those selected columns.
  
  # Apply function rowSums to the given list of vectors
  result <- sapply(mat_gt_ct, function(bcs) Matrix::rowSums(counts(monocle_obj)[, bcs, drop = FALSE]))
  
  #########
  #metadata
  #########
  meta <- data.frame(
    cell = gsub("\\s+$", "", names(mat_gt_ct)),
    genotype = "NA",
    sample = "NA",
    celltype = "NA",
    tissue = tissue
  )
  
  for (gene in unique(monocle_obj@colData$genotype)) {
    meta$genotype[grepl(gene, meta$cell)] <- gene
  }
  
  for (sam in unique(monocle_obj@colData$sample)) {
    meta$sample[grepl(sam, meta$cell)] <- sam
  }
  
  for (ct in unique(monocle_obj@colData$celltype)) {
    meta$celltype[grepl(ct, meta$cell)] <- ct
  }
  
  return(list(result = result, meta = meta))
}

# load Monocle Objects
mobjs <- list()
tissue<-c("ex.vivo","in.vivo")
for(tissuex in tissue){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}

tissue_n<- names(mobjs)
tissues<-mobjs[tissue_n]

#annotations from 
annotations<- readRDS(inDir1("ProjMonocle_celltypes.RDS"))

# Process ex.vivo data
Output_ex <- process_monocle_data(tissues$ex.vivo, annotations,"ex.vivo")
df_counts_ex<-Output_ex$result
meta_ex <- Output_ex$meta

# Process in.vivo data
Output_in <- process_monocle_data(tissues$in.vivo, annotations,"in.vivo")
df_counts_in <-Output_in$result
meta_in <- Output_in$meta

# Merge and write combined counts
df_counts <- merge(df_counts_ex, df_counts_in, by = "row.names")
write.tsv(df_counts,out("combined_in_ex_counts.tsv"))

# Combine metadata
combined_meta <- rbind(meta_ex, meta_in)
colnames(combined_meta) <- c("cell", "genotype", "sample", "celltype", "tissue")
rownames(combined_meta)<-combined_meta$cell
write.table(combined_meta, file = out("metadata.tsv"), sep = "\t", row.names = FALSE)

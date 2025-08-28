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
out <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide/")
#inDir1<- dirout_load("/SCRNA_10_collect_UMAPs")
InDir <- dirout_load("Ag_SCRNA_02_01_Integration/")
InDir1 <- dirout("Ag_SCRNA_05_01_UMAPs_and_celltypes/")
#############
#Functions----
#############
#monocle_obj<-tissues$in.vivo
process_monocle_data <- function(monocle_obj, annotations, tissue){
  # Extract annotations
  #monocle_obj<-tissues$ex.vivo
  annot <- annotations[rn %in% colnames(monocle_obj)][, c("rn", "functional.cluster")]
  colnames(annot) <- c("rn", "celltype")
  
  # Add cell type to colData
  monocle_obj@colData["rn"] <- rownames(monocle_obj@colData)
  monocle_obj <- monocle_obj[, annot$rn]
  stopifnot(all(annot$rn == colnames(monocle_obj)))
  monocle_obj@colData["celltype"] <- annot$celltype
  
  # Create genotype column
  monocle_obj@colData["genotype"] <- gsub("_.+$", "", monocle_obj@colData$guide)
  monocle_obj@colData$guide
  # Split metadata based on celltype, genotype-guide, and sample
  # also split based on sample because when doing DE analysis, otherwise you only get one sample per condition.
  # Split metadata based on celltype, genotype, and sample
  # also split based on sample because when doing DE analysis, otherwise you only get one sample per condition.
  mat_gt_ct <- with(monocle_obj@colData, split(rn, paste(guide, celltype, sample)))
  
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
    tissue = tissue,
    guide ="NA",
    mixscape_global ="NA"
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
  for (gd in unique(monocle.obj@colData$guide)){
    meta$guide[grepl(gd, meta$cell)]<-gd
  }
  for (mix in unique(monocle.obj@colData$mixscape_class.global)){
    meta$mixscape_global[grepl(mix, meta$cell)]<-mix
  }
  
  return(list(result = result, meta = meta))
}
#################
#
##############

# load Monocle Objects
mobjs <- list()
tissue <-c("ex.vivo","in.vivo")#,"leukemia")
for(tissuex in tissue){
  (base::load(InDir(paste0(tissuex,"/soupx/MonocleObject.RData"))))
  mobjs[[tissuex]] <- monocle.obj
}

tissue_n <- names(mobjs)
tissues <- mobjs[tissue_n]

#annotations from invivo projections



annotations <- readRDS(InDir1("ProjVivo_celltypes.RDS"))
unique(annotations$sample)
tissues$ex.vivo@colData
# Process ex.vivo data
Output_ex <- process_monocle_data(tissues$ex.vivo, annotations,"ex.vivo")
df_counts_ex<-Output_ex$result
meta_ex <- Output_ex$meta

# Process in.vivo data
Output_in <- process_monocle_data(tissues$in.vivo, annotations,"in.vivo")
df_counts_in <- Output_in$result
meta_in <- Output_in$meta
unique(meta_in$celltype)
unique(meta_ex$celltype)
#############################

#For ex.vivo_with_Mye and in.vivo
df_counts <- merge(df_counts_ex, df_counts_in, by = "row.names")
write.tsv(df_counts,out("combined_in_ex.vivo_with_Mye_counts_guide.tsv"))

# Combine metadata
combined_meta <- rbind(meta_ex, meta_in)
colnames(combined_meta) <- c("cell", "genotype", "sample", "celltype", "tissue","guide","mixscape_global")
rownames(combined_meta) <- combined_meta$cell
meta <- combined_meta
meta$rowname <- rownames(meta)



# Correct the celltype
meta <- meta %>%
  mutate(
    # Check for discrepancies based on rowname and correct celltype
    celltype = case_when(
      grepl("GMP \\(early\\)", rowname) & celltype != "GMP.early" ~ "GMP.early", 
      grepl("GMP \\(late\\)", rowname) & celltype != "GMP.late" ~ "GMP.late",
      grepl("Gran\\. P", rowname) & celltype != "Gran.P" ~ "Gran.P",
      grepl("MEP \\(G1\\)" , rowname) & celltype != "MEP.G1"  ~ "MEP.G1" ,
      grepl("MEP \\(pert\\.\\)" , rowname) & celltype != "MEP.pert."  ~ "MEP.pert." ,
      grepl("MEP \\(S\\)"  , rowname) & celltype != "MEP.S"   ~ "MEP.S" ,
      grepl("MEP \\(early\\)"  , rowname) & celltype != "MEP.early" ~ "MEP.early" ,
      grepl("Imm. B-cell"  , rowname) & celltype == "Imm.B.cell"  ~ "Imm.B.cell", 
      TRUE ~ celltype
    )
  )

write.table(meta, file = out("metadata_guide_with Mye.tsv"), sep = "\t", row.names = F)


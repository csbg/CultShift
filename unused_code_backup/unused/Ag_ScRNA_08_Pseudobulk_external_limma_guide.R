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
library(biomaRt)
#############
#paths------
#############
out <- dirout("/Ag_ScRNA_08_Pseudobulk_external_limma_guide")
inDir1<- dirout_load("/SCRNA_10_collect_UMAPs")


#############
#Functions----
#############
#monocle_obj<-tissues$in.vivo
process_monocle_data <- function(monocle_obj, tissue){
  # Extract annotations
  #monocle_obj<-tissues$ex.vivo
  if (!"sample" %in% colnames(colData(monocle_obj))) {
    colData(monocle_obj)$sample <- colData(monocle_obj)$library_id
  }
  
  colData(monocle_obj)$celltype <- colData(monocle_obj)$functional.cluster
  # Create genotype column
  
  # Split metadata based on celltype, genotype-guide, and sample
  # also split based on sample because when doing DE analysis, otherwise you only get one sample per condition.
  # Split metadata based on celltype, genotype, and sample
  # also split based on sample because when doing DE analysis, otherwise you only get one sample per condition.
  mat_gt_ct <- with(monocle_obj@colData, split(rn, paste(celltype, sample)))
  
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
    genotype = "WT",
    sample = "NA",
    celltype = "NA",
    tissue = tissue
   
  )
  
  for (sam in unique(monocle_obj@colData$sample)) {
    meta$sample[grepl(sam, meta$cell)] <- sam
  }
  
  for (ct in unique(monocle_obj@colData$celltype)) {
    meta$celltype[grepl(ct, meta$cell)] <- ct
  }
  
  
  return(list(result = result, meta = meta))
}
#################
#
##############

# load Monocle Objects
mobjs <- list()

dataset <-"Anna_et_al"
for (dataset in c("Anna_et_al", "Bet_et_al")){
  
  path <- dirout(paste0("Ag_SCRNA_02_01_Integration/",dataset))
  (base::load(path("MonocleObject.RData")))
  
  mobjs[[dataset]] <- cds
}

tissue_n <- names(mobjs)
tissues <- mobjs[tissue_n]


# Process ex.vivo data
Output_Anna <- process_monocle_data(tissues$Anna_et_al, tissue = "Anna_et_al")
df_counts_Anna <- Output_Anna$result
meta_Anna <- Output_Anna$meta

# Process Bet_et_al data
Output_Bet <- process_monocle_data(tissues$Bet_et_al, tissue = "Bet_et_al")
df_counts_Bet <- Output_Bet$result
meta_Bet <- Output_Bet$meta
unique(meta_Bet$celltype)
head(df_counts_Anna)
#############################
# Step 1: Convert counts matrix to data.frame and preserve Ensembl IDs (without version)

# Connect to Ensembl BioMart (mouse dataset)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Extract Ensembl IDs without version numbers
ensembl_ids <- sub("\\..*", "", rownames(df_counts_Anna))

# Get gene symbols from Ensembl
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)
head(df_counts_Anna)
# Match and rename rownames in df_counts_Anna
# Add a new column to preserve rownames
df_counts_Anna <- as.data.frame(df_counts_Anna)
df_counts_Anna$ensembl_gene_id <- sub("\\..*", "", rownames(df_counts_Anna))
head(df_counts_Anna)
df_counts_Anna <- merge(gene_map, df_counts_Anna, by = "ensembl_gene_id")
df_counts_Anna <- df_counts_Anna[!duplicated(df_counts_Anna$external_gene_name), ]


rownames(df_counts_Anna) <- df_counts_Anna$external_gene_name

# Drop helper columns
df_counts_Anna$ensembl_gene_id <- NULL
df_counts_Anna$external_gene_name <- NULL

# Now df_merged has gene symbols as rownames
#For Anna_et_al and Bet_et_al
df_counts_Anna$gene <- rownames(df_counts_Anna)

df_counts_Bet <- as.data.frame(df_counts_Bet)
df_counts_Bet$gene <- rownames(df_counts_Bet)
df_counts <- merge(df_counts_Anna, df_counts_Bet, by = "gene")
rownames(df_counts) <- df_counts$gene




# Combine metadata
combined_meta <- rbind(meta_Anna, meta_Bet)
colnames(combined_meta) <- c("cell",  "genotype","sample", "celltype", "tissue")
rownames(combined_meta) <- combined_meta$cell
meta <- combined_meta


meta <- meta %>%
  rownames_to_column("rowname") %>%  # Convert row names to a column
  mutate(
    celltype = case_when(
      grepl("GMP \\(early\\)", rowname) & celltype != "GMP.early" ~ "GMP.early", 
      grepl("GMP \\(late\\)", rowname) & celltype != "GMP.late" ~ "GMP.late",
      grepl("Gran\\. P", rowname) & celltype != "Gran.P" ~ "Gran.P",
      grepl("MEP \\(G1\\)", rowname) & celltype != "MEP.G1" ~ "MEP.G1",
      grepl("MEP \\(pert\\.\\)", rowname) & celltype != "MEP.pert." ~ "MEP.pert.",
      grepl("MEP \\(S\\)", rowname) & celltype != "MEP.S" ~ "MEP.S",
      grepl("MEP \\(early\\)", rowname) & celltype != "MEP.early" ~ "MEP.early",
      grepl("Imm. B-cell", rowname) & celltype != "Imm.B.cell" ~ "Imm.B.cell",
      TRUE ~ celltype
    )
  )


clean_names <- function(x) {
  x <- gsub("/$", "", x)        # Remove trailing slash
  x <- trimws(x)                # Remove leading/trailing spaces
  x <- gsub(" ", ".", x)        # Replace spaces with dots
  return(x)
}

colnames(df_counts) <- clean_names(colnames(df_counts))
rownames(meta) <- clean_names(rownames(meta))
clean_names <- function(x) {
  x <- gsub("/$", "", x)        # Remove trailing slash
  x <- trimws(x)                # Remove leading/trailing spaces
  x <- gsub(" ", ".", x)        # Replace spaces with dots
  x <- gsub("\\.+", ".", x)     # Replace multiple dots with single dot
  return(x)
}


colnames(df_counts) <- clean_names(colnames(df_counts))
rownames(meta) <- clean_names(rownames(meta))
#write files
write.table(meta, file = out("metadata_guide_external_data.tsv"), sep = "\t", row.names = T)
df_counts$gene <- NULL
write.table(df_counts,out("combined_Anna_Bet_counts_guide.tsv"),row.names = T)



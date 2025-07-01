source("src/00_init.R")
library(monocle3)
InDir <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")
for (tissue in c("exvivo","invivo")){
  monocle_obj <- read_rds(InDir(paste0(tissue, "_monocle_proj.rds")))
  annotations <- monocle_obj@colData
  annotations <- as.data.frame(annotations)
  write_tsv(annotations,InDir(paste0(tissue,"_","Annotations.tsv")))
}



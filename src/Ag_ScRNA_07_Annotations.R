source("src/00_init.R")
library(monocle3)
InDir <- dirout("Ag_SCRNA_02_01_Integration/")
out <- dirout("/Ag_ScRNA_07_Annotations/")
tissue <-c("ex.vivo_with_Mye","in.vivo")#,"leukemia")
for(tissuex in tissue){
  (base::load(InDir(paste0("in.vivo","/soupx/MonocleObject.RData"))))
  annotations <- monocle.obj@colData
  annotations <- as.data.frame(annotations)
  annotations$rn <- rownames(annotations)
  write.tsv(annotations,out(paste0(tissuex,"_","Annotations.tsv")))
}

t<- read.delim(out("ex.vivo_with_Mye_Annotations.tsv"))

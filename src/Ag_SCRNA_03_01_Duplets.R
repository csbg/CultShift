source("src/00_init.R")
library(purrr)
library(data.table)
require("scds")
require(doMC)
basedir <- "Ag_SCRNA_03_01_Duplets/"
out <- dirout(basedir)
inDir <- dirout_load("Ag_SCRNA_01_01_Seurat")


registerDoMC(cores=8)

# Annotation --------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)
SANN <- SANN[tissue != "leukemia",]
# Integrate data with Monocle
sx <- SANN$sample[1]


walk(SANN$sample, function(sx) {
  fx <- inDir("SeuratObj_", sx, ".RData")
  if (!file.exists(fx)) stop(fx, " seurat object not found")
  
  message("Reading ", sx)
  
  base::load(fx)
  
  scds <- seurat.obj %>%
    as.SingleCellExperiment() %>%
    bcds()
  
  res <- data.table(
    rn = colnames(scds),
    dupletScore = colData(scds)$bcds_score
  )
  
  write.tsv(res, out("Duplet_Scores_", sx, ".tsv"))
})


# Plot
ff <- list.files(out(""), pattern="Duplet_Scores_.*.tsv", full.names = TRUE)
names(ff) <- gsub("^Duplet_Scores_(.*).tsv$", "\\1", basename(ff))
ff <- lapply(ff, fread)
pDT <- rbindlist(ff, idcol = "sample")
ggplot(pDT, aes(x=sample, y=dupletScore)) + geom_violin() + theme_bw(12) + xRot()
ggsave(out("plot_Violin.pdf"), w=20,h=5)

ggplot(pDT, aes(x=dupletScore, group=sample)) + stat_ecdf(alpha=0.25) + theme_bw(12) + xRot()
ggsave(out("plot_ECDF.jpg"), w=5,h=5)


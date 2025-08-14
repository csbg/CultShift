require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)

# renv::snapshot(lockfile = "renv_NF.lock")

source("~/code/resources/RFunctions/Basics.R")
out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Normalized_reads"
################################################################################
dirout_jak <- function(out, ext="", init=TRUE){
  out.dir <- paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/", "/JAKSTAT/", out, "/")
  if(init){
    dir.create(out.dir,showWarnings=FALSE); 
    message("Setting output directory: ", out.dir)
  }
  function(...){
    paste0(out.dir, paste0(...), ext)
  }
}
################################################################################
# Settings (from JAK-STAT paper) ------------------------------------------
NORMALIZE <- TRUE
LMM.REML <- FALSE
LMM.VoomWeights <- FALSE
BPPARAM=BiocParallel::SnowParam(10, "SOCK", progressbar=TRUE)
CUTOFFS <- list(
  DE.ATAC = list(
    model="Lme4",
    AveExpr = 0,
    adj.P.Val = 0.05
  ),
  DE.RNA = list(
    model="Dream",
    AveExpr = 0,
    adj.P.Val = 0.05
  ),
  GSEA.REPEATS=10000
)

# Load data ---------------------------------------------------------------

# download and extract data
if(!file.exists(file.path(out, "Data", "RNA_09_CleanData", "RNA_CountData.RData"))){
  download.file("https://medical-epigenomics.org/papers/jakstat2024/Data.zip", destfile = file.path(out, "Data.zip"))
  unzip(file.path(out, "Data.zip"), exdir = out)
}

# load in vivo data at homeostasis
load(file.path(out, "Data", "RNA_09_CleanData", "RNA_CountData.RData"))
counts.in <- RNA.counts
gmap.in <- RNA.gmap
ann.in <- read_tsv(file.path(out, "Data", "RNA_09_CleanData", "Annotation.tsv"))

# load ex vivo data and subset to untreated
load(file.path(out, "Data", "TREAT_RNA_09_CleanData", "RNA_CountData.RData"))
counts.ex <- RNA.counts
gmap.ex <- RNA.gmap
ann.ex <- read_tsv(file.path(out, "Data", "TREAT_RNA_09_CleanData", "Annotation.tsv"))
ann.ex <- ann.ex %>%
  filter(treatment == "UT")
counts.ex <- counts.ex[,ann.ex$sample_name]


# combine data ------------------------------------------------------------

# combine matrices
stopifnot(all(row.names(counts.ex) == row.names(counts.in)))
gmap <- gmap.in
counts <- cbind(counts.ex, counts.in)

# combine annotations
colnames(ann.in)
col_names = c("sample_name", "experiment_id", "genotype", "cell_type", "treatment")
ann <- rbind(
  ann.in %>%select(any_of(col_names)),
  ann.ex %>% select(any_of(col_names))
) %>% 
  mutate(treatment = ifelse(treatment == "H", "in_vivo", "ex_vivo"))
ann <- data.table(ann)

#write_rds(counts, file = file.path(out, "DEG_Counts.RDS"))
#write_rds(gmap, file = file.path(out, "DEG_GMP.RDS"))
#write_rds(ann, file = file.path(out, "DEG_Annotation.RDS"))


# DEG ---------------------------------------------------------------------

Normalized_reads <- list.files(out)

# List all CSV files in the directory
Normalized_files <- list.files(path = out, pattern = ".csv$", full.names = TRUE)
test<- read_csv("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Normalized_reads/Normalized_readsT8.csv")
# Read the CSV files and store them in a named list
Normalized_reads <- lapply(Normalized_files, read_csv)

# Get the base file names (without directory path and extension)
file_names <- basename(Normalized_files)
file_names <- sub("\\.csv$", "", file_names)

# Assign the file names as the names of the list elements
names(Normalized_reads) <- file_names

# Print the list to check
print(Normalized_reads)



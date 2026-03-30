source("src/00_init.R")

library(monocle3)
library(tidyverse)
out <- dirout("/Ag_ScRNA_21_biolord_create_mtx")

InDir <- dirout("Ag_ScRNA_20_proj_celltype_in_monocle_obj")
# Ex vivo -----------------------------------------------------------------

data_exvivo <- readRDS(InDir("ex.vivo_monocle_proj.rds"))

# colnames(data_exvivo) %>% intersect(colnames(monocle.obj))
# rownames(data_exvivo) %>% setdiff(rownames(monocle.obj))
# Extract counts (sparse dgCMatrix)
mat <- counts(data_exvivo)

head(mat)
# Write to Matrix Market format
writeMM(mat, out("exvivo_monocle_proj.mtx"))


# In vivo -----------------------------------------------------------------

data_invivo <- read_rds(InDir("in.vivo_monocle_proj.rds"))


# colnames(data_invivo) %>% setdiff(colnames(monocle.obj))
# rownames(data_invivo) %>% setdiff(rownames(monocle.obj))
mat <- counts(data_invivo)

# Write to Matrix Market format
writeMM(mat, out("invivo_monocle_proj.mtx"))

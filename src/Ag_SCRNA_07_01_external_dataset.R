# -------------------------------------------------------------------
# Initial Setup
# -------------------------------------------------------------------
source("src/00_init.R")

basedir <- "Ag_SCRNA_07_01_external_dataset/"
out <- dirout(basedir)

library(Seurat)
library(SingleR)
library(celldex)
library(tidyverse)
library(CytoTRACE)
library(doMC)
library(GEOquery)
library(data.table)
library(Matrix)

# Human/Mouse gene mapping and annotation
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = TRUE)
SANN <- fread(PATHS$SCRNA$ANN)
inDir <- dirout_load("Ag_SCRNA_01_01_Seurat")

# -------------------------------------------------------------------
# Prepare Directories
# -------------------------------------------------------------------
out <- function(...) file.path("data", ...)
out
dir.create(out(), recursive = TRUE, showWarnings = FALSE)

Anna.file <- out("Anna.RData")

# -------------------------------------------------------------------
# Load or Download and Process Anna et al Data
# -------------------------------------------------------------------
if (file.exists(Anna.file)) {
  load(Anna.file)
  
} else {
  # Step 1: Download data
  geo_accs <- c("GSM5344484", "GSM5344485")
  for (geo_acc in geo_accs) {
    getGEOSuppFiles(geo_acc, fetch_files = TRUE, baseDir = out(""))
  }
  
  # Step 2: Untar downloaded files
  ff <- list.files(out(""), pattern = "\\.tar\\.gz$", recursive = TRUE, full.names = TRUE)
  untar_dir <- out("untarred")
  dir.create(untar_dir, showWarnings = FALSE)
  
  for (fx in ff) {
    message("Extracting: ", fx)
    untar(fx, exdir = untar_dir)
  }
  
  # Step 3: List all untarred files (for reference/log)
  untarred_files <- list.files(untar_dir, recursive = TRUE, full.names = TRUE)
  cat("Extracted files:\n")
  print(untarred_files)
  
  # Step 4: Load 10X data using Seurat
  sample_dirs <- list(
    WT1 = file.path(untar_dir, "Sca1pos_scRNAseq"),
    WT2 = file.path(untar_dir, "Sca1neg_scRNAseq")
  )
  
  mats <- list()
  for (sample in names(sample_dirs)) {
    message("Reading 10X data for: ", sample)
    mat_list <- Read10X(data.dir = sample_dirs[[sample]])
    
    # Handle multi-modal 10X outputs
    if ("Gene Expression" %in% names(mat_list)) {
      mat <- mat_list[["Gene Expression"]]
    } else {
      mat <- mat_list[[1]]
    }
    
    colnames(mat) <- paste0(sample, "_", colnames(mat))
    mats[[sample]] <- mat
  }
  
  # Step 5: Harmonize gene sets
  all_genes <- unique(unlist(lapply(mats, rownames)))
  mats <- lapply(mats, function(m) {
    missing_genes <- setdiff(all_genes, rownames(m))
    if (length(missing_genes) > 0) {
      padding <- Matrix(0, nrow = length(missing_genes), ncol = ncol(m), sparse = TRUE)
      rownames(padding) <- missing_genes
      m <- rbind(m, padding)
    }
    m[all_genes, , drop = FALSE]
  })
  
  # Step 6: Combine matrices
  AnnaMT <- do.call(cbind, mats)
  
  # Step 7: Load metadata
 # Anna.ann <- fread("")
 # Anna.ann <- Anna.ann[grepl("^WT\\d$", orig.ident), ]
 # Anna.ann[, V1 := gsub("WT6", "WT4", V1)]  # standardize barcode prefix if needed
  
  # Step 8: Match annotation and expression
 # stopifnot(length(setdiff(colnames(AnnaMT), Anna.ann$V1)) == 0)
 # AnnaMT <- AnnaMT[, Anna.ann$V1]
 # stopifnot(all(colnames(AnnaMT) == Anna.ann$V1))
  
  # Step 9: Save for reuse
  save(AnnaMT, file = Anna.file)
}


#GSM2388072
# -------------------------------------------------------------------
# Initial Setup
# -------------------------------------------------------------------
source("src/00_init.R")

basedir <- "Ag_ScRNA_23_external_dataset3/"
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
inDir <- dirout_load("SCRNA_01_01_Seurat")

# -------------------------------------------------------------------
# Prepare Directories
# -------------------------------------------------------------------
out <- function(...) file.path("data", ...)

dir.create(out(), recursive = TRUE, showWarnings = FALSE)

Bet.file <- out("Bet.RData")

# -------------------------------------------------------------------
# Load or Download and Process Bet et al Data
# -------------------------------------------------------------------
if (file.exists(Bet.file)) {
  load(Bet.file)
  
} else {
  # Step 1: Download data
  geo_accs <- c("GSM2388072")
  for (geo_acc in geo_accs) {
    getGEOSuppFiles(geo_acc, fetch_files = TRUE, baseDir = out(""))
  }
  
  # Step 2: Untar downloaded files
  ff <- list.files(out(""))#, pattern = "\\.tar\\.gz$", recursive = TRUE, full.names = TRUE)
    
  # Step 3: List all untarred files (for reference/log)
  files <- list.files(out("GSM2388072"), recursive = TRUE, full.names = TRUE)
  cat("Extracted files:\n")
  print(files)
  
  raw_path <- files[2]
  norm_path <- files[1]
  
  # Load raw UMI matrix (assumed gene-by-cell)
  raw_counts <- fread(raw_path)
  norm_counts <- fread(norm_path)
  
  # Inspect
  dim(raw_counts)
  dim(norm_counts)
  



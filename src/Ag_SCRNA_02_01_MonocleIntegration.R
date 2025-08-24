source("src/00_init.R")
basedir <- "Ag_SCRNA_02_01_Integration/"
out.base <- dirout(basedir)
inDir <- dirout_load("SCRNA_01_01_Seurat")
require("sceasy")



# Annotation --------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)
samples.not.use <- SANN[grepl("\\d$", tissue) | toAnalyse == "NO"]$sample
SANN <- SANN[!sample %in% samples.not.use]

# Filter only ex.vivo samples
SANN_ex <- SANN[tissue == "ex.vivo"]

# Iterate over unique tissues in ex.vivo (if multiple)
for(tx in unique(SANN_ex$tissue)){
  out <- dirout(paste0(basedir, "/ex.vivo_with_Mye", "/soupx/"))
  monocle.file <- out("MonocleObject.RData")
  
  monocle.obj.list <- list()
  
  for(sx in SANN_ex[tissue == tx]$sample){
    fx <- inDir("SeuratObj_", sx, ".RData")
    print(paste("Reading", sx))
    
    if(!file.exists(fx)) stop(fx, " seurat object not found")
    base::load(fx)
    
    # Add annotations
    for(x in c("tissue", "markers", "timepoint", "sample", "sample_broad")){
      seurat.obj@meta.data[[x]] <- SANN_ex[sample == sx][[x]]
    }
    
    # Expression data
    mat.use <- seurat.obj@assays$RNA@counts
    stopifnot(!any(duplicated(row.names(mat.use))))
    if(!"GFP" %in% row.names(mat.use)){
      x <- matrix(0, nrow = 2, ncol = ncol(mat.use))
      row.names(x) <- c("GFP", "BFP")
      mat.use <- rbind(mat.use, x)
    }
    
    monocle.obj.list[[sx]] <- new_cell_data_set(
      expression_data = mat.use,
      cell_metadata = seurat.obj@meta.data
    )
    
    # Store CITE-seq data
    if(sx == "DM_CITEseq-2_NA_NM_1"){
      citeseq.MT <- additional.info.x
      save(citeseq.MT, file=out.base("CITESEQ_Antibodies.RData"))
    }
  }
  
  # Check row counts
  stopifnot(length(unique(sapply(monocle.obj.list, nrow))) == 1)
  
  # Combine and process
  monocle.obj <- combine_cds(cds_list = monocle.obj.list, cell_names_unique = FALSE)
  monocle.obj <- preprocess_cds(monocle.obj, verbose = TRUE) %>%
    reduce_dimension(preprocess_method = "PCA", verbose = TRUE)
  
  set.seed(42)
  monocle.obj <- align_cds(monocle.obj,
                           alignment_group = "sample_broad",
                           verbose = TRUE)
  monocle.obj <- reduce_dimension(monocle.obj,
                                  reduction_method = "UMAP",
                                  preprocess_method = "Aligned",
                                  verbose = TRUE)
  
  set.seed(12121)
  
  monocle.obj <- cluster_cells(monocle.obj)
  # Create directory if missing
  dir.create(dirname(monocle.file), recursive = TRUE, showWarnings = FALSE)
  
  # Save object
  dir.create(dirname(monocle.file), recursive = TRUE, showWarnings = FALSE)
  
  save(monocle.obj, file=monocle.file)
}
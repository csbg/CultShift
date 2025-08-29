source("src/00_init.R")
library(monocle3)
out <- dirout("/Ag_ScRNA_17_proj_celltype_in_monocle_obj")
inDir1<- dirout_load("/SCRNA_10_collect_UMAPs")
inDir <- dirout_load("SCRNA_02_01_Integration/soupx/")

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
  monocle_obj@colData["celltype_projection"] <- annot$celltype
  
  # Create genotype column
  monocle_obj@colData["genotype"] <- gsub("_.+$", "", monocle_obj@colData$guide)
  
  # Split metadata based on celltype, genotype-guide, and sample
  # also split based on sample because when doing DE analysis, otherwise you only get one sample per condition.
  # Split metadata based on celltype, genotype, and sample
  # also split based on sample because when doing DE analysis, otherwise you only get one sample per condition.
    
  #monocle_obj %>% write_rds(paste0(tissue,"_monocle_proj.rds"))
  
  
  # Save the processed Monocle object
  saveRDS(monocle_obj, file = out(paste0(tissue, "_monocle_proj.rds")))
  return(monocle_obj)  # Return the processed object
}

#############
# Load Annotations ----
#############
annotations <- readRDS(inDir1("ProjVivo_celltypes.RDS"))

#############
# Process Monocle Objects ----
#############
tissues <- c("ex.vivo", "in.vivo", "leukemia")

##################
mobjs <- lapply(tissues, function(tissuex) {
  tryCatch({
    # Construct file path
    file_path <- PATHS$SCRNA$MONOCLE.DIR(paste0(tissuex, "/soupx/"))
    message("Checking file path for tissue: ", tissuex, " -> ", file_path)
    
    # Check if file exists
    if (!file.exists(file_path)) {
      stop("File not found for tissue: ", tissuex, " -> ", file_path)
    }
    
    # Load the file
    base::load(file_path)
    message("Loaded file for tissue: ", tissuex)
    
    # Check if monocle.obj exists
    if (!exists("monocle.obj")) {
      stop("Monocle object not found in loaded file for tissue: ", tissuex)
    }
    
    # Process and save the Monocle object
    process_monocle_data(monocle.obj, annotations, tissuex)
  }, error = function(e) {
    message("Error processing tissue: ", tissuex, "\n", e)
    NULL  # Return NULL in case of error
  })
})

# Filter out unsuccessful tissues and assign names
mobjs <- Filter(Negate(is.null), mobjs)
names(mobjs) <- tissues[seq_along(mobjs)]

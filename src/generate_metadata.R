source("src/00_init.R")
basedir <- dirout("generate_metadata")
InDir <- dirout("Ag_SCRNA_02_01_Integration/")
InDir1 <- dirout("Ag_SCRNA_05_01_UMAPs_and_celltypes")
InDir2
Celltype_Annotations <- read_rds(InDir1("ProjVivo_celltypes.RDS"))


# --- Load Annotations ----------------------------------------------------------
SANN <- fread("/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis//Ag_SCRNA_01_01_Seurat/SampleAnnotation.tsv")

# --- Load Monocle objects ------------------------------------------------------
#must not use the functional.cluster already present in monocle.object
#use the celltypes obtained from projection to in vivo
mobjs <- list(
  ex.vivo = NULL,
  in.vivo = NULL
)

paths <- c(
  ex.vivo = InDir("ex.vivo/soupx/MonocleObject.RData"),
  in.vivo = InDir("in.vivo/soupx/MonocleObject.RData")
)

for (nm in names(paths)) {
  message("Loading: ", paths[nm])
  base::load(paths[nm])
  mobjs[[nm]] <- monocle.obj
}


# Extract metadata from Monocle objects and convert to data.frame
metadata_list <- lapply(mobjs, function(obj) {
  if (!is.null(obj)) {
    meta <- as.data.frame(pData(obj))  # convert DFrame to data.frame
    meta$rn <- rownames(meta)          # add rownames for merging
    return(meta)
  } else {
    return(NULL)
  }
})

# Remove NULL entries
metadata_list <- metadata_list[!sapply(metadata_list, is.null)]

# Combine metadata from both objects into one data.frame
combined_metadata <- bind_rows(metadata_list, .id = "tissue") 
combined_metadata <- combined_metadata%>%
  dplyr::select(-functional.cluster)
# 'tissue' column will indicate ex.vivo or in.vivo

# Make sure Celltype_Annotations has 'rn' column
if (!"rn" %in% colnames(Celltype_Annotations)) {
  Celltype_Annotations <- as.data.frame(Celltype_Annotations)
  Celltype_Annotations$rn <- rownames(Celltype_Annotations)
}
head(combined_metadata)
# Merge metadata with cell type annotations by 'rn'
final_annotated_metadata <- combined_metadata %>%
  left_join(Celltype_Annotations, by = "rn")

# Check result


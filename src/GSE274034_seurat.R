library(Seurat)

base_path <- "/vscratch/aarathy/external_datasets/murine_HSC/GSE274034/"
datasets <- c("LTHSC", "Progenitor", "ReplicateA", "ReplicateB")
seurat_objects <- list()

for (ds in datasets) {
  data_dir <- file.path(base_path, ds)
  
  # Read 10X data
  counts <- Read10X(data.dir = data_dir)
  
  # If counts is a list (multiple types), take only "Gene Expression"
  if (is.list(counts)) {
    counts <- counts[["Gene Expression"]]
  }
  
  # Create Seurat object
  seurat_objects[[ds]] <- CreateSeuratObject(counts = counts, project = ds)
  
  cat("✅ Loaded dataset:", ds, "\n")
}

seurat_objects$LTHSC@meta.data

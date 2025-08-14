source("src/00_init.R")
library(Matrix)
library(monocle3)
require("sceasy")
basedir <- "SCRNA_02_01_Integration/"
out <- dirout(paste0(basedir, "/Anna_et_al/"))
out.base <- dirout(basedir)
inDir <- dirout_load("Ag_ScRNA_23_external_dataset")
dir <- dirout("Ag_ScRNA_23_external_dataset/untarred")

##############


# Function to load a 10X dataset and return a CDS
make_cds <- function(path, sample_name) {
  mat <- readMM(dir(paste0("/",path, "matrix.mtx.gz")))
  barcodes <- readLines(gzfile(dir(paste0("/",path, "barcodes.tsv.gz"))))
  features <- read.delim(gzfile(dir(paste0("/",path, "features.tsv.gz"))), header = FALSE)
  
  rownames(mat) <- features$V1
  colnames(mat) <- paste(sample_name, barcodes, sep = "_")  # Ensure unique cell names
  
  gene_metadata <- data.frame(gene_short_name = features$V2)
  rownames(gene_metadata) <- features$V1
  
  cell_metadata <- data.frame(barcode = barcodes, sample = sample_name)
  rownames(cell_metadata) <- colnames(mat)
  
  new_cell_data_set(mat,
                    cell_metadata = cell_metadata,
                    gene_metadata = gene_metadata)
}

# Paths to datasets and sample labels
paths <- c("Sca1neg_scRNAseq/",
           "Sca1pos_scRNAseq/")
samples <- c("Sca1neg", "Sca1pos")

# Load and label each dataset
cds_list <- mapply(make_cds, paths, samples, SIMPLIFY = FALSE)

# Combine datasets
cds <- combine_cds(cds_list, keep_all_genes = TRUE)

#Remove cells with zero total counts
cds <- cds[, Matrix::colSums(SingleCellExperiment::counts(cds)) != 0]
# Estimate size factors
cds <- estimate_size_factors(cds)  # 
# Proceed to preprocessing
cds <- preprocess_cds(cds)

#UMAP
cds <- reduce_dimension(cds,preprocess_method = "PCA", verbose = TRUE)
set.seed(42)
cds <- align_cds(cds, alignment_group = "sample", 
                         verbose = TRUE
)
cds <- reduce_dimension(cds,
                                reduction_method = "UMAP",
                                preprocess_method = "Aligned",
                                verbose = TRUE)

plot_cells(cds, color_cells_by = "sample")
# Clustering
set.seed(12121)
cds = cluster_cells(cds)

monocle.file <- out("MonocleObject.RData")
# Store full dataset
save(cds, file=monocle.file)



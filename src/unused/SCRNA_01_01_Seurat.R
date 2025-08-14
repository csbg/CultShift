source("src/00_init.R")
library(SoupX)
library(Matrix)
library(Seurat)
library(ggplot2)
library(data.table)
require("sceasy")

out <- dirout("Ag_SCRNA_01_01_Seurat/")
# Read cellranger analysis results --------------------------------------------
SANN <- fread("metadata/annotation.tsv", fill=TRUE, sep="\t", header = TRUE)
# Remove empty rows
SANN <- SANN[!New_Name %in% c("", "PENDING")]
# Remove bad columns
SANN <- SANN[, sapply(SANN, function(col) sum(!is.na(col)) > 0), with=F]
SANN <- SANN[,-c("Guides", "Comments", "Technology", "atCIMA"),with=F]
# clean markers
SANN[,markers := tolower(gsub("\\,", "", markers))]
SANN[markers %in% c("noMarker", ""),markers := "unsorted"]
table(SANN$markers)
# clean tissues
table(SANN$tissue)
SANN[,tissue := gsub("^exVivo2$", "exVivo", tissue)]
SANN[,tissue := gsub("^exVivo1$", "exVivo", tissue)]
SANN[,tissue := gsub("V", " v", tissue)]
SANN[,tissue := gsub(" ", ".", tissue)]
table(SANN$tissue)
# clean time-points
SANN[is.na(timepoint),timepoint := "d0"]
table(SANN$timepoint)
# add sample number
SANN$s <- paste0("s", 1:nrow(SANN))
#SANN[, sample_new := gsub(" ", "", paste(tissue, markers, timepoint, s, sep="_"))]
SANN[, sample_old := sample]
SANN[, sample_new := New_Name]
SANN[, sample := sample_new]
SANN[, sample_broad := paste0(tissue, "_", gsub("^.+(OP.).+$", "\\1", sample), "_", timepoint)]
SANN[!grepl("OP\\d", sample_broad), sample_broad := sample_new]
#SANN[grepl("CITESEQ", sample_old), sample_broad := sample_old]
SANN[sample %in% SANN[,.N, by="sample"][N > 1]$sample]
stopifnot(!any(duplicated(SANN$sample)))
SANN$New_Name <- NULL


# Match samples to those in the folder ------------------------------------

# based on md5sums
ff.old <- list.files("/media/AGFORTELNY/PROJECTS/TfCf/Data/" ,pattern="metrics_summary.csv", recursive = TRUE, full.names = TRUE)
ff.new <- c()
md5sums <- list()
for(fx in c(ff.old, ff.new)){
  md5sums[[fx]] <- system(paste("md5sum", fx), intern = TRUE)
}
md5sumsDT <- setNames(data.table(do.call(rbind, strsplit(unlist(md5sums), "\\s")))[,c(1,3), with=F], c("md5sum", "file.NF"))
md5sumsDT <- merge(md5sumsDT, SANN[,c("md5sum","md5sumFile"), with=F], by="md5sum", all=TRUE)
md5sumsDT <- md5sumsDT[md5sum != ""]
md5sumsDT[, file.NF := gsub("/media/AGFORTELNY/PROJECTS/TfCf/", "", file.NF)]
md5sumsDT[md5sumFile %in% md5sumsDT[!is.na(md5sumFile),.N, by="md5sumFile"][N > 1]$md5sumFile]
stopifnot(nrow(md5sumsDT[!is.na(md5sumFile),.N, by="md5sumFile"][N > 1]) == 0)
write.tsv(md5sumsDT, out("Files.tsv"))

# Final listing of samples
ff <- list.files
#SANN[is.na(sample_found) & New_Name %in% ff, sample_found := New_Name]
SANN <- SANN[!is.na(sample_found)]
SANN$md5sumFound <- "NA"
for(i in 1:nrow(SANN)){("/media/AGFORTELNY/PROJECTS/TfCf/Data/")
ff[!ff %in% SANN$sample_new]
SANN[sample_new %in% ff, sample_found := sample_new]
  #if(SANN[i]$md5sumFile == "") next
  md5 <- system(paste0("md5sum ", "/media/AGFORTELNY/PROJECTS/TfCf/Data/", "/", SANN[i]$sample_new, "/outs/metrics_summary.csv"), intern = TRUE)
  SANN[i]$md5sumFound <- md5
}
SANN[, md5sumFound := gsub(" .+$", "", md5sumFound)]
stopifnot(nrow(SANN[md5sum != md5sumFound][!grepl("CITEseq", sample, ignore.case = TRUE)]) == 0)
write.tsv(SANN, out("SampleAnnotation.tsv"))


process_with_soupx <- function(seurat.obj, filtered.data, raw.data) {
  
  #filtered_cells <- colnames(filtered.gx)
  #raw_cells <- colnames(raw.gx)
  #common_cells <- intersect(filtered_cells, raw_cells)
  
  #table(is.table(colnames(filtered.gx) %in% colnames(raw.gx)))
  # Subset both datasets to include only common cells
  #filtered.gx <- filtered.gx[, colnames(filtered.gx) %in% common_cells]
  #raw.gx <- raw.gx[, colnames(raw.gx) %in% common_cells]
  sc <- tryCatch({
    SoupChannel(raw.gx, filtered.gx, exprsType = "UMI")
  }, error = function(e) {
    stop("Error creating SoupChannel object: ", e$message)
  })
  # Create Seurat object with the updated filtered data
  seurat.obj <- CreateSeuratObject(counts = filtered.gx)
  
  #sc = SoupChannel(toc,toc,calcSoupProfile=FALSE)
  #And add manually
  #rowSums = Matrix::rowSums
  # soupProf = data.frame(row.names = rownames(toc),est=rowSums(toc)/sum(toc),
  #                       counts=rowSums(toc))
  #sc = setSoupProfile(sc,soupProf)
  
  # Cluster cells using Seurat
  seurat.obj <- SCTransform(seurat.obj, verbose = FALSE)
  seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:30, verbose = FALSE)
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:30, verbose = FALSE)
  seurat.obj <- FindClusters(seurat.obj, verbose = TRUE)
  
  meta <- seurat.obj@meta.data
  head(sc)
  umap <- seurat.obj@reductions$umap@cell.embeddings
  
  # Set clusters and dimensionality reduction in SoupChannel
  print(head(rownames(seurat.obj@meta.data)))
  print(head(colnames(sc$toc)))
  sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc <- setDR(sc, umap)
  
  # Estimate contamination
  sc <- autoEstCont(sc,forceAccept = T)
  pdf(out(paste0(dsx, "SOUPX_plot.pdf")))
  sc <- autoEstCont(sc, forceAccept = TRUE)
  dev.off()
  corrected.data <- adjustCounts(sc)
  
  # Update Seurat object
  if (nrow(corrected.data) > 0 && ncol(corrected.data) > 0) {
    #seurat.obj[["original_counts"]] <- seurat.obj@assays$RNA@counts
    seurat.obj@assays$RNA@counts <- corrected.data
    message("Updated Seurat object with cleaned data")
  } else {
    warning("No valid cleaned data returned from SoupX")
  }
  
  return(seurat.obj)
}



# Processing each dataset
additional.info.x <- list()
# Read data, Seurat processing
dsx <- SANN$sample_found[40]
for(dsx in unique(SANN[tissue %in% c("in.vivo","ex.vivo"),]$sample_found)[30:40]) {
  
  # File
  dsx.md5sum <- SANN[sample_found == dsx]$md5sumFound[1]
  # already processed?
  # if(file.exists(dsx.file)){
  #   # Check using md5sum
  #   if(file.exists(dsx.md5sum.file)){
  #     if(dsx.md5sum == fread(dsx.md5sum.file, header = F)$V1[1]){
  #       message(dsx, " already processed and new data is identical")
  #       next
  #     }
  #   }
    
    # check using md5sum from data
    
    
  
  # File paths
  filtered.path <- paste("/media/AGFORTELNY/PROJECTS/TfCf/Data/", dsx, "outs", "filtered_feature_bc_matrix.h5", sep = "/")
  raw.path <- paste("/media/AGFORTELNY/PROJECTS/TfCf/Data/", dsx, "outs", "raw_feature_bc_matrix.h5", sep = "/")
  path.guides <- paste("/media/AGFORTELNY/PROJECTS/TfCf/Data/", dsx, "outs", "crispr_analysis", "protospacer_calls_per_cell.csv", sep = "/")
  
  print("------------------------------")
  print(dsx)
  print("---------------")
  
  # Check if files exist
  if(!file.exists(filtered.path) || !file.exists(raw.path)) {
    message("Filtered or raw file for ", dsx, " not found")
    next
  }
  
  # Read in filtered and raw data
  filtered.data <- Read10X_h5(filtered.path)[["Gene Expression"]]
  
  raw.data <- Read10X_h5(raw.path)[["Gene Expression"]]
  
  filtered.gx <- NA
  raw.gx <- NA
  
  if (class(filtered.data) == "dgCMatrix") {
    filtered.gx <- filtered.data
  }
  
  if (class(raw.data) == "dgCMatrix") {
    raw.gx <- raw.data
  }
  
  if (is.list(filtered.data) && "Gene Expression" %in% names(filtered.data)) {
    filtered.gx <- filtered.data[["Gene Expression"]]
  }
  
  if (is.list(raw.data) && "Gene Expression" %in% names(raw.data)) {
    raw.gx <- raw.data[["Gene Expression"]]
  }
  
  # # Load additional guide data
  additional.info.x <- list()
  if(file.exists(path.guides)) {
    additional.info.x[["CRISPR_Cellranger"]] <- fread(path.guides)
    if(nrow(additional.info.x[["CRISPR_Cellranger"]]) == 0) additional.info.x[["CRISPR_Cellranger"]] <- NULL
  }
  # 
  # Create Seurat object
  if(class(filtered.gx) != "dgCMatrix") {
    message("Filtered data for ", dsx, " not in the right format")
  }
  
  if(dsx == "DM_Test1_NM_6d_1"){
    matrix_dir = paste(PATHS$LOCATIONS$DATA, "ECCITE1_citeseq_combined//umi_count/", sep="/")
    barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
    features.path <- paste0(matrix_dir, "features.tsv.gz")
    matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
    colnames(mat) = paste0(barcode.names$V1, "-1")
    rownames(mat) = gsub("\\-.+$", "", feature.names$V1)
    mat <- mat[row.names(mat) != "unmapped",]
    additional.info.x <- list("CRISPR.Guide.Capture" = as(mat, "dgCMatrix"))
  }
  
  # Get cellranger guide assignments
  if(file.exists(path.guides)){
    additional.info.x[["CRISPR_Cellranger"]] <- fread(path.guides)
    if(nrow(additional.info.x[["CRISPR_Cellranger"]])  == 0) additional.info.x[["CRISPR_Cellranger"]] <- NULL
  }
  
  # Create Seurat object
  if(class(filtered.gx) != "dgCMatrix"){
    message("Data for ", dsx, " not in the right format")
  }
  seurat.obj <- CreateSeuratObject(counts = filtered.gx, project = dsx)

  
  # Filter and normalize object
  write.tsv(data.table(seurat.obj@meta.data, keep.rownames = TRUE), out(paste0("AnnotationOriginal_",dsx,".tsv")))
  seurat.obj <- process_with_soupx(seurat.obj, filtered.gx, raw.gx)
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
  seurat.obj$md5sum <- dsx.md5sum
  #seurat.obj <- seurat.obj
  # set cutoffs
  cutoffs <- list(
    nFeature_RNA = max(500, quantile(seurat.obj@meta.data$nFeature_RNA, probs=0.90)/5),
    nCount_RNA = max(1000, quantile(seurat.obj@meta.data$nCount_RNA, probs=0.90)/5),
    percent.mt = 10
  )
  # Plots
  for(qcx in names(cutoffs)){
    ggplot(data.frame(Value=seurat.obj@meta.data[[qcx]]), aes(x=Value)) + stat_density() + scale_x_log10() +
      theme_bw(12) +
      geom_vline(xintercept = cutoffs[[qcx]]) +
      xlab(qcx) +
      ggtitle(paste(dsx))
    ggsave(out("QC_", qcx, "_", dsx, ".pdf"), w=6,h=4)
  }
  
  # Subset and process
  seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > cutoffs$nFeature_RNA & nCount_RNA > cutoffs$nCount_RNA & percent.mt < cutoffs$percent.mt)
  seurat.obj <- NormalizeData(seurat.obj, verbose = FALSE)
  seurat.obj <- CellCycleScoring(seurat.obj, s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes,set.ident = TRUE)
  
  # Add additional info to metadata
  tx <- names(additional.info.x)[1]
  for(tx in setdiff(names(additional.info.x), "CRISPR_Cellranger")){
    txn <- make.names(tx)
    print(tx)
    
    # Get dataset
    x <- additional.info.x[[tx]]
    if(class(x) != "dgCMatrix") next
    
    # Get features for each cells
    if(txn == "CRISPR.Guide.Capture"){
      colcnts <- colSums(x)
      colcnts.min <- min(5, median(colcnts[colcnts != 0]))
      # if the median nUMI per cell is lower than 5, then the cutoff is that median. Else it is 5
      x.max <- apply(x, 2, max)
      x.sum <- apply(x, 2, sum)
      x <- as.matrix(x[,x.max/x.sum > 0.75 & x.sum >= colcnts.min])
      ii <- apply(x, 2, function(col) which(col == max(col)))
    } else {
      cutoff <- 5
      x <- as.matrix(x[,apply(x, 2, max) >= cutoff,drop=F]) # only keep cells with any guide counted above cutoff
      ii <- apply(x, 2, function(col) which(col >= cutoff)) # for each cell (column) get the rows that are above the cutoff
    }
    
    # label features and combine them
    if(length(ii) == 0) next
    labelsx <- sapply(ii, function(ix) paste(row.names(x)[ix], collapse = ","))
    seurat.obj@meta.data[[txn]] <- labelsx[row.names(seurat.obj@meta.data)]
    print(table(seurat.obj@meta.data[[txn]]))
  }
  
  # Cellranger guide counts
  if("CRISPR_Cellranger" %in% names(additional.info.x)){
    x <- additional.info.x$CRISPR_Cellranger
    # Get the index of the top guide
    if(any(grepl("\\|", x$feature_call))){
      x$top.guide <- as.numeric(sapply(strsplit(x$num_umis, "\\|"), function(i){
        i <- as.numeric(i)
        #if(sum(i) < 5) return(NA)
        i <- which(i/sum(i) > 0.75) # which guides has more than 75% of reads?
        if(length(i) == 0) NA else paste(i, collapse = "|") # collapse just in case, checked in the next step
      }))
      stopifnot(!any(grepl("\\|", x$top.guide)))
      # Get the name of the guide from the index above
      x$guide <- sapply(1:nrow(x), function(i){
        xx <- x[i]
        strsplit(xx$feature_call, "\\|")[[1]][xx$top.guide]
      })
    } else {
      x$guide <- x$feature_call
    }
    
    seurat.obj@meta.data[["CRISPR_Cellranger"]] <- setNames(x$guide, x$cell_barcode)[row.names(seurat.obj@meta.data)]
  }
  
  # Mixscape
  guides.use <- intersect(c("CRISPR_Cellranger", "CRISPR.Guide.Capture"), colnames(seurat.obj@meta.data))[1]
  if(!is.na(guides.use)){
    orig.guides <- seurat.obj@meta.data[[guides.use]]
    guides.clean <- orig.guides
    guides.clean[grepl(",", orig.guides)] <- NA
    guides.clean[grepl("^NTC_", orig.guides)] <- "NTC"
    guides.clean[guides.clean %in% names(which(table(guides.clean) < 5))] <- NA
    
    if(sum(!is.na(guides.clean)) > 5 & "NTC" %in% guides.clean){
      seurat.obj@meta.data$guide <- guides.clean
      
      # Mixscape
      eccite <- subset(seurat.obj, cells=row.names(seurat.obj@meta.data)[!is.na(seurat.obj$guide)])
      eccite <- subset(eccite, features=names(which(Matrix::rowSums(eccite@assays$RNA@counts) > 10)))
      eccite <- FindVariableFeatures(eccite, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      #eccite <- subset(eccite, features=eccite@assays$RNA@var.features)
      eccite <- ScaleData(eccite, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
      tryCatch(expr = {
        eccite <- RunPCA(eccite, npcs = 30, verbose = FALSE)
      }, error=function(e){
        message("RunPCA failed, running with only 10 PCs next")
        eccite <- RunPCA(eccite, npcs = 10, verbose = FALSE)
      })
      
      # CalcPerturbSig
      for(nn in seq(from = 20,to = 5, by = -5)){
        tryCatch({
          eccite <- CalcPerturbSig(
            object = eccite, assay = "RNA", slot = "data",
            gd.class ="guide",nt.cell.class = "NTC",
            reduction = "pca",ndims = 20, num.neighbors = nn,
            new.assay.name = "PRTB")
          last
        }, error=function(e){
          message(e)
        })
      }
      
      # RunMixscape
      eccite <- ScaleData(object = eccite, assay = "PRTB", do.scale = F, do.center = T)
      eccite <- RunMixscape(
        object = eccite,assay = "PRTB",slot = "scale.data",
        labels = "guide",nt.class.name = "NTC",
        min.de.genes = 5,
        iter.num = 10,
        de.assay = "RNA",
        verbose = TRUE,
        prtb.type = "KO")
      for(cx in c("mixscape_class", "mixscape_class_p_ko", "mixscape_class.global")){
        seurat.obj@meta.data[[cx]] <- eccite@meta.data[row.names(seurat.obj@meta.data),cx]
      }
    }
  }
  
  # Remove metabolic guides from Nisha's side project
  rm.guides <- c("Mat2a", "Nampt", "Ogdh", "Sdha", "Fh1", "Scd1")
  cx <- "guide"
  for(cx in colnames(seurat.obj@meta.data)){
    message(cx)
    gx <- rm.guides[1]
    
    for(gx in rm.guides){
      print(sum(idx <- grepl(gx, seurat.obj@meta.data[[cx]])))
      if(sum(idx) > 0){
        stopifnot(any(is.na(seurat.obj@meta.data[[cx]])))
        seurat.obj@meta.data[[cx]][idx] <- NA
      }
    }
  }
  
  # Export data
  seurat.obj <- SCRNA.DietSeurat(sobj = seurat.obj)
  # Save the updated Seurat object
  save(seurat.obj,additional.info.x, file = out(paste0("SeuratObj_", dsx, ".RData")))
  write(seurat.obj$md5sum[1], out(paste0("Md5sum_", dsx, ".txt")))
  
}

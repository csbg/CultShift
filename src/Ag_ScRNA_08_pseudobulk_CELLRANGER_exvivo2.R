
#############
source("src/00_init.R")
library("dplyr")
library("stringr")
library("tidyverse")
#############
#paths------
#############
out <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide_exvivo2")
InDir <- dirout("Invivo_Invitro_DAVID_new")


sobjs <- list(
  ex.vivo = NULL,
  in.vivo = NULL
)

paths <- c(
  in.vivo = InDir("Invivo4_5/rna_integrated_annotated_invivo.rds"),
  ex.vivo = InDir("aarathy_Invitro1_20260310/rna_integrated_annotated_proj.rds")
)

for (nm in names(paths)) {
  message("Loading: ", paths[nm])
  
  sobjs[[nm]] <- read_rds(paths[nm])
}
files <- list.files(
  "/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis/Invivo_Invitro_DAVID_new/",
  pattern = "Invitro.*protospacer_calls_per_cell.csv",
  full.names = TRUE
)


meta <- sobjs$ex.vivo@meta.data
# initialize metadata columns
meta$crispr_num_features <- NA
meta$crispr_feature_call <- NA
meta$crispr_num_umis <- NA
meta$crispr_cell_barcode <- NA

for (f in files) {
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$cell_barcode <- trimws(df$cell_barcode)
  
  # sample name from file
  sample_name <- sub("_protospacer_calls_per_cell.csv$", "", basename(f))
  head(df)
  # subset metadata for this sample
  cells <- rownames(meta)[meta$orig.ident == sample_name]
  meta_subset <- meta[cells, , drop = FALSE]  # ensure subset is same length as cells
  
  # subset CSV to only barcodes in this sample
  df_sub <- df[df$cell_barcode %in% cells, , drop = FALSE]
  
  # match rows
  match_idx <- match(rownames(meta_subset), df_sub$cell_barcode)
  
  # only valid matches
  valid <- !is.na(match_idx)
  
  if(sum(valid) == 0) {
    warning(paste0("No matching barcodes found for sample: ", sample_name))
    next
  }
  
  # assign values safely
  meta_subset$crispr_num_features[valid] <- df_sub$num_features[match_idx[valid]]
  meta_subset$crispr_feature_call[valid] <- df_sub$feature_call[match_idx[valid]]
  meta_subset$crispr_num_umis[valid] <- df_sub$num_umis[match_idx[valid]]
  meta_subset$crispr_cell_barcode[valid] <- df_sub$cell_barcode[match_idx[valid]]
  
  # write back to full metadata
  meta[cells, ] <- meta_subset
}

unique(sobjs$ex.vivo@meta.data$orig.ident)

cells <- rownames(meta)
df <- read.csv(InDir("Invitro1_LSK2_protospacer_calls_per_cell.csv"), stringsAsFactors = FALSE)
df$cell_barcode <- trimws(df$cell_barcode)

# see intersection
intersect(cells, df$cell_barcode)
# assign updated metadata to Seurat object

for (f in files) {
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$cell_barcode <- trimws(df$cell_barcode)
  
  # sample name from file
  sample_name <- sub("_protospacer_calls_per_cell.csv$", "", basename(f))
  
  # subset metadata for this sample
  cells <- rownames(meta)[meta$orig.ident == sample_name]
  meta_subset <- meta[cells, , drop = FALSE]  # ensure subset is same length as cells
  
  # subset CSV to only barcodes in this sample
  df_sub <- df[df$cell_barcode %in% cells, , drop = FALSE]
  
  # match rows
  match_idx <- match(rownames(meta_subset), df_sub$cell_barcode)
  
  # # only valid matches
  # valid <- !is.na(match_idx)
  # 
  # if(sum(valid) == 0) {
  #   warning(paste0("No matching barcodes found for sample: ", sample_name))
  #   next
  # }
  
  # =========================
  # Count cells with only 1 guide
  # =========================
  num_unique_guides <- sum(df_sub$num_features == 1)
  total_cells <- nrow(df_sub)
  
  cat(sprintf("Sample %s: %d/%d cells have a unique guide (num_features==1)\n",
              sample_name, num_unique_guides, total_cells))
}
####################
unique(sobjs$in.vivo@meta.data$celltype)
unique(sobjs$ex.vivo@meta.data$celltype)

unique(sobjs$ex.vivo@meta.data$orig.ident)
unique(sobjs$in.vivo@meta.data$orig.ident)
sobjs$in.vivo[["RNA"]]
process_seurat_data <- function(
    Seurat_obj,
    celltype_col = "celltype",
    tissue = "in.vivo",
    min_cells = 50
) {
  
  meta <- Seurat_obj@meta.data
  
  # Add barcode column
  meta$rn <- rownames(meta)
  
  # 🔥 Ensure required columns exist
  required_cols <- c(celltype_col, "orig.ident", "guide", "genotype")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # -------------------------
  # Create grouping dataframe
  # -------------------------
  grouping_df <- data.frame(
    cell = meta$rn,
    celltype = meta[[celltype_col]],
    sample = meta$orig.ident,
    guide = meta$guide,
    genotype = meta$genotype,
    stringsAsFactors = FALSE
  )
  
  # Unique group ID (safe, no parsing later)
  grouping_df$group_id <- paste(
    grouping_df$celltype,
    grouping_df$sample,
    grouping_df$guide,
    sep = "||"
  )
  
  # -------------------------
  # Split cells
  # -------------------------
  mat_gt_ct <- split(grouping_df$cell, grouping_df$group_id)
  
  # Filter groups
  cell_counts <- lengths(mat_gt_ct)
  keep <- cell_counts >= min_cells
  mat_gt_ct <- mat_gt_ct[keep]
  cell_counts <- cell_counts[keep]
  
  message("Keeping ", length(mat_gt_ct), " pseudobulk groups with ≥ ", min_cells, " cells")
  
  # -------------------------
  # RNA counts (Seurat v5 safe)
  # -------------------------
  rna_assay <- Seurat_obj[["RNA"]]
  
  if (!("counts" %in% names(rna_assay@layers))) {
    stop("RNA assay has no counts layer.")
  }
  
  counts_layer <- rna_assay@layers$counts
  
  # Gene + cell names
  gene_names <- rownames(Seurat_obj@assays$RNA@features)
  cell_names <- rownames(meta)
  
  rownames(counts_layer) <- gene_names
  colnames(counts_layer) <- cell_names
  
  # -------------------------
  # Pseudobulk aggregation
  # -------------------------
  result <- sapply(mat_gt_ct, function(bcs) {
    idx <- intersect(bcs, cell_names)
    Matrix::rowSums(counts_layer[, idx, drop = FALSE])
  })
  
  rownames(result) <- gene_names
  
  # -------------------------
  # Build metadata properly
  # -------------------------
  group_info <- do.call(rbind, strsplit(names(mat_gt_ct), "\\|\\|"))
  colnames(group_info) <- c("celltype", "sample", "guide")
  
  meta_out <- data.frame(
    group_id = names(mat_gt_ct),
    celltype = group_info[, "celltype"],
    sample = group_info[, "sample"],
    guide = group_info[, "guide"],
    tissue = tissue,
    n_cells = cell_counts,
    stringsAsFactors = FALSE
  )
  
  # 🔥 Add genotype (map from original metadata)
  genotype_map <- unique(meta[, c("guide", "genotype")])
  meta_out <- merge(meta_out, genotype_map, by = "guide", all.x = TRUE)
  
  rownames(meta_out) <- meta_out$group_id
  
  return(list(result = result, meta = meta_out))
}

in.vivo_pseudobulk <- process_seurat_data(
  Seurat_obj = sobjs$in.vivo,
  celltype_col = "celltype",
  tissue = "in.vivo",
  min_cells = 50
)
ex.vivo_pseudobulk <- process_seurat_data(
  Seurat_obj = sobjs$ex.vivo,
  celltype_col = "celltype",
  tissue = "ex.vivo",
  min_cells = 50
)
invivo_counts <- in.vivo_pseudobulk$result
invivo_meta <- in.vivo_pseudobulk$meta

exvivo_counts <- ex.vivo_pseudobulk$result
exvivo_meta <- ex.vivo_pseudobulk$meta


correct_celltype_names <- function(names_vec) {
  names_vec <-  gsub("\\|\\|", "_",names_vec)
  names_vec <- gsub("GMP \\(early\\)", "GMP.early", names_vec)
  names_vec <- gsub("GMP \\(late\\)", "GMP.late", names_vec)
  names_vec <- gsub("Gran\\. P", "Gran.P", names_vec)
  names_vec <- gsub("MEP \\(G1\\)", "MEP.G1", names_vec)
  names_vec <- gsub("MEP \\(pert\\.\\)", "MEP.pert.", names_vec)
  names_vec <- gsub("MEP \\(S\\)", "MEP.S", names_vec)
  names_vec <- gsub("MEP \\(early\\)", "MEP.early", names_vec)
  names_vec <- gsub("Imm\\. B-cell", "Imm.B.cell", names_vec)
  
}
# ex vivo
colnames(exvivo_counts) <- correct_celltype_names(colnames(exvivo_counts))

# in vivo
colnames(invivo_counts) <- correct_celltype_names(colnames(invivo_counts))
#meta
exvivo_meta$celltype <- correct_celltype_names(exvivo_meta$celltype)
invivo_meta$celltype <- correct_celltype_names(invivo_meta$celltype)
exvivo_meta$group_id <- correct_celltype_names(exvivo_meta$group_id)
invivo_meta$group_id <- correct_celltype_names(invivo_meta$group_id)
# ex vivo
rownames(exvivo_meta) <- correct_celltype_names(rownames(exvivo_meta))

# in vivo
colnames(invivo_meta) <- correct_celltype_names(colnames(invivo_meta))
head(exvivo_meta)
head(exvivo_counts)
###############
saveRDS(exvivo_counts,out("ex.vivo_counts.rds"))
saveRDS(invivo_counts,out("in.vivo_counts.rds"))
saveRDS(exvivo_meta,out("ex.vivo_meta.rds"))
saveRDS(invivo_meta,out("in.vivo_meta.rds"))

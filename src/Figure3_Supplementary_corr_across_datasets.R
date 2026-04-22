source("src/00_init.R")
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
# Quick heatmap
library(ComplexHeatmap)
library(circlize)
library(fgsea)
library(msigdbr)
library(readxl)
library(openxlsx)

library(Hmisc)

source("src/Ag_Optimized_theme_fig.R")
# -----------------------------
# Parent folder containing all datasets
# -----------------------------
parent_folder <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis/"
out <- dirout("Figure3_supplementary_logFC_comparison_v8")
# -----------------------------
# List of dataset folder names (relative to parent)
# -----------------------------
dataset_names <- c(
  "GSE169368_invitro_organoids_vs_invivo_mouse",
  "GSE125213_human_CB1",
  "GSE110968_CB2_external_dataset",
  "GSE149244vsGSE125211_LT-HSC",
  "GSE139184_pancreatic_cancer",
  "JAK_STAT_splenic_M_T",
  "Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide",
  "Glioblastoma_limmaRes"
  
)

# Full paths to dataset folders
dataset_folders <- file.path(parent_folder, dataset_names)

# -----------------------------
# Assemble all ranklists
# -----------------------------
assembled_ranklists <- map_dfr(dataset_folders, function(ds_folder){
  
  dataset_name <- basename(ds_folder)  # Just the folder name
  ranklist_folder <- file.path(ds_folder, "Ranklists")  # Ranklists subfolder
  
  # Skip if Ranklists folder does not exist
  if(!dir.exists(ranklist_folder)){
    message("No Ranklists folder in ", ds_folder)
    return(NULL)
  }
  
  # List all *_logFC_ranklist.rds files
  rds_files <- list.files(ranklist_folder, pattern = "_logFC_ranklist\\.rds$", full.names = TRUE)
  
  if(length(rds_files) == 0){
    message("No RDS files found in ", ranklist_folder)
    return(NULL)
  }
  
  # Read each file and add dataset + contrast info
  map_dfr(rds_files, function(f){
    contrast_name <- gsub("_logFC_ranklist\\.rds$", "", basename(f))
    ranks <- readRDS(f)
    
    tibble(
      gene = toupper(names(ranks)),
      logFC = as.numeric(ranks),
      dataset = dataset_name,
      contrast = contrast_name
    )
  })
})

# -----------------------------
# Save assembled table
# -----------------------------
assembled_dir <- file.path(parent_folder, "Ranklists_assembled")
dir.create(assembled_dir, showWarnings = FALSE)

saveRDS(assembled_ranklists, file = file.path(assembled_dir, "all_ranklists.rds"))
write.csv(assembled_ranklists, file = file.path(assembled_dir, "all_ranklists.csv"), row.names = FALSE)

# -----------------------------
# Pivot to wide format for correlation
# -----------------------------
rank_wide <- assembled_ranklists %>%
  unite("dataset_contrast", dataset, contrast, sep = "_") %>%
  pivot_wider(names_from = dataset_contrast, values_from = logFC)

# Convert to matrix for correlation
logFC_mat <- rank_wide %>% column_to_rownames("gene") %>% as.matrix()


colnames(logFC_mat) <- c("Intestinal organoids Mm",
                         "Cord blood 2days ex vivo Hs(1)",
                         "Cord blood 4days ex vivo Hs(1)",
                         "Cord blood 2days ex vivo Hs(2)",
                         "Cord blood 4days ex vivo Hs(2)",
              
                         "LT-HSC ex-vivo Mm",
                         "PDAC ex vivo-2D Hs",
                         "PDAC ex vivo-3D Hs",
                         "PDAC ex vivo-organoid Hs",
                         "PDAC ex vivo-suspension Hs",
                         "PDAC Mm-Xenograft-Hs",
                         "Splenic M ex vivo 20h Mm",
                         "Splenic T-cells ex vivo 20h Mm",
                         "Eo.Ba Hematopoietic ex vivo",
                         "GMP Hematopoietic ex vivo",
                         "Gran. Hematopoietic ex vivo",
                         "Gran.P Hematopoietic ex vivo",
                         "HSC Hematopoietic ex vivo",
                         "MEP Hematopoietic ex vivo",
                         "MkP Hematopoietic ex vivo",
                         "Mono Hematopoietic ex vivo",
                         "Glioblastoma ex vivo"
)
write.csv(logFC_mat, out("logFC_datasets.csv"))
logFC_ranks <- read_csv(out("logFC_datasets.csv"))
colnames(logFC_ranks)[1] <- "Genes"

ann <- as.data.frame(logFC_ranks)

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, sheetName = "logFCs")

# Write data
writeData(wb, sheet = "logFCs", ann, rowNames = FALSE)

# Freeze first row
freezePane(wb, sheet = "logFCs", firstRow = TRUE, firstCol = FALSE)

# Bold header
headerStyle <- createStyle(textDecoration = "bold")
addStyle(wb, sheet = "logFCs", headerStyle,
         rows = 1,
         cols = 1:ncol(ann),
         gridExpand = TRUE)

# Save workbook (IMPORTANT — missing in your code)
saveWorkbook(
  wb,
  file = file.path(out(), "Supplementary_Table3_logFC_across_datasets.xlsx"),
  overwrite = TRUE
)
# Compute Pearson correlation across contrasts
cor_matrix <- cor(logFC_mat, use = "pairwise.complete.obs", method = "pearson")

# Save correlation matrix
write.csv(cor_matrix, file = file.path(assembled_dir, "logFC_correlation_matrix.csv"))



# Define color function for correlation
col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#4C889C", "white", "#D0154E")
)

# Create a tmp folder inside your project
tmpdir <- file.path(parent_folder, "tmp")
dir.create(tmpdir, showWarnings = FALSE)

# Tell R to use it
Sys.setenv(TMPDIR = tmpdir)


library(grid)  # for gpar


pdf(
  out("Supp.Fig.3B.logFC_cor_across_datasets.pdf"),
  width  = 10 / 2.54,
  height = 7 / 2.54,
  compress = FALSE
)

Heatmap(
  cor_matrix,
  name = "Spearman\ncorrelation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  # Row and column label styles
  row_names_gp = gpar(fontsize = 5, fontface = "plain"),
  column_names_gp = gpar(fontsize = 5, fontface = "plain"),
  column_names_rot = 45,  # <-- rotate x-axis labels
  
  # Gaps between cells
  row_gap = unit(0.1, "mm"),
  column_gap = unit(0.1, "mm"),
  
  # Remove dendrograms
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  
  # Legend settings
  heatmap_legend_param = list(
    at = c(-1, -0.5, 0, 0.5, 1),
    legend_height = unit(2, "cm"),
    legend_width = unit(1, "mm"),
    labels_gp = gpar(fontsize = 5, fontface = "plain"),
    title_gp = gpar(fontsize = 5, fontface = "plain")
  ),
  
  # Shrink each cell
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(
      x, y,
      width = width * 0.008,
      height = height * 0.008,
      gp = gpar(fill = fill, col = NA)
    )
  }
)

dev.off()




##############


ranked_lists <- assembled_ranklists %>%
  unite("coef", dataset, contrast, sep = "_") %>%
  group_by(coef) %>%
  summarise(
    ranks = list({
      r <- logFC
      names(r) <- gene
      r <- r[!is.na(r)]
      r <- r[!duplicated(names(r))]
      sort(r, decreasing = TRUE)
    }),
    .groups = "drop"
  )
pathways <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  split(x = .$gene_symbol, f = .$gs_name)

coef_map <- c(
  "GSE169368_invitro_organoids_vs_invivo_mouse_tissueex.vivo" = "Intestinal organoids Mm",
  "GSE125213_human_CB1_tissueex.vivo_2d" = "Cord blood 2days ex vivo Hs(1)",
  "GSE125213_human_CB1_tissueex.vivo_4d" = "Cord blood 4days ex vivo Hs(1)",
  "GSE110968_CB2_external_dataset_tissueex.vivo_2days" = "Cord blood 2days ex vivo Hs(2)",
  "GSE110968_CB2_external_dataset_tissueex.vivo_4days" = "Cord blood 4days ex vivo Hs(2)",
  "GSE149244vsGSE125211_LT-HSC_tissueex.vivo" = "LT-HSC ex vivo Mm",
  "GSE139184_pancreatic_cancer_tissueex.vivo_2D" = "PDAC ex vivo-2D Hs",
  "GSE139184_pancreatic_cancer_tissueex.vivo_3D" = "PDAC ex vivo-3D Hs",
  "GSE139184_pancreatic_cancer_tissueex.vivo_organoid" = "PDAC ex vivo-organoid Hs",
  "GSE139184_pancreatic_cancer_tissueex.vivo_suspension" = "PDAC ex vivo-suspension Hs",
  "GSE139184_pancreatic_cancer_tissuein.vivo_xenograft" = "PDAC Mm-Xenograft-Hs",
  "JAK_STAT_splenic_M_T_ex-vivoMacrophage" = "Splenic M ex vivo 20h Mm",
  "JAK_STAT_splenic_M_T_ex-vivoT-cells" = "Splenic T-cells ex vivo 20h Mm",
  "Glioblastoma_limmaRes_tissueex.vivo" = "Glioblastoma GL261 ex vivo"
)


ranked_lists <- ranked_lists %>%
  mutate(
    coef_pretty = coef_map[coef]
  )
ranked_lists <- ranked_lists %>%
  mutate(coef = coef_map[coef]) %>%
  dplyr::select(coef, ranks) %>%
  na.omit()

#fgsea
# fgsea_results <- ranked_lists %>%
#   mutate(
#     fgsea = map(
#       ranks,
#       ~ fgsea(
#         pathways = pathways,
#         stats = .x,
#         minSize = 1,
#         maxSize = 500,
#         nPermSimple = 100000
#       )
#     )
#   ) %>%
#   dplyr::select(coef, fgsea) %>%
#   unnest(fgsea)
# fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, paste, collapse = ",")
# write.csv(fgsea_results, out("fgsea.csv"), row.names = FALSE)

Fgsea <- read_csv(out("fgsea.csv"))
colnames(Fgsea)[1] <- "Dataset"
colnames(Fgsea) <- c("Dataset","pathway","pval" ,"padj",  "log2err" ,
                     "ES",
                     "NES",         "size",
                     "leadingEdge")
ann <- as.data.frame(Fgsea[,c("Dataset" , "pathway",  "padj", "NES",  "leadingEdge")])

# Create workbook
wb <- createWorkbook()

# Add worksheet
addWorksheet(wb, sheetName = "Fgsea")

# Write data
writeData(wb, sheet = "Fgsea", ann, rowNames = FALSE)

# Freeze first row
freezePane(wb, sheet = "Fgsea", firstRow = TRUE, firstCol = FALSE)

# Bold header
headerStyle <- createStyle(textDecoration = "bold")
addStyle(wb, sheet = "Fgsea", headerStyle,
         rows = 1,
         cols = 1:ncol(ann),
         gridExpand = TRUE)

# Save workbook (IMPORTANT — missing in your code)
saveWorkbook(
  wb,
  file = file.path(out(), "Supplementary_Table4_enrichment_across_datasets.xlsx"),
  overwrite = TRUE
)

fgsea_results <- Fgsea
terms <- fgsea_results %>%
  filter(padj < 0.05, abs(NES)>2)%>%
  pull(pathway)

pDT <- fgsea_results %>%
  filter(pathway %in% terms) %>%
  mutate(
    Dataset = factor(Dataset, levels = rev(unique(Dataset))),
    pathway = factor(pathway, levels = unique(pathway))
  )



plot <- ggplot(
  pDT,
  aes(
    x = Dataset,
    y = pathway,
    color = pmin(pmax(NES, -2), 2),
    size  = pmin(5, -log10(padj))
  )
) +
  geom_point() +
  
  scale_color_gradient2(
    low  = "#4C889C",
    mid  = "white",
    high = "#D0154E",
    name = TeX("NES")
  ) +
  
  scale_size_continuous(
    range = c(0, 1.5),
    name  = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  
  labs(
    y = "Dataset / Condition",
    x = "Pathways",
    title = "Enriched pathways"
  ) +
  
  optimized_theme_fig() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "bottom"
  )

plot
ggsave(out("Sup.Fig.3C.fgsea_across_datasets.pdf"),plot = plot, w = 12,
       h = 10 , units = "cm")
####################################3


# -----------------------------
# Correlation + stats
# -----------------------------
cor_res <- rcorr(logFC_mat, type = "pearson")

cor_matrix <- cor_res$r
p_matrix   <- cor_res$P

# Overlap (number of shared genes)
n_matrix <- outer(
  1:ncol(logFC_mat),
  1:ncol(logFC_mat),
  Vectorize(function(i, j) {
    sum(!is.na(logFC_mat[, i]) & !is.na(logFC_mat[, j]))
  })
)

colnames(n_matrix) <- colnames(logFC_mat)
rownames(n_matrix) <- colnames(logFC_mat)

# -----------------------------
# Tidy table for supplement
# -----------------------------
cor_table <- reshape2::melt(cor_matrix, 
                            varnames = c("dataset1", "dataset2"), value.name = "r") %>%
  mutate(
    p_value = as.vector(p_matrix),
    n_genes = as.vector(n_matrix)
  ) %>%
  filter(dataset1 != dataset2)
# Save
write.csv(
  cor_table,
  file = file.path(out("Supplementary_Table_correlation_stats.csv")),
  row.names = FALSE
)

#################
#Sup.Fig3B
cor_plot_df <- cor_table %>%
  filter(as.character(dataset1) < as.character(dataset2))

p <- cor_plot_df %>%
  mutate(
    bin = cut(
      n_genes,
      breaks = c(0, 5000, 10000, 15000, 20000, Inf),
      labels = c("0–5k", "5k–10k", "10k–15k", "15k–20k", ">20k"),
      include.lowest = TRUE
    )
  ) %>%
  ggplot(aes(x = bin, y = r)) +
  geom_boxplot(size = 0.2) +   # 👈 thinner lines here
  optimized_theme_fig() +
  labs(x = "Gene overlap (k-scale)", y = "Correlation (r)")

ggsave(out("Supp.Fig3C_corr_vs_overlap.pdf"), p,
       width = 4, height = 3, units = "cm")

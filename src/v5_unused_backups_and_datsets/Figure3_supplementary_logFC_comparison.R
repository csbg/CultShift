
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
# Quick heatmap
library(ComplexHeatmap)
library(circlize)
source("src/Ag_Optimized_theme_fig.R")
# -----------------------------
# Parent folder containing all datasets
# -----------------------------
parent_folder <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis/"

# -----------------------------
# List of dataset folder names (relative to parent)
# -----------------------------
dataset_names <- c(
  "GSE169368_invitro_organoids_vs_invivo_mouse",
  "GSE125213_human_CB1",
  "GSE110968_CB2_external_dataset",
  "GSE149244vsGSE125211_LT-HSC",
  "GSE139184_pancreatic_cancer"
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
colnames(logFC_mat)
colnames(logFC_mat) <- c("GSE169368_intestine_organoids_Mm",
                          "GSE125213_CB1_2d_ex_Hs",
                          "GSE125213_CB1_4d_ex_Hs",
                          "GSE110968_CB2_2d_ex_Hs",
                          "GSE110968_CB2_4d_ex_Hs",
                          "GSE149244vsGSE125211_LT-HSC_ex_Mm",
                          "GSE139184_PDAC_ex_2D_Mm",
                          "GSE139184_PDAC_ex_3D_Mm",
                          "GSE139184_PDAC_organoid_Mm",
                          "GSE139184_PDAC_ex_suspension_Mm",
                          "GSE139184_PDAC_xenograft_Mm"
)
# Compute Pearson correlation across contrasts
cor_matrix <- cor(logFC_mat, use = "pairwise.complete.obs", method = "pearson")

# Save correlation matrix
write.csv(cor_matrix, file = file.path(assembled_dir, "logFC_correlation_matrix.csv"))



# Define color function for correlation
col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "white", "#B2182B")
)


# Draw heatmap
Heatmap(
  cor_matrix,
  name = "Spearman\ncorrelation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_dend = F,
  show_row_dend = F,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  heatmap_legend_param = list(
    at = c(-1, -0.5, 0, 0.5, 1),
    legend_height = unit(5, "cm")
  )
)
##############
library(fgsea)
library(msigdbr)

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
  "GSE169368_invitro_organoids_vs_invivo_mouse_tissueex.vivo" =
    "GSE169368_intestine_organoids_Mm",
  
  "GSE125213_human_CB1_tissueex.vivo_2d" =
    "GSE125213_CB1_2d_ex_Hs",
  
  "GSE125213_human_CB1_tissueex.vivo_4d" =
    "GSE125213_CB1_4d_ex_Hs",
  
  "GSE110968_CB2_external_dataset_tissueex.vivo_2days" =
    "GSE110968_CB2_2d_ex_Hs",
  
  "GSE110968_CB2_external_dataset_tissueex.vivo_4days" =
    "GSE110968_CB2_4d_ex_Hs",
  
  "GSE149244vsGSE125211_LT-HSC_tissueex.vivo" =
    "GSE149244vsGSE125211_LT-HSC_ex_Mm",
  
  "GSE139184_pancreatic_cancer_tissueex.vivo_2D" =
    "GSE139184_PDAC_ex_2D_Mm",
  
  "GSE139184_pancreatic_cancer_tissueex.vivo_3D" =
    "GSE139184_PDAC_ex_3D_Mm",
  
  "GSE139184_pancreatic_cancer_tissueex.vivo_organoid" =
    "GSE139184_PDAC_organoid_Mm",
  
  "GSE139184_pancreatic_cancer_tissueex.vivo_suspension" =
    "GSE139184_PDAC_ex_suspension_Mm",
  
  "GSE139184_pancreatic_cancer_tissuein.vivo_xenograft" =
    "GSE139184_PDAC_xenograft_Mm"
)
ranked_lists <- ranked_lists %>%
  mutate(
    coef_pretty = coef_map[coef]
  )
ranked_lists <- ranked_lists %>%
  mutate(coef = coef_map[coef]) %>%
  dplyr::select(coef, ranks)

#fgsea
fgsea_results <- ranked_lists %>%
  mutate(
    fgsea = map(
      ranks,
      ~ fgsea(
        pathways = pathways,
        stats = .x,
        minSize = 1,
        maxSize = 500
      )
    )
  ) %>%
  dplyr::select(coef, fgsea) %>%
  unnest(fgsea)
terms <- fgsea_results %>%
  filter(padj < 0.05, abs(NES)>2)%>%
  pull(pathway)
library(ggplot2)
library(forcats)
pDT <- fgsea_results %>%
  filter(pathway %in% terms) %>%
  mutate(
    coef = factor(coef, levels = rev(unique(coef))),
    pathway = factor(pathway, levels = unique(pathway))
  )

library(latex2exp)
library(ggplot2)

plot <- ggplot(
  pDT,
  aes(
    x = coef,
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
    range = c(0, 1.8),
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

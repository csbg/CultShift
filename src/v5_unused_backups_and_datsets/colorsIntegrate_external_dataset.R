##############################
# Modular GEO RNA-seq pipeline (clean + updated)
# Author: Aarathy (optimized by GPT-5)
##############################

# ---------------------------
# Load environment & packages
# ---------------------------
source("src/00_init.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")

library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
library(GEOquery)
library(R.utils)
library(msigdbr)
library(fgsea)
library(latex2exp)
library(Biobase)
library(circlize)
library(patchwork)
library(stringr)
library(biomaRt)
library(scales)

# ---------------------------
# Input / Output directories
# ---------------------------
InDir <- dirout("Ag_top_pathway_genes")
out   <- dirout("Ag_external_datasets_generalized")

# ---------------------------
# Load and prepare gene sets
# ---------------------------
genesets <- fread(InDir("goi_logFC_perturb_seq.tsv")) %>%
  select(genes, logFC, adj.P.Val, pathway, celltype, group)

get_genes_for_pathway <- function(pattern, genesets) {
  unique(genesets[grepl(pattern, pathway, ignore.case = TRUE)]$genes)
}

ISGs         <- get_genes_for_pathway("ISGs", genesets)
ISG_core     <- get_genes_for_pathway("ISG_core", genesets)
mTORC1       <- get_genes_for_pathway("mTORC1", genesets)
Cholesterol  <- get_genes_for_pathway("Cholesterol", genesets)
Housekeeping <- c("Tbp", "Ubc")

mouse_genes <- list(
  ISGs        = ISGs,
  ISG_core    = ISG_core,
  mTORC1      = mTORC1,
  Cholesterol = Cholesterol,
  Housekeeping = Housekeeping
)

# ---------------------------
# Convert mouse → human (using local map)
# ---------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = TRUE)
hm <- unique(hm.map[Human.gene.name != "", .(Mouse = Gene.name, Human = Human.gene.name)])

map_mouse_to_human_keep_all <- function(mouse_genes_list, hm) {
  lapply(mouse_genes_list, function(genes) {
    human_genes <- hm$Human[match(genes, hm$Mouse)]
    unmapped <- genes[is.na(human_genes)]
    human_genes[is.na(human_genes)] <- unmapped
    toupper(unique(human_genes))
  })
}
human_genesets <- map_mouse_to_human_keep_all(mouse_genes, hm)

# ---------------------------
# Dataset definitions
# ---------------------------
datasets <- list(
  AML_primary_vs_cell_line = dirout("GSE229032_AML_primary_vs_cell_line")("AML_NTCs_dataVoom.rds"),
  CB1_human                = dirout("GSE125213_human_CB1")("CB1_NTCs_dataVoom.rds"),
  CB2_human                = dirout("GSE110968_CB2_external_dataset")("CB_NTCs_dataVoom.rds"),
  PDAC_human               = dirout("GSE139184_pancreatic_cancer")("GSE139184_PDAC.rds"),
  Organoids_mouse          = dirout("GSE169368_invitro_organoids_vs_invivo_mouse")("GSE169368_mouse_organoids.rds"),
  LT_HSC_mouse             = dirout("GSE149244vsGSE125211_LT-HSC")("LT-HSC_mouse.rds.rds")
)

species_map <- list(
  AML_primary_vs_cell_line = "human",
  CB1_human                = "human",
  CB2_human                = "human",
  PDAC_human               = "human",
  Organoids_mouse          = "mouse",
  LT_HSC_mouse             = "mouse"
)

# ---------------------------
# Compute z-scores per dataset
# ---------------------------
compute_dataset_zscores <- function(data, dataset_name, species, gene_sets, selected_sets) {
  selected_gene_sets <- gene_sets[selected_sets]
  if (species == "human") {
    data <- data %>% mutate(genes = toupper(genes))
    selected_gene_sets <- lapply(selected_gene_sets, toupper)
  }
  
  genes_of_interest <- unique(unlist(selected_gene_sets))
  data <- data %>% filter(genes %in% genes_of_interest)
  
  data <- data %>%
    mutate(geneset = case_when(
      genes %in% selected_gene_sets$ISG_core    ~ "ISG_core",
      genes %in% selected_gene_sets$ISGs        ~ "ISGs",
      genes %in% selected_gene_sets$mTORC1      ~ "mTORC1",
      genes %in% selected_gene_sets$Cholesterol ~ "Cholesterol",
      genes %in% selected_gene_sets$Housekeeping ~ "Housekeeping",
      TRUE ~ "Other"
    )) %>%
    group_by(genes) %>%
    mutate(zscore = (Expression - mean(Expression, na.rm = TRUE)) / sd(Expression, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(dataset = dataset_name, organism = species)
  return(data)
}

selected_sets <- c("ISGs", "ISG_core", "mTORC1", "Housekeeping", "Cholesterol")

all_zscore_data <- bind_rows(lapply(names(datasets), function(ds) {
  data <- read_rds(datasets[[ds]])
  species <- species_map[[ds]]
  compute_dataset_zscores(data, ds, species, mouse_genes, selected_sets)
}))

# ---------------------------
# Refine tissue categories
# ---------------------------
tissue_levels <- sort(unique(all_zscore_data$tissue))
in_vivo_levels <- grep("^in\\.vivo", tissue_levels, value = TRUE, ignore.case = TRUE)
ex_vivo_levels <- grep("^ex\\.vivo", tissue_levels, value = TRUE, ignore.case = TRUE)
other_levels   <- setdiff(tissue_levels, c(in_vivo_levels, ex_vivo_levels))
# ---------------------------
# Generate visually distinct color palette for tissue subtypes
# ---------------------------

library(colorspace)
library(scales)

# Base colors
base_colors <- c(
  "ex.vivo" = "#6a3d9aff",  # purple–blue base
  "in.vivo" = "#d38d5fff"   # orange–yellow base
)

# Get unique tissue levels
tissue_levels <- sort(unique(all_zscore_data$tissue))

# Separate them by family
in_vivo_levels <- grep("^in\\.vivo", tissue_levels, value = TRUE, ignore.case = TRUE)
ex_vivo_levels <- grep("^ex\\.vivo", tissue_levels, value = TRUE, ignore.case = TRUE)
other_levels   <- setdiff(tissue_levels, c(in_vivo_levels, ex_vivo_levels))

# Function: generate evenly spaced distinct colors in a hue range
generate_distinct_colors <- function(base, n, hue_shift = 20) {
  if (n == 1) return(base)
  base_hcl <- coords(as(hex2RGB(base), "polarLUV"))
  base_hue <- base_hcl[,"H"]
  hues <- seq(base_hue - hue_shift, base_hue + hue_shift, length.out = n)
  hcl(h = hues, c = 80, l = 65)
}

# Distinct shades for each group
in_vivo_colors <- generate_distinct_colors(base_colors["in.vivo"], length(in_vivo_levels), hue_shift = 40)
ex_vivo_colors <- generate_distinct_colors(base_colors["ex.vivo"], length(ex_vivo_levels), hue_shift = 60)

# Optional fallback colors for others using a valid colorspace palette
if(length(other_levels) > 0){
  other_colors <- qualitative_hcl(length(other_levels), palette = "Dark 3")
  names(other_colors) <- other_levels
} else {
  other_colors <- c()
}

# Assign names
names(in_vivo_colors) <- in_vivo_levels
names(ex_vivo_colors) <- ex_vivo_levels

# Combine all
condition_colors_detailed <- c(in_vivo_colors, ex_vivo_colors, other_colors)

# Apply to data
all_zscore_data <- all_zscore_data %>%
  mutate(tissue_label = factor(tissue, levels = tissue_levels))

# ---------------------------
# Box + Density plots
# ---------------------------
for(gs in unique(all_zscore_data$geneset)) {
  gs_data <- all_zscore_data %>% filter(geneset == gs)
  
  for(ds in unique(gs_data$dataset)) {
    ds_df <- gs_data %>% filter(dataset == ds)
    author_label <- unique(ds_df$author)
    if (length(author_label) == 0) author_label <- ds
    author_label <- stringr::str_wrap(author_label, 35)
    
    # sample order
    sample_order_median <- ds_df %>%
      group_by(tissue_label, sample) %>%
      summarise(median_z = median(zscore, na.rm = TRUE), .groups = "drop") %>%
      arrange(factor(tissue_label, levels = tissue_levels), desc(median_z)) %>%
      pull(sample)
    
    ds_df <- ds_df %>% mutate(sample_label = factor(sample, levels = sample_order_median))
    
    # boxplot
    p_box <- ggplot(ds_df, aes(x = sample_label, y = zscore, color = tissue_label)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5) +
      scale_color_manual(values = condition_colors_detailed, name = "Culture subtype") +
      optimized_theme_fig() +
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "bottom") +
      labs(y = "Z-score", x = NULL, caption = author_label)
    
    # density
    p_density <- ggplot(ds_df, aes(x = zscore, fill = tissue_label)) +
      geom_density(alpha = 0.6, adjust = 1.5) +
      coord_flip(expand = FALSE) +
      scale_fill_manual(values = condition_colors_detailed) +
      optimized_theme_fig() +
      theme(axis.text = element_blank(), legend.position = "none")
    
    combined_plot <- p_box + p_density +
      plot_layout(widths = c(4, 1), guides = "collect")
    
    ggsave(out(paste0("box_density_", ds, "_", gs, ".pdf")),
           combined_plot, width = 6, height = 2, units = "in")
  }
}

# ---------------------------
# Heatmaps
# ---------------------------
for(gs in unique(all_zscore_data$geneset)) {
  gs_data <- all_zscore_data %>% filter(geneset == gs)
  
  for(ds in unique(gs_data$dataset)) {
    ds_df <- gs_data %>% filter(dataset == ds)
    
    mat <- ds_df %>%
      dplyr::select(genes, sample, zscore) %>%
      pivot_wider(names_from = sample, values_from = zscore) %>%
      column_to_rownames("genes") %>%
      as.matrix()
    
    mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
    mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
    if (nrow(mat) < 3 || ncol(mat) < 3) next
    
    ann_data <- ds_df %>%
      distinct(sample, tissue_label) %>%
      filter(sample %in% colnames(mat)) %>%
      arrange(match(sample, colnames(mat)))
    
    ha <- HeatmapAnnotation(
      Condition = ann_data$tissue_label,
      height = unit(0.2, "cm"),
      col = list(Condition = condition_colors_detailed)
    )
    
    legend_list <- list(
      Legend(title = "Culture Type",
             labels = c("in vivo", "ex vivo"),
             legend_gp = gpar(fill = c(base_colors["in.vivo"], base_colors["ex.vivo"])))
    )
    
    pdf(out(paste0("Heatmap_", ds, "_", gs, ".pdf")),
        width = max(4, ncol(mat) * 0.2), height = max(4, nrow(mat) * 0.2))
    print(
      Heatmap(
        mat, name = "Z-score",
        col = colorRamp2(c(-2, 0, 2), c("#4C889C", "white", "#D0154E")),
        top_annotation = ha,
        cluster_rows = TRUE, cluster_columns = FALSE,
        show_column_names = FALSE, show_row_names = FALSE
      ),
      annotation_legend_list = legend_list,
      merge_legends = TRUE
    )
    dev.off()
  }
}

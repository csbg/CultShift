##############################
# Modular GEO RNA-seq pipeline
##############################
#Aarathy
###############
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
InDir <- dirout("Ag_top_pathway_genes")

out <- dirout("Ag_external_datasets_generalized")
##############
genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))
genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

# Original gene sets
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)
Housekeeping <- c("Tbp", "Ubc")
#
# Convert to uppercase for human
genesets <- list(
  ISGs        = ISGs,
  ISG_core    = ISG_core,
  mTORC1      = mTORC1,
  Cholesterol = Cholesterol,
  Housekeeping = Housekeeping
)
# Example mouse genesets
mouse_genes <- list(
  ISGs        = ISGs,
  ISG_core    = ISG_core,
  mTORC1      = mTORC1,
  Cholesterol = Cholesterol,
  Housekeeping = Housekeeping
)
ENRICHR <- dirout(paste0("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/ENRICHR"))
ENRICHR.DBS <-union(ENRICHR.DBS,
                    c("GO_Biological_Process_2021",
                      "TRRUST_Transcription_Factors_2019",
                      "Reactome_2022",
                      "GO_Molecular_Function_2023",
                      "GO_Biological_Process_2023",
                      "CellMarker_2024"))
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
# # save(enr.terms, file=out("Genesets_Human.RData"))
# Convert to mouse --------------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
# Map mouse genes to human using 'hm' table
# Map mouse genes to human, keep unmapped genes
map_mouse_to_human_keep_all <- function(mouse_genes_list, hm) {
  lapply(mouse_genes_list, function(genes) {
    human_genes <- hm$Human[match(genes, hm$Mouse)]
    
    # For unmapped genes, use original mouse gene
    unmapped_genes <- genes[is.na(human_genes)]
    human_genes[is.na(human_genes)] <- unmapped_genes
    
    # Uppercase and unique
    toupper(unique(human_genes))
  })
}

# Apply the function
human_genesets <- map_mouse_to_human_keep_all(mouse_genes, hm)




# humnan--------------------------
InDir1 <- dirout("GSE229032_AML_primary_vs_cell_line")

"AML_NTCs_dataVoom.rds"

InDir2 <- dirout("GSE125213_human_CB1")
"CB1_NTCs_dataVoom.rds"


InDir3 <- dirout("GSE110968_CB2_external_dataset")
"CB_NTCs_dataVoom.rds"


InDir4 <- dirout("GSE169368_invitro_organoids_vs_invivo_mouse")
"GSE169368_mouse_organoids.rds"

InDir5 <- dirout("GSE149244vsGSE125211_LT-HSC")
"LT-HSC_mouse.rds.rds"

InDir6 <- dirout("GSE139184_pancreatic_cancer")

#################
datasets <- list(
  AML_primary_vs_cell_line = InDir1("AML_NTCs_dataVoom.rds"),
  CB1_human                = InDir2("CB1_NTCs_dataVoom.rds"),
  CB2_human                = InDir3("CB_NTCs_dataVoom.rds"),
  PDAC_human               = InDir6("GSE139184_PDAC.rds"),
  Organoids_mouse          = InDir4("GSE169368_mouse_organoids.rds"),
  LT_HSC_mouse             = InDir5("LT-HSC_mouse.rds.rds")
)

# ---- Define which species each dataset is ----
species_map <- list(
  AML_primary_vs_cell_line = "human",
  CB1_human                = "human",
  CB2_human                = "human",
  PDAC_human               = "human",
  Organoids_mouse          = "mouse",
  LT_HSC_mouse             = "mouse"
)

data_list <- lapply(datasets, read_rds)

##############
# ===========================
# STEP 1: Precompute z-scores per dataset
# ===========================

compute_dataset_zscores <- function(data, dataset_name, species, gene_sets, selected_sets) {
  # 1️⃣ Select only the gene sets we care about
  selected_gene_sets <- gene_sets[selected_sets]
  
  # 2️⃣ Adjust for human datasets only
  if (species == "human") {
    data <- data %>% mutate(genes = toupper(genes))
    selected_gene_sets <- lapply(selected_gene_sets, toupper)
  }
  
  # 3️⃣ Get all genes of interest
  genes_of_interest <- unique(unlist(selected_gene_sets))
  
  # 4️⃣ Filter data for genes of interest
  data <- data %>% filter(genes %in% genes_of_interest)
  
  # 5️⃣ Assign each gene to its geneset
  data <- data %>%
    mutate(geneset = case_when(
      genes %in% selected_gene_sets$ISG_core    ~ "ISG_core",
      genes %in% selected_gene_sets$ISGs        ~ "ISGs",
      genes %in% selected_gene_sets$mTORC1      ~ "mTORC1",
      genes %in% selected_gene_sets$Cholesterol ~ "Cholesterol",
      genes %in% selected_gene_sets$Housekeeping ~ "Housekeeping",
      TRUE ~ "Other"
    ))
  
  # 6️⃣ Compute z-scores per gene
  data_z <- data %>%
    group_by(genes) %>%
    mutate(zscore = (Expression - mean(Expression, na.rm = TRUE)) / sd(Expression, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      dataset = dataset_name,
      organism = species
    )
  
  return(data_z)
}

# ===========================
# STEP 2: Loop over all datasets and combine the zscores calculated within datasets
# ===========================

selected_sets <- c("ISGs", "ISG_core", "mTORC1", "Housekeeping", "Cholesterol")

all_zscore_data <- bind_rows(
  lapply(names(datasets), function(ds) {
    data <- read_rds(datasets[[ds]])
    species <- species_map[[ds]]
    compute_dataset_zscores(data, ds, species, genesets, selected_sets)
  })
)
unique(all_zscore_data$tissue)
unique(all_zscore_data$organism)
unique(all_zscore_data$dataset)
unique(all_zscore_data$geneset)
#####################
#plots- box and density------------
# Make sure tissue is a factor with consistent levels
all_zscore_data <- all_zscore_data %>%
  mutate(
    #tissue = factor(tissue, levels = c("in.vivo", "ex.vivo")),
    dataset = factor(dataset, levels = unique(dataset)),
    geneset = factor(geneset, levels = c("ISGs", "ISG_core", "mTORC1", "Cholesterol", "Housekeeping"))
  )
all_zscore_data <- all_zscore_data %>%
  mutate(
    tissue_label = case_when(
      grepl("primary|in\\.vivo", as.character(tissue), ignore.case = TRUE) ~ "in.vivo",
      grepl("ex|CL|culture", as.character(tissue), ignore.case = TRUE)     ~ "ex.vivo",
      TRUE ~ "other"
    ),
    tissue_label = factor(tissue_label, levels = c("in.vivo", "ex.vivo"))
  )
unique(all_zscore_data$tissue)
# Plot
# Define colors
condition_colors <- c(
  "ex.vivo" = "#6a3d9aff", 
  "in.vivo" = "#d38d5fff"
)

# Make tissue a factor (order important)
all_zscore_data <- all_zscore_data %>%
  mutate(
    tissue_label = factor(tissue_label, levels = c("in.vivo", "ex.vivo")),
    dataset = factor(dataset, levels = unique(dataset))
  )

# List of genesets to plot
geneset_levels <- c("ISGs", "ISG_core", "mTORC1", "Housekeeping", "Cholesterol")

plots_list <- list()

# -----------------------------
# PARAMETERS
# -----------------------------
max_width_cm <- 18
#max_height_cm <- 5
density_panel_cm <- 1.5
box_width_cm <- 0.25
padding_cm <- 1  # extra space for legend/margins

max_width_in <- max_width_cm /3
#plot_height_in <- max_height_cm / 3
density_panel_in <- density_panel_cm / 3
box_width_in <- box_width_cm / 3
padding_in <- padding_cm / 3
fixed_plot_height_in <- 1.667
condition_colors <- c("ex.vivo"="#6a3d9aff","in.vivo"="#d38d5fff")

# -----------------------------
# LOOP OVER GENESSETS & DATASETS
# -----------------------------

col_fun <- colorRamp2(
  breaks = c(-2, 0, 2),
  colors = c("#4C889C",
             "white",
             "#D0154E")  # blue-white-red, colorblind-safe
)

for(gs in unique(all_zscore_data$geneset)){
  
  gs_data <- all_zscore_data %>% filter(geneset == gs)
  
  for(ds in unique(gs_data$dataset)){
    
    ds_df <- gs_data %>% filter(dataset == ds)

    author_label <- ds_df %>%
      pull(author) %>%
      unique()
    if (length(author_label) == 0) author_label <- ds
    author_label <- stringr::str_wrap(author_label, width = 35)
    # -----------------------------
    # Compute median z-score per sample within tissue
    # -----------------------------
    sample_order_median <- ds_df %>%
      group_by(tissue_label, sample) %>%
      summarise(median_z = median(zscore, na.rm = TRUE), .groups = "drop") %>%
      arrange(factor(tissue_label, levels = c("in.vivo","ex.vivo","other")),
              desc(median_z)) %>%
      pull(sample)
    
    ds_df <- ds_df %>%
      mutate(sample_label = factor(sample, levels = sample_order_median)) %>%
      arrange(factor(tissue_label, levels = c("in.vivo","ex.vivo","other")), sample_label)
    
    n_samples <- length(unique(ds_df$sample))
    n_genes <- length(unique(ds_df$genes))
    
    # -----------------------------
    # Compute total width (adjust if exceeds max)
    # -----------------------------
    total_width_in <- n_samples*box_width_in + density_panel_in + padding_in
    if(total_width_in > max_width_in){
      scale_factor <- max_width_in / total_width_in
      box_width_in <- box_width_in * scale_factor
      total_width_in <- max_width_in
    }
    
    # Dashed line index
    boundary_index <- sum(ds_df$tissue_label=="in.vivo")
    if(boundary_index==0 || boundary_index==n_samples) boundary_index <- NA
    
    # -----------------------------
    # Boxplot
    # -----------------------------
    p_box <- ggplot(ds_df, aes(x = sample_label, y = zscore, color = tissue_label)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5, position = position_dodge(width = 0.8)) 
    # +
    #   geom_jitter(aes(color = tissue_label), width = 0.15, size = 0.2, alpha = 0.2)      # geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.5, position = position_dodge(width = 0.8)) +
      # geom_jitter(aes(color = tissue_label), width = 0.15, size = 0.2, alpha = 0.2)
      # 
    if(!is.na(boundary_index)){
      p_box <- p_box + geom_vline(xintercept = boundary_index + 0.2, linetype = "dashed", color = "grey50")
    }
    
    p_box <- p_box +
      scale_color_manual(values = condition_colors, name = "culture model") +
      scale_x_discrete(drop = FALSE) +
      labs(x = NULL, y = "Z-score") +
      optimized_theme_fig() +
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "bottom")
    
    # -----------------------------
    # Density panel
    # -----------------------------
    p_density <- ggplot(ds_df, aes(x = zscore, fill = tissue_label)) +
      geom_density(alpha = 0.6, adjust = 1.5) +
      coord_flip(expand = FALSE) +
      scale_fill_manual(values = condition_colors) +
      labs(x = NULL, y = NULL) +
      optimized_theme_fig() +
      theme(
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0,0,0,0)
      )
    
    # -----------------------------
    # Combine box + density
    # -----------------------------
    # combined_plot <- p_box + p_density +
    #   plot_layout(widths = c(n_samples*box_width_in, density_panel_in), guides = "collect") &
    #   theme(legend.position = "bottom")
    # Remove internal spacing
    p_box <- p_box + theme(plot.margin = margin(0, 0, 0, 0))
    p_density <- p_density + theme(plot.margin = margin(0, 0, 0, 0))
    
    # Combine and add a small outer margin (e.g., 5 mm all around)
    # combined_plot <- p_box + p_density +
    #   plot_layout(
    #     widths = c(n_samples * box_width_in, density_panel_in),
    #     guides = "collect"
    #   ) &
    #   theme(
    #     legend.position = "bottom",
    #     plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "mm")  # ← small outer margin
    #   )+
    #   plot_annotation(
    #     caption = author_label  # acts as shared x-axis label
    #   )
    
    combined_plot <- p_box + p_density +
      plot_layout(
        widths = c(n_samples * box_width_in, density_panel_in),
        guides = "collect"
      ) +
      plot_annotation(
        caption = author_label,
        theme = theme(
          plot.caption = element_text(size = 5, hjust = 0.5)  # center caption
        )
      ) &
      theme(
        legend.position = "bottom",
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "mm")
      )
    
    
    # -----------------------------
    # Save box/density plot
    # -----------------------------
    ggsave(
      filename = out(paste0("box_density_", ds, "_", gs, ".pdf")),
      plot = combined_plot,
      width = total_width_in,
      height = fixed_plot_height_in,
      units = "in",
      limitsize = FALSE
    )
    
    
    
  }
}
#
#heatmap--------------
#
library(ComplexHeatmap)
library(circlize)
library(grid)

for(gs in unique(all_zscore_data$geneset)) {
  gs_data <- all_zscore_data %>% filter(geneset == gs)
  
  for(ds in unique(gs_data$dataset)) {
    ds_df <- gs_data %>% filter(dataset == ds)
    
    # pivot matrix
    mat <- ds_df %>%
      dplyr::select(genes, sample, zscore) %>%
      pivot_wider(names_from = sample, values_from = zscore) %>%
      column_to_rownames("genes") %>%
      as.matrix()
    
    # drop all-NA rows/cols
    mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
    mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
    
    # skip tiny matrices
    if(nrow(mat) < 3 || ncol(mat) < 3) next
    
    # sample ordering
    sample_order <- ds_df %>%
      group_by(tissue_label, sample) %>%
      summarise(median_z = median(zscore, na.rm = TRUE), .groups = "drop") %>%
      arrange(factor(tissue_label, levels = c("in.vivo","ex.vivo","other")), desc(median_z)) %>%
      pull(sample)
    sample_order <- intersect(sample_order, colnames(mat))
    mat <- mat[, sample_order, drop = FALSE]
    
    # Cap matrix values at -2/+2
    mat[mat > 2] <- 2
    mat[mat < -2] <- -2
    
    # annotation
    ann_data <- ds_df %>% distinct(sample, tissue_label) %>%
      filter(sample %in% colnames(mat)) %>%
      arrange(match(sample, colnames(mat)))
    
    ha <- HeatmapAnnotation(
      Condition = ann_data$tissue_label,
      height = unit(0.05, "cm"),  # <--- must be a unit object
      col = list(Condition = condition_colors)
    )
    
    
    # fixed color scale from -2 to 2
    col_fun <- colorRamp2(c(-2, 0, 2), c("#4C889C", "white", "#D0154E"))
    
    # numeric PDF width/height (inches)
    pdf_file <- out(paste0("Heatmap_", ds, "_", gs, ".pdf"))
    pdf(pdf_file, width = max(4, ncol(mat) * 0.2), height = max(4, nrow(mat) * 0.2))
    
    # Force evaluation of Heatmap
    print(
      Heatmap(
        mat,
        name = "Z-score",
        col = col_fun,
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = F,
        #row_names_gp = gpar(fontsize = 6),
        width = unit(0.2, "cm") * ncol(mat),
        height = unit(0.2, "cm") * nrow(mat),
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 6),
          labels_gp = gpar(fontsize = 5),
          legend_height = unit(1, "cm"),
          legend_width = unit(0.1, "cm")
    )))
    
    dev.off()
    message("Heatmap generated for ", ds, " / ", gs)
  }
}


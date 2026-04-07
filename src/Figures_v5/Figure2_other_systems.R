##############################
# Modular GEO RNA-seq pipeline
##############################
#Aarathy
###############
source("src/00_init.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_function_save_zscore_table.R")
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
InDir1 <- dirout("Fig.Glio")
out <- dirout("Figure2_other_systems_v8")

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
ROS <- get_genes_for_pathway("ROS", genesets)
ISGs_glioblastoma <- read_rds(InDir1("ISG_glioblastoma.rds"))
#InDir1()
# Convert to uppercase for human
genesets <- list(
  ISGs        = ISGs,
  ISG_core    = ISG_core,
  mTORC1      = mTORC1,
  ROS         = ROS,
  Cholesterol = Cholesterol,
  Housekeeping = Housekeeping,
  ISGs_glioblastoma = ISGs_glioblastoma
)
# Example mouse genesets
mouse_genes <- list(
  ISGs        = ISGs,
  ISG_core    = ISG_core,
  mTORC1      = mTORC1,
  ROS         = ROS,
  Cholesterol = Cholesterol,
  Housekeeping = Housekeeping,
  ISGs_glioblastoma = ISGs_glioblastoma
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

InDir7 <- dirout("Glioblastoma_limmaRes")
InDir8 <- dirout("GSE74516_humon_brain_tumor")
InDir9 <- dirout("Figure4")

#################
datasets <- list(
  AML_primary_vs_cell_line = InDir1("AML_NTCs_dataVoom.rds"),
  CB1_human                = InDir2("CB1_NTCs_dataVoom.rds"),
  CB2_human                = InDir3("CB_NTCs_dataVoom.rds"),
  PDAC_human               = InDir6("GSE139184_PDAC.rds"),
  Organoids_mouse          = InDir4("GSE169368_mouse_organoids.rds"),
  LT_HSC_mouse             = InDir5("LT-HSC_mouse.rds.rds"),
  Spleen_M                 = InDir9("Spleen_M.rds"),
  Spleen_T8                = InDir9("Spleen_T8.rds"),
  Glioblastoma             = InDir7("Glioblastoma_NTCs.rds"),
  Glioblastoma3565         = InDir8("GSE74516_GBM_primary.rds")
)

# ---- Define which species each dataset is ----
species_map <- list(
  AML_primary_vs_cell_line = "human",
  CB1_human                = "human",
  CB2_human                = "human",
  PDAC_human               = "human",
  Organoids_mouse          = "mouse",
  LT_HSC_mouse             = "mouse",
  Glioblastoma             = "mouse",
  Glioblastoma3565         = "human",
  Spleen_T8                = "mouse",
  Spleen_M                = "mouse"
  
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
      genes %in% selected_gene_sets$ROS         ~ "ROS",
      genes %in% selected_gene_sets$Cholesterol ~ "Cholesterol",
      genes %in% selected_gene_sets$Housekeeping ~ "Housekeeping",
      genes %in% selected_gene_sets$ISGs_glioblastoma ~ "ISGs_glioblastoma",
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
  # Write Excel file with one sheet per geneset
  write_zscore_tables(data_z,paste0(out(dataset_name, "_geneset_tables.xlsx")))
  
  return(data_z)
}

# ===========================
# STEP 2: Loop over all datasets and combine the zscores calculated within datasets
# ===========================

selected_sets <- c("ISGs", "ISG_core", "mTORC1",
                   "ROS","Housekeeping", "Cholesterol", "ISGs_glioblastoma")

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
    geneset = factor(geneset, levels = c("ISGs", "ISG_core", "mTORC1", "Cholesterol",
                                         "ROS","Housekeeping","ISGs_glioblastoma"))
  )

all_zscore_data <- all_zscore_data %>%
  mutate(
    tissue_label = case_when(
      grepl("primary|in\\.vivo|xenograft", as.character(tissue), ignore.case = TRUE) ~ "in.vivo",
      grepl("ex\\.vivo|CL|culture|organoid", as.character(tissue), ignore.case = TRUE)     ~ "ex.vivo",
      TRUE ~ "other"
    ),
    tissue_label = factor(tissue_label, levels = c("in.vivo", "ex.vivo"))
  )


# Make tissue a factor (order important)
all_zscore_data <- all_zscore_data %>%
  mutate(
    tissue_label = factor(tissue_label, levels = c("in.vivo", "ex.vivo")),
    dataset = factor(dataset, levels = unique(dataset))
  )

# List of genesets to plot
geneset_levels <- c("ISGs", "ISG_core", "mTORC1", "Housekeeping", "Cholesterol",
                    "ROS","ISGs_glioblastoma")

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
dataset_color_map[[ds]] <- c("ex.vivo"="#6a3d9aff","in.vivo"="#d38d5fff")
##############


# Define your manual palettes
in_vivo_palette <- c("#D38D5F", "#E5971E", "#F9C319", "#D4E062")
ex_vivo_palette <- c("#4A2C69", "#613488", "#882A73", "#AE77AF")

get_dataset_colors <- function(df) {
  tissue_levels <- unique(df$tissue)
  
  # Identify groups
  in_vivo_levels <- tissue_levels[grepl("in\\.vivo|primary|xenograft",
                                        tissue_levels, ignore.case = TRUE)]
  ex_vivo_levels <- tissue_levels[grepl("ex\\.vivo|CL|culture|organoid",
                                        tissue_levels, ignore.case = TRUE)]
  
  # Assign colors manually, recycling if needed
  in_vivo_colors <- setNames(
    in_vivo_palette[seq_len(min(length(in_vivo_levels), length(in_vivo_palette)))],
    in_vivo_levels[seq_len(min(length(in_vivo_levels), length(in_vivo_palette)))]
  )
  
  # If more in.vivo samples than colors, recycle the palette
  if(length(in_vivo_levels) > length(in_vivo_palette)) {
    extra <- length(in_vivo_levels) - length(in_vivo_palette)
    in_vivo_colors <- c(in_vivo_colors, setNames(in_vivo_palette[seq_len(extra)], in_vivo_levels[(length(in_vivo_palette)+1):length(in_vivo_levels)]))
  }
  
  # Same for ex.vivo
  ex_vivo_colors <- setNames(
    ex_vivo_palette[seq_len(min(length(ex_vivo_levels), length(ex_vivo_palette)))],
    ex_vivo_levels[seq_len(min(length(ex_vivo_levels), length(ex_vivo_palette)))]
  )
  
  if(length(ex_vivo_levels) > length(ex_vivo_palette)) {
    extra <- length(ex_vivo_levels) - length(ex_vivo_palette)
    ex_vivo_colors <- c(ex_vivo_colors, setNames(ex_vivo_palette[seq_len(extra)], ex_vivo_levels[(length(ex_vivo_palette)+1):length(ex_vivo_levels)]))
  }
  
  # Combine
  c(in_vivo_colors, ex_vivo_colors)
}

# Apply per dataset
dataset_color_map <- all_zscore_data %>%
  group_by(dataset) %>%
  summarise(colors = list(get_dataset_colors(cur_data_all()))) %>%
  deframe()

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
    # ---- Assign exact precomputed colors to each tissue ----
    if (!ds %in% names(dataset_color_map)) {
      stop(paste("Dataset", ds, "not found in dataset_color_map"))
    }
    
    # Map the exact dataset colors into ds_df
    ds_df <- ds_df %>%
      mutate(color = dataset_color_map[[ds]][as.character(tissue)])  # add explicit color hex
    
    # Verify
    print(unique(ds_df[, c("tissue", "color")]))
    

    author_label <- ds_df %>%
      pull(author) %>%
      unique()
    if (length(author_label) == 0) author_label <- ds
    author_label <- stringr::str_wrap(author_label, width = 35)
    # -----------------------------
    # Compute median z-score per sample within tissue
    # -----------------------------
    tissue_levels_all <- unique(ds_df$tissue)
    
    in_vivo_levels <- sort(tissue_levels_all[grepl("in\\.vivo|primary|xenograft",
                                                   tissue_levels_all, ignore.case = TRUE)])
    ex_vivo_levels <- sort(tissue_levels_all[grepl("ex\\.vivo|CL|culture",
                                                   tissue_levels_all, ignore.case = TRUE)])
    other_levels   <- setdiff(tissue_levels_all, c(in_vivo_levels, ex_vivo_levels))
    
    # Combine, ensuring unique levels and correct order
    tissue_levels_ordered <- unique(c(in_vivo_levels, ex_vivo_levels, other_levels))
    
    sample_order_median <- ds_df %>%
      group_by(tissue, sample) %>%
      summarise(median_z = median(zscore, na.rm = TRUE), .groups = "drop") %>%
      arrange(factor(tissue, levels = tissue_levels_ordered),
              desc(median_z)) %>%
      pull(sample)
    ds_df$tissue
    ds_df <- ds_df %>%
      mutate(sample_label = factor(sample, levels = sample_order_median)) %>%
      arrange(factor(tissue, levels = tissue_levels_ordered), sample_label)
    
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
    # BOX
    p_box <- ggplot(ds_df, aes(x = sample_label, y = zscore)) +
      geom_boxplot(aes(color = tissue), outlier.shape = NA, alpha = 0.8, width = 0.5) +
      labs(x = NULL, y = "Z-score") +
      optimized_theme_fig() +
      scale_fill_manual(values = dataset_color_map[[ds]] ) +
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none")  # legend is unnecessary now
    
    if(!is.na(boundary_index)){
      p_box <- p_box + geom_vline(xintercept = boundary_index + 0.2, linetype = "dashed", color = "grey50")
    }
    
    p_box <- p_box +
      scale_color_manual(values = dataset_color_map[[ds]], name = "culture model") +
      scale_x_discrete(drop = FALSE) +
      labs(x = NULL, y = "Z-score") +
      optimized_theme_fig() +
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "bottom")
    
    # -----------------------------
    # Density panel
    # -----------------------------
    # DENSITY
    p_density <- ggplot(ds_df, aes(x = zscore)) +
      geom_density(aes(fill = tissue), alpha = 0.6, adjust = 1.5) +
      coord_flip(expand = FALSE) +
      scale_fill_manual(values = dataset_color_map[[ds]] ) +
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
# library(ComplexHeatmap)
# library(circlize)
# library(grid)
# 
# for(gs in unique(all_zscore_data$geneset)) {
#   gs_data <- all_zscore_data %>% filter(geneset == gs)
#   
#   for(ds in unique(gs_data$dataset)) {
#     ds_df <- gs_data %>% filter(dataset == ds)
#     
#     # pivot matrix
#     mat <- ds_df %>%
#       dplyr::select(genes, sample, zscore) %>%
#       pivot_wider(names_from = sample, values_from = zscore) %>%
#       column_to_rownames("genes") %>%
#       as.matrix()
#     
#     # drop all-NA rows/cols
#     mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
#     mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
#     
#     # skip tiny matrices
#     if(nrow(mat) < 3 || ncol(mat) < 3) next
#     
#     tissue_levels_all <- unique(ds_df$tissue)
#     
#     in_vivo_levels <- sort(tissue_levels_all[grepl("in\\.vivo|primary|xenograft", tissue_levels_all, ignore.case = TRUE)])
#     ex_vivo_levels <- sort(tissue_levels_all[grepl("ex\\.vivo|CL|culture|organoid", tissue_levels_all, ignore.case = TRUE)])
#     other_levels   <- setdiff(tissue_levels_all, c(in_vivo_levels, ex_vivo_levels))
#     
#     # Combine, ensuring unique levels and correct order
#     tissue_levels_ordered <- unique(c(in_vivo_levels, ex_vivo_levels, other_levels))
#     
#     sample_order_median <- ds_df %>%
#       group_by(tissue, sample) %>%
#       summarise(median_z = median(zscore, na.rm = TRUE), .groups = "drop") %>%
#       arrange(factor(tissue, levels = tissue_levels_ordered),
#               desc(median_z)) %>%
#       pull(sample)
#     sample_order <- intersect(sample_order_median, colnames(mat))
#     mat <- mat[, sample_order, drop = FALSE]
#     
#     # Cap matrix values at -2/+2
#     mat[mat > 2] <- 2
#     mat[mat < -2] <- -2
#     
#     # annotation
#     ann_data <- ds_df %>% distinct(sample, tissue) %>%
#       filter(sample %in% colnames(mat)) %>%
#       arrange(match(sample, colnames(mat)))
#     
#     ha <- HeatmapAnnotation(
#       Condition = ann_data$tissue,
#       height = unit(0.05, "cm"),  # <--- must be a unit object
#       col = list(Condition = dataset_color_map[[ds]])
#     )
#     
#     
#     # fixed color scale from -2 to 2
#     col_fun <- colorRamp2(c(-2, 0, 2), c("#4C889C", "white", "#D0154E"))
#     
#     # numeric PDF width/height (inches)
#     pdf_file <- out(paste0("Heatmap_", ds, "_", gs, ".pdf"))
#     pdf(pdf_file, width = max(4, ncol(mat) * 0.2), height = max(4, nrow(mat) * 0.2))
#     
#     # Calculate plot size safely
#     heatmap_width_cm  <- max(0.2 * ncol(mat), 2)   # min 2 cm
#     heatmap_height_cm <- max(0.2 * nrow(mat), 2)   # min 2 cm
#     
#     print(
#       Heatmap(
#         mat,
#         name = "Z-score",
#         col = col_fun,
#         top_annotation = ha,
#         cluster_rows = F,
#         cluster_columns = FALSE,
#         show_column_names = FALSE,
#         show_row_names = T,
#         width = unit(heatmap_width_cm, "cm"),
#         height = unit(heatmap_height_cm, "cm"),
#         heatmap_legend_param = list(
#           title_gp = gpar(fontsize = 6),
#           labels_gp = gpar(fontsize = 5),
#           legend_height = unit(2, "cm"),
#           legend_width = unit(0.2, "cm")
#         )
#       )
#     )
#     
#     
#     dev.off()
#     message("Heatmap generated for ", ds, " / ", gs)
#   }
# }
# ds <- "LT_HSC_mouse"
# gs <- "ISG_core"
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
    ds_df$tissue <- factor(ds_df$tissue, levels = unique(ds_df$tissue))
    unique(ds_df$tissue)
    tissue_levels_all <- unique(ds_df$tissue)
    
    in_vivo_levels <- sort(tissue_levels_all[grepl("in\\.vivo|primary|xenograft", tissue_levels_all, ignore.case = TRUE)])
    ex_vivo_levels <- sort(tissue_levels_all[grepl("ex\\.vivo|CL|culture|organoid", tissue_levels_all, ignore.case = TRUE)])
    #other_levels   <- setdiff(tissue_levels_all, c(in_vivo_levels, ex_vivo_levels))
    
    # Combine, ensuring unique levels and correct order
    tissue_levels_ordered <- unique(c(in_vivo_levels, ex_vivo_levels))
    
    # Enforce tissue-level ordering: in.vivo → ex.vivo
    ds_df <- ds_df %>%
      mutate(
        tissue_group = case_when(
          grepl("in\\.vivo|primary|xenograft", tissue, ignore.case = TRUE) ~ "in.vivo",
          grepl("ex\\.vivo|CL|culture|organoid", tissue, ignore.case = TRUE) ~ "ex.vivo",
          TRUE ~ "other"
        ),
        tissue_group = factor(tissue_group, levels = c("in.vivo", "ex.vivo", "other"))
      )
    
    # Define sample order purely by tissue group
    sample_order <- ds_df %>%
      arrange(tissue_group, tissue, sample) %>%
      distinct(sample) %>%
      pull(sample)
    
    # Keep only samples that exist in the matrix
    sample_order <- intersect(sample_order, colnames(mat))
    mat <- mat[, sample_order, drop = FALSE]
    
    # Match annotation order
    ann_data <- ds_df %>%
      distinct(sample, tissue) %>%
      filter(sample %in% colnames(mat)) %>%
      arrange(match(sample, colnames(mat)))
    
    mat <- mat[, sample_order, drop = FALSE]
    
    # Cap matrix values at -2/+2
    mat[mat > 2] <- 2
    mat[mat < -2] <- -2
    
    # annotation
    
    
    # fixed color scale from -2 to 2
    col_fun <- colorRamp2(c(-2, 0, 2), c("#4C889C", "white", "#D0154E"))
    
    # numeric PDF width/height (inches)
    pdf_file <- out(paste0("Heatmap_", ds, "_", gs, "WITHLABEL.pdf"))
    pdf(pdf_file, width = max(7, ncol(mat) * 0.1), height = max(7, nrow(mat) * 0.1))
    
    # Calculate plot size safely
    heatmap_width_cm  <- max(0.1 * ncol(mat), 2)   # min 2 cm
    heatmap_height_cm <- max(0.1 * nrow(mat), 2)   # min 2 cm
    
  
    ha <- HeatmapAnnotation(
      Condition = ann_data$tissue,
      show_legend = TRUE,
      height = unit(0.2, "cm"),               # bar thickness
      col = list(Condition = dataset_color_map[[ds]]),
      simple_anno_size = unit(0.2, "cm"),     # bar tile thickness
      annotation_name_gp = gpar(fontsize = 5),# annotation name font (on top of bar, optional)
      annotation_legend_param = list(
        title = NULL,                          # remove legend title if you want
        labels_gp = gpar(fontsize = 5, fontfamily = "sans"),  # smaller text
        legend_height = unit(1, "cm"),         # smaller overall height
        legend_width  = unit(0.1, "cm")        # smaller color squares
      )
    )
    
    heatmap_legend_param = list(
      title_gp  = gpar(fontsize = 6),                     # keep title as is
      labels_gp = gpar(fontsize = 5, fontfamily = "sans"), # legend labels
      legend_height = unit(1.5, "cm"),
      legend_width  = unit(0.15, "cm")
    )
    
    print(
      Heatmap(
        mat,
        name = "Z-score",
        col = col_fun,
        top_annotation = ha,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 5, fontfamily = "sans"), # <-- font size 5
        row_names_side = "left" ,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        width = unit(0.25 * ncol(mat), "cm"),  # slightly wider to fit labels
        height = unit(0.2 * nrow(mat), "cm"),  # keep height same
        heatmap_legend_param = list(
          title_gp  = gpar(fontsize = 6),
          labels_gp = gpar(fontsize = 5, fontfamily = "sans"),
          legend_height = unit(1.5, "cm"),
          legend_width  = unit(0.15, "cm")
        ),
        cell_fun = function(j, i, x, y, w, h, fill) {
          grid.rect(x = x, y = y, width = min(w, h), height = min(w, h),
                    gp = gpar(col = NA, fill = fill))
        }
      )
    )
    
    
    
    dev.off()
    message("Heatmap generated for ", ds, " / ", gs)
  }
}
##########
#density plot alone-----------
#########

library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

for(gs in unique(all_zscore_data$geneset)){
  
  gs_data <- all_zscore_data %>% filter(geneset == gs)
  
  for(ds in unique(gs_data$dataset)){
    
    ds_df <- gs_data %>% filter(dataset == ds)
    
    # ---- Assign exact precomputed colors to each tissue ----
    if (!ds %in% names(dataset_color_map)) {
      stop(paste("Dataset", ds, "not found in dataset_color_map"))
    }
    
    ds_df <- ds_df %>%
      mutate(color = dataset_color_map[[ds]][as.character(tissue)])
    
    # Author label
    author_label <- ds_df %>% pull(author) %>% unique()
    if(length(author_label)==0) author_label <- ds
    author_label <- str_wrap(author_label, width = 35)
    
    # ---- Sample order within tissues ----
    tissue_levels_all <- unique(ds_df$tissue)
    in_vivo_levels <- sort(tissue_levels_all[grepl("in\\.vivo|primary|xenograft", tissue_levels_all, ignore.case=TRUE)])
    ex_vivo_levels <- sort(tissue_levels_all[grepl("ex\\.vivo|CL|culture", tissue_levels_all, ignore.case=TRUE)])
    other_levels   <- setdiff(tissue_levels_all, c(in_vivo_levels, ex_vivo_levels))
    tissue_levels_ordered <- unique(c(in_vivo_levels, ex_vivo_levels, other_levels))
    
    sample_order_median <- ds_df %>%
      group_by(tissue, sample) %>%
      summarise(median_z = median(zscore, na.rm=TRUE), .groups="drop") %>%
      arrange(factor(tissue, levels = tissue_levels_ordered),
              desc(median_z)) %>%
      pull(sample)
    
    ds_df <- ds_df %>%
      mutate(sample_label = factor(sample, levels = sample_order_median)) %>%
      arrange(factor(tissue, levels=tissue_levels_ordered), sample_label)
    
    # ---- Median per tissue ----
    median_df <- ds_df %>%
      group_by(tissue) %>%
      summarise(median_z = median(zscore, na.rm=TRUE), .groups="drop")
    
    
    
    # Your density panel
    p_density_panel <- ggplot(ds_df, aes(x = zscore)) +
      geom_density(aes(color = tissue), alpha=0.6, adjust=1.5) +
      geom_vline(
        data = median_df,
        aes(xintercept = median_z, color = tissue),
        linetype = "dashed",
        size = 0.4,
        show.legend = FALSE
      ) +
      scale_color_manual(values = dataset_color_map[[ds]]) +
      labs(x="Z-score", y="Density") +
      optimized_theme_fig() +
      theme(
        legend.position = "bottom",
        plot.margin = margin(0,0,0,0, "mm")
      )
    
    # Empty padding around the plot to fix width
    padding_left <- ggplot() + theme_void()
    padding_right <- ggplot() + theme_void()
    empty_space <- ggplot() + theme_void()
    
    # Combine horizontally to center panel
    panel_centered <- padding_left + p_density_panel + padding_right + 
      plot_layout(widths = c(2, 3, 2)) # 4 units for plot panel, 2 units on sides
    
    # Stack vertically: panel on top, space below
    final_plot <- panel_centered / empty_space +
      plot_layout(heights = c(2, 9)) +  # small panel, big space below
      plot_annotation(
        caption = author_label,
        theme = theme(plot.caption = element_text(size=6, hjust=0.5))
      )
    
    # Save PDF
    ggsave(
      filename = out(paste0("density_only_", ds, "_", gs, ".pdf")),
      plot = final_plot,
      width = 8,   # total PDF width
      height = 12, # total PDF height
      units = "cm",
      limitsize = FALSE
    )
    
    
    # ---- Ensure directory exists ----
    dir.create(dirname(out("dummy")), recursive = TRUE, showWarnings = FALSE)
    
    # ---- Save PDF ----
    ggsave(
      filename = out(paste0("density_only_", ds, "_", gs, ".pdf")),
      plot = final_plot,
      width = 8,   # total PDF width in cm
      height = 12, # total PDF height in cm
      units = "cm",
      limitsize = FALSE
    )
    
  }
}

##########
#density plot alone FILLED------
#########

library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

library(patchwork)
library(cowplot)

for(gs in unique(all_zscore_data$geneset)){
  
  gs_data <- all_zscore_data %>% filter(geneset == gs)
  
  for(ds in unique(gs_data$dataset)){
    
    ds_df <- gs_data %>% filter(dataset == ds)
    
    # ---- Assign exact precomputed colors to each tissue ----
    if (!ds %in% names(dataset_color_map)) {
      stop(paste("Dataset", ds, "not found in dataset_color_map"))
    }
    
    ds_df <- ds_df %>%
      mutate(color = dataset_color_map[[ds]][as.character(tissue)])
    
    # Author label
    author_label <- ds_df %>% pull(author) %>% unique()
    if(length(author_label)==0) author_label <- ds
    author_label <- str_wrap(author_label, width = 35)
    
    # ---- Sample order within tissues ----
    tissue_levels_all <- unique(ds_df$tissue)
    in_vivo_levels <- sort(tissue_levels_all[grepl("in\\.vivo|primary|xenograft", tissue_levels_all, ignore.case=TRUE)])
    ex_vivo_levels <- sort(tissue_levels_all[grepl("ex\\.vivo|CL|culture", tissue_levels_all, ignore.case=TRUE)])
    other_levels   <- setdiff(tissue_levels_all, c(in_vivo_levels, ex_vivo_levels))
    tissue_levels_ordered <- unique(c(in_vivo_levels, ex_vivo_levels, other_levels))
    
    sample_order_median <- ds_df %>%
      group_by(tissue, sample) %>%
      summarise(median_z = median(zscore, na.rm=TRUE), .groups="drop") %>%
      arrange(factor(tissue, levels = tissue_levels_ordered),
              desc(median_z)) %>%
      pull(sample)
    
    ds_df <- ds_df %>%
      mutate(sample_label = factor(sample, levels = sample_order_median)) %>%
      arrange(factor(tissue, levels=tissue_levels_ordered), sample_label)
    
    # ---- Median per tissue ----
    median_df <- ds_df %>%
      group_by(tissue) %>%
      summarise(median_z = median(zscore, na.rm=TRUE), .groups="drop")
    
    # ---- Density panel (original) ----
    p_density_panel <- ggplot(ds_df, aes(x = zscore)) +
      geom_density(aes(fill = tissue, color = tissue), alpha=0.6, adjust=1.5) +
      geom_vline(
        data = median_df,
        aes(xintercept = median_z, color = tissue),
        linetype = "dashed",
        size = 0.4,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = dataset_color_map[[ds]]) +
      scale_color_manual(values = dataset_color_map[[ds]]) +
      labs(x="Z-score", y="Density") +
      optimized_theme_fig() +
      theme(
        legend.position = "right",
        legend.box = "vertical", 
        legend.box.spacing = unit(0, "mm"),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(0.7, "mm"),
        panel.grid.major =  element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0,0,0,0, "mm"),
        axis.text.x  = element_text(angle = 0)
      )
    
    # ---- FIX: Extract legend so it does not resize the plot ----
    legend_density <- cowplot::get_legend(
      p_density_panel + theme(legend.position = "right")
    )
    
    p_density_nolegend <- p_density_panel +
      theme(legend.position = "none")
    
    # ---- Padding panels ----
    padding_left <- ggplot() + theme_void()
    padding_right <- ggplot() + theme_void()
    empty_space <- ggplot() + theme_void()
    
    # ---- Fixed-width plot panel ----
    panel_centered <- padding_left + p_density_nolegend + padding_right + 
      plot_layout(widths = c(2, 1.5, 2))
    
    # ---- Add legend BELOW while keeping plot width constant ----
    final_plot <- panel_centered / wrap_elements(full = legend_density) +
      plot_layout(heights = c(1, 3)) +   # REDUCED HEIGHT BALANCE
      plot_annotation(
        caption = author_label,
        theme = theme(plot.caption = element_text(size=6, hjust=0.5))
      )
    
    # ---- Save PDF WITH REDUCED HEIGHT ----
    ggsave(
      filename = out(paste0("density_only_FILLED", ds, "_", gs, ".pdf")),
      plot = final_plot,
      width = 12,
      height = 8,     # REDUCED HEIGHT
      units = "cm",
      limitsize = FALSE
    )
    
  }
}


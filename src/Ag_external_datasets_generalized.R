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

# Convert to uppercase for human
human_genesets <- list(
  ISGs        = toupper(ISGs),
  ISG_core    = toupper(ISG_core),
  mTORC1      = toupper(mTORC1),
  Cholesterol = toupper(Cholesterol),
  Housekeeping = toupper(Housekeeping)
)

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

InDir6 <-""
#################
datasets <- list(
  AML_primary_vs_cell_line = InDir1("AML_NTCs_dataVoom.rds"),
  CB1_human                = InDir2("CB1_NTCs_dataVoom.rds"),
  CB2_human                = InDir3("CB_NTCs_dataVoom.rds"),
  Organoids_mouse          = InDir4("GSE169368_mouse_organoids.rds"),
  LT_HSC_mouse             = InDir5("LT-HSC_mouse.rds.rds")
)

# ---- Define which species each dataset is ----
species_map <- list(
  AML_primary_vs_cell_line = "human",
  CB1_human                = "human",
  CB2_human                = "human",
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
    geneset = factor(geneset, levels = c("ISGs", "ISG_core", "mTORC1", "Housekeeping"))
  )
all_zscore_data <- all_zscore_data %>%
  mutate(
    tissue_label = case_when(
      grepl("primary|in\\.vivo", as.character(tissue), ignore.case = TRUE) ~ "in.vivo",
      grepl("ex|CL|culture", as.character(tissue), ignore.case = TRUE)     ~ "ex.vivo",
      TRUE ~ "other"
    ),
    tissue_label = factor(tissue_label, levels = c("in.vivo", "ex.vivo", "other"))
  )

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
max_height_cm <- 8
density_panel_cm <- 2
box_width_cm <- 0.25
padding_cm <- 1  # extra space for legend/margins

max_width_in <- max_width_cm / 2.54
plot_height_in <- max_height_cm / 2.54
density_panel_in <- density_panel_cm / 2.54
box_width_in <- box_width_cm / 2.54
padding_in <- padding_cm / 2.54

condition_colors <- c("ex.vivo"="#6a3d9aff","in.vivo"="#d38d5fff")

# -----------------------------
# LOOP OVER GENESSETS & DATASETS
# -----------------------------
gs <- "ISGs"
ds <- "LT_HSC_mouse"
col_fun <- colorRamp2(
  breaks = c(-2, 0, 2),
  colors = c("#4C889C",
             "white",
             "#D0154E")  # blue-white-red, colorblind-safe
)

for(gs in unique(all_zscore_data$geneset)){
  
  gs_data <- all_zscore_data %>% filter(geneset==gs)
  
  for(ds in unique(gs_data$dataset)){
    
    ds_df <- gs_data %>% filter(dataset==ds) %>%
      arrange(factor(tissue_label, levels=c("in.vivo","ex.vivo","other")), sample) %>%
      mutate(sample_label=factor(sample, levels=unique(sample)))
    
    n_samples <- length(unique(ds_df$sample))
    
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
    p_box <- ggplot(ds_df, aes(x=sample_label, y=zscore, color=tissue_label)) +
      geom_boxplot(outlier.shape=NA, alpha=0.8, width=0.6, position=position_dodge(width=0.8)) +
      geom_jitter(aes(color=tissue_label), width=0.15, size=0.2, alpha=0.2)
    
    if(!is.na(boundary_index)){
      p_box <- p_box + geom_vline(xintercept=boundary_index+0.5, linetype="dashed", color="grey50")
    }
    
    p_box <- p_box +
      scale_color_manual(values=condition_colors, name="Condition") +
      scale_x_discrete(drop=FALSE) +
      labs(x=NULL, y="Z-score") +
      optimized_theme_fig() +
      theme(axis.text.x=element_blank(),
            legend.position="bottom")
    
    # -----------------------------
    # Density panel
    # -----------------------------
    p_density <- ggplot(ds_df, aes(x=zscore, fill=tissue_label)) +
      geom_density(alpha=0.6, adjust=1.5) +
      coord_flip(expand=FALSE) +
      scale_fill_manual(values=condition_colors) +
      labs(x=NULL, y=NULL) +
      optimized_theme_fig() +
      theme(
        legend.position="none",
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        plot.margin=margin(0,0,0,0)
      )
    
    # -----------------------------
    # Combine box + density
    # -----------------------------
    combined_plot <- p_box + p_density +
      plot_layout(widths=c(n_samples*box_width_in, density_panel_in), guides="collect") &
      theme(legend.position="bottom")
    
    # -----------------------------
    # Save
    # -----------------------------
    ggsave(
      filename = out(paste0("box_density_", ds, "_", gs, ".png")),
      plot = combined_plot,
      width = total_width_in,
      height = plot_height_in,
      units="in",
      limitsize=FALSE
    )
    
    message("Generating heatmap for ", ds, " / ", gs)
    
    if (nrow(ds_df) > 1) {
      
      # Prepare z-score matrix
      mat <- ds_df %>%
        dplyr::select(genes, sample, zscore) %>%
        pivot_wider(names_from = sample, values_from = zscore) %>%
        column_to_rownames("genes") %>%
        as.matrix()
      
      # Order columns by tissue_label
      sample_order <- ds_df %>%
        distinct(sample, tissue_label) %>%
        arrange(factor(tissue_label, levels = c("in.vivo","ex.vivo","other"))) %>%
        pull(sample)
      mat <- mat[, sample_order, drop = FALSE]
      
      # Column annotation
      ha <- HeatmapAnnotation(
        Condition = ds_df %>%
          distinct(sample, tissue_label) %>%
          arrange(factor(tissue_label, levels = c("in.vivo","ex.vivo","other"))) %>%
          pull(tissue_label),
        col = list(Condition = condition_colors),
        annotation_legend_param = list(title = "Condition")
      )
      
      # Heatmap size
      n_genes <- nrow(mat)
      n_samples <- ncol(mat)
      cell_size_in <- 0.15
      heatmap_width <- unit(n_samples * cell_size_in, "in")
      heatmap_height <- unit(n_genes * cell_size_in, "in")
      
      #---------------------------------------------
      # 🔸 Draw Heatmap (Fixed)
      #---------------------------------------------
      ds_name <- ds
      gs_name <- gs
      pdf(
        out(paste0("Heatmap_", ds_name, "_", gs_name, ".pdf")),
        width = convertWidth(heatmap_width, "in", valueOnly = TRUE),
        height = convertHeight(heatmap_height, "in", valueOnly = TRUE)
      )
      pdf(
        out(paste0("Heatmap_", ds_name, "_", gs_name, ".pdf")),
        width = 20,
        height = 20
      )
      Heatmap(
        mat,
        name = "Z-score",
        col = col_fun,
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 6),
        width = heatmap_width,
        height = heatmap_height,
        cell_fun = function(j, i, x, y, w, h, fill) {
          grid.rect(x = x, y = y, width = unit(cell_size_in, "in"),
                    height = unit(cell_size_in, "in"),
                    gp = gpar(col = NA, fill = fill))
        }
      )
      
     
      dev.off()
      
  }
}
}


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
library(biomaRt)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
library(GEOquery)
library(R.utils)
library(msigdbr)
library(fgsea)
library(readxl)
library(latex2exp)
library(Biobase)
InDir <- dirout("Ag_top_pathway_genes")
InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
InDir3 <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
out <- dirout("GSE169368_invitro_organoids_vs_invivo_mouse")
# -----------------------------
# USER PARAMETERS
# -----------------------------

##############
genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))
genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

# gene subsets
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)

# -----------------------------
# Download GEO data
# -----------------------------
# GEO bulk RNA-seq counts file (gzipped)
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE169368&format=file&file=GSE169368%5FProcessed%5FNormalizedCountsByDESeq2%5FbulkRNAseq%2Ecsv%2Egz"
dest <- "GSE169368_Processed_NormalizedCountsByDESeq2_bulkRNAseq.csv.gz"
download.file(url, dest, mode = "wb")



# fread automatically detects headers unless it’s confused by quotes, so:
tbl <- read.csv(dest, header = TRUE)
rownames(tbl) <- make.unique(tbl$X)

############
gse <- getGEO("GSE169368", GSEMatrix = TRUE)
# Access the GPL21103 platform
metadata <- gse[["GSE169368_series_matrix.txt.gz"]]

# Extract phenoData from gse[[1]] (AnnotatedDataFrame)
adf <- phenoData(gse[["GSE169368_series_matrix.txt.gz"]])
# Use pData() from Biobase explicitly
metadata <- Biobase::pData(adf)
rownames(metadata) <- metadata$title
metadata <- metadata %>%
  filter(`age:ch1` == "young")%>%
  filter(`age:ch1` == "young")

metadata$tissue <- ifelse(metadata$`sample type:ch1` == "Organoids","ex.vivo",
                          "in.vivo")

# Keep only overlapping samples
# Original column names
tbl <- tbl[-1]
# Original metadata rownames
meta_rows <- rownames(metadata)

# Replace female/male
meta_rows_fixed <- meta_rows %>%
  gsub("female", "F", .) %>%
  gsub("male", "M", .) %>%
  gsub("_replication", "_rep", .) %>%
  gsub("_", "_", .)  # keep underscores intact

# Assign back
rownames(metadata) <- meta_rows_fixed
metadata$sample <- rownames(metadata)
# Check
head(rownames(metadata))


samples_keep <- intersect(rownames(metadata), colnames(tbl))
metadata <- metadata[samples_keep,]
metadata$tissue <- fct_relevel(metadata$tissue,ref = "in.vivo")
counts <- tbl[,samples_keep]

############################
design <- model.matrix(~tissue, data = metadata)
dge <- DGEList(counts)
dge <- calcNormFactors(dge, method = "TMM")
#design <- model.matrix(as.formula(design_formula), metadata)
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr,]
# -----------------------------

# voom transformation
dataVoom <- voom(dge, design, plot = TRUE)


# -----------------------------
# limma fit
# -----------------------------
limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)

# Extract results
limmaRes <- map_dfr(colnames(coef(limmaFit)), function(coefx) {
  topTable(limmaFit, coef = coefx, number = Inf) %>%
    rownames_to_column("genes") %>%
    filter(coefx != "(Intercept)") %>%
    mutate(coef = coefx,
           group = case_when(
             logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
             logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
             TRUE ~ "n.s"
           ))
})
unique(limmaRes$coef)
#
# -----------------------------
#saving D.E table----------------
coef_names <- setdiff(colnames(coef(limmaFit)), "(Intercept)")

limmaRes <- map_dfr(coef_names, function(coefx) {
  
  topTable(limmaFit, coef = coefx, number = Inf, sort.by = "none") %>%
    rownames_to_column("gene") %>%
    mutate(
      celltype = coefx,   # 👈 important for your export function
      group = case_when(
        logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
        logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
        TRUE ~ "n.s"
      )
    ) %>%
    dplyr::select(
      gene,
      celltype,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B,
      group
    )
})
export_by_celltype <- function(df, 
                               output_dir, 
                               output_file, 
                               sheet_columns = NULL, 
                               freeze_first_row = TRUE) {
  
  if(!is.null(sheet_columns)){
    df <- df[, intersect(sheet_columns, colnames(df)), drop = FALSE]
  }
  
  df <- as.data.frame(df)
  
  # 👇 KEY CHANGE
  cell_types <- unique(df$celltype)
  split_sheets <- length(cell_types) > 1
  
  wb <- createWorkbook()
  
  if(split_sheets){
    
    for(ct in cell_types){
      ann_ct <- df %>% filter(celltype == ct)
      
      if("adj.P.Val" %in% colnames(ann_ct)){
        ann_ct$adj.P.Val <- format_padj(ann_ct$adj.P.Val)
      }
      ann_ct <- format_numbers(ann_ct)
      
      sheet_name <- substr(gsub("[\\/:*?\\[\\]]", "_", ct), 1, 31)
      
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, ann_ct, rowNames = FALSE)
      
      freezePane(wb, sheet = sheet_name, firstRow = freeze_first_row)
      
      headerStyle <- createStyle(textDecoration = "bold")
      addStyle(wb, sheet = sheet_name, headerStyle,
               rows = 1, cols = 1:ncol(ann_ct), gridExpand = TRUE)
    }
    
  } else {
    # 👇 SINGLE SHEET MODE
    
    ann <- df
    
    if("adj.P.Val" %in% colnames(ann)){
      ann$adj.P.Val <- format_padj(ann$adj.P.Val)
    }
    ann <- format_numbers(ann)
    
    addWorksheet(wb, sheetName = "Results")
    writeData(wb, sheet = "Results", ann, rowNames = FALSE)
    
    freezePane(wb, sheet = "Results", firstRow = freeze_first_row)
    
    headerStyle <- createStyle(textDecoration = "bold")
    addStyle(wb, sheet = "Results", headerStyle,
             rows = 1, cols = 1:ncol(ann), gridExpand = TRUE)
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  saveWorkbook(wb, file = file.path(output_dir, output_file), overwrite = TRUE)
}
export_by_celltype(
  df = limmaRes,
  output_dir = out("DE_tables"),
  output_file = "Supplementary_Table4_DE_intestinal_organoids_limma.xlsx",
  sheet_columns = colnames(limmaRes)
)

###################################
#ranklist
coef_names <- setdiff(colnames(coef(limmaFit)), "(Intercept)")


# Directory for rank lists
dir.create(out("Ranklists"), showWarnings = FALSE)

ranklists <- map(coef_names, function(coefx) {
  
  df <- topTable(
    limmaFit,
    coef = coefx,
    number = Inf,
    sort.by = "none"
  ) %>%
    rownames_to_column("gene")
  
  # Named vector: logFC
  ranks <- df$logFC
  names(ranks) <- df$gene
  
  # Remove NA + duplicates
  ranks <- ranks[!is.na(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  
  # Sort decreasing (important for GSEA)
  ranks <- sort(ranks, decreasing = TRUE)
  
  # Save as RDS (best for R)
  saveRDS(
    ranks,
    file = file.path(out("Ranklists"), paste0(coefx, "_logFC_ranklist.rds"))
  )
  
  # Save as TXT (universal)
  write.table(
    data.frame(gene = names(ranks), logFC = ranks),
    file = file.path(out("Ranklists"), paste0(coefx, "_logFC_ranklist.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  ranks
})

names(ranklists) <- coef_names

#fgsea----------------

longer_dataVoom <-  counts %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(metadata, by ="sample")
longer_dataVoom$tissue <- factor(longer_dataVoom$tissue, 
                                 levels = c("in.vivo" ,"ex.vivo"))
longer_dataVoom$organ <- "intestine"
longer_dataVoom$author <- "Lu et al., Scientific Reports, 2021"

longer_dataVoom %>%
  dplyr::select(c("genes","tissue","sample","Expression","author","organ"))%>%
  write_rds(out("GSE169368_mouse_organoids.rds"))

################
prefix <- "mouse_organoids_vs_in_vivo"
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% ISGs) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  dplyr::select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }



# Make tissue a factor first
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, 
                                        levels = c("in.vivo" ,"ex.vivo"))

# Reorder samples so ex.vivo samples come first
longer_dataVoom_zscore$sample <- factor(
  longer_dataVoom_zscore$sample,
  levels = longer_dataVoom_zscore %>%
    arrange(tissue, sample) %>%
    pull(sample) %>%
    unique()
)




# Filter only genes in your ISGs list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% ISGs)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "Average ISG z-score per sample (colored by tissue)",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ggplot(plot_df, aes(x = tissue, y = zscore, fill = tissue)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
#   labs(
#     title = paste0("Z-score of ISGs across tissues"),
#     x = "Tissue",
#     y = "Z-score"
#   ) +
#   facet_wrap(~ genes, scales = "free_y") +  # one subplot per gene
#   optimized_theme_fig() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#########################
plot_df$sample <- factor(plot_df$sample, levels = unique(plot_df$sample))

# Plot heatmap
ggplot(plot_df, aes(x = sample, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Z-score of ISGs across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")


ggsave(out("per_gene_per_sample_heatmap_ISGs.png"))
########################################
# Define your gene sets (character vectors)
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)

longer_dataVoom <-  counts %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(metadata, by ="sample")
longer_dataVoom$tissue <- factor(longer_dataVoom$tissue, 
                                 levels = c("in.vivo" ,"ex.vivo"))
unique(longer_dataVoom$source_name_ch1)
longer_dataVoom %>%
 dplyr::select(c("genes",
                 "Expression",
                 "sample",
                 "tissue",
                 )

write_rds(out("GSE169368_longer_dataVoom.rds"))
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% ISG_core) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  dplyr::select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }
# Filter only genes in your mTORC1 list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% ISG_core)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "Average ISG_core z-score per sample (colored by tissue)",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ggplot(plot_df, aes(x = tissue, y = zscore, fill = tissue)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
#   labs(
#     title = paste0("Z-score of ISGs across tissues"),
#     x = "Tissue",
#     y = "Z-score"
#   ) +
#   facet_wrap(~ genes, scales = "free_y") +  # one subplot per gene
#   optimized_theme_fig() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#########################
plot_df$sample <- factor(plot_df$sample, levels = unique(plot_df$sample))

# Plot heatmap
ggplot(plot_df, aes(x = sample, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Z-score of ISGs across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")


ggsave(out("per_gene_per_sample_heatmap_ISG_core.png"))
##################

library(ggplot2)
library(patchwork)

# Boxplot
p_box <- ggplot(plot_df, aes(x = reorder(sample, zscore, FUN = median), y = zscore, fill = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1, alpha = 0.8) +
  labs(x = "Sample", y = "Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Density plot (flipped to match boxplot y-axis)
p_density <- ggplot(plot_df, aes(x = zscore, fill = tissue)) +
  geom_density(alpha = 0.6) +
  coord_flip(expand = FALSE) +  # remove extra space
  labs(x = NULL, y = NULL) +
  optimized_theme_fig() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0, 0, 0, 0)  # remove gap
  )

# Combine and collect the legend on the bottom
combined <- p_box + p_density + 
  plot_layout(widths = c(4, 1), guides = "collect") & 
  theme(legend.position = "bottom")

combined
ggsave(out("combined_plot.png"))
###################
#mTORC1
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% mTORC1) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  dplyr::select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }
# Filter only genes in your mTORC1 list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% mTORC1)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "Average mTORC1 z-score per sample (colored by tissue)",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ggplot(plot_df, aes(x = tissue, y = zscore, fill = tissue)) +
#   geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
#   labs(
#     title = paste0("Z-score of ISGs across tissues"),
#     x = "Tissue",
#     y = "Z-score"
#   ) +
#   facet_wrap(~ genes, scales = "free_y") +  # one subplot per gene
#   optimized_theme_fig() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#########################
plot_df$sample <- factor(plot_df$sample, levels = unique(plot_df$sample))

# Plot heatmap
ggplot(plot_df, aes(x = sample, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Z-score of mTORC1 genes across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")


ggsave(out("per_gene_per_sample_heatmap_mTORC1.png"))

longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% Cholesterol) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  dplyr::select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }
# Filter only genes in your mTORC1 list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% Cholesterol)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "Average Cholesterol z-score per sample (colored by tissue)",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#########################
plot_df$sample <- factor(plot_df$sample, levels = unique(plot_df$sample))

# Plot heatmap
ggplot(plot_df, aes(x = sample, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Z-score of cholesterol genes across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")


ggsave(out("per_gene_per_sample_heatmap_cholesterol.png"))
##############
Housekeeping <-  c("Gapdh","Tbp","Ubc","Actb","Pgk1")
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% Housekeeping) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  dplyr::select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }



# Make tissue a factor first
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, 
                                        levels = c("in.vivo" ,"ex.vivo"))

# Reorder samples so ex.vivo samples come first
longer_dataVoom_zscore$sample <- factor(
  longer_dataVoom_zscore$sample,
  levels = longer_dataVoom_zscore %>%
    arrange(tissue, sample) %>%
    pull(sample) %>%
    unique()
)




# Filter only genes in your ISGs list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% Housekeeping)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "housekeeping",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(ggplot2)
library(patchwork)

# Boxplot
p_box <- ggplot(plot_df, aes(x = reorder(sample, zscore, FUN = median), y = zscore, fill = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = tissue), width = 0.1, size = 1, alpha = 0.8) +
  labs(x = "Sample", y = "Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Density plot (flipped to match boxplot y-axis)
p_density <- ggplot(plot_df, aes(x = zscore, fill = tissue)) +
  geom_density(alpha = 0.6) +
  coord_flip(expand = FALSE) +  # remove extra space
  labs(x = NULL, y = NULL) +
  optimized_theme_fig() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(0, 0, 0, 0)  # remove gap
  )

# Combine and collect the legend on the bottom
combined <- p_box + p_density + 
  plot_layout(widths = c(4, 1), guides = "collect") & 
  theme(legend.position = "bottom")

combined
ggsave(out("combined_plot_housekeeping.png"))

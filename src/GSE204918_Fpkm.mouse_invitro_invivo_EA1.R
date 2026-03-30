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

out <- dirout("GSE204918_mouse_invitro_invivo")



# -----------------------------
# USER PARAMETERS
# -----------------------------
geo_id <- "GSE204918"           # GEO dataset
count_file <- "raw_counts_GRCh38.tsv.gz"  # raw counts filename in GEO
annot_file <- "annot_GRCh38.tsv.gz"      # gene annotation filename in GEO
design_formula <- "~tissue"     # design formula for limma/voom
tissue_levels <- c("in.vivo","ex.vivo")  # optional factor ordering
prefix <- "example"
series_matrix <- paste(geo_id,"_series_matrix.txt.gz")

##############
genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))
genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

#genesets$genes <- toupper(genesets$genes)
#

# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

# gene subsets
# gene subsets
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)

# -----------------------------
# Download GEO data
# -----------------------------
# DOWNLOAD GEO SUPPLEMENTARY FILES
# -----------------------------

# Metadata
gse <- getGEO(geo_id, GSEMatrix = TRUE)
metadata <- gse$GSE204918_series_matrix.txt.gz


# Extract phenoData from gse[[1]] (AnnotatedDataFrame)
adf <- phenoData(gse[[1]])

# Use pData() from Biobase explicitly
metadata <- Biobase::pData(adf)
#modify metadata
metadata$tissue <- metadata$`condition:ch1` %>%
  gsub("in vitro","ex.vivo",.) %>%
  gsub("in vivo","in.vivo",.)
metadata$tissue <- factor(metadata$tissue, levels = tissue_levels)
metadata <- metadata %>% select(title, geo_accession, tissue)

# Keep only overlapping samples
samples_keep <- intersect(rownames(metadata), colnames(tbl))
metadata <- metadata[samples_keep,]
counts <- tbl[,samples_keep]

# -----------------------------
# DGEList + normalization
# -----------------------------
dge <- DGEList(counts)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(as.formula(design_formula), metadata)
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr,]

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

# -----------------------------
# Pivot to long format
# -----------------------------
long_data <- as.data.frame(dataVoom$E) %>%
  rownames_to_column("genes") %>%
  pivot_longer(cols = -genes, names_to = "sample", values_to = "Expression") %>%
  inner_join(metadata, by = c("sample"))

# -----------------------------
# Z-score calculation for gene sets
# -----------------------------
compute_zscore <- function(df, genes) {
  df %>%
    filter(genes %in% genes) %>%
    group_by(genes) %>%
    mutate(
      mean_expr = mean(Expression, na.rm = TRUE),
      sd_expr   = sd(Expression, na.rm = TRUE),
      zscore    = (Expression - mean_expr)/sd_expr
    ) %>%
    ungroup()
}
#
gene_sets <- gene_sets <- list(
  "ISGs" = ISGs,
  "mTORC1" = mTORC1,
 " Cholesterol" = Cholesterol
)
# Apply to all gene sets
long_data_z <- map2_dfr(gene_sets, names(gene_sets), function(gs, pathway_name) {
  compute_zscore(long_data, gs) %>%
    mutate(pathway = pathway_name)
})

# -----------------------------
# Plots
# -----------------------------
# 1. Heatmap of genes x samples
ggplot(long_data_z, aes(x = title, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  facet_grid(pathway ~ tissue , scales = "free_x") +
  labs(title = "Gene z-scores by sample and tissue", x = "Sample", y = "Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Boxplot of z-scores per tissue
ggplot(long_data_z, aes(x = tissue, y = zscore, fill = tissue)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  facet_wrap(~ genes, scales = "free_y") +
  labs(title = "Gene z-score distribution per tissue", x = "Tissue", y = "Z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
########
longer_dataVoom <-  dataVoom$E %>%
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
                                 levels = c("in.vivo" ,"ex.vivo_2d", "ex.vivo_4d"))


prefix <- "CB1_human"
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% ISGs) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }



# Make tissue a factor first
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, 
                                        levels = c("in.vivo" ,"ex.vivo_2d", "ex.vivo_4d"))

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


ggplot(plot_df, aes(x = tissue, y = zscore, fill = tissue)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  labs(
    title = paste0("Z-score of ISGs across tissues"),
    x = "Tissue",
    y = "Z-score"
  ) +
  facet_wrap(~ genes, scales = "free_y") +  # one subplot per gene
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
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

longer_dataVoom <-  dataVoom$E %>%
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
                                 levels = c("in.vivo" ,"ex.vivo_2d", "ex.vivo_4d"))
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% mTORC1) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }
# Filter only genes in your mTORC1 list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% mTORC1)


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


ggsave(out("per_gene_per_sample_heatmap_mTORC1.png"))

###################
#cholesterol

longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% Cholesterol) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }
# Filter only genes in your mTORC1 list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% Cholesterol)


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


ggsave(out("per_gene_per_sample_heatmap_cholesterol.png"))




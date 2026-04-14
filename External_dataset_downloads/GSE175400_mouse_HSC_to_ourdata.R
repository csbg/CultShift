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
InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
InDir3 <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
out <- dirout("GSE175400_mouse_HSC_to_our_data")
# -----------------------------
# USER PARAMETERS
# -----------------------------
geo_id <- "GSE207740"           # GEO dataset


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
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE207740&format=file&file=GSE207740%5Fall%5Fcounts%5Fmetadata%2Etxt%2Egz"

# Save locally
dest <- "GSE207740_all_counts_metadata.txt.gz"
download.file(url, dest, mode = "wb")

# fread automatically detects headers unless it’s confused by quotes, so:
tbl <- read.delim(dest, header = TRUE, quote = "\"", fill = TRUE)
unique_genes <- make.unique(tbl$gene_name)

rownames(tbl) <- tbl$Geneid


# 3️⃣ Check first few columns
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207740/matrix/GSE207740_series_matrix.txt.gz"
dest <- "GSE207740_series_matrix.txt.gz"
download.file(url, dest, mode = "wb")

gse <- getGEO("GSE207740", GSEMatrix = TRUE)
# Access the GPL21103 platform
metadata_1 <- gse[["GSE207740_series_matrix.txt.gz"]]

# Extract phenoData from gse[[1]] (AnnotatedDataFrame)
adf <- phenoData(gse[["GSE207740_series_matrix.txt.gz"]])
# Use pData() from Biobase explicitly
metadata_1 <- Biobase::pData(adf)
rownames(metadata_1) <- colnames(tbl)

# Keep only overlapping samples
samples_keep <- intersect(rownames(metadata_1), colnames(tbl))
metadata_1 <- metadata_1[samples_keep,]
counts_1 <- tbl[,samples_keep]

metadata_1$tissue <- "ex.vivo"

metadata_1 <- metadata_1[c("tissue")]
metadata_1$sample <- rownames(metadata_1)
#############
#merge_our_data
counts <- read.delim(InDir3("combined_in_ex.vivo_with_Mye_counts_guide.tsv"), row.names = 1)
# Remove the first column (it’s now in rownames)
NTC_metadata <- read_rds(InDir2("NTC_meta.rds"))
NTC_metadata <- NTC_metadata[c("tissue","sample", "celltype")]
NTC_metadata <- NTC_metadata %>%
  filter(celltype == "HSC")

NTC_metadata$sample <- rownames(NTC_metadata)
NTC_metadata <- NTC_metadata %>%
  filter(tissue == "in.vivo")
NTC_metadata <- NTC_metadata[c("tissue","sample")]
NTC_counts <- counts[,rownames(NTC_metadata)]


metadata <- rbind(metadata_1, NTC_metadata)
both <- intersect(rownames(NTC_counts), rownames(counts))

counts_all <- cbind(NTC_counts,counts_1[both,])

counts_all <- counts_all %>% na.omit()
# -----------------------------
# DGEList + normalization
# -----------------------------
design <- model.matrix(~tissue, data = metadata)
dge <- DGEList(counts_all)
dge <- calcNormFactors(dge, method = "TMM")
#design <- model.matrix(as.formula(design_formula), metadata)
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr,]
# -----------------------------
# Align metadata and counts
# -----------------------------
# Ensure sample names match
all(colnames(counts_all) %in% rownames(metadata))
all(rownames(metadata) %in% colnames(counts_all))


metadata$tissue <- factor(metadata$tissue,
                          levels = c("in.vivo","ex.vivo"))
counts_all <- counts_all %>% na.omit()
# Rebuild design
design <- model.matrix(~tissue, data = metadata)

# DGEList + normalization
dge <- DGEList(counts_all)
dge <- calcNormFactors(dge, method = "TMM")

keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr, ]

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
                                 levels = c("in.vivo" ,"ex.vivo"))


prefix <- "mouseHSC_external_and_our"
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


ggplot(plot_df, aes(x = reorder(sample, zscore, FUN = median), y = zscore, fill = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Boxplot per sample, transparent
  geom_jitter(width = 0.2, size = 1, alpha = 0.8) +  # Jitter points for each observation
  labs(title = "ISG z-score distribution per sample",
       x = "Sample",
       y = "Z-score") +
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
                                 levels = c("in.vivo" ,"ex.vivo"))
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% ISG_core) %>%             # keep only your gene set
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




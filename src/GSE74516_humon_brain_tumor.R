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

out <- dirout("GSE74516_humon_brain_tumor")



# -----------------------------
# USER PARAMETERS
# -----------------------------
geo_id <- "GSE74516"           # GEO dataset
design_formula <- "~tissue"     # design formula for limma/voom
out <- dirout("GSE74516_humon_brain_tumor")


##############
genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))
genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

genesets$genes <- toupper(genesets$genes)
#

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
message("Downloading data...")

urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE74516&format=file&file=GSE74516_raw_counts_GRCh38.p13_NCBI.tsv.gz"
path <- paste(urld, "acc=GSE74516", "file=GSE74516_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&")
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
# Gene annotation
path <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
gmap <- as.data.frame(fread(path))


# Filter low-count genes
rownames(tbl) <- gmap$Symbol[match(rownames(tbl), gmap$GeneID)]
keep <- rowSums(tbl >= 10) >= 2
tbl <- tbl[keep,]
colnames(tbl)
# Metadata
gse <- getGEO("GSE74516", GSEMatrix = TRUE)
metadata <- gse$GSE74516_series_matrix.txt.gz


# Extract phenoData from gse[[1]] (AnnotatedDataFrame)
adf <- phenoData(gse[[1]])

# Use pData() from Biobase explicitly
metadata <- Biobase::pData(adf)
colnames(metadata)
#modify metadata
metadata$sample <- metadata$title %>%
  gsub("IN528ic Rep1","IN528_xeno_rep1",.) %>%
  gsub("IN528ic Rep2","IN528_xeno_rep2",.) %>%
  gsub("IN528culture Rep1","IN528_culture_rep1",.)%>%
  gsub("IN528culture Rep2","IN528_culture_rep2",.)%>%
  gsub("3565_IC_rep1","3565_xeno_rep1",.)%>%
  gsub("3565_IC_rep2","3565_xeno_rep2",.)%>%
  gsub("3565_IC_rep3","3565_xeno_rep3",.)
metadata$tissue <- ifelse(grepl("xeno",metadata$sample),"in.vivo","ex.vivo")
metadata <- metadata %>% dplyr::select(sample, geo_accession, tissue)
GBM3565 <- grep("3565", metadata$sample, value = T)
#GBM3565 <- grep("3565", metadata$tissue, value = T)
metadata <- metadata %>%
  filter(sample %in% GBM3565)
# metadata_3565 <- metadata %>%
#   filter(tissue %in% GBM3565)
# Keep only overlapping samples
samples_keep <- intersect(rownames(metadata), colnames(tbl))
metadata <- metadata[samples_keep,]
counts <- tbl[,samples_keep]
colnames(counts) <- metadata$sample[match(colnames(counts), metadata$geo_accession)]
rownames(metadata) <- metadata$sample
metadata$tissue <- factor(metadata$tissue, 
                          levels = c("in.vivo","ex.vivo"))
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
# Plots
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

#save---------
longer_dataVoom$organ <- "Glioblastoma_primary"
longer_dataVoom$author <- "Miller TE et al., Nature, 2017"

longer_dataVoom %>%
  dplyr::select(c("genes","tissue","sample","Expression","author","organ"))%>%
  write_rds(out("GSE74516_GBM_primary.rds"))





prefix <-"human_brain_tumor"
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
  filter(genes %in% ISGs)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "Average ISGs z-score per sample (colored by tissue)",
       x = "Sample",
       y = "Mean Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


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
ggsave(out("combined_plot_ISG.png"))
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

longer_dataVoom$organ <- "Glioblastoma_primary"
longer_dataVoom$author <- "Miller TE et al., Nature, 2017"

longer_dataVoom %>%
  dplyr::select(c("genes","tissue","sample","Expression","author","organ"))%>%
  write_rds(out("glioblastoma_primary_dataVoom.rds"))


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
#####################
ggplot(limmaRes,aes(
  x = logFC,
  y = -log10(adj.P.Val)))+
  geom_point()+
  facet_grid(cols = vars(coef))

#############
ranked_lists <- limmaRes %>%
  
  # Aggregate duplicates: take mean t-statistic per gene per coef
  group_by(coef, genes) %>%
  summarise(logFC = mean(logFC), .groups = "drop") %>%
  
  # Create named vector of ranks per coef
  group_by(coef) %>%
  summarise(
    ranks = list(setNames(logFC, genes)),
    .groups = "drop"
  )


pathways <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)


fgsea_results <- ranked_lists %>%
  mutate(
    fgsea = purrr::map(ranks, ~ fgsea(pathways = pathways,
                                      stats = .))   # no nperm
  ) %>%
  dplyr::select(coef, fgsea) %>%
  unnest(fgsea)

write_rds(fgsea_results, out("FGSEA_3565.rds"))
terms <- fgsea_results %>%
  filter(padj < 0.05) %>%
  pull(pathway)%>%
  unique()
fgsea_plot <- fgsea_results %>%
  filter(pathway %in% terms)
#####
# Number of pathways to plot
n_pathways <- length(unique(fgsea_plot$pathway))

# Number of coefficients (if multiple coefs are plotted)
n_coefs <- length(unique(fgsea_plot$coef))

# Define dynamic size (tweak multipliers as needed)
height <- max(4, n_pathways * 0.3)   # 0.3 inches per pathway, min 4 in
width  <- max(6, n_coefs * 2)        # 2 inches per coef, min 6 in

ggplot(fgsea_plot, aes(x = reorder(pathway, NES), y = coef,
                       fill = NES, size = -log10(padj))) +
  geom_point(shape = 21, color = "black") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(
    x = "Pathway",
    y = "Coefficient",
    fill = "NES",
    size = "-log10(padj)",
    title = "FGSEA by Coefficient"
  ) +
  optimized_theme_fig()

#######

##############
genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))


genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

genesets$genes <- toupper(genesets$genes)
#
limmaRes$celltype <- limmaRes$coef
# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

# gene subsets
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)

#
limmaRes_ISGs <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "ISGs")
limmaRes_ISG_core <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISG_core) %>%
  mutate(pathway = "ISGs")
limmaRes_Cholesterol <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "Cholesterol")
limmaRes_mTORC1 <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "mTORC1")


# Combine them into a named list
pathway_list <- list(
  ISGs = ISGs,
  ISG_core = ISG_core,
  Cholesterol = Cholesterol,
  mTORC1 = mTORC1
)

# Iterate over each pathway and gene set, build filtered results, and bind them together
limmaRes_all <- map_dfr(names(pathway_list), function(pathway_name) {
  limmaRes %>%
    dplyr::select(genes, logFC, adj.P.Val, celltype, group) %>%
    filter(genes %in% pathway_list[[pathway_name]]) %>%
    mutate(pathway = pathway_name)
})

# Optionally combine with genesets if you want one big table
limmaRes_all <- bind_rows(genesets, limmaRes_all)


# Function to generate the dot plot for a given geneset
plot_geneset <- function(data, geneset_name) {
  ggplot(data %>% filter(pathway == geneset_name), 
         aes(x = celltype, y = genes,
             color = pmin(3, pmax(-3, logFC)),
             size = pmin(5, -log10(adj.P.Val)))) +
    geom_point() +
    scale_color_gradient2(
      low = "#4C889C",
      mid = "white",
      high = "#D0154E",
      name = TeX("log_{2}(FC)")
    ) +
    scale_size_continuous(
      range = c(0, 1.8),
      name = TeX("$-\\log_{10}(p_{adj})$")
    ) +
    labs(
      title = paste0("Geneset: ", geneset_name),
      y = "Genes",
      x = "Cell Type"
    ) +
    #facet_grid(rows = vars(pathway), scales = "free", space = "free") +
    optimized_theme_fig()
}

# Make three separate plots
ISG <- plot_geneset(limmaRes_all, "ISGs")
ISG_core <- plot_geneset(limmaRes_all, "ISG_core")
mTORC1 <- plot_geneset(limmaRes_all, "mTORC1")
Cholesterol <- plot_geneset(limmaRes_all, "Cholesterol")

# Save each plot as a PNG
ggsave(out("ISG.png"), ISG, width = 8, height = 6, dpi = 300)
ggsave(out("ISG_core.png"), ISG_core, width = 8, height = 6, dpi = 300)
ggsave(out("mTORC1.png"), mTORC1, width = 8, height = 6, dpi = 300)
ggsave(out("Cholesterol.png"), Cholesterol, width = 8, height = 6, dpi = 300)


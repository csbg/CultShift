source("src/00_init.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)
library(tibble)
library(data.table)
library(msigdbr)
library(org.Hs.eg.db) 
library(fgsea)# for human; use org.Mm.eg.db for mouse

InDir <- dirout("Ag_top_pathway_genes")
##*****
##
# -----------------------------
# USER PARAMETERS
# -----------------------------
geo_id <- "GSE110968"           # GEO dataset
count_file <- "raw_counts_GRCh38.tsv.gz"  # raw counts filename in GEO
annot_file <- "annot_GRCh38.tsv.gz"      # gene annotation filename in GEO
design_formula <- "~tissue"     # design formula for limma/voom

prefix <- "example"
out <- dirout("GSE110968_CB2_external_dataset")
# Metadata
gse <- getGEO("GSE110968", GSEMatrix = TRUE)
metadata <- gse$GSE110968_series_matrix.txt.gz
# Extract phenoData from gse[[1]] (AnnotatedDataFrame)
adf <- phenoData(gse[[1]])

# Use pData() from Biobase explicitly
metadata <- Biobase::pData(adf)

##############
genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))
genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

genesets$genes <- toupper(genesets$genes)


# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

# gene subsets
ISGs <- toupper(get_genes_for_pathway("ISGs", genesets))
ISG_core <- toupper(get_genes_for_pathway("ISG_core", genesets))
mTORC1 <- toupper(get_genes_for_pathway("mTORC1", genesets))
Cholesterol <- toupper(get_genes_for_pathway("Cholesterol", genesets))
########****
out <- dirout("GSE110968_CB2_external_dataset")
counts <- fread(out("../../GSE110968_CB2_raw_counts_GRCh38.p13_NCBI.tsv/GSE110968_raw_counts_GRCh38.p13_NCBI.tsv"))
colnames(counts)
colnames(counts) <- c("GENEID",
                      "uncultured1",
                      "uncultured2",
                      "uncultured3",
                      "cultured_2days1",
                      "cultured_2days2",
                      "cultured_2days3",
                      "ex1",
                      "ex2",
                      "ex3",
                      "cultured_4days1",
                      "cultured_4days2",
                      "cultured_4days3",
                      "ex4",
                      "ex5",
                      "ex6")

#change ID to names

# convert counts to data.table if not already
counts <- as.data.table(counts)

# Keep original IDs as a column
counts$ORIG_ID <- counts$GENEID

# Map Entrez IDs to symbols
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(counts$GENEID),
  keytype = "ENTREZID",
  columns = c("SYMBOL"))

# Merge
counts <- merge(gene_map, counts, by.x = "ENTREZID", by.y = "GENEID", all.y = TRUE)

# Replace missing SYMBOLs with original IDs
counts$SYMBOL[is.na(counts$SYMBOL) | counts$SYMBOL == ""] <- 
  as.character(counts$ORIG_ID[is.na(counts$SYMBOL) | counts$SYMBOL == ""])

# Make SYMBOLs unique
counts$SYMBOL <- make.unique(as.character(counts$SYMBOL))

# Set rownames
rownames(counts) <- counts$SYMBOL

# Drop helper columns
counts$SYMBOL <- NULL
counts$ENTREZID <- NULL
counts$ORIG_ID <- NULL



head(counts)

counts <- counts[,c(
  "uncultured1",
  "uncultured2",
  "uncultured3",
  "cultured_2days1",
  "cultured_2days2",
  "cultured_2days3",
  "cultured_4days1",
  "cultured_4days2",
  "cultured_4days3"
)]
metadata <- data.frame(
  sample = colnames(counts)
) %>%
  mutate(
    tissue = ifelse(grepl("uncul", sample), "in.vivo",
                    ifelse(grepl("2days", sample), "ex.vivo_2days", "ex.vivo_4days"))
  )

rownames(metadata) <- metadata$sample

# Make sure factor levels match actual values
metadata$tissue <- factor(metadata$tissue, 
                          levels = c("in.vivo","ex.vivo_2days","ex.vivo_4days"))

# Double check order
stopifnot(all(colnames(counts) == rownames(metadata)))

# Create design matrix
design <- model.matrix(~tissue, data = metadata)
rownames(design) <- rownames(metadata)

##################################
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")
keep <- filterByExpr(d0,design)
d <- d0[keep,]


###############################################################################
#setting the model
###############################################################################


# Normalization
dataVoom <- voom(d, design, plot = T)

dataVoom %>% write_rds(out("dataVoom_perCTex.vivovsin.vivo.rds"))
########################
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
coef_names <- setdiff(colnames(coef(limmaFit)), "(Intercept)")


# Directory for rank lists
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

# -----------------------------
# Pivot to long format
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
                                 levels = c("in.vivo","ex.vivo_2days","ex.vivo_4days"))
longer_dataVoom$organ <- "cord-blood"
longer_dataVoom$author <- "Papa L et al., Experimental Hematology, 2023, Papa L et al., Blood Advances, 2018"

longer_dataVoom %>%
  dplyr::select(c("genes","tissue","sample","Expression","author","organ"))%>%
  write_rds(out("CB_NTCs_dataVoom.rds"))

################
prefix <- "CB_ntc"
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% ISGs
  ) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  dplyr::select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }



# # Make tissue a factor first
# longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, 
#                                         levels = c("in.vivo","ex.vivo_2days","ex.vivo_4days"))
###############
#plots controls
#Housekeeping <-  c("Gapdh","Tbp","Ubc","Actb","Pgk1")
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% toupper(ISGs)) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  dplyr::select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }



# # Make tissue a factor first
# longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, 
#                                         levels = c("in.vivo","ex.vivo_2days","ex.vivo_4days"))

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
  filter(genes %in% toupper(ISGs))


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "ISG core",
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
ggsave(out("combined_plot_ISG.png"))
#####################
ggplot(limmaRes,aes(
  x = logFC,
  y = -log10(adj.P.Val)))+
  geom_point()+
  facet_grid(cols = vars(coef))

#############
ranked_lists <- limmaRes %>%
  
  # Aggregate duplicates: take mean t-statistic per gene per coef
  group_by(coef, ensg) %>%
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
fgsea_results$dataset <-"glioblastoma"
fgsea_results$organ <-"glioblastoma"
fgsea_results$author <- "Liu SJ et al., Genome Biology, 2024"
fgsea_results %>% write_rds(out("NTC_fgsea_glioblastoma.rds"))
write_rds(fgsea_results, out("FGSEA_CB_2.rds"))
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
ggsave(out("fgsea_CB1.png"))
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


# Define your gene sets (character vectors)
ISGs <- get_genes_for_pathway("ISGs", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)

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
###################################
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
                                 levels = c("in.vivo" ,"ex.vivo_2ays", "ex.vivo_4days"))


prefix <- "CB2_human"
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
                                        levels = c("in.vivo" ,"ex.vivo_2days", "ex.vivo_4days"))

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
                                 levels = c("in.vivo" ,"ex.vivo_2days", "ex.vivo_4days"))
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




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

out <- dirout("GSE229032_AML_primary_vs_cell_line")




urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE229032", "file=GSE229032_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&")
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
#

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]
counts <- tbl
#
path <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
gmap <- as.data.frame(fread(path))
# Ensure your gmap has unique GeneID values
rownames(counts) <- gmap$Symbol[match(rownames(counts), gmap$GeneID)]
#metadata
gse <- getGEO("GSE229032", GSEMatrix = TRUE)
metadata <- gse$GSE229032_series_matrix.txt.gz


# Extract phenoData from gse[[1]] (AnnotatedDataFrame)
adf <- phenoData(gse[[1]])

# Use pData() from Biobase explicitly
metadata <- Biobase::pData(adf)
metadata$source_name_ch1 <- gsub("Acute myeloid leukemia cell line","AML_CL",
                                 metadata$source_name_ch1)

metadata$source_name_ch1 <- gsub("AML Primary sample","AML_primary",
                                 metadata$source_name_ch1)
metadata$tissue <- metadata$source_name_ch1
metadata <- metadata %>%
  dplyr::select(c("title","geo_accession","tissue"))
colnames(counts)
rownames(metadata)
list <- rownames(metadata)[rownames(metadata)%in%colnames(counts)]
metadata <- metadata[list,]
counts <- counts[,list]
metadata$tissue <- factor(metadata$tissue,
                          levels = c("AML_primary","AML_CL"))

design <- model.matrix(~tissue, data=metadata)

#########
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")
keep <- filterByExpr(d0,design)
d0 <- d0[keep,]

# voom transformation
dataVoom <- voom(d0, design, plot=TRUE)
head(dataVoom$E)

limmaFit <- lmFit(dataVoom, design)
limmaFit <- eBayes(limmaFit)
colnames(coef(limmaFit))
#

# Initialize list to store top table results
limmaRes <- list() # start an empty list
for(coefx in colnames(coef(limmaFit))){ # run a loop for each coefficient
  print(coefx)
  # topTable returns the statistics of our genes. We then store the result of each coefficient in a list.
  limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) %>%
    rownames_to_column("genes")
}
limmaRes <- bind_rows(limmaRes, .id = "coef") # bind_rows combines the results and stores the name of the coefficient in the column "coef"
limmaRes <- filter(limmaRes, coef != "(Intercept)") # then we keep all results except for the intercept
limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))


#######
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
fgsea_results$dataset <-"AML"
fgsea_results$organ <-"AML"
fgsea_results$author <- "Paudel BB et al., Blood Advances, 2024"
fgsea_results %>% write_rds(out("NTC_fgsea_AML.rds"))

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
limmaRes$celltype <- "AML_CL"
# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

# gene subsets
ISGs <- get_genes_for_pathway("ISGs", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)
ISG_core <- get_genes_for_pathway("ISG_core", genesets)
#
limmaRes_ISGs <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "ISG_core")
limmaRes_Cholesterol <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "Cholesterol")
limmaRes_mTORC1 <- limmaRes %>%
  dplyr::select(c("genes","logFC","adj.P.Val","celltype","group")) %>%
  filter(genes %in% ISGs) %>%
  mutate(pathway = "mTORC1")


# Define your gene sets (character vectors)
ISGs <- get_genes_for_pathway("ISG_core", genesets)
Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)
mTORC1 <- get_genes_for_pathway("mTORC1", genesets)

# Combine them into a named list
pathway_list <- list(
  ISG_core = ISGs,
  Cholesterol = Cholesterol,
  mTORC1 = mTORC1
)

# Iterate over each pathway and gene set, build filtered results, and bind them together
limmaRes_all <- map_dfr(names(pathway_list), function(pathway_name) {
  limmaRes %>%
    select(genes, logFC, adj.P.Val, celltype, group) %>%
    filter(genes %in% pathway_list[[pathway_name]]) %>%
    mutate(pathway = pathway_name)
})

# Optionally combine with genesets if you want one big table
limmaRes_all <- bind_rows(genesets, limmaRes_all)



ggplot(limmaRes_all, aes(x = celltype, y = genes,
                                    color = pmin(3, pmax(-3, logFC)),
                                    size = pmin(5,-log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name=TeX("log_{2}(FC)"))+
  geom_point() +
  scale_size_continuous(
    range = c(0, 1.8),
    #breaks = c(0,2,5,10),
    #limits = c(0, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  labs(title = "Downregulation of ISGs and upregulation
       of growth/metabolic genes in ex vivo ",
       y = "Genes",
       x = "Cell Type") +
  facet_grid(rows = vars(pathway), scales = "free", space = "free") +
  # coord_flip() +
  optimized_theme_fig()

#colnames(genesets) <- c("genes","logFC_pert","adj.P.Val_pert","pathway","celltype")
################################################################################


################################################################################
metadata$sample <- metadata$geo_accession
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
longer_dataVoom$sample <- longer_dataVoom$title
longer_dataVoom$tissue <- factor(longer_dataVoom_zscore$tissue, levels = c("AML_primary", "AML_CL"))
longer_dataVoom$organ <- "AML"
longer_dataVoom$author <- "Paudel BB et al., Blood Advances, 2024"

longer_dataVoom %>%
  dplyr::select(c("genes","tissue","sample","Expression","author","organ"))%>%
  write_rds(out("AML_NTCs_dataVoom.rds"))



longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% ISGs) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  dplyr::select(tissue,  genes, sample, Expression, zscore, title) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }


prefix <-"AML"
# Make tissue a factor first
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, levels = c("AML_primary", "AML_CL"))

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
  group_by(title, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(title, mean_z), y = mean_z, fill = tissue)) +
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
ggplot(plot_df, aes(x = title, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4C889C", mid = "white", high = "#D0154E", midpoint = 0) +
  labs(title = "Z-score of ISGs across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(out("per_gene_per_sample_heatmap_ISGs.png"))


##########
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
  select(tissue,  genes, sample, Expression, zscore, title) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }



# Make tissue a factor first
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, levels = c("AML_primary", "AML_CL"))

# Reorder samples so ex.vivo samples come first
longer_dataVoom_zscore$sample <- factor(
  longer_dataVoom_zscore$sample,
  levels = longer_dataVoom_zscore %>%
    arrange(tissue, sample) %>%
    pull(sample) %>%
    unique()
)




# Filter only genes in your Cholesterol list
plot_df <- longer_dataVoom_zscore %>%
  filter(genes %in% Cholesterol)


summary_df <- plot_df %>%
  group_by(title, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(title, mean_z), y = mean_z, fill = tissue)) +
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
    title = paste0("Z-score of Cholesterol across tissues"),
    x = "Tissue",
    y = "Z-score"
  ) +
  facet_wrap(~ genes, scales = "free_y") +  # one subplot per gene
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
#########################
plot_df$sample <- factor(plot_df$sample, levels = unique(plot_df$sample))

# Plot heatmap
ggplot(plot_df, aes(x = title, y = genes, fill = zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Z-score of Cholesterol across samples",
       x = "Sample",
       y = "Gene") +
  optimized_theme_fig() +
  facet_grid(cols = vars(tissue), scales = "free")
theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(out("per_gene_per_sample_heatmap_Cholesterol.png"))
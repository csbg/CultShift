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

out <- dirout("GSE112786_neuroendocrine_prostate_cancer")
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE112786", "file=GSE112786_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
counts <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
path <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
gmap <- as.data.frame(fread(path))
# Ensure your gmap has unique GeneID values
rownames(counts) <- gmap$Symbol[match(rownames(counts), gmap$GeneID)]

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( counts >= 10 ) >= 2
counts <- counts[keep, ]

# log transform raw counts
# instead of raw counts can display vst(as.matrix(counts)) i.e. variance stabilized counts
dat <- log10(counts + 1)

# box-and-whisker plot
par(mar=c(7,4,2,1))
boxplot(dat, boxwex=0.7, notch=T, main="GSE112786", ylab="lg(cnt + 1)", outline=F, las=2)

#
colnames(counts) <- c("PDOX1","PDOX2","PDOX3",
                   "organoid2D1","organoid1","organoid2D2",
                   "organoid2","organoid3","organoid4")

# get sample names
samples <- colnames(counts)

# tissue: DMSO samples are ex vivo, others in vivo
tissue <- ifelse(grepl("PDOX", samples), "xenograft", 
                 ifelse(grepl("organoid2D", samples), "organoid2D","organoid"))
tissue <- factor(tissue, levels = c("xenograft", "organoid2D", "organoid"))
# celltype: extract prefix before first underscore
# for DMSO samples, it seems CB1/CB2/CB3 info is after first underscore, so we want the part before "_CB"

# create metadata dataframe
metadata <- data.frame(
  sample = samples,
  tissue = tissue
  
)

design <- model.matrix(~tissue, data=metadata)

#########
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")
keep <- filterByExpr(d0,design)
d0 <- d0[keep,]

#########
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")
keep <- filterByExpr(d0,design)
d0 <- d0[keep,]

# voom transformation
dataVoom <- voom(d0, design, plot=TRUE)
head(dataVoom$E)
#######

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



genes_fig_1 <- read_rds(InDir_genes("genes_fig1.rds"))
colnames(genes_fig_1) <- c("genes","pathway")
ISG_fig_1 <- genes_fig_1[genes_fig_1$pathway =="ISG core",]$genes
prefix <- "ISG_fig_1"
gene_set <- toupper(ISG_fig_1)
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% gene_set) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }

longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% gene_set) %>%             # keep only your gene set
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
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, levels = c("xenograft", "organoid2D",
                                                                                  "organoid"))

# Reorder samples so ex.vivo samples come first
longer_dataVoom_zscore$sample <- factor(
  longer_dataVoom_zscore$sample,
  levels = longer_dataVoom_zscore %>%
    arrange(tissue, sample) %>%
    pull(sample) %>%
    unique()
)


ggplot(longer_dataVoom_zscore, aes(x = sample, y = zscore, color = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # box per sample
  #geom_jitter(width = 0.2, alpha = 0.7, color = "blue") +  # points per gene
  labs(title = paste0("Z-score of ", prefix, " genes per sample"),
       x = "Sample",
       y = "Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
##########
#cholesterol
genes_fig_1 <- read_rds(InDir_genes("genes_fig1.rds"))
colnames(genes_fig_1) <- c("genes","pathway")
cholesterol <- genes_fig_1[genes_fig_1$pathway =="mTORC1_or_Cholesterol",]$genes
cholesterol <- toupper(cholesterol)
prefix <- "ISG_fig_1"
gene_set <- cholesterol
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% gene_set) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
  select(tissue,  genes, sample, Expression, zscore) %>%
  { write_rds(., out(paste0("zscore_plot_external_", prefix, ".rds"))); . }

longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% gene_set) %>%             # keep only your gene set
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
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, levels = c("xenograft", "organoid2D",
                                                                                  "organoid"))


# Reorder samples so ex.vivo samples come first
longer_dataVoom_zscore$sample <- factor(
  longer_dataVoom_zscore$sample,
  levels = longer_dataVoom_zscore %>%
    arrange(tissue, sample) %>%
    pull(sample) %>%
    unique()
)


ggplot(longer_dataVoom_zscore, aes(x = sample, y = zscore, color = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # box per sample
  #geom_jitter(width = 0.2, alpha = 0.7, color = "blue") +  # points per gene
  labs(title = paste0("Z-score of ", prefix, " genes per sample"),
       x = "Sample",
       y = "Z-score") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

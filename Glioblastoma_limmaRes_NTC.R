##############################
# Modular GEO RNA-seq pipeline
##############################
#Aarathy
###############
source("src/00_init.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_enrichR_mouse_genes.R")
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
library(Matrix)
source("src/Ag_ScRNA_11_invivo_exvivo_KO_limma_function.R")
InDir <- dirout("Ag_top_pathway_genes")

InDir1 <- dirout("Glioblastoma")
###########
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

#############
out <- dirout("Glioblastoma_limmaRes")
meta <- read_rds(InDir1("metadata.rds"))
counts <- read_rds(InDir1("counts.rds"))
meta$genotype <- gsub("non-targeting","NTC",meta$genotype)

#subsetting
meta$tissue <- ifelse(grepl("invitro", meta$condition), "ex.vivo",
                      ifelse(grepl("preinf", meta$condition), "in.vivo",
                             "CED"))
meta <- meta %>%
  filter(tissue != "CED")%>%
  filter(RT_status =="noRT")%>%
  filter(genotype == "NTC")

meta$sample <- rownames(meta)

#########
metadata <- meta
metadata$tissue <- factor(metadata$tissue)
metadata$tissue <- relevel(metadata$tissue, ref = "in.vivo")

counts <- counts[, rownames(metadata)]
############################
# -----------------------------
# DGEList + normalization
# -----------------------------
dge <- DGEList(counts)
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(as.formula(~tissue), metadata)
keep_expr <- filterByExpr(dge, design)
dge <- dge[keep_expr,]

# voom transformation------------
dataVoom <- voom(dge, design, plot = TRUE)

# -----------------------------
# limma Res------------------
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

up <- limmaRes %>%
  filter(group == "up")%>%
  pull(genes)%>%
  unique()
enrich_results_up <- enrichr(up, ENRICHR.DBS[3])
enrich_up <- enrich_results_up$MSigDB_Hallmark_2020
down <- limmaRes %>%
  filter(group == "down")%>%
  pull(genes)%>%
  unique()
enrich_results_down <- enrichr(down, ENRICHR.DBS[3])
enrich_down <- enrich_results_down$MSigDB_Hallmark_2020
interferon_ov <- enrich_down%>%
  filter(Term %in% c("Interferon Gamma Response", "Interferon Alpha Response"))%>%
  pull(Genes)%>%
  strsplit(";") %>%           # split both strings by ';'
  unlist() %>%
  trimws() %>%                # remove extra spaces
  unique() %>%                # remove duplicates
  paste(collapse = ",") 
#######
#fgsea-----
gsea.res <- data.table()

for (dbx in names(enr.terms)) {
    
    stats <- with(limmaRes, setNames(logFC, nm = genes))
    
    if (any(is.na(stats))) {
      next
    }
    
    fgsea_output <- fgsea(pathways = enr.terms[[dbx]], stats = stats)
    
    if (length(fgsea_output) > 0) {
      gsea.res <- rbind(gsea.res, data.table(fgsea_output,  db = dbx))
    }
  }
gsea.res$dataset <-"glioblastoma"
gsea.res$organ <-"glioblastoma"
gsea.res$author <- "Liu SJ et al., Genome Biology, 2024"
gsea.res %>% write_rds(out("NTC_fgsea_glioblastoma.rds"))

# -----------------------------
# Plots
# -----------------------------
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

# -----------------------------
colnames(longer_dataVoom)
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
longer_dataVoom$organ <- "brain_glioblastoma"
longer_dataVoom$author <- "Liu et al., Genome Biology, 2024"

longer_dataVoom %>%
  dplyr::select(c("genes","tissue","sample","Expression","author","organ"))%>%
  write_rds(out("Glioblastoma_NTCs.rds"))

################
prefix <- "glioblastoma_ntc"
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% interferon_ov
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



# Make tissue a factor first
longer_dataVoom_zscore$tissue <- factor(longer_dataVoom_zscore$tissue, 
                                        levels = c("in.vivo" ,"ex.vivo"))
###############
#plots controls
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
ggsave(out("combined_plot.png"))
##################

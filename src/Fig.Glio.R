##############################
# Modular GEO RNA-seq pipeline
##############################
#Aarathy
###############
source("src/00_init.R")
#source("src/Ag_enrichR_mouse_genes.R")
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
library(Matrix)
library(circlize)
library(grid)

InDir1 <- dirout("Ag_top_pathway_genes")
InDir <- dirout("Glioblastoma_limmaRes.3way.mod")
out <- dirout("Fig.Glio_Fig5")
limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))
dataVoom <- read_rds(InDir("Glioblastoma_dataVoom.rds"))
unique(limmaRes$coef)
head(limmaRes)
#correlation----------------
#Effect of culture on Radiotherapy
RT_effect_ex.vs.in <- limmaRes %>%
  filter(coef %in% c("tissueex.vivo.RT_statusRT",
                     "RT_statusRT","tissue.ex.vivo","RT_status_ex.vivo"))

# Subset the coefficients of interest
rt_effects_ntc <- limmaRes %>%
  filter(coef %in% c("tissueex.vivo", "tissueex.vivo.RT_statusRT"))

# Spread to wide format for correlation calculation
rt_wide <- rt_effects_ntc %>%
  dplyr::select(ensg, coef, logFC) %>%   # adjust 'gene' column name if different
  tidyr::pivot_wider(names_from = coef, values_from = logFC)


cor_test <- cor.test(
  abs(rt_wide$tissueex.vivo),
  abs(rt_wide$tissueex.vivo.RT_statusRT),
  use = "pairwise.complete.obs",
  method = "pearson"
)

# Print results
cor_test


ggplot(rt_wide,
       aes(x = pmin(abs(tissueex.vivo), 5),
           y = pmin(abs(tissueex.vivo.RT_statusRT), 5))) +
  geom_hex(bins = 80) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b",
                      trans = "log10") +
  geom_smooth(method = "lm", se = FALSE,
              color = "#e41a1c", size = 1, linetype = "dotted") +
  labs(
    x = "abs(culture effect)",
    y = "abs(RT interaction effect_NTC)",
    title = paste0(
      "Pearson r = ", round(cor_test$estimate, 3),
      ", p = ", signif(cor_test$p.value))
  ) +
  optimized_theme_fig()


ggsave(out("correlation_plot_RT_vs_culture.pdf"),w = 4.5, h = 3.5, units = "cm")
# Extract correlation value
cor_val <- cor_test$estimate

# Create a small data frame for plotting
df_cor <- data.frame(
  metric = "Culture vs RT interaction",
  correlation = cor_val
)

ggplot(df_cor, aes(x = metric, y = correlation, fill = correlation)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, color = "black") +
  geom_text(
    aes(label = round(correlation, 3)),
    vjust = ifelse(cor_val > 0, -0.5, 1.2),
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Correlation\n(–1 to 1)"
  ) +
  labs(
    x = NULL,
    y = "Pearson correlation",
    title = "Correlation between culture effect and RT interaction effect"
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave(
  out("barplot_correlation_RT_vs_culture.pdf"),
  w = 4, h = 3, units = "cm"
)

################
#heatmap------------------
limma_results <- limmaRes%>%
  filter(coef == "tissueex.vivo.RT_statusRT")
get_top_genes <- function(
                          limma_results,
                          selected_coef,
                          logFC_threshold = 1,
                          pval_threshold = 0.05, 
                          top_n = 10) {
  # Filter limma results for the pathway genes
  filtered_genes <- limma_results %>%
    filter(coef == selected_coef) %>%
    filter(group != "n.s") %>%
    filter( adj.P.Val < pval_threshold) %>%
    arrange(adj.P.Val) %>%
    #mutate(abs_logFC = abs(logFC)) %>%
    arrange(abs(logFC)) %>%
    slice_head(n = top_n) %>%
    pull(ensg)%>%unique()
  
  
  
  return(filtered_genes)
}

top_interaction <- get_top_genes(
  limmaRes,
  selected_coef = "tissueex.vivo.RT_statusRT",
  logFC_threshold = 1,
  pval_threshold = 0.05, 
  top_n = Inf)
top_NTC <- get_top_genes(
  limmaRes,
  selected_coef = "tissueex.vivo",
  logFC_threshold = 1,
  pval_threshold = 0.05, 
  top_n = Inf)
# Get intersected gene IDs
common_genes <- intersect(top_interaction, top_NTC)

# Filter limma results for interaction effect only, among common genes
top20_interaction <- limmaRes %>%
  filter(coef == "tissueex.vivo.RT_statusRT",
         ensg %in% common_genes) %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%   # sort by padj (ascending) then logFC (descending)
  slice_head(n = 20)%>%
  pull(ensg)

# View results

dataVoom <- read_rds(InDir("Glioblastoma_mouse.rds.rds"))

dataVoom_selected <- dataVoom%>%
  filter(genes %in% top20_interaction)


# Ensure proper ordering of samples based on tissue and RT status
# You can adapt this if your sample naming differs
dataVoom_selected <- dataVoom_selected %>%
  mutate(
    RT_status = case_when(
      grepl("noRT", sample, ignore.case = TRUE) ~ "noRT",
      grepl("RT", sample, ignore.case = TRUE) ~ "RT"
    ),
    condition = case_when(
      grepl("ex\\.vivo", tissue) ~ paste0("ex.vivo_", RT_status),
      grepl("in\\.vivo", tissue) ~ paste0("in.vivo_", RT_status)
    )
  )


# Pivot to wide format (genes as rows, samples as columns)
expr_mat <- dataVoom_selected %>%
  dplyr::select(genes, sample, Expression) %>%
  pivot_wider(names_from = sample, values_from = Expression) %>%
  column_to_rownames("genes") %>%
  as.matrix()

# Optional: scale expression per gene (z-score)
expr_scaled <- t(scale(t(expr_mat)))

# Define color palette
col_fun <- colorRamp2(
  breaks = c(-2, 0, 2),
  colors = c(low = "#4C889C",
             mid = "white",
             high = "#D0154E")
)
cn <- colnames(expr_scaled)


# Create annotation for sample groups
sample_anno <- data.frame(
  Condition = dplyr::case_when(
    grepl("48hit_RT", cn, ignore.case = TRUE)       ~ "ex.vivo_RT",
    grepl("48hit_noRT", cn, ignore.case = TRUE)     ~ "ex.vivo_noRT",
    grepl("_RT_preinf", cn, ignore.case = TRUE)     ~ "in.vivo_RT",
    grepl("_noRT_preinf", cn, ignore.case = TRUE)   ~ "in.vivo_noRT",
    TRUE                                             ~ NA_character_ ))
rownames(sample_anno) <- colnames(expr_scaled)
sample_order <- c("ex.vivo_noRT", "ex.vivo_RT", "in.vivo_noRT", "in.vivo_RT")

# Reorder columns by condition
col_order <- order(sample_anno$Condition)
expr_scaled_ordered <- expr_scaled[, col_order]
sample_anno_ordered <- sample_anno[col_order, , drop = FALSE]

# Define colors for each condition
cond_colors <- c(
  "ex.vivo_noRT" = "#6a3d9aff",  # vivid blue
  "ex.vivo_RT"   = "#976CAB",  # purple
  "in.vivo_noRT" = "#d38d5fff",  # cyan-blue
  "in.vivo_RT"   = "#D2B370"   # light lavender-blue
)

# Create top annotation
ha <- HeatmapAnnotation(
  Condition = sample_anno_ordered$Condition,
  col = list(Condition = cond_colors),
  height = unit(0.2, "cm"),               # bar thickness
  simple_anno_size = unit(0.2, "cm"),     # bar tile thickness
  
  # 🔹 Make top annotation text match your heatmap labels
  annotation_name_gp = gpar(fontsize = 5, fontfamily = "sans", fontface = "plain"),
  
  annotation_legend_param = list(
    title_gp  = gpar(fontsize = 5, fontfamily = "sans", fontface = "plain"),
    labels_gp = gpar(fontsize = 5, fontfamily = "sans", fontface = "plain")
))
# Assume expr_scaled_ordered is your data matrix
n_rows <- nrow(expr_scaled_ordered)
n_cols <- ncol(expr_scaled_ordered)

# Set exact size per cell (for example, 4 mm each)
cell_size <- unit(2, "mm")

# Compute total heatmap dimensions
heatmap_width  <- cell_size * n_cols
heatmap_height <- cell_size * n_rows

pdf(out("Heatmap_NTC.pdf"), 
    width = convertWidth(heatmap_width, "cm", valueOnly = TRUE),
    height = convertHeight(heatmap_height, "cm", valueOnly = TRUE))

Heatmap(
  expr_scaled_ordered,
  name = "Z-score",
  col = col_fun,
  top_annotation = ha,
  cluster_rows = TRUE,
  show_row_dend = F,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  
  # 🔹 Consistent font for both row names and legend labels
  row_names_gp = gpar(fontsize = 5, fontfamily = "sans", fontface = "plain"),
  
  
  
  width = heatmap_width,
  height = heatmap_height,
  
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 5, fontfamily = "sans", fontface = "plain"),
    labels_gp = gpar(fontsize = 5, fontfamily = "sans", fontface = "plain")
  ),
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x = x, y = y, width = cell_size, height = cell_size,
              gp = gpar(col = NA, fill = fill))
  }
)


dev.off()
#box plot
# Define outline colors by tissue
outline_colors <- c(
  "ex.vivo_noRT" = "#6a3d9aff",
  "ex.vivo_RT" = "#976cabff",      # violet
  "in.vivo_noRT" = "#d38d5fff",
  "in.vivo_RT"  = "#d2b370ff"      # yellow
)

# Define nicer x-axis labels
nice_labels <- c(
  "ex.vivo_noRT" = "ex vivo noRT",
  "ex.vivo_RT"   = "ex vivo RT",
  "in.vivo_noRT" = "in vivo noRT",
  "in.vivo_RT"   = "in vivo RT"
)

# Plot
ggplot(dataVoom_selected[dataVoom_selected$genes=="Ifi47",], 
       aes(x = condition, y = Expression)) +
  
  # Draw boxes (outline)
  geom_boxplot(aes(color = condition),
               fill = "white",
               position = position_dodge(width = 0.8),
               size = 0.2,
               outlier.shape = NA) +
  
  scale_color_manual(values = outline_colors) +
  scale_x_discrete(labels = nice_labels) +   # <-- nicer labels
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  
  labs(
    x = "Condition",
    y = "Expression",
    title = "Ifi47 Expression"
  ) +
  optimized_theme_fig()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Save
ggsave(out("Expression_Ifi47.pdf"), w = 2.5, h = 4, units = "cm")

# Enrichment analysis---------------
gsea.res <- read_rds(InDir("gsea.res.ntc.rds"))
unique(gsea.res$pathway)
write_rds(gsea.res, out("Glioblastoma_noRT.gsea.rds"))
head(gsea.res)
interferon_pathways <- gsea.res[pathway %in% c("Interferon Alpha Response", "Interferon Gamma Response")]

interferon_genes <- unlist(lapply(interferon_pathways$leadingEdge, as.character))
interferon_genes <- unlist(strsplit(interferon_genes, split = ","))
ISG_glioblastoma <- interferon_genes
ISG_glioblastoma %>% write_rds(out("ISG_glioblastoma.rds"))

gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge,
                                      function(vec) paste(vec[1:10], collapse = ","))
dbx<-"MSigDB_Hallmark_2020"
pDT <- gsea.res[db == dbx]
## Splitting the task to handle both ends of the NES spectrum-positive and negative
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=4)]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=4)]$pathway)

# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(c(pw.display.pos, pw.display.neg))
pDT <- pDT[pathway %in% pw.display]

if (nrow(pDT) > 0){
  # Step 1: Aggregate NES values across all celltypes (mean of NES per pathway)
  pDT_agg <- pDT %>%
    group_by(pathway) %>%
    summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
    arrange(desc(average_NES))  # Ordering pathways by the average NES, highest first
  
  # Step 2: Create a factor for pathway that reflects the aggregated NES order
  pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)
  
  
}
enrichment <- ggplot(pDT, aes(y=pathway, x= pmin(10, -log10(padj)), fill = NES )) +
  
  scale_fill_gradient2(low = "#4C889C",
                       mid = "white",
                       high = "#D0154E",
                       name=TeX("NES"))+
  #name=TeX("log_{2}(FC)"))+
  geom_col() +
  # scale_size_continuous(
  #   range = c(0, 1.8),
  #   #limits = c(0, 5),
  #   name=TeX("$-\\log_{10}(p_{adj})$"))+
  # 
  
  #xRot() +
  #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
  labs(y = "Pathways",
       x = "-log10(Padj)",
       title = "Enriched pathways")+
  
  #coord_flip()+
  
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,
  ),
  legend.position = "right", legend.direction = "vertical",
  legend.justification = "bottom")+
  optimized_theme_fig()+
  theme(panel.grid.major = element_blank())

ggsave(out("fgsea.ntc.pdf"), w = 7, h = 5, units = "cm")
##################
unique(limmaRes$coef)




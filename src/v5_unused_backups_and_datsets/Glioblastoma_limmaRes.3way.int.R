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
library(Matrix)
source("src/Ag_ScRNA_11_invivo_exvivo_KO_limma_function.R")
InDir1 <- dirout("Ag_top_pathway_genes")

InDir <- dirout("Glioblastoma")

out <- dirout("Glioblastoma_limmaRes")
meta <- read_rds(InDir("metadata.rds"))
counts <- read_rds(InDir("counts.rds"))
meta$genotype <- gsub("non-targeting","NTC",meta$genotype)
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
unique(meta$condition)
#subsetting
meta$tissue <- ifelse(grepl("invitro", meta$condition), "ex.vivo",
                      ifelse(grepl("preinf", meta$condition), "in.vivo",
                             "CED"))
meta <- meta %>%
  filter(tissue != "CED")
#meta$genotype <-   paste(meta$genotype, meta$RT_status, sep = "_")

counts <- counts[, rownames(meta)]
#########
meta$tissue <- factor(meta$tissue)
meta$genotype <- factor(meta$genotype)
meta$RT_status <- factor(meta$RT_status)

meta$genotype <- relevel(meta$genotype, ref = "NTC")
meta$tissue <- relevel(meta$tissue, ref = "in.vivo")
meta$RT_status <- relevel(meta$RT_status, ref = "noRT")

head(meta)
metadata <- meta
metadata$sample <- rownames(metadata)
####################
model_formula <- "~ tissue * genotype * RT_status"
design <- model.matrix(~ tissue * genotype * RT_status, data = meta)
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")
keep <- filterByExpr(d0,design)
d0 <- d0[keep,]
cleanDev(); pdf(out("Design.pdf"),w=30,h=30)
pheatmap(design)
dev.off()


# voom transformation
dataVoom <- voom(d0, design, plot=TRUE)
cleanDev(); pdf(out("Voom_", "_Before.pdf"), w=6,h=6)
voom(d0, design, plot=TRUE)
dev.off()

head(dataVoom$E)
colnames(design) <- make.names(colnames(design))
limmaFit <- lmFit(dataVoom, design)
colnames(design)
non_estimable_coefs <-  nonEstimable(design)

# Remove non-estimable coefficients from the model matrix

estimable_coefs <- setdiff(colnames(design), non_estimable_coefs)
design_estimable <- design[, estimable_coefs, drop = FALSE]
# Re-fit model using estimable coefficients only
limmaFit_estimable <- lmFit(dataVoom, design_estimable)
limmaFit_estimable <- eBayes(limmaFit_estimable)
colnames(coef(limmaFit_estimable))

# Generate valid contrasts
genotypes_estimable <- grep("^genotype", colnames(design_estimable), value = TRUE)

genotypes <- grep("^genotype", colnames(design_estimable), value = TRUE)
genotypes <- genotypes[!grepl("RT_statusRT", genotypes)] # remove RT terms

contrast_flat_expr <- list()

for (i in seq_along(genotypes)) {
  g <- genotypes[i]
  
  # Build potential column names
  noRT <- paste0("tissueex.vivo.", g)
  RT_1 <- paste0("tissueex.vivo.", g)
  RT_2 <- paste0("tissueex.vivo.", g, ".RT_statusRT")
  
  # For noRT: only if the column exists
  if (noRT %in% colnames(design_estimable)) {
    contrast_flat_expr[[paste0(g, "_exvivo_vs_invivo_noRT")]] <- noRT
  }
  
  # For RT: only if both exist
  if (all(c(RT_1, RT_2) %in% colnames(design_estimable))) {
    contrast_flat_expr[[paste0(g, "_exvivo_vs_invivo_RT")]] <- paste(RT_1, RT_2, sep = " + ")
  }
}

# Add RT_status.ex.vivo contrast if both columns exist
if (all(c("tissueex.vivo.RT_statusRT", "RT_statusRT") %in% colnames(design_estimable))) {
  contrast_flat_expr[["RT_status_ex.vivo"]] <- "tissueex.vivo.RT_statusRT + RT_statusRT"
}

# Finally make the contrast matrix
contrast_matrix <- makeContrasts(contrasts = contrast_flat_expr, levels = design_estimable)


# 2️⃣ Make contrast matrix
if (length(contrast_flat_expr) > 0) {
  contrast_matrix <- makeContrasts(contrasts = contrast_flat_expr, levels = design_estimable)
  
  # 3️⃣ Fit contrasts
  limmaFit.contrast <- contrasts.fit(limmaFit_estimable, contrast_matrix)
  limmaFit.contrast <- eBayes(limmaFit.contrast)
  
  # 4️⃣ Extract results for each contrast
  limmaRes.contrast <- map_dfr(colnames(contrast_matrix), function(contrast_name) {
    topTable(limmaFit.contrast, coef = contrast_name, number = Inf) %>%
      rownames_to_column("ensg") %>%
      mutate(coef = contrast_name,
             group = case_when(
               logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
               logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
               TRUE ~ "n.s")) %>%
      dplyr::select(-contains("Intercept"))
  }, .id = "coefficient")
  
  # 5️⃣ Extract results for all main model coefficients
  limmaRes <- map_dfr(colnames(limmaFit_estimable$coef), function(coef_name) {
    topTable(limmaFit_estimable, coef = coef_name, number = Inf) %>%
      rownames_to_column("ensg") %>%
      mutate(coef = coef_name,
             group = case_when(
               logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
               logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
               TRUE ~ "n.s")) %>%
      dplyr::select(-contains("Intercept"))
  }, .id = "coefficient")
  
  # 6️⃣ Combine contrast results with main results
  limmaRes <- rbind(limmaRes.contrast, limmaRes)
  
} else {
  cat("No valid contrasts found, proceeding with the main results only.\n")
}
# 2️⃣ Make contrast matrix
if (length(contrast_flat_expr) > 0) {
  contrast_matrix <- makeContrasts(contrasts = contrast_flat_expr, levels = design_estimable)
  
  # 3️⃣ Fit contrasts
  limmaFit.contrast <- contrasts.fit(limmaFit_estimable, contrast_matrix)
  limmaFit.contrast <- eBayes(limmaFit.contrast)
  
  # 4️⃣ Extract results for each contrast
  limmaRes.contrast <- map_dfr(colnames(contrast_matrix), function(contrast_name) {
    topTable(limmaFit.contrast, coef = contrast_name, number = Inf) %>%
      rownames_to_column("ensg") %>%
      mutate(coef = contrast_name,
             group = case_when(
               logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
               logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
               TRUE ~ "n.s")) %>%
      dplyr::select(-contains("Intercept"))
  }, .id = "coefficient")
  
  # 5️⃣ Extract results for all main model coefficients
  limmaRes <- map_dfr(colnames(limmaFit_estimable$coef), function(coef_name) {
    topTable(limmaFit_estimable, coef = coef_name, number = Inf) %>%
      rownames_to_column("ensg") %>%
      mutate(coef = coef_name,
             group = case_when(
               logFC >= 1 & adj.P.Val <= 0.05 ~ "up",
               logFC <= -1 & adj.P.Val <= 0.05 ~ "down",
               TRUE ~ "n.s")) %>%
      dplyr::select(-contains("Intercept"))
  }, .id = "coefficient")
  
  # 6️⃣ Combine contrast results with main results
  limmaRes <- rbind(limmaRes.contrast, limmaRes)
  
} else {
  cat("No valid contrasts found, proceeding with the main results only.\n")
}
unique(limmaRes$coef)

limmaRes %>%
  write_rds(out("limmaRes_threeway.rds"))

#Effect of culture on Radiotherapy
RT_effect_ex.vsin <- limmaRes %>%
  filter(coef == "tissueex.vivo.RT_statusRT")
##############################################
#plots controls

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
Housekeeping <-  c("Gapdh","Tbp","Ubc","Actb","Pgk1")
prefix <- "glioblastoma"
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% Housekeeping) %>%             # keep only your gene set
  group_by(genes) %>%                          # calculate z-score per gene across all samples
  mutate(
    mean_expr_gene = mean(Expression, na.rm = TRUE),
    sd_expr_gene   = sd(Expression, na.rm = TRUE),
    zscore         = (Expression - mean_expr_gene) / sd_expr_gene
  ) %>%
  ungroup() %>%
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
ggsave(out("combined_plot_housekeeping.pdf"), w= 30, h =10)

##############################
# Modular GEO RNA-seq pipeline
##############################
#Aarathy
###############
source("src/00_init.R")
source("src/Ag_enrichR_mouse_genes.R")
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

out <- dirout("Glioblastoma_limmaRes.3way.mod")
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


#########
meta$tissue <- factor(meta$tissue)
meta$genotype <- factor(meta$genotype)
meta$RT_status <- factor(meta$RT_status)

meta$genotype <- relevel(meta$genotype, ref = "NTC")
meta$tissue <- relevel(meta$tissue, ref = "in.vivo")
meta$RT_status <- relevel(meta$RT_status, ref = "noRT")


metadata <- meta
metadata$sample <- rownames(metadata)

# Assuming your metadata is in a dataframe called `metadata`
# and contains columns: genotype, RT_status, tissue

# Step 1: summarize which conditions exist per genotype
geno_summary <- metadata %>%
  group_by(genotype, tissue, RT_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = c(tissue, RT_status),
    values_from = n,
    values_fill = 0
  )

# Step 2: keep only genotypes with both RT + noRT in both tissues
valid_genotypes <- geno_summary %>%
  filter(
    in.vivo_noRT  > 0,
    in.vivo_RT    > 0,
    ex.vivo_noRT  > 0,
    ex.vivo_RT    > 0
  ) %>%
  pull(genotype)

# Step 3: subset metadata to only those genotypes
metadata <- metadata %>%
  filter(genotype %in% valid_genotypes)

# Step 4 (optional): check how many genotypes remain
message("✅ Genotypes with RT and noRT in both in.vivo and ex.vivo:")
print(valid_genotypes)
counts <- counts[,rownames(metadata)]

all(colnames(counts)%in% rownames(metadata))
all(rownames(metadata)%in%colnames(counts))
####################
model_formula <- "~ tissue * genotype * RT_status"
design <- model.matrix(~ tissue * genotype * RT_status, data = metadata)
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
  Tex.GenX.noRT <- paste0("tissueex.vivo.", g)
  Tex.GenX.RT   <- paste0("tissueex.vivo.", g, ".RT_statusRT")
  GenX.RT       <- paste0(g, ".RT_statusRT")
  GenX          <- paste0(g)
  
  # invivo_noRT
  if (GenX %in% colnames(design_estimable)) {
    cname <- paste0(g, "_invivo_noRT")
    contrast_flat_expr[[cname]] <- GenX
    message("Added contrast: ", cname)
  }
  
  # exvivo_noRT
  if (all(c(GenX, Tex.GenX.noRT) %in% colnames(design_estimable))) {
    cname <- paste0(g, "_exvivo_noRT")
    contrast_flat_expr[[cname]] <- paste(GenX, Tex.GenX.noRT, sep = "+")
    message("Added contrast: ", cname)
  }
  
  # invivo_RT
  if (all(c(GenX, GenX.RT) %in% colnames(design_estimable))) {
    cname <- paste0(g, "_invivo_RT")
    contrast_flat_expr[[cname]] <- paste(GenX, GenX.RT, sep = "+")
    message("Added contrast: ", cname)
  }
  
  # exvivo_RT
  if (all(c(GenX, Tex.GenX.noRT, Tex.GenX.RT, GenX.RT) %in% colnames(design_estimable))) {
    cname <- paste0(g, "_exvivo_RT")
    contrast_flat_expr[[cname]] <- paste(GenX, GenX.RT, Tex.GenX.noRT, Tex.GenX.RT, sep = "+")
    message("Added contrast: ", cname)
  }
}



# Add RT_status.ex.vivo contrast if both columns exist
if (all(c("tissueex.vivo.RT_statusRT", "RT_statusRT") %in% colnames(design_estimable))) {
  contrast_flat_expr[["RT_status_ex.vivo"]] <- "tissueex.vivo.RT_statusRT + RT_statusRT"
}

# Finally make the contrast matrix
contrast_matrix <- makeContrasts(contrasts = contrast_flat_expr, levels = design_estimable)

unique(coef(limmaFit.contrast))
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
# Create a named vector of renames based on your genotypes
coef_rename <- c()


# Step 1 — find genotypes that appear in all 4 contrast types
# Define your required conditions
required_patterns <- c(
  "invivo_noRT", "exvivo_noRT",
  "invivo_RT", "exvivo_RT"
)

genotypes_alone <- gsub("genotype","", genotypes)
# Step 3 — Build a mapping table for renaming
coef_rename <- c()

for (g in genotypes_alone) {
  invivo_noRT <- paste0("genotype", g)
  exvivo_noRT <- paste0("genotype", g, "+tissueex.vivo.genotype", g)
  invivo_RT   <- paste0("genotype", g, "+genotype", g, ".RT_statusRT")
  exvivo_RT   <- paste0("genotype", g, "+genotype", g, 
                        ".RT_statusRT+tissueex.vivo.genotype",
                        g, "+tissueex.vivo.genotype", g, ".RT_statusRT")
  
  coef_rename[invivo_noRT] <- paste0(g, "_invivo_noRT")
  coef_rename[exvivo_noRT] <- paste0(g, "_exvivo_noRT")
  coef_rename[invivo_RT]   <- paste0(g, "_invivo_RT")
  coef_rename[exvivo_RT]   <- paste0(g, "_exvivo_RT")
}

# Step 4 — rename coefficients robustly (exact string match)
limmaRes$coef <- plyr::mapvalues(
  limmaRes$coef,
  from = names(coef_rename),
  to = coef_rename,
  warn_missing = FALSE
)

# Step 5 — verify
unique(limmaRes$coef)

limmaRes$coef <- gsub("exvivo","ex.vivo",limmaRes$coef)
limmaRes$coef <- gsub("invivo","in.vivo",limmaRes$coef)
limmaRes %>%
  write_rds(out("limmaRes_threeway.rds"))

#Effect of culture on Radiotherapy
RT_effect_ex.vs.in <- limmaRes %>%
  filter(coef %in% c("tissueex.vivo.RT_statusRT",
         "RT_statusRT","tissue.ex.vivo","RT_status_ex.vivo"))



# Subset the coefficients of interest
rt_effects <- limmaRes %>%
  filter(coef %in% c("tissueex.vivo", "tissueex.vivo.RT_statusRT"))

# Spread to wide format for correlation calculation
rt_wide <- rt_effects %>%
  dplyr::select(ensg, coef, logFC) %>%   # adjust 'gene' column name if different
  tidyr::pivot_wider(names_from = coef, values_from = logFC)

# Calculate correlation
correlation <- cor(abs(rt_wide$tissueex.vivo), abs(rt_wide$tissueex.vivo.RT_statusRT), use = "pairwise.complete.obs")

correlation
#####
gsea.res <- data.table()

for (dbx in names(enr.terms)) {
    subset_limmaRes_NTC <- limmaRes[limmaRes$coef == "tissueex.vivo", ]
    stats <- with(subset_limmaRes_NTC, setNames(logFC, nm = ensg))
    
    if (any(is.na(stats))) {
      next
    }
    
    fgsea_output <- fgsea(pathways = enr.terms[[dbx]], stats = stats)
    
    if (length(fgsea_output) > 0) {
      gsea.res <- rbind(gsea.res, data.table(fgsea_output,  db = dbx))
    }
  }
###
gsea.res %>% write_rds(out("gsea.res.ntc.rds"))
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge,
                                      function(vec) paste(vec[1:10], collapse = ","))
dbx<-"MSigDB_Hallmark_2020"
pDT <- gsea.res[db == dbx]
## Splitting the task to handle both ends of the NES spectrum-positive and negative
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=4),
                                                       by=c("celltype")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=4),
                                                      by=c("celltype")]$pathway)

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

# Step 3: Plot with the new pathway order (highest NES first)
Fig1C <- ggplot(pDT, aes(y=db, x=pathway, color = pmin(pmax(NES, -2), 2), size=pmin(5, -log10(padj)))) +
  
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name=TeX("NES"))+
  #name=TeX("log_{2}(FC)"))+
  geom_point() +
  scale_size_continuous(
    range = c(0, 1.8),
    #limits = c(0, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  
  
  #xRot() +
  #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
  labs(y = "Cell type",
       x = "Pathways",
       title = "Enriched pathways")+
  
  #coord_flip()+
  optimized_theme_fig()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,
  ),
  legend.position = "right", legend.direction = "vertical",
  legend.justification = "bottom")

Fig1C

#######
# Filter for the interferon pathways
interferon_genes <- gsea.res %>%
  filter(pathway %in% c("Interferon Alpha Response", "Interferon Gamma Response")) %>%
  pull(leadingEdge)%>%
  unlist() %>%
  trimws() %>%                                      # remove whitespace
  unique() 
  
head(interferon_genes)

##############################################
#plots controls
head(metadata)
longer_dataVoom <-  dataVoom$E %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(metadata, by ="sample") %>%
  filter(genotype == "NTC")

longer_dataVoom$organ <- "glioblastoma_liu"
longer_dataVoom$author <- "Liu S. J. et al., Genome Biology, 2024"

longer_dataVoom %>%
  dplyr::select(c("genes","tissue","sample","Expression","author","organ"))%>%
  write_rds(out("Glioblastoma_mouse.rds.rds"))



prefix <- "glioblastoma"
longer_dataVoom_zscore <- longer_dataVoom %>%
  filter(genes %in% ISGs) %>%             # keep only your gene set
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
  filter(genes %in% ISGs)


summary_df <- plot_df %>%
  group_by(sample, tissue) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop")

ggplot(summary_df, aes(x = reorder(sample, mean_z), y = mean_z, fill = tissue)) +
  geom_col() +
  labs(title = "interferon",
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
ggsave(out("combined_plot_interferon.pdf"), w= 30, h =10)
intersect(ISGs, interferon_genes)
intersect(ISG_core, interferon_genes)

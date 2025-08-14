source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
library(tidyverse)
library(enrichR)
library(purrr)
library("scales")
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library("ggridges")

#
InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
InDir3 <-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
InDir4 <- dirout("Figure2_paper")
InDir5 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
basedir <- dirout("Figure4_Supplementary")
#
limmaRes <- read_rds(InDir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))%>%
  mutate(genes = ensg)
limmaRes_NTC <- read_rds(InDir2("limma_perCTex.vivovsin.vivo.rds"))
merged_data <- limmaRes %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)
meta <- fread(InDir5("meta_cleaned.tsv")) # Read data
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   

meta <- meta[, -1, drop = FALSE] 
colnames(meta) <- gsub("rowname","sample1", colnames(meta))
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop")%>%
  mutate(coef = genotype)

replicates_per_ko <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(
    valid_ko = any(valid_ko),
    total_in_vivo = sum(in.vivo, na.rm = TRUE),
    total_ex_vivo = sum(ex.vivo, na.rm = TRUE)
    , .groups = "drop") %>%
  mutate(coef = genotype)

selected_KOs <- meta %>%
  group_by(genotype, tissue, celltype) %>%                  # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample1), .groups = 'drop') %>% # Count distinct samples for each group
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>% # Spread tissue to separate columns (in.vivo and ex.vivo)
  group_by(genotype) %>%                                    # Regroup by genotype
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) %>%              # Keep genotypes that have at least one celltype with 3+ samples in both tissues
  pull(genotype) %>% unique()

adj_p_cutoff <- 0.05
logfc_cutoff <- 1

summary_df <- limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")


count_threshold = 10
coefficients  <-  summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)%>%
  pull(coef)%>%
  unique()
correlation_deg <- read_rds(InDir4("correlation_deg.rds"))


koi <- Reduce(intersect, list(selected_KOs,  coefficients)) #only valid_ko

#Supplementary-------------------3.1
# Step 1: Calculate observed correlations with p-values
observed_correlations <- merged_data %>%
  filter(coef %in% koi) %>%
  merge(ko_flags, by = c("coef", "celltype")) %>%
  filter(valid_ko)%>%
  group_by(coef, celltype) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,
    .groups = 'drop'
  ) %>%
  mutate(type = "Observed")  # Label as observed correlations

# Step 2: Shuffle values and calculate correlations again
set.seed(42)  # Set seed for reproducibility
shuffled_correlations <- merged_data %>%
  merge(ko_flags, by = c("coef", "celltype")) %>%
  filter(valid_ko)%>%
  filter(coef %in% koi) %>%
  group_by(coef, celltype) %>%
  mutate(shuffled_logFC_KO = sample(logFC_KO)) %>%  # Shuffle logFC_KO within each group
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(shuffled_logFC_KO), method = "pearson"),
    .groups = 'drop'
  ) %>%
  mutate(type = "Shuffled")  # Label as shuffled correlations

# Step 3: Combine observed and shuffled correlations
combined_correlations <- bind_rows(
  observed_correlations %>% dplyr::select(coef, celltype, cor_abs, type),
  shuffled_correlations %>% dplyr::select(coef, celltype, cor_abs, type)
)

Sup.Fig.4A <- ggplot(combined_correlations, aes(x = cor_abs, y = celltype, fill = type)) +
  geom_density_ridges(alpha = 0.4, scale = 1) +
  labs(
    title = "Correlation of culture effects to observed vs 
shuffled interaction effects across KOs",
    x = "Correlation",
    y = "Cell type",
    fill = "Correlation Type"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Observed" = "blue", "Shuffled" = "gray"))+
  optimized_theme_fig()


ggsave(basedir("Sup.Fig.4A.pdf"),plot = Sup.Fig.4A,
       w=8,h=4, units = "cm")


###
#Supplementary-----------------
#
InDir7  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide")
gsea.res <- read_rds(InDir7("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
gsea.res$coef <- gsub("interaction","",gsea.res$coef )
# Step 2: Summarize to find KOs with at least one valid cell type
valid_ko_summary <- ko_flags %>%
  group_by(genotype) %>%
  summarize(has_valid_celltype = any(valid_ko), .groups = 'drop')


# Step 3: Filter GSEA results based on valid KOs

db = "MSigDB_Hallmark_2020"
create_gsea_plot <- function(db) {
  
  # Filter the data for the given database
  pDT <- gsea.res %>%
    filter(coef %in% koi) %>%
    filter(db == !!db) %>%  # Correctly reference the current database
    left_join(ko_flags, by = c("coef", "celltype")) %>%
    left_join(summary_df, by = c("coef", "celltype")) %>%
    filter(Count > 5) %>%
    filter(valid_ko)%>%
    dplyr::select(-c("Count", "Regulation"))
  # Merge with KO flags
  
  # Step 3 continued: Keep only valid KOs for the specific cell type
  pDT <- pDT %>% filter(valid_ko == TRUE, padj < 0.05)
  
  # Select the pathways for plotting (both positive and negative)
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n = 7), by = c("coef", "celltype", "pathway")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n = 7), by = c("coef", "celltype", "pathway")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(union(pw.display.pos, pw.display.neg))
  # Remove duplicate rows
  pDT <- pDT %>% distinct()
  
  # Filter pDT to include only selected pathways
  # **Convert list-columns to character format**
  dat <- pDT %>%
    mutate(across(where(is.list), ~sapply(., toString)))  # Converts lists to comma-separated strings
  
  # Save the filtered data table
  write.table(dat, file = basedir(paste0("fgsea_", db, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # Create the plot for the current database
  fig <- ggplot(pDT, aes(x = coef, y = pathway, color = NES, size = pmin(5, -log10(padj)))) +
    geom_point() + 
    scale_color_gradient2(low = "#4C889C", mid = "white", high = "#D0154E", name = TeX("NES")) +
    geom_point(data = pDT[padj < 0.05], shape = 1) +
    scale_size_continuous(range = c(0, 2), limits = c(0, 5), name = TeX("$-\\log_{10}(p_{adj})$")) +
    theme_bw() +
    xRot() +
    labs(x = "KOs") +
    facet_grid(cols = vars(celltype), scales = "free", space = "free") +  # Create facets for each cell type
    optimized_theme_fig()+
    theme(strip.text.x = element_text(angle = 90))
  #strip.text.y = element_text(angle = 0)) +
  
  # Save the plot for the current database
  ggsave(basedir("fig2.2_supplementary_fgsea", db, "_3rep.pdf"), fig,
         w =11, h =12, units = "cm")
  
  
  
}

# Example usage: 
# Create the plot for "MSigDB_Hallmark_2020"
create_gsea_plot("MSigDB_Hallmark_2020")
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
# You can also loop through the databases to generate plots for all of them
lapply(databases, create_gsea_plot)
##

#only top 10 genes per pathway
combined_results %>%
  filter(pathway == "TNF") %>%
  pull(ensg) %>%
  unique()

#Sup.Fig.4B--------------
# Filter for Smc3 and HSC
scatter_data <- merged_data %>%
  filter(coef %in% c("Brd9","Wdr82","Cebpa"), celltype == "GMP")

# Plot
Sup.Fig.4Ba <- ggplot(scatter_data, aes(x = abs(logFC_NTC), y = abs(logFC_KO))) +
  geom_hex(bins = 50) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b", name = "Gene count") + 
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c") +
  labs(
    
    x = expression(atop(atop("Absolute interaction effect", paste(log[2](FC))), "Eo.Ba")),
    y = expression(atop("Absolute culture effect", 
                        paste( log[2](FC))))
    
  ) +
  facet_grid(cols = vars(coef))+
  optimized_theme_fig()
ggsave(basedir("Sup.Fig.4Ba.pdf"),plot = Sup.Fig.4Ba,
       w=12,h=5, units = "cm")
#########
scatter_data_2 <- merged_data %>%
  filter(coef %in% c("Brd9","Wdr82","Smc3"), celltype == "HSC")

# Plot
Sup.Fig.4Bb <- ggplot(scatter_data_2, aes(x = abs(logFC_NTC), y = abs(logFC_KO))) +
  geom_hex(bins = 50) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b", name = "Gene count") + 
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c") +
  labs(
    x = expression(atop(atop("Absolute interaction effect", paste(log[2](FC))), "HSC")),
    y = expression(atop("Absolute culture effect", 
                        paste( log[2](FC))))
    
  ) +
  facet_grid(cols = vars(coef))+
  optimized_theme_fig()
ggsave(basedir("fig4Bb.pdf"),plot = Sup.Fig.4Bb,
       w=12,h=5, units = "cm")

second_row <- (Sup.Fig.4Ba / Sup.Fig.4Bb) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Predicted vs actual ko effects",
    theme = theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.title = element_text(hjust = 0.5, size = 7, face = "bold", color = "black")
    )
  )

ggsave(basedir("Sup.Fig.4B.pdf"),plot = second_row,
       w=12,h=11, units = "cm")
#########


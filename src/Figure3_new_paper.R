source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_Optimized_theme.R")
source("src/Ag_top_genes_per_pathway.R")
library(tidyverse)
library(enrichR)
library(purrr)
library("scales")
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
#directories ------
base<-"Figure3_new_paper"
basedir <- dirout("Figure3_new_paper")
#Data and function
Indir1 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/")
InDir5 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
InDir3 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
InDir4 <- dirout("Figure2_paper")
InDir6 <- dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_per_pathway_fgsea_in.vivo")
#
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
enr.terms <- enrichrGetGenesets(databases)

# Convert to mouse --------------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl){
  dbl <- lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})
#

#Fig3.1------------(DEGS)
InDir2 <- dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide/")


#gsea.res <- read_rds(InDir2("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
# mutate(coef = gsub("interaction","",coef))

limmaRes <- read_rds(InDir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
################################################################################
#gene_based analysis

# Step 1: Filter for significant genes (adj.P.Val < 0.05 and abs(logFC) > 1)
limmaRes_significant <- limmaRes %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)  # Only significantly altered genes


meta <- fread(InDir5("meta_cleaned.tsv")) # Read data
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   

meta <- meta[, -1, drop = FALSE] 
colnames(meta) <- gsub("rowname","sample1", colnames(meta))
# Check if there are at least 2 distinct samples per tissue for each genotype and celltype
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

KO_list <- correlation_deg %>% filter(num_degs >= 10 & correlation < 0.5) %>%
  pull(genotype)%>%
  unique()
koi <- Reduce(intersect, list(selected_KOs,  KO_list, coefficients)) #only valid_ko


data <- summary_df %>%
  filter(coef %in% koi) %>%
  #group_by(celltype, coef) %>%
  # Filter for counts > 10 and KOs of interest only
  filter(Count >= 10) %>%
  inner_join(ko_flags, by = c("celltype", "coef")) %>%
  # Filter only valid KOs for each cell type
  filter(valid_ko)
data %>% write_rds(basedir("kos_fig3.rds"))
# Plot
#Fig3.1 DEG interaction logFC ---------------------------------------------
#fig ----
# Fig3.1 <- ggplot(data,aes(
#   x = coef,
#   y = ifelse(Regulation == "Downregulated", -log10(Count), log10(Count)),
#   fill = Regulation
# )) +
#   geom_col() +
#   # Custom colors for upregulated and downregulated bars
#   scale_fill_manual(values = c("Upregulated" = "#D0154E", "Downregulated" = "#4C889C")) +
#   # Facet by celltype with free space for flexibility in cell widths
#   facet_grid(cols = vars(celltype), space = "free", scales = "free") +
#   labs(
#     #title = "No. of Genes with interaction effect of culture condition in KOs",
#     x = NULL,
#     y = "log10(Number of Genes)"
#   ) +
#   # Custom theme with no legend if not needed
#   optimized_theme_fig() + 
#   theme(legend.position = "right",
#         axis.text.x = element_blank(),
#         strip.text.x = element_text(angle = 90,
#                                     hjust = 0)
#   )

summary_total <- data %>%
  group_by(celltype, coef, genotype, valid_ko) %>%
  summarise(Total_Regulated = sum(Count), .groups = "drop")

Fig3.1 <- ggplot(summary_total,aes(
  x = coef,
  y = log10(Total_Regulated)
)) +
  geom_col(fill = "#b3b3b3ff") +
  # Custom colors for upregulated and downregulated bars
  #scale_fill_manual(values = c("Upregulated" = "#D0154E", "Downregulated" = "#4C889C")) +
  # Facet by celltype with free space for flexibility in cell widths
  facet_grid(cols = vars(celltype), space = "free", scales = "free") +
  labs(
    #title = "No. of Genes with interaction effect of culture condition in KOs",
    x = NULL,
    y = TeX("$\\log_{10}\\; (No.of\\; genes\\; with\\; interaction\\; effects)$")
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
        axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90,
                                    hjust = 0))


# Display the plot
Fig3.1

ggsave(basedir("Fig3.1_3rep.pdf"),plot=Fig3.1,
       w=10,h=4, units = "cm")

#

#atleast 3 rep in both tissues, 
#Fig3.1---------------------------(Correllation)----------
InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
InDir3 <-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")

limmaRes_int <- read_rds(InDir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))%>%
  mutate(genes = ensg)
limmaRes_NTC <- read_rds(InDir2("limma_perCTex.vivovsin.vivo.rds"))
merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)




# Filter for Smc3 and HSC
scatter_data <- merged_data %>%
  filter(coef %in% c("Brd9","Wdr82","Cebpa"), celltype == "GMP")

# Plot
scatter_example <- ggplot(scatter_data, aes(x = abs(logFC_NTC), y = abs(logFC_KO))) +
  geom_hex(bins = 50) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") + 
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c") +
  labs(
    #title = "Correlation of logFC: NTC vs KO (Brd9, HSC)",
    x = TeX("Absolute $\\log_{2}(\\FC)\\; of\\; interaction\\; effect$"),
    y = TeX("Absolute $\\log_{2}(\\FC)\\; of\\; culture\\; effect$")
  ) +
  facet_grid(cols = vars(coef))+
  optimized_theme_fig()
ggsave(basedir("fig4.2_example_scatter_3rep.pdf"),plot = scatter_example,
       w=12,h=5, units = "cm")
#########
scatter_data_2 <- merged_data %>%
  filter(coef %in% c("Brd9","Wdr82","Smc3"), celltype == "HSC")

# Plot
scatter_example_2 <- ggplot(scatter_data_2, aes(x = abs(logFC_NTC), y = abs(logFC_KO))) +
  geom_hex(bins = 50) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") + 
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c") +
  labs(
    #title = "Correlation of logFC: NTC vs KO (Brd9, HSC)",
    x = TeX("Absolute $\\log_{2}(\\FC)\\; of\\; interaction\\; effect$"),
    y = TeX("Absolute $\\log_{2}(\\FC)\\; of\\; culture\\; effect$")
  ) +
  facet_grid(cols = vars(coef))+
  optimized_theme_fig()
ggsave(basedir("fig4.2_example_scatter_3rep_2.pdf"),plot = scatter_example_2,
       w=12,h=5, units = "cm")
# Step 1: Calculate correlation with p-value for each KO and celltype
correlation_results <- merged_data %>%
  inner_join(ko_flags, by = c("coef","celltype")) %>%
  filter(valid_ko) %>%
  group_by(coef, celltype) %>%
  filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results <- correlation_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(coef %in% koi)

# Step 3: Add significance labels based on p-value thresholds
correlation_results <- correlation_results %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )
# Step 4: Plot the correlation results with significance annotations
Fig3.2a <- ggplot(correlation_results, aes(x = coef, y = cor_abs)) +
  geom_bar(stat = "identity", position = "dodge", fill ="#b3b3b3ff") +
  geom_text(aes(label = significance), 
            vjust = -0.4, 
            color = "black",
            size = 2.5) + 
  
  labs(
    #title = TeX("$Correlation of $\\log_{2}$(\\FC)$(ex-vivo vs in-vivo in NTC) vs KO-tissue interaction"),
    
    x = "KO",
    y = TeX("$Correlation\\;of\\;interaction\\;and\\;culture\\;effects$")
    ,
    fill = "Cell Type"
  ) +
  facet_grid(cols = vars(celltype), scales = "free_x", space = "free_x") +
  optimized_theme_fig()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        strip.text.x = element_blank()) 
# combine with enrichment to ntc
#Fig3.3
enrichment_ntc_in.vivo <- read_rds(InDir6("enrichment_to_NTC_genes.rds"))
combined <- data %>%
  dplyr::select(coef,celltype)%>%
  distinct()%>%
  left_join(enrichment_ntc_in.vivo, by = c("coef","celltype"))

# Step 1: Combine the correlation data and enrichment data
combined <- combined %>%
  mutate(significance_en = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))
combined <- combined %>%
  mutate(log2.odds.ratio =log2(odds.ratio))

combined <- combined %>%
  mutate(log2.odds.ratio = case_when(
    coef == "Smc3" ~ 0,
    TRUE ~ odds.ratio
  ))



Fig3.2b <- ggplot(combined, aes(x = coef, y = pmin(log2.odds.ratio, 5))) +  # Capping at 5
  geom_bar(stat = "identity", position = position_dodge(), fill ="#b3b3b3ff") +
  geom_text(aes(label = significance_en), vjust = -0.4, size = 2.5) +
  facet_grid(cols = vars(celltype), scales = "free_x", space = "free_x") +
  labs(
    x = "KO",
    y = TeX("$\\log_{2}\\; (Odds.ratio) (y-axis\\capped\\at\\5)$"),
    title = ""  # Add capping info to title
  ) +
  scale_y_continuous(limits = c(-1, 5)) +  # Set y-axis limit (just for visualization)
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text.x = element_blank()
  )





fig3.1_2 <-  Fig3.2a / Fig3.1 /  Fig3.2b
#fig3.1_2 <- Fig3.1 / Fig3.2 + plot_layout(heights = c(1, 1.3))
fig3.1_2
ggsave(basedir("Fig3.1_2.pdf"),fig3.1_2, w = 15,h = 10, units = "cm")
#######################
#density-------------
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

fig4.1_supplementary <- ggplot(combined_correlations, aes(x = cor_abs, y = celltype, fill = type)) +
  geom_density_ridges(alpha = 0.4, scale = 1) +
  labs(
    title = "Observed vs Shuffled Correlations by Cell Type",
    x = "Absolute Correlation (abs logFC)",
    y = "Cell Type",
    fill = "Correlation Type"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Observed" = "blue", "Shuffled" = "gray"))+optimized_theme_fig()


ggsave(basedir("fig4.1_supplementary_density_ridge_3rep.pdf"),plot = fig4.1_supplementary,
       w=8,h=4, units = "cm")



#significance
# Perform statistical testing for each KO and tissue separately

#fig3.2-------------
InDir5 <- dirout("Figure1")
# grouping gene sets specifically for ISG and cholesterol
genes_fig_1 <- read_rds(InDir5("genes_fig1.rds"))
colnames(genes_fig_1) <- c("ensg","pathway")
genes_fig_1$pathway <- gsub("mTORC1_or_Cholesterol","mTORC1/Cholesterol", genes_fig_1$pathway)

filtered_genes <- limmaRes %>%
  inner_join(summary_df,by =c("coef","celltype"))%>%
  filter(Count >10)%>%
  #filter(group != "n.s")%>%
  filter(ensg %in% genes_fig_1$ensg, coef %in% koi) %>%
  group_by(celltype, coef) %>%
  merge(genes_fig_1, by = "ensg") %>%
  left_join(ko_flags, by = c("coef" = "genotype", "celltype")) %>%  # Merge with KO flags per cell type
  filter(valid_ko == TRUE)  # Keep only valid KOs for the specific cell type
unique(filtered_genes$pathway)
# Recode pathways for better labeling


# fig---------
Fig3.2 <- ggplot(filtered_genes %>%
                   filter(pathway == "ISG core"), aes(x = coef, y = ensg,
                                     color = pmin(2, pmax(-2, logFC)) ,
                                     size = pmin(3, -log10(adj.P.Val))
)) +  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(1,3),
    #limits = c(0,5),
    #breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(#title = "ISG core genes",
       x = "KOs",
       y = "Genes")+
  facet_grid(cols = vars(celltype), rows = vars(pathway), scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()+theme(
    legend.position = "right",
    strip.text.x = element_text(angle = 90, hjust = 0)
  )

Fig3.2

#paper-----------
#paper--------------
ggsave(
  filename = basedir("Fig3.2_ISG_core.pdf"),
  plot = Fig3.2,
  width = 11,
  height = 8.5,
  units = "cm"
)
ggsave(
  filename = basedir("Fig3.2_without_legend.pdf"),
  plot = Fig3.1+theme(
    legend.position = "none"
  ),
  width = 17,
  height = 7,
  units = "cm"
)
#
#fig 3.2---------------
#

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
#################################################################################
# # Pathways of interest (you can adjust these as per your needs)
# get_top_genes <- function(pathway_name,
#                           pathway_genes,
#                           limma_results,
#                           logFC_threshold = 1,
#                           pval_threshold = 0.05,
#                           top_n = 5) {
#   # Filter limma results for the pathway genes
#   
#   filtered_genes <- limma_results %>%
#     filter(group != "n.s") %>%  # Skip non-significant genes
#     filter(adj.P.Val < pval_threshold, abs(logFC) > logFC_threshold) %>%
#     filter(toupper(ensg) %in% toupper(pathway_genes)) %>%  # Match pathway genes
#     group_by(celltype, coef) %>%
#     arrange(adj.P.Val) %>% 
#     arrange(desc(abs(logFC)))%>%# Sort by adjusted p-value
#     slice_head(n = top_n) %>%  # Get the top 5 genes per group (celltype)
#     pull(ensg) %>% 
#     unique()
#   
#   # Filter the main table to return data for plotting
#   pathway_plot <- limma_results %>%
#     filter(toupper(ensg) %in% toupper(filtered_genes), coef %in% koi)%>%
#     mutate(pathway = pathway_name)
#   
#   return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
# }
# 
# pathway_list <- list(
#   "EMT" = enr.terms$MSigDB_Hallmark_2020$`Epithelial Mesenchymal Transition`,
#   "ROS" = enr.terms$MSigDB_Hallmark_2020$`Reactive Oxygen Species Pathway`,
#   #"E2F targets" = enr.terms$MSigDB_Hallmark_2020$`E2F Targets`,
#   "TNF" = enr.terms$MSigDB_Hallmark_2020$`TNF-alpha Signaling via NF-kB`,
# # Glycolysis = enr.terms$MSigDB_Hallmark_2020$`Glycolysis`,
# "Protein_loc" = Reduce(union, list(
#   enr.terms$MSigDB_Hallmark_2020$`Myc Targets V1`,
#   enr.terms$GO_Biological_Process_2023$`Protein Import (GO:0017038)`,
#   enr.terms$GO_Biological_Process_2023$`Protein Insertion Into ER Membrane (GO:0045048)`,
#   enr.terms$GO_Biological_Process_2023$`Protein Insertion Into Membrane (GO:0051205)`,
#   enr.terms$GO_Biological_Process_2021$`RNA biosynthetic process (GO:0032774)`,
#   enr.terms$MSigDB_Hallmark_2020$`Protein Secretion`
# )),
# "Glycolysis" = Reduce(union, list(enr.terms$MSigDB_Hallmark_2020$`Glycolysis`,
# "NADH_metabolism" = enr.terms$GO_Biological_Process_2023$`NADH Dehydrogenase Complex Assembly (GO:0010257)`,
# enr.terms$GO_Biological_Process_2023$`NADH Metabolic Process (GO:0006734)`)))
# 
# # Loop through each pathway to extract and plot the top genes
# # by default, top 5 genes per celltype per KO are obtained 
# pathway_results <- lapply(names(pathway_list), function(pathway_name) {
#   pathway_genes <- pathway_list[[pathway_name]]
#   get_top_genes(pathway_name, pathway_genes, limmaRes)
# })
# 
# 
# # Combine results from all pathways into one dataframe
# pathway_plot_data <- bind_rows(lapply(pathway_results, function(x) x$pathway_plot))
# 
# filtered_results <- combined_results %>%
#   group_by(pathway, celltype, coef) %>%
#   arrange(adj.P.Val, desc(abs(logFC))) %>%  # Sort by p-value, then effect size
#   slice_head(n = 5) %>%  # Take the top 5 genes
#   ungroup()%>%
#   dplyr::select(ensg, pathway)%>%
#   distinct()
# pathway_plot_data <- limmaRes %>%
#   filter(ensg %in% filtered_results$ensg)%>%
#   filter(coef %in% koi)
# pathway_plot_data <- left_join(filtered_results,pathway_plot_data, by = "ensg") %>%
#   left_join(ko_flags[,c("coef","celltype","valid_ko")], by = c("coef", "celltype"))%>%
#   filter(valid_ko)
# unique(pathway_plot_data$ensg)  
# 
# 
# # Plot the results
# Fig3.3.2_supplementary <- ggplot(pathway_plot_data, aes(x = coef, y = ensg,
#                                                        color = pmin(2, pmax(-2, logFC)),  # Clamp logFC between -2 and 2
#                                                        size = pmin(5, -log10(adj.P.Val))  # Set size by adjusted p-value
# )) +
#   geom_point() +  # Create points
#   scale_color_gradient2(
#     low = "#4C889C",
#     mid = "white",
#     high = "#D0154E",
#     name =TeX("$\\log_{2}\\; (FC)$")
#   ) +
#   scale_size_continuous(
#     range = c(0,2),
#     limits = c(0,3),
#     
#     name =TeX("$-\\log_{10}(p_{adj})$")
#   )+
#   labs(title = "Differentially expressed genesets",
#        x = "KOs",
#        y = "Genes")+
#   facet_grid(rows = vars(pathway), cols = vars(celltype),
#              scales = "free", space = "free") +
#   theme_bw() +
#   optimized_theme_fig()+theme(
#     axis.text = element_text(size=5)
#   )
# #+coord_flip()
# Fig3.3.2_supplementary
# 
# # Save the plot
# ggsave(basedir(paste0("Fig3_supplementary.pdf")), plot = Fig3.3.2_supplementary,
#        width =18,
#        height = 24.7, units = "cm" )
#
#######
#SUPP FIG 3.2from Figure1_Supplementary
################################################################################
#only top 10 genes per pathway
combined_results %>%
  filter(pathway == "TNF") %>%
  pull(ensg) %>%
  unique()


#
#load libraries and functions--------------------------------------------------
#
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
library(ggridges)
library(ggpubr)
#directories ------
#
base <- "Figure2_paper"
basedir <- dirout("Figure2_paper")

########################
ENRICHR <- dirout(paste0("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment_guide/ENRICHR"))
ENRICHR.DBS <-union(ENRICHR.DBS,
                    c("GO_Biological_Process_2021",
                      "TRRUST_Transcription_Factors_2019",
                      "Reactome_2022",
                      "GO_Molecular_Function_2023",
                      "GO_Biological_Process_2023",
                      "CellMarker_2024"))
enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
# # save(enr.terms, file=out("Genesets_Human.RData"))
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

ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)
IFN_genes = union(enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`,
                              enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`) 

IFN_genes <- union(ISG_core, IFN_genes)
IFN_genes <-IFN_genes[IFN_genes != ""]
#Data and function
InDir5 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
InDir3 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
limmaRes <- read_rds(InDir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))

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
replicates_per_ko %>% write_rds(basedir("replicates_per_ko.rds"))
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
correlation_deg <- read_rds(basedir("correlation_deg.rds"))
KO_list <- correlation_deg %>% filter(num_degs >= 10 & correlation < 0.5) %>%
  pull(genotype)%>%
  unique()
koi <- Reduce(intersect, list(selected_KOs,  KO_list, coefficients)) #only valid_ko
###########

#fig2.1-------------
InDir <- dirout("FIG_cluster_enrichment/")
fish <- read_rds(InDir("detailed_filtered_cluster_enrichment.rds"))
# Count occurrences of ex vivo and in vivo for each genotype
genotype_counts <- fish %>%
  group_by(mixscape_class, tissue) %>%
  summarise(n = n(), log2OR_positive = sum(log2OR != 0), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = c(n, log2OR_positive), values_fill = 0)

# Filter genotypes where:
# 1. ex vivo AND in vivo occur >1 times
# 2. At least one ex vivo and one in vivo have log2OR > 0
valid_genotypes <- genotype_counts %>%
  filter(n_ex.vivo > 1 & n_in.vivo > 1 & 
           log2OR_positive_ex.vivo != 0 & 
         log2OR_positive_in.vivo != 0) %>%
  pull(mixscape_class)

# Keep only the valid genotypes in the original dataset
fish <- fish %>%
  filter(mixscape_class %in% valid_genotypes)
# 
# fish$Clusters <- factor(fish$Clusters,
#                        levels = c("HSC", "MEP.early", "MkP",
#                                                    "GMP", "Gran.P", "Gran.",
#                                                    "Mono", "Eo.Ba"))
fig2.1 <- ggplot(fish, aes(y=Clusters, x=tissue, size=sig.frac, color=log2OR_cap)) + 
  themeNF(rotate=TRUE) +
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (Odds\\.Ratio)$")
    ) + 
  scale_size_continuous(name="sign.frac.", range=c(0,1.9)) +
  geom_point() +
  facet_grid(cols = vars(mixscape_class))+
  #geom_point(shape=1, color="lightgrey") +
  xlab("Gene") + ylab("Cell type")+
  optimized_theme_fig()+
  theme(strip.text.x = element_text(angle = 90,
                                    hjust = 0),
        legend.position = "right",
        panel.spacing = unit(0.05, "lines"))
# Display the combined plot
#paper----
fig2.1
ggsave(basedir("fig2.1_3rep.pdf"),plot =fig2.1,
       w=18.1, height = 5, 
       units = "cm")

#fig2.2------------
#fig ------------
# Read in the correlation matrix and DEG data
Indir1 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/")
merged_logFC <- read_rds(Indir1("in.vivo_ex.vivo_logFC.rds"))
correlation_matrix_data <- read_rds(Indir1("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(Indir1("DEGs_per_tissue.rds"))%>%
  select(-correlation)

# Prepare correlation data
correlation_matrix_no_na <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)

# Perform hierarchical clustering
hc_cols <- hclust(dist(t(correlation_matrix_no_na)), method = "ward.D2")
column_order <- colnames(correlation_matrix_data)[hc_cols$order]

# Reshape correlation data
correlation <- correlation_matrix_data %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(cols = colnames(correlation_matrix_data),
               names_to = "genotype",
               values_to = "correlation",
               values_drop_na = FALSE)


# Prepare DEG data
deg_plot_data <- deg_plot_data %>%
  na.omit() %>%
  group_by(genotype, celltype) %>%
  filter(num_degs == max(num_degs, na.rm = TRUE)) %>%
  dplyr::select(-condition) %>%
  unique()

# Merge DEG data with correlation data
correlation_deg <- inner_join(deg_plot_data, correlation,
                              by = c("celltype", "genotype"))
correlation_deg %>% write_rds(basedir("correlation_deg.rds"))
# Merge with KO flags to include valid KO status
correlation_deg_flagged <- correlation_deg %>%
  inner_join(ko_flags, by = c("genotype", "celltype")) %>%
  filter(valid_ko)%>%
  na.omit()%>%
  filter(valid_ko)
# Update the genotype factor levels for plotting
correlation_deg_flagged$genotype <- factor(correlation_deg_flagged$genotype,
                                           levels = column_order)

#fig-----
# Create the plot
Fig2.2 <- correlation_deg_flagged %>%
 #filter(genotype %in% koi) %>%
  ggplot() +
  geom_point(aes(
    x = genotype,
    y = celltype,
    size = pmin(3,log10(num_degs)),
    fill = correlation  # Set transparency based on KO validity
  ),
  shape = 21,           # Use shape 21 to enable fill and color
  color = "black",       # Black outline
  stroke = 0.5          # Adjust the width of the outline
  ) +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E"
  ) +
  scale_size_continuous(
    range = c(0,2.5),
    #limits = c(0,2.5),
    breaks = c(1,2,3),
    name =TeX("$\\log_{10}\\; (\\No.\\; of \\;DEGs)$")
    
  )+
  labs(x = "KOs",
       y = "Celltype") +
  optimized_theme_fig()+
  theme(
    
    legend.justification = "right"
  )

Fig2.2
#combine with example
example <- merged_logFC %>%
  filter(genotype == "Wdr82",
         celltype == "Eo.Ba") %>%
  ggplot(aes(x = logFC_ex.vivo, y = logFC_in.vivo)) +
  geom_hex(bins = 50) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") + 
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c") +
  #scale_fill_viridis_c() +
  #scale_fill_viridis_c() +  # Optional: nice continuous color scale
  
  labs(
    title = "logFC ex vivo vs in vivo",
    x = "logFC (Ex Vivo)",
    y = "logFC (In Vivo)",
    fill = "Gene Count"
  )+ optimized_theme_fig()

# Save the version with the legend
Fig2.2 <- example + Fig2.2 + plot_layout(widths = c(1, 5))
#paper--------------
ggsave(
  filename = basedir("Fig2.2_with_legend_3rep_06.04.25.pdf"),
  plot = Fig2.2,
  width = 18,
  height = 5,
  units = "cm"
)

# Save the version without the legend
ggsave(
  filename = basedir("Fig2.2_no_legend__3rep.pdf"),
  plot = Fig2.2 ,
  width = 16,
  height = 4.5,
  units = "cm"
)

koi_t <- as.data.frame(koi,"koi")
write.table(koi_t,basedir("genotypes.tsv"), row.names = F)
################################################################################



#Supplementary-------------------2.1
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
##############
#Fig2.3 current

InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
meta <- fread(InDir2("meta_cleaned.tsv")) # Read data
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   

meta <- meta[, -1, drop = FALSE] 
meta$sample1 <- rownames(meta)

limmaRes <- read_rds(InDir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
# Check if there are at least 3 distinct samples per tissue for each genotype and celltype
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop")%>%
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

koi <- Reduce(intersect, list(selected_KOs))

#
dataVoom_Eo.Ba <- read_rds(InDir3("Eo.Ba_dataVoom.rds"))
dataVoom_Mono <- read_rds(InDir3("Mono_dataVoom.rds"))
dataVoom_MkP <- read_rds(InDir3("MkP_dataVoom.rds"))
dataVoom_GMP <- read_rds(InDir3("GMP_dataVoom.rds"))
dataVoom_HSC <- read_rds(InDir3("HSC_dataVoom.rds"))
dataVoom_MEP.early <- read_rds(InDir3("MEP.early_dataVoom.rds"))
dataVoom_Gran. <- read_rds(InDir3("Gran._dataVoom.rds"))
dataVoom_Gran.P <- read_rds(InDir3("Gran.P_dataVoom.rds"))
KO <- koi[1]
ct <- unique(meta$celltype)[1]


dat.list <-list()
for (KO in koi){
  list_of_genes <- c("Oas2","Gbp3","Tnfaip6",
                     "Oas3","Irf7","Gvin1",
                     "Msmo1","Mthfd2","Idi1","Ccnd1")
  for (ct in unique(meta$celltype)) {
    # Get the dataVoom object corresponding to the current cell type
    dataVoom_ct <- get(paste0("dataVoom_", ct))
    #CAPITALIZE
    
    # Check if goi exists in the row names of dataVoom_ct$E
    if (any(rownames(dataVoom_ct$E) %in% unique(list_of_genes))){
      for (goi in unique(list_of_genes)) {
        # Proceed only if goi exists in the row names of dataVoom_ct$E
        if (goi %in% rownames(dataVoom_ct$E)) {
          # Subset the metadata and E values for the current gene and cell type
          gene_data <- meta[names(dataVoom_ct$E[goi,]),] %>%
            mutate(E = dataVoom_ct$E[goi,]) %>%
            rownames_to_column("samples") %>%
            filter(genotype %in% c(KO, "NTC")) %>%
            mutate(scaled_E = scale(E)) %>%
            mutate(gene = goi)%>%
            mutate(celltype=ct)%>%
            mutate(comparison=KO)
          
          # Store the gene data in the list
          dat.list[[paste0(ct, "_", goi,KO)]] <- gene_data
        }
      }
    }
  }
}

goi_exp <- bind_rows(dat.list,.id = "celltype_gene_genotype")

# Define the gene and cell type
goi <- "Gbp3"
ct <- "Eo.Ba"

# Define the KOs of interest
kos <- c("Brd9", "Wdr82", "Rcor1")

# Manually define the effect labels for each KO
effect_labels <- c("Brd9" = "Opposite trend", 
                   "Wdr82" = "Exaggerated trend", 
                   "Rcor1" = "similar trend")


# Generate the plots using lapply() instead of a for loop and store plots by KO name
#fig-----
stat_tests_all <- list()  # Initialize a list to store results

# Loop over all KOs
for (KO in kos) {
  
  # Filter the gene expression data for this KO
  filtered_data <- goi_exp %>%
    filter(gene == goi, celltype == ct, comparison == KO)
  
  filtered_data$genotype <- factor(filtered_data$genotype,
                                   levels = c("NTC", setdiff(filtered_data$genotype, "NTC")))
  
  filtered_data$E <- as.numeric(filtered_data$E)
  
  if (nrow(filtered_data) > 0) {
    for (tissue_type in unique(filtered_data$tissue)) {
      filtered_tissue <- filtered_data %>% filter(tissue == tissue_type)
      
      if (nrow(filtered_tissue) > 1 && length(unique(filtered_tissue$genotype)) > 1) {
        
        E_values <- filtered_tissue$E
        
        normality_p <- shapiro.test(E_values)$p.value
        test_method <- ifelse(normality_p > 0.05, "t.test", "wilcox.test")
        
        # Safe comparison using tryCatch
        stat_test <- tryCatch({
          compare_means(
            E ~ genotype,
            data = filtered_tissue,
            method = test_method
          ) %>%
            mutate(tissue = tissue_type, KO = KO)
        }, error = function(e) {
          message(paste("compare_means failed for", KO, "in", tissue_type, ":", e$message))
          NULL
        })
        
        if (!is.null(stat_test)) {
          stat_tests_all[[paste(KO, tissue_type, sep = "_")]] <- stat_test
        }
      }
    }
  } else {
    message(paste("No data available for gene:", goi, "in celltype:", ct, "for KO:", KO))
  }
}

# Combine all the results into a single data frame
stat_tests_combined <- bind_rows(stat_tests_all)

# Add significance stars based on p-values
stat_tests_combined <- stat_tests_combined %>%
  mutate(significance = case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))


# Save the table to a CSV file
write.csv(stat_tests_combined, "statistical_tests_results.csv", row.names = FALSE)


plots <- lapply(kos, function(KO) {
  
  # Filter limma and expression data
  filtered_limma <- limmaRes %>%
    filter(ensg == goi, coef == KO, celltype == ct)
  
  effect_label <- effect_labels[KO]
  
  filtered_data <- goi_exp %>%
    filter(gene == goi, celltype == ct, comparison == KO)
  filtered_data$genotype <- factor(filtered_data$genotype,
                                   levels = c("NTC", setdiff(filtered_data$genotype, "NTC")))
  
  if (nrow(filtered_data) > 0) {
    
    # Get significance annotations for this KO
    stat_subset <- stat_tests_combined %>%
      filter(KO == !!KO) %>%
      select(tissue, significance)
    
    # Merge significance with a default y-position for label placement
    annotation_data <- filtered_data %>%
      group_by(tissue) %>%
      summarize(y_pos = max(E, na.rm = TRUE) * 1.05) %>%  # Place label above the max value
      left_join(stat_subset, by = "tissue")
    
    # Build the plot
    p <- ggplot(filtered_data, aes(x = genotype, y = E, fill = tissue)) + 
      geom_boxplot(
        outlier.shape = NA,
        position = position_dodge(width = 0.8),
        color = "black",
        size = 0.2
      ) +
      geom_jitter(
        position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
        alpha = 0.5
      ) +
      facet_grid(cols = vars(tissue), scales = "free") +
      scale_fill_manual(values = c("ex.vivo" = "#6F987B", "in.vivo" = "#764BAB"),
                        
                        name = "Culture model") +
      labs(title = NULL,
           y = "scaled Exp") +
      xlab(paste0(KO, " KO ", "(", effect_label, ")")) +
      theme(legend.position = "none") +
      optimized_theme_fig() +
      geom_text(data = annotation_data,
                aes(x = 1.5, y = y_pos, label = significance),  # x = 1.5 centers between 2 groups
                inherit.aes = FALSE,
                size = 4)
    
    return(p)
    
  } else {
    message(paste("No data available for gene:", goi, "in celltype:", ct, "for KO:", KO))
    return(NULL)
  }
})

names(plots) <- kos
# Assign legend positions to individual plots
brd9 <- plots[["Brd9"]] + theme(legend.position = "none")  # Legend at bottom for Brd9

wdr82 <- plots[["Wdr82"]] + theme(legend.position = "none")   # No legend for Wdr82
Rcor1 <- plots[["Rcor1"]] + theme(legend.position = "none")   # No legend for Rcor1

# Combine plots using patchwork
Fig.2.3 <- brd9 + wdr82 + Rcor1 + 
  plot_layout(ncol = 3,guides= "collect") + 
  theme(
    legend.position = "right"
  )  # Collects legends into one

# Print the combined plot to check before saving
Fig.2.3

# Save the version with the legend
ggsave(
  filename = basedir(paste0(goi,"Fig.2.3_with_legend_new",".pdf")),
  plot = Fig.2.3,
  width = 18,
  height = 6 ,
  units = "cm"
)


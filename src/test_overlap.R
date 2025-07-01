#figure 1
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggtext)

##################################################################################
base<-"Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/"
inp<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
InDir1 <-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
out <- "Figure1"
outdir <- dirout("Figure1")
source("src/Ag_Optimized_theme_fig.R")
limmaRes_int <- read_rds(InDir1("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))%>%
  mutate(genes = ensg)
limmaRes_NTC <- read_rds(inp("limma_perCTex.vivovsin.vivo.rds"))
head(limmaRes_int)
head(limmaRes_NTC)
# Assuming the `coef` column identifies each KO type in `limmaRes_int`.
# Join the datasets by matching rows as necessary (e.g., gene or ID column).
ko_flags <- meta %>%
  group_by(genotype, tissue, celltype) %>%
  summarize(num_sample = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>%
  mutate(valid_ko = if_else(in.vivo >= 3 & ex.vivo >= 3, TRUE, FALSE)) %>%
  dplyr::select(genotype, celltype, valid_ko)%>%
  mutate(coef = genoty)

merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)
head(merged_data)
selected_KOs <- meta %>%
  group_by(genotype, tissue, celltype) %>%                  # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample1), .groups = 'drop') %>% # Count distinct samples for each group
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>% # Spread tissue to separate columns (in.vivo and ex.vivo)
  group_by(genotype) %>%                                    # Regroup by genotype
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) #%>%              # Keep genotypes that have at least one celltype with 3+ samples in both tissues
  #pull(genotype) %>% unique()
# Function to calculate correlation for each KO
correlation_results <- merged_data %>%
  group_by(coef,celltype) %>%  # Group by each KO type
  summarize(
    #cor_raw = cor(logFC_NTC, logFC_KO, method = "pearson"),
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson")
  )

# Ensure `coef` column is named appropriately to represent KO identifiers

# Step 1: Reshape the data, handling duplicates by averaging them
abs_logFC_data <- merged_data %>%
  dplyr::select(genes, coef, logFC_KO, logFC_NTC, celltype) %>%
  pivot_wider(names_from = coef, values_from = logFC_KO, values_fn = mean) %>%  # Handle duplicates by averaging
  mutate(logFC_NTC = abs(logFC_NTC)) %>%
  mutate(across(where(is.numeric), abs))  # Apply abs() only to numeric columns

# Check the output to ensure it's correctly reshaped and numeric
str(abs_logFC_data)

# Step 2: Calculate the correlation matrix for the KO and NTC columns, grouped by celltype
cor_matrix <- abs_logFC_data %>%
  group_by(celltype) %>%
  dplyr::select(-genes) %>%
  cor(method = "pearson", use = "complete.obs")  # Calculate correlation within each celltype

# Step 3: Convert correlation matrix into a format suitable for ggplot2
cor_data <- as.data.frame(as.table(cor_matrix))
names(cor_data) <- c("KO1", "KO2", "Correlation")

# Step 4: Plot the heatmap
ggplot(cor_data, aes(x = KO1, y = KO2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    title = "Pairwise Correlation Heatmap of Absolute logFC for KOs and NTC",
    x = "KO",
    y = "KO"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ celltype, scales = "free")  # Facet by celltype for separate heatmaps

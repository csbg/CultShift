#limma_NTC
#Comparison:ex.vivo-in.vivo

###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)
library(ComplexHeatmap)
##################################################################################
inDir<-dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
out<-"Ag_ScRNA_09_pseudobulk_per_celltype_correlation_of_AveGene_exp/"
outdir<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_correlation_of_AveGene_exp")
source("src/Ag_Optimized_theme.R")
##################################################################################
#load data
########################
#metadata
meta<-read.delim(inDir("metadata_guide.tsv"),row.names=1)
counts <- read.delim(inDir("combined_in_ex_counts_guide.tsv"), row.names = 1)
#meta data
unique(meta$celltype)
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell","MEP","Imm. B-cell")
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                    "Hba-a2","Pim1","Fabp5","Fdps","Cd9")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]
#only select the genotypes present in both tissue conditions
meta<-meta[meta$genotype %in% meta[meta$tissue=="ex.vivo",]$genotype,]

# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)-]", ".", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta<-meta%>%filter(!grepl("NA",rownames(meta)))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))
#selecting only NTC
NTC_counts<-counts[,grep("NTC",colnames(counts),value = T)]
NTC_meta<-meta[grep("NTC",rownames(meta),value = T),]
basedir<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
NTC_meta%>%write_rds(basedir("NTC_meta.rds"))
#counts<-counts[,rownames(NTC_meta)]
counts<-counts[!(rownames(counts) %in% genes_to_exclude),rownames(NTC_meta)]
stopifnot(all(colnames(counts)==rownames(NTC_meta)))
meta<-NULL

##################################
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")threshold <- 30



threshold <- 30
drop <- which(apply(edgeR::cpm(d0), 1, max) < threshold)
d <- d0[-drop,] 

cpm_values_before <- cpm(d0)

# Calculate row sums of CPM values before filtering
row_sums_cpm_before <- rowSums(cpm_values_before)

# Convert row sums to a data frame for ggplot2
row_sums_cpm_before_df <- data.frame(RowSumsCPM = row_sums_cpm_before)

# Plot the density of row sums of CPM values before filtering
p_before <- ggplot(row_sums_cpm_before_df, aes(x = log10(RowSumsCPM + 1))) +
  geom_density(fill = "blue", alpha = 0.4) +
  labs(
    title = "Density Plot of Row Sums of CPM Values (Log Scale) Before Filtering",
    x = "Row Sums of CPM (Log Scale)",
    y = "Density"
  ) +
  theme_minimal()

# Save the plot before filtering with default filename
ggsave(basedir(paste0("counts_density_before_filtering.pdf")), plot = p_before)

# Report the number of genes remaining after filtering
cat("Number of genes remaining after threshold filtering:",nrow(d$counts), "\n")

# Calculate CPM for the filtered DGEList
cpm_values_after <- cpm(d)

# Calculate row sums of CPM values after filtering
row_sums_cpm_after <- rowSums(cpm_values_after)

# Convert row sums to a data frame for ggplot2
row_sums_cpm_after_df <- data.frame(RowSumsCPM = row_sums_cpm_after)

# Plot the density of row sums of CPM values after filtering
p_after <- ggplot(row_sums_cpm_after_df, aes(x = log10(RowSumsCPM + 1))) +
  geom_density(fill = "blue", alpha = 0.4) +
  labs(
    title = "Density Plot of Row Sums of CPM Values (Log Scale) After Filtering",
    x = "Row Sums of CPM (Log Scale)",
    y = "Density"
  ) +
  theme_minimal()

# Save the plot after filtering with threshold in filename
ggsave(filename = basedir(paste0("counts_density_after_filtering_threshold_", threshold, ".pdf")), plot = p_after)
###############################################################################
#setting the model
###############################################################################
group <- interaction(NTC_meta$celltype, NTC_meta$tissue)

mm <- model.matrix(~0 + group)
rownames(mm) <- rownames(NTC_meta)

# Normalization
dataVoom <- voom(d, mm)

# Extract normalized expression values (log-transformed)
norm_data <- dataVoom$E  # This contains the normalized expression values

# Extract relevant metadata to create separate subsets for ex.vivo and in.vivo
NTC_meta$celltype <- factor(NTC_meta$celltype)
NTC_meta$tissue <- factor(NTC_meta$tissue)

# List of unique celltypes
celltypes <- unique(NTC_meta$celltype)

# Initialize a list to store correlations

correlations <- list()  # Initialize a list to store correlations for each cell type

for (celltype in unique(NTC_meta$celltypes)) {
  # Subset the data for this celltype
  celltype_data <- norm_data[, NTC_meta$celltype == celltype]
  celltype_meta <- NTC_meta[NTC_meta$celltype == celltype, ]
  
  # Subset by tissue type (ex.vivo vs in.vivo)
  ex_vivo_data <- celltype_data[, celltype_meta$tissue == "ex.vivo"]
  in_vivo_data <- celltype_data[, celltype_meta$tissue == "in.vivo"]
  
  # Calculate average expression across all samples within each condition
  ex_vivo_avg <- rowMeans(ex_vivo_data, na.rm = TRUE)
  in_vivo_avg <- rowMeans(in_vivo_data, na.rm = TRUE)
  
  # Ensure the dimensions are correct (number of genes must match)
  stopifnot(length(ex_vivo_avg) == length(in_vivo_avg))
  
  # Calculate correlation (Spearman or Pearson)
  cor_value <- cor(ex_vivo_avg, in_vivo_avg, method = "pearson")
  
  # Store the correlation value in the list
  correlations[[celltype]] <- cor_value
}

# Print or save the results
correlations_df <- data.frame(Celltype = names(correlations), Correlation = unlist(correlations))
print(correlations_df)
write.csv(correlations_df, outdir("correlations_per_celltype.csv"), row.names = FALSE)
# Generate a heatmap

ggplot(correlations_df, aes(x = Celltype, y = Correlation)) +
  geom_bar(stat = "identity") +
  
  labs(
    title = "Gene Expression Correlation (Ex Vivo vs In Vivo) by Celltype",
    x = "Celltype",
    y = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.position = "none"
  )
ggsave(outdir("correlation_percelltype_NTC.png"))

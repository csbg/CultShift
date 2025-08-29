###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)

##################################################################################
inDir<-dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
base<-"Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/"
basedir<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
##################################################################################
#load data
########################
#metadata
meta<-read.delim(inDir("metadata_guide.tsv"),row.names=1)
counts <- read.delim(inDir("combined_in_ex_counts_guide.tsv"), row.names = 1)

#meta data
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell","Gran.","MEP")
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
inDir1<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
NTC_meta%>%write_rds(inDir1("NTC_meta.rds"))
counts<-counts[!(rownames(counts) %in% genes_to_exclude),rownames(NTC_meta)]

stopifnot(all(colnames(counts)==rownames(NTC_meta)))
meta<-NULL

##################################
#'d0' is your DGEList object
QC<-dirout(paste0(base,"/QC"))
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

cat("Number of genes remaining after filterByExpr filtering:", nrow(d$counts), "\n")

# Assuming 'd0' is your DGEList object created with your counts table
# Step 1: Calculate CPM values
cpm_values <- cpm(d0)

# Define the specific CPM cutoff values and minimum number of samples
cpm_thresholds <- c(0,1, 5, 10,20,30, 50, 80, 100)
min_samples <- 2

# Initialize a vector to store the proportions
proportions <- numeric(length(cpm_thresholds))

# Loop through each CPM threshold
for (i in seq_along(cpm_thresholds)) {
  threshold <- cpm_thresholds[i]
  
  # Create a logical matrix where CPM > threshold
  logical_matrix <- cpm_values > threshold
  
  # Count how many samples have CPM > threshold for each gene
  gene_expression_count <- rowSums(logical_matrix)
  
  # Determine if each gene has CPM > threshold in at least 'min_samples' samples
  genes_meeting_criteria <- gene_expression_count >= min_samples
  
  # Calculate the proportion of genes that meet the criteria
  proportions[i] <- mean(genes_meeting_criteria)
}

# Print the results
plot_data <- data.frame(CPM_Threshold = cpm_thresholds, Proportion_Genes_Retained = proportions)

ggplot(plot_data, aes(x = CPM_Threshold, y = Proportion_Genes_Retained)) +
  geom_line() +
  geom_point() +
  labs(title = "Proportion of Genes Retained vs. CPM Threshold",
       x = "CPM Threshold",
       y = "Proportion of Genes Retained") +
  theme_minimal()
ggsave(QC("genes_retained_vs_threshold_filter.pdf"))
################################################################################
# Calculate total CPM for each gene across all samples
cpm_values <- cpm(d0)
total_cpm_per_gene <- rowSums(cpm_values)

# Filter out zero read genes
filtered_cpm_values <- cpm_values[!zero_read_genes, ]
cpm_thresholds <- c(0, 1, 5, 10, 20,30, 50, 80, 100)
min_samples <- 2

# Initialize vectors to store the number of genes and proportions
num_genes_retained <- numeric(length(cpm_thresholds))

# Loop through each CPM threshold
for (i in seq_along(cpm_thresholds)) {
  threshold <- cpm_thresholds[i]
  
  # Create a logical matrix where CPM > threshold
  logical_matrix <- filtered_cpm_values >= threshold
  
  # Count how many samples have CPM > threshold for each gene
  gene_expression_count <- rowSums(logical_matrix)
  
  # Determine if each gene has CPM > threshold in at least 'min_samples' samples
  genes_meeting_criteria <- gene_expression_count >= min_samples
  
  # Calculate the number of genes that meet the criteria
  num_genes_retained[i] <- sum(genes_meeting_criteria)
}

# Print the results
results <- data.frame(CPM_Threshold = cpm_thresholds, Number_of_Genes_Retained = num_genes_retained)

# Plotting the number of genes retained at each CPM threshold
ggplot(results, aes(x = factor(CPM_Threshold), y = Number_of_Genes_Retained)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Number of Genes Retained vs. CPM Threshold (>=threshold)",
       x = "CPM Threshold",
       y = "Number of Genes Retained") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(QC("number_of_genes_vs_threshold.pdf"))

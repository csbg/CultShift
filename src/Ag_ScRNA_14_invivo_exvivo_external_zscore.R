##################################################
#***Ag_ScRNA_13_invivo_izzo_limma********#
##################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#04-02-25
#Aarathy


########################################################################
source("src/00_init.R")
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(tidyverse)
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")

# Input/Output directories-----------

InDir3 <- dirout("Ag_ScRNA_08_Pseudobulk_external_limma_guide/")
InDir5 <- dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
base <- "Ag_ScRNA_14_invivo_exvivo_external_zscore/"
basedir <- dirout("Ag_ScRNA_14_invivo_exvivo_external_zscore/")

########################################################################
#load data and clean metadata
########################################################################
#metadata from in-vivo ex-vivo
meta <- fread(InDir_NTC("meta_cleaned.tsv"))
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]

meta <- meta[grep("NTC",rownames(meta),value = T),]
#selecting only NTC

meta <- meta[, -1, drop = FALSE]
meta <- meta[,c("sample","celltype","tissue")]

#metadata fromizzo
meta_izzo <- fread(InDir5("izzo_metadata.tsv"))
meta_izzo <- as.data.frame(meta_izzo)               # Convert to dataframe (optional)
 
colnames(meta_izzo) <- c("sample","replicate", "celltype","tissue")
rownames(meta_izzo) <- meta_izzo$sample
meta_izzo <- meta_izzo[,c("sample","celltype","tissue")]

#meta from external
meta_ext <- fread(InDir3("metadata_guide_external_data.tsv"))
meta_ext <- as.data.frame(meta_ext)  
rownames(meta_ext) <- meta_ext$cell
head(meta_ext)
clean_names <- function(x) {
  x <- gsub("/$", "", x)        # Remove trailing slash
  x <- trimws(x)                # Remove leading/trailing spaces
  x <- gsub(" ", ".", x)        # Replace spaces with dots
  return(x)
}


rownames(meta_ext) <- clean_names(rownames(meta_ext))
clean_names <- function(x) {
  x <- gsub("/$", "", x)        # Remove trailing slash
  x <- trimws(x)                # Remove leading/trailing spaces
  x <- gsub(" ", ".", x)        # Replace spaces with dots
  x <- gsub("\\.+", ".", x)     # Replace multiple dots with single dot
  return(x)
}
rownames(meta_ext) <- clean_names(rownames(meta_ext))

meta_ext <- meta_ext[,c("sample","celltype","tissue")]
meta_ext$sample <- gsub("/$", "",meta_ext$sample)

meta <- rbind(meta,meta_izzo)
meta <- rbind(meta,meta_ext)
# Replace spaces, parentheses, and hyphens with periods
rownames(meta) <- gsub("[\\s\\(\\)-]", ".", rownames(meta))
# Replace only spaces with periods
rownames(meta) <- gsub(" ", ".", rownames(meta))


rownames(meta) <- gsub("[\\s\\(\\)\\-]", ".", rownames(meta))
# filtering steps
# celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell",
#                           "MEP","MEP.pert.","MEP.S", "MEP.G1","GMP.late",
#                           "GMP.early")
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                      "Hba-a2","Pim1","Fabp5","Fdps","Cd9")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]


# Standardize specific terms in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))
#rownames(meta) <- gsub("IMP1|IMP2", "GMP", rownames(meta))

# Replace terms across all columns
meta <- meta %>%
  rownames_to_column(var = "cell_id") %>%  # Convert rownames to a proper column
  mutate(across(everything(), ~ gsub("Eo/Ba", "Eo.Ba", .))) %>%
  #mutate(across(everything(), ~ gsub("IMP1|IMP2", "GMP", .))) %>%
  mutate(celltype = ifelse(celltype %in% c("Eo", "Ba"), "Eo.Ba", celltype)) %>%
  filter(!celltype %in% c("B-cell", "Ery", "T-cell-Cd3d+", "E/B", "Imm. B-cell")) %>% #"Neu"
  filter(!grepl("T\\.cell", cell_id)) %>%  # Now filtering works
  column_to_rownames(var = "cell_id") %>%
  filter(tissue != "Bet_et_al")

# Check unique cell types in ex vivo condition


##############

#counts
counts_david <- read.delim(InDir5("combined_in_ex.vivo_with_Mye_counts_guide.tsv"), row.names = 1)
counts_izzo <- read.table(InDir5("izzo_counts.tsv"), row.names = 1)
counts_ext <- read.table(InDir3("combined_Anna_Bet_counts_guide.tsv"), row.names = 1)
colnames(counts_ext)
# Clean column names
clean_colnames <- function(x) {
  gsub("(WT[0-9]).*", "\\1",
       gsub("[\\ \\(\\)-]", ".",
            gsub("Eo/Ba", "Eo.Ba", x)))
}

colnames(counts_izzo) <- clean_colnames(colnames(counts_izzo))
colnames(counts_ext)  <- clean_colnames(colnames(counts_ext))

# Subset and merge count data
unique(meta$tissue)
counts_izzo <- as.matrix(counts_izzo[, rownames(meta[meta$tissue == "izzo", ])])
NTC_counts <- counts_david[rownames(counts_izzo), rownames(meta[meta$tissue %in% c("in.vivo","ex.vivo"),])]
counts <- cbind(NTC_counts, counts_izzo)
counts <- cbind(counts,counts_ext[rownames(counts),])
grep("MEP",unique(meta$celltype), value = T)
# Define a cleanup function to apply to both
clean_sample_names <- function(x) {
  x %>%
    gsub("Eo/Ba", "Eo.Ba", .) %>%
    gsub("ba.al", "basal", .) %>%
    gsub("_\\.cRNA\\.eq", "_scRNAseq", .) %>%
    gsub("_", ".", .) %>%
    gsub("\\.+", ".", .) %>%             # collapse multiple dots
    gsub("\\.$", "", .) %>%              # remove trailing dot
    trimws()
}

# Clean both names
colnames(counts) <- clean_sample_names(colnames(counts))
rownames(meta)   <- clean_sample_names(rownames(meta))
rownames(meta) <- gsub("Sca1po\\.", "Sca1pos.", rownames(meta))
meta$celltype[meta$celltype %in% c("MEP.pert.", "MEP.S", "MEP (S)", "MEP.G1","MEP")] <- "MEP"

counts <- counts[!rownames(counts) %in% genes_to_exclude,rownames(meta)]
# Ensure column names match row names of meta
stopifnot(all(colnames(counts) == rownames(meta)))

# Remove any NA values
counts <- na.omit(counts)

# Create DGEList and normalize
d <- DGEList(counts)
d <- calcNormFactors(d)

# Filter out low-expression genes
cutoff <- 30
d <- d[apply(cpm(d), 1, max) >= cutoff, ]

# Define groups and design matrix
group <- interaction(meta$celltype, meta$tissue)
mm <- model.matrix(~0 + group, data = meta)
rownames(mm) <- rownames(meta)

# Apply voom transformation
dataVoom <- voom(d, mm)

# Save and inspect results
write_rds(dataVoom$E, basedir("dataVoom_in_vivo_ex.vivo_vs_external.rds"))
dataVoom <- NULL
datavoom <- read_rds(basedir("dataVoom_in_vivo_ex.vivo_vs_external.rds"))
# set rownames to column sample
meta$sample <- rownames(meta)
# modify dataVoom
longer_dataVoom <-  datavoom %>%
  as.data.frame() %>%
  rownames_to_column("genes") %>%
  as_tibble() %>%
  pivot_longer(
    cols = -genes,     # Keep 'genes' as the identifier column
    names_to = "sample",  # Create a new column for previous column names
    values_to = "Expression"  # Create a new column for values
  )%>%
  inner_join(meta, by ="sample")%>%
  mutate(tissue_celltype =paste0(tissue,"_",celltype))

# Step 2: Pivot to wide format to have genes as rows and samples as columns
wide_data_ex_in <- longer_dataVoom %>%
  #filter(tissue != "izzo") %>%
  dplyr::select(genes, sample, Expression) %>%
  pivot_wider(names_from = sample, values_from = Expression) %>%
  column_to_rownames("genes")

# Step 3: Calculate pairwise correlation between samples for each gene
# Now, we'll compute the correlation matrix across samples for each gene

correlation_matrix <- cor(wide_data_ex_in, use = "complete.obs")

diag(correlation_matrix) <- NA
# Convert correlation to distance (1 - correlation)
dist_matrix <- as.dist(1 - correlation_matrix)

# Perform hierarchical clustering
row_clust <- hclust(dist_matrix, method = "ward.D2")


pdf(file = basedir("corheatmap_ex_in_external_celltype.pdf")
)
Heatmap(correlation_matrix,cluster_rows = row_clust, cluster_columns = row_clust)
dev.off()



#Selection of genes

genes_fig_1 <- read_rds(InDir_genes("genes_fig1.rds"))
colnames(genes_fig_1) <- c("genes","pathway")
ISG_fig_1 <- genes_fig_1[genes_fig_1$pathway =="ISG core",]$genes

# Step 2: Calculate Z-score within Each Tissue
ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)


# Step 1: Calculate the z-score per gene, per celltype, across samples (tissue)
longer_dataVoom_zscore <- longer_dataVoom %>%
  group_by(celltype, genes) %>%
  mutate(
    mean_expr = mean(Expression, na.rm = TRUE),         # Calculate mean for each gene per celltype
    sd_expr = sd(Expression, na.rm = TRUE),             # Calculate standard deviation for each gene per celltype
    zscore = (Expression - mean_expr) / sd_expr         # Compute the z-score
  ) %>%
  ungroup()%>%  # Remove grouping after calculation
  filter(genes %in% ISG_fig_1)%>%
  filter(celltype != "B.cell") %>%
  filter(!(genes %in% c("Stat1","Ube2l6"))) %>%
  dplyr::select(tissue, celltype,tissue_celltype, genes, sample, Expression, zscore)  # Keep relevant columns
longer_dataVoom_zscore %>% write_rds(basedir("zscore_plot_external.rds"))
# View the resulting dataframe with z-scores

unique(longer_dataVoom_zscore$tissue)
# Create the plot
ggplot(longer_dataVoom_zscore, aes(x = sample, y = genes , fill = pmin(2,zscore))) +
  geom_tile(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) +  # Scatter plot with jitter for better visualization
  facet_grid(cols = vars(tissue_celltype), scales = "free", space = "free") +  # Separate by tissue, free_x ensures that each tissue has its own x-axis range
  labs(title = "Z-score of Gene Expression Across Tissues", 
       x = "Genes", 
       y = "Z-score") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  optimized_theme_fig()+
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(angle = 90,hjust = 0))+
  theme(panel.spacing = unit(0.2, "lines")) # Rotate x-axis labels for readability
         # Remove minor gridlines

ggsave(basedir(paste0("Z-score_Distribution_", "_per_Tissue_line.pdf")), w=18,
       h=8, units = "cm")
####################
generate_zscore_plots <- function(gene_set) {
  # Filter zscore_tissue for the given gene set
  zscore_tissue_filtered <- gene_mean %>%
    group_by(tissue_celltype) %>%
    mutate(
      mean_tissue = mean(mean_expr, na.rm = TRUE),  # Mean expression for the group
      sd_tissue = sd(mean_expr, na.rm = TRUE),      # Standard deviation for the group
      zscore = (mean_expr - mean_tissue) / sd_tissue # Z-score for each gene
    ) %>%
    ungroup() %>%  # Remove grouping after calculation
    filter(genes %in% gene_set) %>%
    filter(celltype != "B.cell")
  
  # Generate plot 1: Z-score distribution with median line
  zscore_tissue_filtered %>%
    group_by(tissue, celltype) %>%
    mutate(median_zscore = median(zscore, na.rm = TRUE)) %>%
    ungroup() %>%
    ggplot(aes(x = tissue, y = zscore, group = genes)) +
    geom_line(color = "grey", alpha = 0.7) +  # Grey lines connecting the same gene across tissues
    geom_point(aes(color = tissue), size = 2) +  # Points colored by tissue
    # Add colored median line
    geom_line(aes(y = median_zscore, group = celltype), color = "black", size = 1) + 
    theme_minimal() +
    facet_grid(rows = vars(celltype)) +
    labs(title = paste("Z-score Distribution of Genes in", "per Tissue"),
         x = "Tissue",
         y = "Z-score") +
    optimized_theme_fig() +
    ggsave(basedir(paste0("Z-score_Distribution_", "_per_Tissue_line.pdf")))
  
  # Generate plot 2: Boxplot for Z-score distribution across genes
  zscore_tissue_filtered %>%
    ggplot(aes(x = tissue, y = zscore, color = tissue)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # Boxplot for each gene's zscore per tissue
    facet_wrap(~ genes, scales = "free_y") +  # Facet by gene, separate scales for each gene
    theme_minimal() +
    labs(title = paste("Z-score Distribution of Genes", "Across Tissues"),
         x = "Tissue",
         y = "Z-score") +
    optimized_theme_fig() +
    ggsave(basedir(paste0("per_gene_Z-score_Distribution_", "_Across_Tissues.pdf")))
  #generate plot 3
  # --- Plot 2: Heatmap ---
  zscore_tissue_filtered %>%
    ggplot(aes(x = tissue, y = genes, fill = zscore)) +
    geom_tile() +  # Heatmap for z-scores
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +  # Custom color scale
    facet_grid(cols = vars(celltype)) +  # Facet by celltype
    theme_minimal() +
    labs(title = paste("Gene Z-scores "),
         x = "Tissue",
         y = "Gene",
         fill = "Z-score") +
    optimized_theme_fig() +
    ggsave(basedir(paste0("Z-score_Distribution_", "_per_Tissue_heatmap.pdf")))
  #Boxplot
  zscore_tissue_filtered %>%
    ggplot(aes(x = tissue, y = zscore, color = tissue)) +
    geom_boxplot() +  # Heatmap for z-scores
    scale_color_manual(values=c("#c09fab","#87b1d6ff","#2b0000ff"))+
    facet_grid(cols = vars(celltype)) +  # Facet by celltype
    theme_minimal() +
    labs(title = paste("Gene Z-scores "),
         x = "Tissue",
         y = "Gene",
         fill = "Z-score") +
    optimized_theme_fig() +
    ggsave(basedir(paste0("Z-score_Distribution_", "_per_Tissue_boxplot.pdf")))
}

# Example: Calling the function with ISG_fig_1 gene set
generate_zscore_plots(ISG_fig_1)

#################################

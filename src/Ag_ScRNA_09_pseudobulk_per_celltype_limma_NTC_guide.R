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

##################################################################################
inDir<-dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")
base<-"Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/"
basedir <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
source("src/Ag_Optimized_theme.R")
##################################################################################
#load data
########################
#metadata
meta <- read.delim(inDir("metadata_guide.tsv"))
rownames(meta) <- meta$rowname

meta <- meta %>%
  mutate(
    # Check for discrepancies based on rowname and correct celltype
    celltype = case_when(
      grepl("GMP \\(early\\)", rowname) & celltype != "GMP.early" ~ "GMP.early", 
      grepl("GMP \\(late\\)", rowname) & celltype != "GMP.late" ~ "GMP.late",
      grepl("Gran\\. P", rowname) & celltype != "Gran.P" ~ "Gran.P",
      grepl("MEP \\(G1\\)" , rowname) & celltype != "MEP.G1"  ~ "MEP.G1" ,
      grepl("MEP \\(pert\\.\\)" , rowname) & celltype != "MEP.pert."  ~ "MEP.pert." ,
      grepl("MEP \\(S\\)"  , rowname) & celltype != "MEP.S"   ~ "MEP.S" ,
      grepl("MEP \\(early\\)"  , rowname) & celltype != "MEP.early" ~ "MEP.early" ,
      grepl("Imm. B-cell"  , rowname) & celltype == "Imm\\. B-cell"  ~ "Imm.B.cell", 
      TRUE ~ celltype
    )
  )


#meta data

celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell",
                          "Imm.B.cell","MEP","MEP.pert.","MEP.S", "MEP.G1","GMP.late",
                          "GMP.early","Imm. B-cell")

celltypes_to_exclude %>% write_rds(basedir("celltypes.exclude.rds"))
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                    "Hba-a2","Pim1","Fabp5","Fdps","Cd9")

# samples_28d <- unique(grep("_28d_",meta$sample, value = T))
# s28 <- meta %>%
#   filter(!(celltype %in% celltypes_to_exclude))%>%
#   filter(sample %in% samples_28d)%>%
#   group_by(celltype,genotype) %>%
#   summarise(n = n(), .groups = "drop")
# 28 days not excluded!!
meta <- meta %>%
  filter(!(celltype %in% celltypes_to_exclude))#%>%
  #filter(!(sample %in% samples_28d))


# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)-]", ".", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta <- meta%>%filter(!grepl("NA",rownames(meta)))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))

unique(meta$sample)
#only select the genotypes present in both tissue conditions
genotypes <- unique(meta[meta$tissue=="ex.vivo",]$genotype)
meta <- meta %>% 
  filter(genotype %in% genotypes)
write.table(meta,basedir("meta_cleaned.tsv"))
NTC_meta <- meta[grep("NTC",rownames(meta),value = T),]
NTC_meta %>% write_rds(basedir("NTC_meta.rds"))
#selecting only NTC
write.table(meta,basedir("meta_cleaned.tsv"))
counts <- read.delim(inDir("combined_in_ex_counts_guide.tsv"), row.names = 1)
NTC_counts <- counts[,grep("NTC",colnames(counts),value = T)]

#counts<-counts[,rownames(NTC_meta)]

counts <- counts[!(rownames(counts) %in% genes_to_exclude),rownames(NTC_meta)]
stopifnot(all(colnames(counts)==rownames(NTC_meta)))


##################################
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0,method = "TMM")
threshold <- 30
drop <- which(apply(cpm(d0), 1, max) < threshold)
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
dataVoom %>% write_rds(basedir("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaFit <- lmFit(dataVoom, mm)
limmaFit <- eBayes(limmaFit)

# Generate comparisons based on tissue types
tissue_type1<-"ex.vivo"
tissue_type2<-"in.vivo"
comparisons <- lapply(unique(NTC_meta$celltype), function(celltype) {
  comp1 <- paste0("group", celltype, ".",tissue_type1)
  comp2 <- paste0("group", celltype, ".", tissue_type2)
  list(comp1, comp2, celltype)
})

# Initialize list to store top table results
top_table <- list()

# Perform differential expression analysis for each comparison
for (comp_list in comparisons) {
  group1 <- comp_list[[1]]
  group2 <- comp_list[[2]]
  id <- comp_list[[3]]
  
  # Create contrast matrix
  contrast <- paste0(group1, "-", group2)
  
  # Debugging: Print contrast matrix
  print(contrast)
  
  # Ensure that contrast inputs are accepted
  if (!(all(c(group1, group2) %in% colnames(limmaFit)))) {
    stop("Invalid column names for comparison:", group1, " and ", group2)
  }
  
  # Fit model and perform differential expression analysis
  tmp <- contrasts.fit(limmaFit, makeContrasts(contrast, levels = colnames(coef(limmaFit))))
  tmp <- eBayes(tmp)
  top_table[[id]] <- topTable(tmp, sort.by = "P", n = Inf)
}

# Bind results and add NTC_metadata
ex_in_NTC_per_ct <- bind_rows(top_table, .id = 'celltype') %>%
  mutate(genes = gsub("\\...\\d+", "", rownames(.)))

# Create a column 'group' to group the genes as up, down, or ns
ex_in_NTC_per_ct$group <- ifelse(ex_in_NTC_per_ct$logFC >= 1 & 
                                ex_in_NTC_per_ct$adj.P.Val <= 0.05, "up", 
                              ifelse(ex_in_NTC_per_ct$logFC <= -1 & 
                                       ex_in_NTC_per_ct$adj.P.Val <= 0.05, "down", "n.s"))

ex_in_NTC_per_ct %>% write_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
################################################################################
# Plotting and saving-----------------------------------------------------------
################################################################################
# Generate comparisons based on tissue types
#figure 1

ex_in_NTC_per_ct <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))

top_genes <- ex_in_NTC_per_ct[ex_in_NTC_per_ct$genes %in% c("Idi1","Msmo1","Iigp1","Oas2"),] %>%
  unique()
ggplot() +
  # Hexbin plot for the "others" group
  stat_bin_hex(data = filter(ex_in_NTC_per_ct, group == "n.s"), 
               aes(x = logFC, y = -log10(adj.P.Val), fill = ..count..), 
               bins = 20, color = NA, alpha = 0.7) +
  scale_fill_gradient(low = "lightgrey", high = "black",
                      limits = c(1, 5000), name = "Gene Count") +
  
  # Overlay points for the main groups
  geom_point(data = filter(ex_in_NTC_per_ct, group != "n.s"), 
             aes(x = logFC, y = -log10(adj.P.Val), color = group), 
             alpha = 0.9, size = 2.5) +

  # Add text labels for top genes with ggrepel
  geom_text_repel(
    data = top_genes,
    aes(x = logFC, y = -log10(adj.P.Val), label = genes),
    size = 3, 
    color = "black",
    max.overlaps = Inf,  # Ensure no labels are omitted
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'black',  # Color for the line pointing to the gene
    segment.size = 0.5,       # Thickness of the line
    force = 2,                # Increase the repelling force to reduce overlap
    force_pull = 0.5,         # Pull force for the label towards the point
    min.segment.length = 0,   # Ensure segments are always drawn, even for close points
    arrow = arrow(length = unit(0.02, "npc"), type = "closed", angle = 15)  # Add arrows/lines pointing to the gene
  )+

  # Manually setting colors for groups
  scale_color_manual(values = c(
    "up" = "#D0154E", 
    "down" = "#4C889C"
  ), name = "Group") +
  
  labs(title = "Ex-vivo vs in-vivo differentially expressed genes",
       x = "logFC",
       y = "-log10(adj.P)") +
  facet_wrap(vars(celltype), scales = "free") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

ggsave(basedir(paste0("NTC_D.E_genes_percelltype_pvalue_.pdf")))


# Assuming ex_in_NTC_per_ct is already loaded and 'group' column is created
# Assuming ex_in_NTC_per_ct is already loaded and 'group' column is created
# Filter out 'n.s' group
filtered_data <- ex_in_NTC_per_ct %>%
  filter(group != "n.s")%>%
  filter(abs(logFC)>1)

# Count the number of up and down genes for each cell type
gene_counts <- filtered_data %>%
  group_by(celltype, group) %>%
  summarise(count = n()) %>%
  ungroup()

# Log10-transform the counts, adding 1 to avoid log(0)
gene_counts <- gene_counts %>%
  mutate(log10_count = log10(count + 1))

# Plotting
ggplot(gene_counts, aes(x = celltype, y = log10_count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("down" =  "#4C889C", "up" = "#D0154E")) +
  labs(title = "Number of differentially expressed genes per celltype",
       x = "Cell Type",
       y = "Log10(Number of Genes + 1)",
       fill = "Regulation") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        title = element_text(size=18))+
  theme(axis.text.x = element_blank()) +
  facet_wrap(~ celltype, scales = "free_x")+
  theme(strip.text = element_text(size = 18))
ggsave(basedir(paste0("NTC_D.E_genes_percelltype_padj_logFC.pdf")))
##########################################
# Plot logFC distribution for each coefficient
ggplot(ex_in_NTC_per_ct, aes(x = logFC, fill = celltype)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ celltype, scales = "free") +
  labs(title = "Log-Fold Change (logFC) Distribution for Each Celltype",
       x = "logFC",
       y = "Density") +optimized_theme_fig()
  
  
ggsave(basedir("logFC_distribution_per_celltype_TMM_threshold_30.pdf"))




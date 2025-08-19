
InDir_NTC <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide_Mye/")
InDir_int <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide_Mye/")
InDir_cor <-  dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation_Mye/")
InDir <- dirout("Figure2_Mye")      
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))


#--- Load and Prepare Metadata ---#
meta <- fread(InDir_NTC("meta_cleaned.tsv"))         # Read metadata
meta <- as.data.frame(meta)                       # Optional: convert to dataframe
rownames(meta) <- meta[[1]]                       # Set row names
meta <- meta[, -1, drop = FALSE]                  # Remove first column (now row names)
colnames(meta) <- gsub("rowname", "sample1", colnames(meta))  # Standardize sample name column


#--- Compute valid_ko flags per genotype, celltype ---#
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop") %>%
  mutate(coef = genotype)


#--- Compute replicate counts and valid_ko per KO ---#
replicates_per_ko <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(
    valid_ko = any(valid_ko),
    total_in_vivo = sum(in.vivo, na.rm = TRUE),
    total_ex_vivo = sum(ex.vivo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(coef = genotype)

replicates_per_ko %>% write_rds(InDir_int("replicates_per_ko.rds"))


#--- Select KOs valid in any cell type (≥3 samples per tissue) ---#
selected_KOs <- meta %>%
  group_by(genotype, tissue, celltype) %>%
  summarize(num_sample = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>%
  group_by(genotype) %>%
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) %>%
  pull(genotype) %>%
  unique()


#--- Summarize DEGs per KO/celltype ---#
adj_p_cutoff <- 0.05
logfc_cutoff <- 1

summary_df <- limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Upregulated, Downregulated),
    names_to = "Regulation",
    values_to = "Count"
  )


#--- Filter KOs with enough DEGs ---#
count_threshold = 10

coefficients <- summary_df %>%
  filter(Count != 0) %>%
  filter(Count >= count_threshold) %>%
  pull(coef) %>%
  unique()


#--- Filter KOs with poor correlation and enough DEGs ---#
# correlation_deg <- read_rds(InDir("correlation_deg.rds"))
# 
# KO_list <- correlation_deg %>%
#   filter(num_degs >= 10 & correlation < 0.5) %>%
#   pull(genotype) %>%
#   unique()
# 

#--- Final intersection: only valid and interesting KOs ---#
koi <- Reduce(intersect, list(selected_KOs,  coefficients))

# Step 2: Summarize to find KOs with at least one valid cell type
valid_ko_summary <- ko_flags %>%
  group_by(genotype) %>%
  summarize(has_valid_celltype = any(valid_ko), .groups = 'drop')
###############

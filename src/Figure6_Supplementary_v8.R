source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
library(tidyverse)
library(enrichR)
library(purrr)
library("scales")
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library(readr)

#directories ------

basedir <- dirout("Figure6_supplementary_v8")
## InDir4 <- dirout("Figure2_Mye")

InDir6 <- dirout("Ag_ScRNA_12_fgsea_overlap")



#Sup.Fig.6A---correlation density----------------------------


ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)

#Data and function


#
merged_logFC <- read_rds(InDir_cor("in.vivo_ex.vivo_logFC.rds"))
correlation_matrix_data <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(InDir_cor("DEGs_per_tissue.rds")) %>%
  dplyr::select(-correlation)

# Clean correlation matrix
correlation_matrix_no_na <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)

# Hierarchical clustering for column order
hc_cols <- hclust(dist(t(correlation_matrix_no_na)), method = "ward.D2")
column_order <- colnames(correlation_matrix_data)[hc_cols$order]
write_rds(column_order, basedir("column_order.rds"))

# Reshape correlation data
correlation <- correlation_matrix_data %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(
    cols = colnames(correlation_matrix_data),
    names_to = "genotype",
    values_to = "correlation",
    values_drop_na = FALSE
  )
correlation$dataset <- str_wrap("Perturb-seq(Hematopoietic)", width = 11)
correlation_perturb <- correlation
# Prepare DEG data: select max num_degs per genotype/celltype
deg_plot_data <- deg_plot_data %>%
  na.omit() %>%
  group_by(genotype, celltype) %>%
  filter(num_degs == max(num_degs, na.rm = TRUE)) %>%
  dplyr::select(-condition) %>%
  unique()

# Merge DEG data with correlation data
correlation_deg <- inner_join(deg_plot_data, correlation,
                              by = c("celltype", "genotype"))

# Merge with KO flags to include valid KO status
correlation_deg_flagged <- correlation_deg %>%
  inner_join(ko_flags, by = c("genotype", "celltype")) %>%
  filter(valid_ko) %>%
  na.omit()

# Update genotype factor levels based on hierarchical clustering and mean correlation
ko_order <- correlation_deg_flagged %>%
  group_by(genotype) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_corr)) %>%
  pull(genotype)

correlation_deg_flagged <- correlation_deg_flagged %>%
  mutate(genotype = factor(genotype, levels = rev(ko_order)))


# splenic--------------------
out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"

res <- read.delim(paste0(out,"/DEG_ResultsT8.tsv"))
#,
# read.delim(paste0(out,"/DEG_ResultsM.tsv"))
#)
res$probe <- res$rn
res$rn<-NULL
res <- as.data.frame(res)
gmap <- as.data.frame(read_rds(file = file.path(out, "DEG_GMP.RDS")))
res <- merge(res,gmap[,c("probe","gene")],by="probe")

# Prepare ex.vivo data
ex_vivo_data <- res %>%
  filter(coef != "treatmentex_vivo") %>%
  filter(str_detect(coef, "^ex.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "ex.vivo\\.", "")) %>%
  mutate(genotype = gsub("ex.vivo_genotype", "", coef)) %>%
  dplyr::select("gene", "logFC", "cell_type", "genotype", "adj.P.Val")

# Prepare in.vivo data
in_vivo_data <- res %>%
  filter(!(coef %in% grep("treatmentex_vivo", res$coef, value = T))) %>%
  filter(str_detect(coef, "^genotype")) %>%
  mutate(coef_clean = str_replace(coef, "genotype\\.", "")) %>%
  mutate(genotype = gsub("genotype", "", coef)) %>%
  dplyr::select("gene", "logFC", "cell_type", "genotype", "adj.P.Val")
interaction_data <- res %>%
  # keep only interactions
  filter(str_detect(coef, ":treatmentex_vivo$")) %>%
  
  mutate(
    # Replace "genotype" prefix with "Interaction_"
    coef_clean = str_replace(coef, "^genotype", "Interaction_"),
    
    # Remove the ":treatmentex_vivo" suffix
    coef_clean = str_replace(coef_clean, ":treatmentex_vivo$", ""),
    
    # Extract genotype (everything after Interaction_)
    genotype = str_remove(coef_clean, "^Interaction_")
  ) %>%
  
  dplyr::select(gene, logFC, cell_type, genotype, adj.P.Val, coef_clean)



# Combine datasets
merged_data_spleen <- ex_vivo_data %>%
  inner_join(in_vivo_data, by = c("gene", "cell_type", "genotype"), 
             suffix = c("_ex.vivo", "_in.vivo"))


correlation_results <- merged_data_spleen %>%
  group_by(cell_type, genotype) %>%
  summarise(correlation = cor(logFC_ex.vivo, logFC_in.vivo, use = "complete.obs")) %>%
  ungroup()
correlation_results$dataset <- "Spleen"
correlation_spleen <- correlation_results
#deg
deg_ex_vivo <- ex_vivo_data %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype, cell_type) %>%
  summarise(degs_ex_vivo = n())

deg_in_vivo <- in_vivo_data %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype, cell_type) %>%
  summarise(degs_in_vivo = n())

# Merge the DEG counts
deg_counts <- full_join(deg_ex_vivo, deg_in_vivo, by = c("genotype", "cell_type")) %>% na.omit()
deg_plot_data <- deg_counts %>%
  pivot_longer(cols = starts_with("degs"), names_to = "condition", values_to = "num_degs") %>%
  mutate(condition = ifelse(condition == "degs_ex_vivo", "Ex Vivo", "In Vivo")) %>%
  left_join(correlation_results, by = c("genotype", "cell_type")) %>%
  group_by(cell_type, genotype) %>%
  slice_max(order_by = num_degs, n = 1, with_ties = FALSE) %>%  # Retain only the row with max num_degs
  ungroup()

# Compute mean correlation per genotype across cell types
ko_order <- deg_plot_data %>%
  group_by(genotype) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_corr)) %>%  # descending order
  pull(genotype)

# Update the factor levels for plotting
deg_plot_data <- deg_plot_data %>%
  mutate(genotype = factor(genotype, levels = ko_order))

deg_plot_data %>%
  filter(genotype %in% c("IRF9KO","STAT1KO","STAT2KO","TYK2CMV",
                         "TYK2KE","STAT3VAV", "STAT1BB"))%>%
  mutate(genotype = ifelse(grepl("IRF9KO",genotype),"Irf9KO",ifelse(
    grepl("STAT1KO",genotype),"Stat1KO",ifelse(
      grepl("STAT2KO",genotype),"StatKO",ifelse(
        grepl("TYK2CMV",genotype),"Tyk2CMV",ifelse(
          grepl("TYK2KE",genotype),"Tyk2KE",ifelse(
            grepl("STAT3VAV",genotype),"Stat3Vav","Stat1BB"
          )))))))
deg_plot_data %>%write_rds(basedir("JAK_Stat_in_ex_cor.rds"))

#Glioblastoma------------
InDir <- dirout("Glioblastoma_limmaRes.3way.mod")

limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))
#invivo vs ex.vivo

limmaRes_invivo <- limmaRes %>% filter(coef %in% grep("[A-Za-z0-9]*_in.vivo_noRT",limmaRes$coef, value = T))
invivo_glio <- limmaRes_invivo %>%
  filter(str_detect(coef, "_in.vivo_noRT")) %>%
  mutate(genotype = gsub("_in.vivo_noRT", "", coef)) %>%
  mutate(gene = ensg) %>%
  dplyr::select("gene", "logFC",  "genotype", "adj.P.Val")
limmaRes_exvivo <- limmaRes %>% filter(coef %in% grep("[A-Za-z0-9]*_ex.vivo_noRT",limmaRes$coef, value = T))
exvivo_glio <- limmaRes_exvivo %>%
  filter(str_detect(coef, "_ex.vivo_noRT")) %>%
  mutate(genotype = gsub("_ex.vivo_noRT", "", coef)) %>%
  mutate(gene = ensg) %>%
  dplyr::select("gene", "logFC",  "genotype", "adj.P.Val")

# Combine datasets
merged_data_glio <- exvivo_glio%>%
  inner_join(invivo_glio, by = c("gene", "genotype"), 
             suffix = c("_ex.vivo", "_in.vivo"))


correlation_results <- merged_data_glio %>%
  group_by(genotype) %>%
  summarise(correlation = cor(logFC_ex.vivo, logFC_in.vivo, use = "complete.obs")) %>%
  ungroup()
correlation_results$celltype <- "Glioblastoma"
correlation_results$dataset <-"Glioblastoma"
correlation_glio <- correlation_results
#deg
deg_ex_vivo <- exvivo_glio %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype) %>%
  summarise(degs_ex_vivo = n())

deg_in_vivo <- invivo_glio %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype) %>%
  summarise(degs_in_vivo = n())

# Merge the DEG counts
deg_counts <- full_join(deg_ex_vivo, deg_in_vivo, by = c("genotype")) %>% na.omit()
deg_plot_data <- deg_counts %>%
  pivot_longer(cols = starts_with("degs"), names_to = "condition", values_to = "num_degs") %>%
  mutate(condition = ifelse(condition == "degs_ex_vivo", "Ex Vivo", "In Vivo")) %>%
  left_join(correlation_results, by = c("genotype")) %>%
  group_by(genotype) %>%
  slice_max(order_by = num_degs, n = 1, with_ties = FALSE) %>%  # Retain only the row with max num_degs
  ungroup()

# Compute mean correlation per genotype across cell types
ko_order <- deg_plot_data %>%
  group_by(genotype) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_corr)) %>%  # descending order
  pull(genotype)

# Update the factor levels for plotting
deg_plot_data <- deg_plot_data %>%
  mutate(genotype = factor(genotype, levels = rev(ko_order)))%>%
  mutate(celltype = "glioblastoma")


#Combined correlation plot--------------------------
names(correlation_perturb)
names(correlation_spleen) <- c("celltype",
                               "genotype",
                               "correlation",
                               "dataset")
correlation_all <- rbind(correlation_perturb, correlation_spleen, correlation_glio)
unique(correlation_all$dataset)


condition_colors <- c("#4E92D0" ,
                      "#478E5C",
                      "#AD1F69"
)
ggplot(correlation_all, aes(x = correlation,  color = dataset)) +
  geom_density(alpha = 0.4, size = 0.3) +        # semi-transparent fill
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept  = 0, linetype ="dotted", color ="red")+
  #scale_fill_manual(values = condition_colors)+
  scale_color_manual(values = condition_colors)+
  labs(
    title = "Correlation between in vivo and ex vivo KO effects",
    x = "Correlation",
    y = "Density",
    fill = "Dataset",
    color = "Dataset"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "top"
  )+optimized_theme_fig()+theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank())
ggsave(basedir("Sup.Fig.6A.Correlation_density_in.vivo_ex.vivo.pdf"), w = 7, h = 3, units = "cm")




#bFig.Sup.6B-----
correlation_matrix_data <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(InDir_cor("DEGs_per_tissue.rds")) %>%
  dplyr::select(-correlation)


# Clean correlation matrix
correlation_matrix_no_na <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)

# Hierarchical clustering for column order
hc_cols <- hclust(dist(t(correlation_matrix_no_na)), method = "ward.D2")

# Reshape correlation data
correlation <- correlation_matrix_data %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(
    cols = colnames(correlation_matrix_data),
    names_to = "genotype",
    values_to = "correlation",
    values_drop_na = FALSE
  )
correlation$dataset <- str_wrap("Perturb-seq(Hematopoietic)", width = 11)
correlation_perturb <- correlation
# Prepare DEG data: select max num_degs per genotype/celltype
deg_plot_data <- deg_plot_data %>%
  na.omit() %>%
  group_by(genotype, celltype) %>%
  filter(num_degs == max(num_degs, na.rm = TRUE)) %>%
  dplyr::select(-condition) %>%
  unique()

# Merge DEG data with correlation data
correlation_deg <- inner_join(deg_plot_data, correlation,
                              by = c("celltype", "genotype"))

# Merge with KO flags to include valid KO status
correlation_deg_flagged <- correlation_deg %>%
  inner_join(ko_flags, by = c("genotype", "celltype")) %>%
  filter(valid_ko) %>%
  na.omit()

# Update genotype factor levels based on hierarchical clustering and mean correlation
ko_order <- correlation_deg_flagged %>%
  group_by(genotype) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_corr)) %>%
  pull(genotype)

correlation_deg_flagged <- correlation_deg_flagged %>%
  mutate(genotype = factor(genotype, levels = rev(ko_order)))
mean_corr <- correlation_deg_flagged %>%
  group_by(coef) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(coef = factor(coef, levels = rev(ko_order)))



# Load limma results
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds")) %>%
  mutate(coef = gsub("interaction", "", coef))

# Step 1: Filter for significant genes
limmaRes_significant <- limmaRes %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# summary_df described in KO classification
data <- summary_df %>%
  filter(coef %in% koi) %>%
  filter(Count >= 10) %>%
  inner_join(ko_flags, by = c("celltype", "coef")) %>%
  filter(valid_ko)

# Load enrichment results
enrichment_ntc_in.vivo <- read_rds(InDir6("enrichment_to_NTC_genes.rds"))
enrichment_ntc_in.vivo <- enrichment_ntc_in.vivo %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  mutate(genotype = coef)
enrichment_ntc_in.vivo %>%
  dplyr::select(genotype,celltype,sig_genes, ref_genes,overlap, p.value, padj,odds.ratio) %>%
  write_rds(basedir("hematopoietic_enrichment_overlap.rds"))
# Combine KO/celltype with enrichment results
combined <- data %>%
  dplyr::select(coef, celltype) %>%
  distinct() %>%
  left_join(enrichment_ntc_in.vivo, by = c("coef", "celltype"))

# Significance stars
combined <- combined %>%
  mutate(significance_en = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01 ~ "**",
    p.value < 0.05 ~ "*",
    TRUE ~ ""
  ))

combined <- combined %>%
  mutate(log2.odds.ratio = log2(odds.ratio))

# Special case
combined <- combined %>%
  mutate(log2.odds.ratio = case_when(
    coef == "Smc3" ~ 0,
    TRUE ~ log2.odds.ratio
  ))

# Filter: keep only if overlap > 5
combined <- combined %>%
  mutate(log2.or.filtered = ifelse(overlap > 5, pmin(log2.odds.ratio, 7), NA))

# Cap -log10(p) at 10
combined <- combined %>%
  mutate(minuslog10p = pmin(-log10(padj), 10))

# === NEW: Apply identical KO ordering to BOTH datasets ===
ko_levels <- rev(ko_order)

combined$coef  <- factor(combined$coef,  levels = ko_levels)
mean_corr$coef <- factor(mean_corr$coef, levels = ko_levels)

# optional: drop unused KOs in heatmap
mean_corr$coef <- droplevels(mean_corr$coef)


# --- Enrichment dot plot ---
p <- ggplot(
  combined,
  aes(
    x = coef,
    y = celltype,
    fill = log2.or.filtered,
    size = minuslog10p
  )
) +
  geom_point(shape = 21, colour = "black", stroke = 0.2) +
  scale_fill_gradient2(
    low  = "#CCE5FF",
    mid  = "white",
    high = "#FF9933",
    midpoint = 0,
    na.value = "grey80",
    name = "log2(OR)"
  ) +
  scale_size_continuous(range = c(0, 2), name = "-log10(Padj)") +
  labs(
    title = str_wrap("Overlap enrichment (NTC in vivo) per KO and cell type", width = 90),
    x = "KO",
    y = "Cell type"
  ) +
  optimized_theme_fig() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
ko_present <- combined %>%
  droplevels() %>%        # remove unused factor levels
  pull(coef) %>%          # extract factor
  levels()                # get only levels that occur in the data

ko_present

heatmap_corr <- ggplot(mean_corr %>%
                         filter(coef %in% ko_present)
                       , aes(x = coef, y = 1, fill = mean_corr)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    midpoint = 0,
    na.value = "grey80",
    name = "Mean correlation",
    limits = c(-1, 1)
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


final_plot <- heatmap_corr / p + plot_layout(heights = c(0.2, 1))
final_plot

ggsave(basedir("Overlap_enrichment_mean_cor.pdf"),
       plot = final_plot,
       width = 10, height = 6, units = "cm")

#enrichment plot_overlap--glioFig.Sup.9C--------------------------------------
#

InDir <- dirout("Fig.Glio")
data_glio <- read_rds(InDir("overlap_enr_corr_inex.rds"))
data_glio %>%
  mutate(genotype = coef)%>%
  dplyr::select(genotype,sig_genes, ref_genes,overlap, p.value, p.adj,odds.ratio) %>%
  write_rds(basedir("glio_enrichment_overlap.rds"))

# Create bubble plot
Fig_overlap_bubble <- ggplot(data_glio%>%
                               filter(overlap > 5), aes(
                                 x = KO,
                                 y = 1,
                                 size = pmin(-log10(p.adj),5),
                                 fill = log2.odds.ratio
                               )) +
  geom_point(shape = 21, colour = "black", stroke = 0.2) +
  scale_fill_gradient2(
    low  = "#CCE5FF",
    mid  = "white",
    high = "#FF9933",
    midpoint = 0,
    na.value = "grey80",
    name = "log2(OR)"
  ) +
  scale_size_continuous(range = c(0, 2), name = "-log10(Padj)") +
  labs(
    title = str_wrap("Overlap enrichment (NTC in vivo) per KO and cell type", width = 90),
    x = "KO",
    y = "Cell type"
  ) +
  optimized_theme_fig() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
Fig_overlap_bubble

# Heatmap with correlation as fill
heatmap <- ggplot(data_glio%>%
                    filter(overlap > 5), aes(x = KO, y = 1, fill = cor)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    midpoint = 0,
    na.value = "grey80",
    name = "Mean correlation",
    limits = c(-1, 1)
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )



final_plot <- heatmap / Fig_overlap_bubble +
  plot_layout(heights = c(1, 0.6))

ggsave(basedir("Overlap_enrichment_mean_cor_Glio.pdf"),
       plot = final_plot,
       width = 10, height = 4.5, units = "cm")




#enrichment_plot_overlap--jakstatFig.Sup.9D------------
out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"
res <- read.delim(paste0(out,"/DEG_ResultsT8.tsv"))

res$probe <- res$rn
res$rn<-NULL
res <- as.data.frame(res)
gmap <- as.data.frame(read_rds(file = file.path(out, "DEG_GMP.RDS")))
res <- merge(res,gmap[,c("probe","gene")],by="probe")
unique(res$coef)
# Prepare ex.vivo data
ex_vivo_data <- res %>%
  filter(coef != "treatmentex_vivo") %>%
  filter(str_detect(coef, "^ex.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "ex.vivo\\.", "")) %>%
  mutate(genotype = gsub("ex.vivo_genotype", "", coef)) %>%
  dplyr::select("gene", "logFC", "cell_type", "genotype", "adj.P.Val")

# Prepare in.vivo data
in_vivo_data <- res %>%
  filter(!(coef %in% grep("treatmentex_vivo", res$coef, value = T))) %>%
  filter(str_detect(coef, "^genotype")) %>%
  mutate(coef_clean = str_replace(coef, "genotype\\.", "")) %>%
  mutate(genotype = gsub("genotype", "", coef)) %>%
  dplyr::select("gene", "logFC", "cell_type", "genotype", "adj.P.Val")
InDir <- dirout("Ag_ScRNA_28_Correlations_across_datasets")
combined <- read_rds(InDir("JAK_Stat_in_ex_cor.rds"))

res$group <- ifelse(res$logFC >= 1 & 
                      res$adj.P.Val <= 0.05, "up", 
                    ifelse(res$logFC <= -1 & 
                             res$adj.P.Val <= 0.05, "down", "n.s"))
ref_gene_sets <- res %>%
  filter(coef == "treatmentex_vivo")%>%
  filter(group != "n.s") %>%
  summarise(ref_genes = list(unique(gene)), .groups = "drop")

input_gene_sets <- in_vivo_data  %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype) %>%
  summarise(sig_genes = list(unique(gene)), .groups = "drop")


# All genes for background
all_genes <- unique(limmaRes$ensg)

# Join input and reference by celltype
joined_sets <- input_gene_sets 
joined_sets$ref_genes <- ref_gene_sets$ref_genes


run_fisher <- function(sig, ref, bg) {
  # Basic input checks
  if (length(sig) == 0 || length(ref) == 0 || length(bg) == 0) {
    return(summary.frame(overlap = NA, p.value = NA, odds.ratio = NA))
  }
  
  # Ensure unique genes
  sig <- unique(sig)
  ref <- unique(ref)
  bg <- unique(bg)
  
  # Calculate values
  x <- length(intersect(sig, ref))
  K <- length(sig)
  M <- length(ref)
  N <- length(bg)
  
  # Check for any invalid matrix values
  if (any(c(x, K, M, N) < 0) || x > K || x > M || K > N || M > N) {
    return(summary.frame(overlap = NA, p.value = NA, odds.ratio = NA))
  }
  
  # Create contingency table
  mat <- matrix(c(
    x,
    K - x,
    M - x,
    N - K - M + x
  ), nrow = 2)
  
  # Avoid invalid table (e.g., negative cell)
  if (any(mat < 0) || any(!is.finite(mat))) {
    return(summary.frame(overlap = NA, p.value = NA, odds.ratio = NA))
  }
  
  # Run Fisher test
  ft <- fisher.test(mat, alternative = "greater")
  
  # Return results
  data.frame(
    overlap = x,
    p.value = ft$p.value,
    odds.ratio = unname(ft$estimate)
  )
}


# Apply to each row in joined_sets
enrichment_results <- joined_sets %>%
  mutate(result = pmap(list(sig_genes, ref_genes), ~ run_fisher(..1, ..2, all_genes))) %>%
  unnest(result)

# View the result
enrichment_ntc_in.vivo <- enrichment_results %>%
  mutate(padj = p.adjust(p.value, method = "BH"))
combined <- combined %>%
  distinct()%>%
  left_join(enrichment_ntc_in.vivo, by = c("genotype"))
combined <- combined %>%filter(genotype %in% c("IRF9KO","STAT1KO","STAT2KO","TYK2CMV",
                                               "TYK2KE","STAT3VAV", "STAT1BB"))%>%
  mutate(genotype = ifelse(grepl("IRF9KO",genotype),"Irf9KO",ifelse(
    grepl("STAT1KO",genotype),"Stat1KO",ifelse(
      grepl("STAT2KO",genotype),"StatKO",ifelse(
        grepl("TYK2CMV",genotype),"Tyk2CMV",ifelse(
          grepl("TYK2KE",genotype),"Tyk2KE",ifelse(
            grepl("STAT3VAV",genotype),"Stat3Vav","Stat1BB"
          )))))))
ko_order <- combined %>%
  group_by(genotype) %>%
  arrange(correlation) %>%
  pull(genotype)
combined$genotype <- factor(combined$genotype, levels = ko_order)
combined %>% 
  dplyr::select(genotype,cell_type,sig_genes, ref_genes,overlap, p.value, padj,odds.ratio) %>%
  write_rds(basedir("spleen_enrichment_overlap.rds"))
#
Fig_overlap_bubble <- ggplot(combined%>%
                               filter(overlap > 5), aes(
                                 x = genotype,
                                 y = 1,
                                 size = pmin(-log10(padj),5),
                                 fill = log2(odds.ratio)
                               )) +
  geom_point(shape = 21, colour = "black", stroke = 0.2) +
  scale_fill_gradient2(
    low  = "#CCE5FF",
    mid  = "white",
    high = "#FF9933",
    midpoint = 0,
    na.value = "grey80",
    name = "log2(OR)"
  ) +
  scale_size_continuous(range = c(0, 2), 
                        limits =c(0,5),
                        name = "-log10(Padj)") +
  labs(
    title = str_wrap("Overlap enrichment (NTC in vivo) per KO and cell type", width = 90),
    x = "KO",
    y = "Cell type"
  ) +
  optimized_theme_fig() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
Fig_overlap_bubble

# Heatmap with correlation as fill
heatmap <- ggplot(combined%>%
                    filter(overlap > 5), aes(x = genotype, y = 1, fill = correlation)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    midpoint = 0,
    na.value = "grey80",
    name = "Mean correlation",
    limits = c(-1, 1)
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )



final_plot <- heatmap / Fig_overlap_bubble +
  plot_layout(heights = c(1, 0.6))

ggsave(basedir("Overlap_enrichment_mean_cor_Jakstat.pdf"),
       plot = final_plot,
       width = 5.2, height = 4.5, units = "cm")


#####################################
#

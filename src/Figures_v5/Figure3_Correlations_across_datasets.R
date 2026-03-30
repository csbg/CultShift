#
#load libraries and functions--------------------------------------------------
#
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
source("src/Ag_ko_classification.R")

library("scales")
library(tidyverse)
library(enrichR)
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library(ggridges)
library(ggpubr)
#directories ------
selected_KOs
#
base <- "Figure3_Correlations_across_datasets"
basedir <- dirout("Figure3_Correlations_across_datasets")

ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)

#Data and function

#Fig3A------------
#fig ------------
# Read in the correlation matrix and DEG data
# Load required data

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


#combine with example
Fig3A_example <- merged_logFC %>%
  filter(genotype == "Wdr82",
         celltype == "Eo.Ba") %>%
  ggplot(aes(x = logFC_ex.vivo, y = logFC_in.vivo)) +
  geom_hex(bins = 50) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") + 
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c") +
  #scale_fill_viridis_c() +
  #scale_fill_viridis_c() +  # Optional: nice continuous color scale
  
  labs(
    #title = "logFC ex vivo vs in vivo",
    x = "logFC (Ex vivo)",
    y = "logFC (In vivo)",
    
    fill = "Gene Count"
  )+ optimized_theme_fig()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(
  filename = basedir("Fig3A_example.pdf"),
  plot = Fig3A_example,
  width = 4.5,
  height = 2.5,
  units = "cm"
)


# Save the version with the legend
#fig-----
# Create the plot
Fig3B <- correlation_deg_flagged %>%
  filter(genotype %in% koi) %>%
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
    high = "#D0154E",
    midpoint = 0,       # center the gradient at 0
    limits = c(-1, 1), 
    name = expression("Pearson's correlation")
  ) +
  scale_size_continuous(
    range = c(0,2.5),
    #limits = c(0,2.5),
    breaks = c(1,2,3),
    name = expression(atop("No. of genes", log[10](n)))
  )+
  labs(x = "KOs",
       y = "Cell type",
       title = 
         "In vivo vs ex vivo KO effect concordance is gene-dependent \n
         Correlation of KO effects between in vivo and ex vivo (per cell type)"
  )+
  optimized_theme_fig()+
  theme(
    legend.position = "bottom",
    legend.justification = "right"
  )

Fig3B
n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$celltype))
ggsave(
  filename = basedir("Fig3A_perturb.pdf"),
  plot = Fig3B,
  width = n_col * 0.28,    # PDF becomes wider
  height = n_row * 0.8,
  units = "cm"
)


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

#merged_interaction_NTC

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
# Plot with genotypes ordered by aggregated correlation
Fig3C <- deg_plot_data %>%
  ggplot() +
  geom_point(aes(
    x = genotype,
    y = cell_type,
    size = pmin(3, log10(num_degs)),
    fill = correlation
  ),
  shape = 21,
  color = "black",
  stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "#4C889C",   # blue for negative correlation
    mid = "white",      # zero correlation
    high = "#D0154E",  # red for positive correlation
    midpoint = 0,       # center the gradient at 0
    limits = c(-1, 1),  # fixed scale from -1 to 1
    name = "Pearson's correlation"
  ) +
  scale_size_continuous(
    range = c(0, 2.5),
    breaks = c(1, 2, 3),
    name = expression(atop("No. of genes", log[10](n)))
  ) +
  labs(
    x = "KOs",
    y = "Cell type",
    title = str_wrap("Discordent KO effects ex vivo vs in vivo for JAK-STAT component KOs", width = 40)
  ) +
  optimized_theme_fig() +
  theme(
    legend.position = "bottom",
    legend.justification = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Wrap plot with extra space on the right for the legend
p <- ggdraw(Fig3C) + 
  theme(plot.margin = margin(t = 2, r = 50, b = 10, l = 70))

n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$cell_type))
# Save larger PDF; plot itself remains the same size
ggsave(
  filename = basedir("Fig3C_splenic_cells.pdf"),
  plot = p,
  width =  n_col,   # bigger PDF width to include legend
  height = 4.5 * n_row,   # keep plot height small
  units = "cm"
)







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

# Plot with genotypes ordered by aggregated correlation
Fig3D <- deg_plot_data %>%
  
  ggplot() +
  geom_point(aes(
    x = genotype,
    y = celltype,
    size = pmin(3, log10(num_degs)),
    fill = correlation
  ),
  shape = 21,
  color = "black",
  stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "#4C889C",   # blue for negative correlation
    mid = "white",      # zero correlation
    high = "#D0154E",  # red for positive correlation
    midpoint = 0,       # center the gradient at 0
    limits = c(-1, 1),  # fixed scale from -1 to 1
    name = "Pearson's correlation"
  ) +
  scale_size_continuous(
    range = c(0, 2.5),
    breaks = c(1, 2, 3),
    name = expression(atop("No. of genes", log[10](n)))
  ) +
  labs(
    x = "KOs",
    y = "Cell type",
    title = str_wrap("In vivo vs ex vivo KO effect concordance is gene-dependent", width = 70)
  ) +
  optimized_theme_fig() +
  theme(
    legend.position = "bottom",
    legend.justification = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Wrap plot with extra space on the right for the legend
p <- ggdraw(Fig3D) + 
  theme(plot.margin = margin(t = 2, r = 50, b = 10, l = 70))


n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$celltype))
# Save larger PDF; plot itself remains the same size
ggsave(
  filename = basedir("Fig3D_glioblastoma_cells.pdf"),
  plot = p,
  width =  0.45* n_col,   # bigger PDF width to include legend
  height = 4.5 * n_row,   # keep plot height small
  units = "cm"
)
###################


koi_t <- as.data.frame(koi,"koi")
write.table(koi_t,basedir("genotypes.tsv"), row.names = F)
################################################################################
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
ggsave(basedir("Correlation_density_in.vivo_ex.vivo.pdf"), w = 7, h = 3, units = "cm")
##############
#Fig3C current

InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")


limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
limmaRes_all <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))
#consistent
# function
sign <- limmaRes_NTC %>%
  filter(group == "n.s") %>%
  pull(genes)%>%
  unique()

get_consistent_genes <- function(limma_df, coef, pval_thresh = 0.01, logfc_thresh = 1) {
  coefs_in <- paste0("in.vivo",coef)
  coefs_ex <- paste0("ex.vivo",coef)
  coefs_int <- paste0("interaction",coef)
  genes_consistent <- limma_df %>%
    filter(coef %in% c(coefs_in,coefs_ex), group != "n.s",coef != coefs_int)%>%
    filter(adj.P.Val < pval_thresh, abs(logFC) > logfc_thresh) %>%
    group_by(ensg, celltype) %>%
    filter(n() == 2) %>%            # gene must be in both coefs per cell type
    arrange(desc(logFC)) %>%
    ungroup() %>%
    distinct(ensg)
  
  return(genes_consistent$ensg)
}
genes_Chd4 <- get_consistent_genes(limmaRes_all, coef = "Chd4")
genes_Kmt2d <- get_consistent_genes(limmaRes_all, coef = "Kmt2d")
genes_Cbx3 <- get_consistent_genes(limmaRes_all, coef = "Cbx3")



#
dataVoom_Eo.Ba <- read_rds(InDir_int("Eo.Ba_dataVoom.rds"))
dataVoom_Mono <- read_rds(InDir_int("Mono_dataVoom.rds"))
dataVoom_MkP <- read_rds(InDir_int("MkP_dataVoom.rds"))
dataVoom_GMP <- read_rds(InDir_int("GMP_dataVoom.rds"))
dataVoom_HSC <- read_rds(InDir_int("HSC_dataVoom.rds"))
dataVoom_MEP.early <- read_rds(InDir_int("MEP.early_dataVoom.rds"))
dataVoom_Gran. <- read_rds(InDir_int("Gran._dataVoom.rds"))
dataVoom_Gran.P <- read_rds(InDir_int("Gran.P_dataVoom.rds"))
KO <- koi[1]
ct <- unique(meta$celltype)[1]



dat.list <-list()
non_affected <- c("Chd4","Prmt5")
for (KO in c(selected_KOs,non_affected)){
  list_of_genes <- c("Oas2","Gbp3","Tnfaip6","Rasl2","Pgam2","Slc4a1", "Klk1",
                     "Gng3","Cebpa","Rab44",
                     "Oas3","Irf7","Gvin1","Ifit1","Myc", "Fxyd1",
                     "Msmo1","Idi1","Myc",
                     "Dppa5a","Rbakdn","Slc4a1","Aqp1","Myo1b",
                     "Atp7b",
                     "Rps27l","Rps2","Pop5","Myc","Bcl2",
                     "Stat5")
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

goi_exp %>% write_rds(basedir("expression.rds"))
goi_exp_only <- goi_exp

limmaRes_all$comparison <- gsub("^(ex\\.vivo|in\\.vivo|interaction)", "", limmaRes_all$coef)

limmaRes_all <- limmaRes_all %>%
  mutate(tissue = str_extract(limmaRes_all$coef, "^(ex\\.vivo|in\\.vivo)"),
                              gene = ensg)
  
goi_exp_limma <- merge(goi_exp, limmaRes_all, by = c("celltype", "comparison", "tissue", "gene"))

analyze_kos <- function(goi, ct, kos, effect_labels, goi_exp_limma, geneset) {
  
  # Step 1: Subset to the relevant gene + celltype
  filtered_data <- goi_exp_limma %>%
    filter(gene == goi, celltype == ct, comparison %in% kos)
  
  if (nrow(filtered_data) == 0) {
    message(paste("No data available for", goi, "in", ct))
    return(NULL)
  }
  
  # Step 2: Add significance from limma
  filtered_data <- filtered_data %>%
    mutate(significance = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE              ~ "ns"
    ))
  
  # y-position for significance labels per tissue & KO
  filtered_data <- filtered_data %>%
    group_by(comparison, tissue) %>%
    mutate(y_pos = max(E, na.rm = TRUE) * 1.1) %>%
    ungroup()
  
  # Step 3: Generate plots for each KO
  plots <- lapply(kos, function(KO) {
    subset_data <- filtered_data %>%
      filter(comparison == KO)
    
    if (nrow(subset_data) == 0) {
      message(paste("No data for", goi, "in", ct, "KO:", KO))
      return(NULL)
    }
    
    effect_label <- effect_labels[KO]
    
    p <- ggplot(subset_data, aes(x = genotype, y = E, color = tissue)) + 
      geom_boxplot(aes(color = tissue),
                   outlier.shape = NA,
                   position = position_dodge(width = 0.8),
                   size = 0.2) +
      # geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      #             alpha = 0.5) +
      facet_grid(
        cols = vars(tissue),
        scales = "free",
        labeller = labeller(tissue = c("ex.vivo" = "ex vivo", "in.vivo" = "in vivo"))
      ) +
      scale_color_manual(
        values = c("ex.vivo" = "#6a3d9aff", "in.vivo" = "#d38d5fff"),
        name = expression("Culture model")
      ) +
      labs(
        title = bquote(atop(.(paste0(goi, ": ", geneset)), .(ct))),
        y = "Expression") +
      xlab(paste0(KO, " KO (", effect_label, ")")) +
      theme(legend.position = "none") +
      optimized_theme_fig() +
      theme(panel.grid = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_text(
        data = subset_data %>% distinct(tissue, comparison, significance, y_pos),
        aes(x = 1.5, y = y_pos, label = significance),
        inherit.aes = FALSE,
        size = 2.5
      )
    
    return(p)
  })
  
  names(plots) <- kos
  
  return(list(
    # stat_tests = filtered_data %>% 
    #   select(gene, celltype, comparison, tissue, coefficient, logFC, P.Value, adj.P.Val, significance),
    plots = plots
  ))
}

goi_exp_limma <- goi_exp_limma %>%
  mutate(genotype = factor(genotype, levels = c("NTC", setdiff(unique(genotype), "NTC"))))

run_and_extract <- function(goi, ct, kos, effect_labels, geneset, goi_exp_limma) {
  result <- analyze_kos(
    goi = goi,
    ct = ct,
    kos = kos,
    effect_labels = effect_labels,
    goi_exp_limma = goi_exp_limma,
    geneset = geneset
  )
  
  # list(
  #   stat = result$stat_tests,
  #   plots = lapply(names(result$plots), function(koname) {
  #     result$plots[[koname]] + theme(legend.position = "none")
  #   }) %>% setNames(names(result$plots))
  
}

Klk1_Cbx3 <- run_and_extract(
  goi = "Klk1",
  ct = "Mono",
  kos = c("Cbx3"),
  effect_labels = c("Cbx3" = "Consistent trend"),
  geneset = "Serine protease",
  goi_exp_limma = goi_exp_limma
)
Gng3_Chd4 <- run_and_extract(
  goi = "Gng3", ct = "Mono", kos = c("Chd4"),
  effect_labels = c("Chd4" = "Consistent trend"),
  geneset = "G protein signaling",
  goi_exp_limma = goi_exp_limma
)

Ifit1_Brd9 <- run_and_extract(
  goi = "Ifit1", ct = "Eo.Ba", kos = c("Brd9"),
  effect_labels = c("Brd9" = "Opposite trend"),
  geneset = "ISG",
  goi_exp_limma = goi_exp_limma
)

# Run for all genes of interest
genes_to_run <- list(
  list(goi="Ifit1", ct="Eo.Ba", kos=c("Brd9"), effect_labels=c("Brd9"="Opposite trend"), geneset="ISG"),
  list(goi="Ifit1", ct="Eo.Ba", kos=c("Rcor1"), effect_labels=c("Rcor1"="De-novo effect"), geneset="ISG"),
  list(goi="Myc", ct="GMP", kos=c("Setdb1"), effect_labels=c("Setdb1"="Effect not captured"), geneset="Growth regulator"),
  list(goi="Atp7b", ct="Mono", kos=c("Cbx3","Brd9"), effect_labels=c("Cbx3"="Consistent trend","Brd9"="Opposite trend"), geneset="Copper homeostasis"),
  list(goi="Myc", ct="GMP", kos=c("Rcor1"), effect_labels=c("Rcor1"="No effect"), geneset="Growth regulator"),
  list(goi="Klk1", ct="Mono", kos=c("Cbx3"), effect_labels=c("Cbx3"="Consistent trend"), geneset="Serine protease"),
  list(goi="Cebpa", ct="GMP", kos=c("Brd9"), effect_labels=c("Brd9"="Effect not captured"), geneset="Myeloid activator")
  #list(goi="Rab44", ct="Eo.Ba", kos=c("Brd9"), effect_labels=c("Brd9"="Consistent trend"), geneset="Consistent trend")
)

results_list <- lapply(genes_to_run, function(x) {
  analyze_kos(x$goi, x$ct, x$kos, x$effect_labels, goi_exp_limma, x$geneset)
})


# Example: Combine first KO plots into a multi-panel figure
Fig.2F <- results_list[[1]]$plots[[1]] + 
  results_list[[7]]$plots[[1]] + 
  results_list[[3]]$plots[[1]] + 
  results_list[[2]]$plots[[1]] + 
  plot_layout(ncol=4, guides="collect") &
  theme(legend.position="right")




Fig.2F <- Fig.2F + 
  plot_annotation(
    title = str_wrap("Gene expression representing consistent and inconsistent KO-effects between experimental models", 80)
  ) &
  optimized_theme_fig(
               )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
Fig.2F
ggsave(
  filename = basedir(paste0("Fig3F",".pdf")),
  plot = Fig.2F,
  width = 18,
  height = 5 ,
  units = "cm"
)

Sup.Fig6B <- results_list[[4]]$plots[[1]] + 
  results_list[[6]]$plots[[1]] + 
  plot_layout(ncol=4, guides="collect") &
  theme(legend.position="right")
Sup.Fig6B <- Sup.Fig6B + 
  plot_annotation(
    title = str_wrap("Gene expression representing consistent KO-effects between experimental models", 80)
  ) &
  optimized_theme_fig(
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

Sup.Fig6B
ggsave(
  filename = basedir(paste0("Sup.Fig6B.pdf")),
  plot = Sup.Fig6B,
  width = 12,
  height = 5.5,
  units = "cm"
)

#####
#celltype distribution
inDir1 <- dirout_load("/Ag_SCRNA_05_01_UMAPs_and_functional.clusters")

xu <- xlab("UMAP 1")
yu <- ylab("UMAP 2")

# SETTINGS ----------------------------------------------------------------

annList <- readRDS(inDir1("ProjVivo_functional.clusters.RDS"))

# Other projections
umap.proj <- list(
  original = readRDS(inDir1("ProjMonocle.RDS")),
  #izzo = readRDS(dirout_load("Ag_SCRNA_UMAPs")("ProjIzzo.RDS")),
  in.vivo = readRDS(inDir1("ProjVivo.RDS")),
  in.vivo.X = readRDS(inDir1("ProjVivoX.RDS"))
)

# SIMPLE SETUP ENDS HERE ---------------------------------------------------------

mobjs <- list()

tissue<-c("ex.vivo","in.vivo")
for(tissuex in tissue){
  (base::load(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/Analysis//Ag_SCRNA_02_01_Integration/",tissuex,"/soupx/MonocleObject.RData")))
  mobjs[[tissuex]] <- monocle.obj
}

# Remove duplets
for(tissuex in names(mobjs)){
  mobjs[[tissuex]] <- mobjs[[tissuex]][,colnames(mobjs[[tissuex]]) %in% annList$rn]
}

# Initialize the guide and tissue columns in annList
annList$guide <- NA
annList$tissue <- NA

# Match and extract guide and tissue information from `ex.vivo`
match_ex_vivo <- match(annList$rn, colnames(mobjs$ex.vivo))
idx <- which(is.na(annList$guide) & !is.na(match_ex_vivo))
annList$guide[idx] <- mobjs$ex.vivo@colData$guide[match_ex_vivo[idx]]
#annList$guide[is.na(annList$guide) &!is.na(match_ex_vivo)] <- mobjs$ex.vivo@colData$guide[match_ex_vivo[(match_ex_vivo)]]
annList$tissue[!is.na(match_ex_vivo)] <- "ex.vivo"  # Assign "ex.vivo" for matched rows

# Match and extract guide and tissue information from `in.vivo`
match_in_vivo <- match(annList$rn, colnames(mobjs$in.vivo))
annList$guide[is.na(annList$guide) & !is.na(match_in_vivo)] <- mobjs$in.vivo@colData$guide[match_in_vivo[!is.na(match_in_vivo)]]
annList$tissue[is.na(annList$tissue) & !is.na(match_in_vivo)] <- "in.vivo"  # Assign "in.vivo" for matched rows

# Filter annList to keep only rows where guide is "NTC"
annList <- annList[annList$guide %in% c("NTC","Brd9_AS_70306", "Brd9_BR_48004", "Brd9_BR_48005"), ]
unique(annList$guide)

annList <- annList %>%
  mutate(
    functional.cluster = case_when(
      functional.cluster %in% c("MEP (G1)", "MEP (pert.)", "MEP (S)") ~ "MEP",  # Grouping all MEP categories to "MEP"
      functional.cluster == "Imm. B-cell" ~ "Imm.B.cell",  # Correct "Imm. B-cell" to "Imm.B.cell"
      is.na(functional.cluster) ~ "Unknown",  # Handle NA values if necessary
      TRUE ~ functional.cluster  # Leave other values unchanged
    )
  )%>%
  mutate(genotype = gsub("_.*","",guide))
unique(annList$genotype)
#annList$tissue <- gsub("ex.vivo", "ex.vivo", annList$tissue)
in.vivo.X <- umap.proj$in.vivo.X
# Filter in.vivo.X based on annList rn to keep only NTC samples

in.vivo.X <- inner_join(in.vivo.X, annList, by = c("rn", "tissue"))

# Generate hexbin object based on filtered in.vivo.X (only NTC samples)
hex.obj <- hexbin::hexbin(x = in.vivo.X$UMAP_1, y = in.vivo.X$UMAP_2, xbins = 100, IDs = TRUE)
in.vivo.X <- cbind(in.vivo.X, data.table(hex.x=hex.obj@xcm, hex.y=hex.obj@ycm, hex.cell=hex.obj@cell)[match(hex.obj@cID, hex.cell),])
pDT <- in.vivo.X
pDT <- pDT[, .N, by = c("hex.x", "hex.y", "functional.cluster")]
pDT[, sum := sum(N), by = c("hex.x", "hex.y")]
pDT[, frac := N / sum]

# Filter clusters with significant fractions
pDT <- pDT[frac > 0.25]

# Merge summary data back with the original dataset
merged_data <- inner_join(in.vivo.X, pDT, by = c("hex.x", "hex.y", "functional.cluster"), all.x = TRUE)


# Check the unique values in functional.cluster to ensure it's working

# Generate cluster labels for significant clusters
pDT.labels <- pDT[frac > 0.25, .(hex.x = median(hex.x), hex.y = median(hex.y)), by = c("functional.cluster")]
pDT.labels %>% write_rds(basedir("pDT.labels.rds"))
# color coding
cluster_colors <- c(
  "Mono" = "#E69F00",      # Orange
  "Eo/Ba" = "#56B4E9",     # Sky Blue
  "GMP" = "#009E73",       # Green
  "MEP (early)" = "#F0E442", # Yellow
  "MkP" = "#CC79A7",       # Pink/Purple
  "Gran. P" = "#0072B2",   # Blue
  "Gran." = "#D55E00",     # Reddish Orange
  "HSC" = "#A020F0",       # Purple
  "GMP (early)" = "#999999",  # Light Gray (unchanged)
  "CLP" = "#D9D9D9",       # Light gray for CLP
  "unclear" = "#B0B0B0",    # Gray for unclear
  "Imm. B-cell" = "#8DA0CB", # Soft Blue
  "MEP" = "#D3D3D3",        # Lighter gray for MEP
  "Ery" = "#A9A9A9",        # Slightly darker gray for Ery
  "Imm.B.cell" = "gray"     # Other
)

unique(merged_data$guide)
merged_data %>% write_rds(basedir("Cross_projected_on_in.vivo.rds"))
# Ensure factor ordering for correct label display
merged_data$functional.cluster <- factor(merged_data$functional.cluster, 
                                         levels = names(cluster_colors))

ggplot(merged_data[tissue != "leukemia"], aes(x = UMAP_1, y = UMAP_2)) + 
  
  geom_point(aes(color = functional.cluster), size = 0.00000001 ) + 
  
  geom_text_repel(data = pDT.labels %>%
                    filter(functional.cluster %in% c("Mono", "Eo/Ba", "GMP", "MEP (early)",
                                                     "MkP", "Gran. P", "Gran.", "HSC"
                    )),
                  aes(x = hex.x, y = hex.y, label = functional.cluster),
                  size = 1,                  # Adjust text size
                  box.padding = 0.21,         # Distance from points
                  point.padding = 0.21,       # Distance from label anchor
                  segment.color = "black",   # Line color
                  segment.size = 0.004,        # Line thickness
                  force = 10,                # Repelling force
                  max.overlaps = Inf) +
  facet_grid(cols = vars(tissue), rows = vars(genotype)) + 
  # Defining color manual scale for clusters
  # Defining color manual scale for clusters
  scale_color_manual(name = "functional.cluster", values = cluster_colors) + 
  
  # Adjusting the alpha scale for fraction
  #scale_alpha_continuous(name = "Fraction", range = c(0, 1)) +  # To use fraction values for transparency
  
  # Adjusting the theme
  optimized_theme_fig() +
  
  # Positioning legends separately
  # Positioning legends separately
  theme(
    axis.text.x = element_text(angle = 0),
    legend.position = "bottom",  # Color legend at the bottom
    legend.box = "horizontal",   # Horizontal alignment of legends
    legend.text = element_text(size = 5),     # Adjust legend text font size
    legend.title = element_blank(), # Remove title for color legend
    legend.spacing = unit(0.5, "cm"),  # Adjust spacing between legends
    legend.key.size = unit(0.5, "lines")  # Adjust size of the legend keys
  ) 
# Correct axis labels (you can define xu and yu separately if needed)
ggsave(outBase("UMAP_InvivoX_NTC.pdf"), w=3.5,h=2.5)


ggplot(merged_data[tissue != "leukemia"], aes(x=UMAP_1, UMAP_2)) + 
  themeNF() +
  stat_binhex(aes(fill=log10(..count..)), bins=100) + 
  scale_fill_gradientn(colors = c("lightgrey", "#1f78b4", "#e31a1c", "#ff7f00")) +
  facet_grid(. ~ tissue) +
  xu + yu+
  optimized_theme_fig()
ggsave(basedir("Sup.Fig.7.UMAP_InvivoX_NTC_Brd9_distribution.pdf"), w = 10,h= 8, units = "cm")
ggsave(basedir("Sup.Fig.7.UMAP_InvivoX_NTC_Brd9_distribution.png"), w = 10,h= 8, units = "cm")
######################
InDir6 <- dirout("Ag_ScRNA_12_fgsea_overlap")

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

combined <- combined %>%
  mutate(
    log2.or.filtered = ifelse(overlap > 5, pmin(log2.odds.ratio, 7), NA)  # keep value only if overlap>5
  )

Fig3E <- ggplot(combined, aes(x = coef, y = log2.or.filtered)) +
  geom_col(color = "darkgrey", fill =NA,  width = 0.6) +
  
  # significance labels (for overlap > 5 only, already handled by filtering)
  geom_text(
    data = combined %>% filter(overlap > 5),
    aes(label = significance_en),
    y = 6.5,
    color = "black",
    size = 1.5
  ) +
  
  # add "NA" label for overlap <= 5
  geom_text(
    data = combined %>% filter(overlap <= 5),
    aes(label = "NA"),
    y = 1.5,   # adjust position
    color = "black",
    size = 1.5
  ) +
  
  facet_grid(cols = vars(celltype), scales = "free", space = "free_x") +
  labs(
    x = "KOs",
    y = TeX("$\\log_{2}\\; (Odds ratio)$"),
    title = "KOs with large overlap to culture effect genes show more discordant effects between ex vivo and in vivo models"
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),  # 👈 center facet labels
    panel.spacing = unit(0.02, "cm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


Fig3E

ggsave(basedir("Fig3E.pdf"),plot=Fig3E,
       w=14.8,h=4, units = "cm")
#

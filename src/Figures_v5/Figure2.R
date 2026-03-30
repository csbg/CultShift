#
#load libraries and functions--------------------------------------------------
#
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
source("src/Ag_ko_classification.R")
source("src/Ag_enrichR_mouse_genes.R")
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
base <- "Figure2"
basedir <- dirout("Figure2")

ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)
IFN_genes = union(enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`,
                              enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`) 

IFN_genes <- union(ISG_core, IFN_genes)
IFN_genes <-IFN_genes[IFN_genes != ""]
#Data and function

#Fig2A------------
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

#fig-----
# Create the plot
Fig2A <- correlation_deg_flagged %>%
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

Fig2A
#combine with example
Fig2A_example <- merged_logFC %>%
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
  filename = basedir("Fig2A_example_wd.pdf"),
  plot = Fig2A_example,
  width = 4.5,
  height = 2.5,
  units = "cm"
)

n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$celltype))
ggsave(
  filename = basedir("Fig2A_perturb.pdf"),
  plot = Fig2A,
  width = n_col * 0.28,    # PDF becomes wider
  height = n_row * 0.8,
  units = "cm"
)

# Save the version with the legend
Fig2A_all <- Fig2A_example + Fig2A + plot_layout(widths = c(0.8, 5))
#paper--------------
ggsave(
  filename = basedir("Fig.2A.pdf"),
  plot = Fig2A,
  width = 18,
  height = 6.5,
  units = "cm"
)

# splenic----------
out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"

res <- rbind(read.delim(paste0(out,"/DEG_ResultsT8.tsv")),
             read.delim(paste0(out,"/DEG_ResultsM.tsv"))
)
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

# Combine datasets
merged_data_spleen <- ex_vivo_data %>%
  inner_join(in_vivo_data, by = c("gene", "cell_type", "genotype"), 
             suffix = c("_ex.vivo", "_in.vivo"))


correlation_results <- merged_data_spleen %>%
  group_by(cell_type, genotype) %>%
  summarise(correlation = cor(logFC_ex.vivo, logFC_in.vivo, use = "complete.obs")) %>%
  ungroup()

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

# Plot with genotypes ordered by aggregated correlation
Fig2B <- deg_plot_data %>%
  filter(genotype %in% c("IRF9KO","STAT1KO","STAT2KO","TYK2CMV"))%>%
  mutate(genotype = ifelse(grepl("IRF9KO",genotype),"Irf9KO",ifelse(
    grepl("STAT1KO",genotype),"Stat1KO",ifelse(
      grepl("STAT2KO",genotype),"StatKO","Tyk2CMV"
    )))) %>%
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
p <- ggdraw(Fig2B) + 
  theme(plot.margin = margin(t = 2, r = 50, b = 10, l = 70))

n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$cell_type))
# Save larger PDF; plot itself remains the same size
ggsave(
  filename = basedir("Fig2B_splenic_cells.pdf"),
  plot = p,
  width = 1 * n_col,   # bigger PDF width to include legend
  height = 2.5 * n_row,   # keep plot height small
  units = "cm"
)
#Glioblastoma-----------------
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
Fig2C <- deg_plot_data %>%
  
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
p <- ggdraw(Fig2C) + 
  theme(plot.margin = margin(t = 2, r = 50, b = 10, l = 70))


n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$celltype))
# Save larger PDF; plot itself remains the same size
ggsave(
  filename = basedir("Fig2C_glioblastoma_cells.pdf"),
  plot = p,
  width =  0.45* n_col,   # bigger PDF width to include legend
  height = 4.5 * n_row,   # keep plot height small
  units = "cm"
)

# ggsave(
#   filename = basedir("row1_for_legend.pdf"),
#   plot = Fig2A_all,
#   width = 18,
#   height = 6.5,
#   units = "cm"
# )

koi_t <- as.data.frame(koi,"koi")
write.table(koi_t,basedir("genotypes.tsv"), row.names = F)
################################################################################

# merged_logFC <- read_rds(InDir_cor("in.vivo_ex.vivo_logFC.rds"))
# correlation_matrix_data <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds"))
# deg_plot_data <- read_rds(InDir_cor("DEGs_per_tissue.rds"))%>%
#   dplyrselect(-correlation)
# 
# # Merge with KO flags to include valid KO status
# correlation_deg_flagged <- correlation_deg %>%
#   inner_join(ko_flags, by = c("genotype", "celltype")) %>%
#   filter(valid_ko)%>%
#   na.omit()%>%
#   filter(valid_ko)
# # Update the genotype factor levels for plotting
# correlation_deg_flagged$genotype <- factor(correlation_deg_flagged$genotype,
#                                            levels = column_order)
# 
# # Calculate mean correlation per KO for ordering
# ko_order <- correlation_deg_flagged %>%
#   group_by(genotype) %>%
#   summarise(mean_corr = mean(correlation, na.rm = TRUE)) %>%
#   arrange(desc(mean_corr)) %>%
#   pull(genotype)
# 
# # Update factor levels of genotype based on mean correlation
# correlation_deg_flagged$genotype <- factor(correlation_deg_flagged$genotype, levels = ko_order)
# 
# 
# #
# Fig2B <- correlation_deg_flagged %>%
#   filter(valid_ko) %>%
#   filter(coef %in% koi)%>%
#   ggplot(aes(
#     x = reorder(genotype, correlation, FUN = mean),
#     y = correlation
#   )) +
#   
#   geom_boxplot(outlier.shape = NA, fill = "white", color = "#213c75ff",
#                width = 0.6,) +
#   geom_jitter(aes(
#     shape = celltype,
#     size = pmin(3, log10(num_degs))
#   ),
#   width = 0.15,
#   alpha = 0.7, color = "#808080ff"
#   ) +
#   scale_shape_manual(
#     values = 0:7,
#     name = expression("Cell type")
#   )+
# # Use distinct shapes (you can customize)
#   scale_size_continuous(
#     range = c(0, 1.8),
#     breaks = c(1, 2, 3),
#     name = expression(atop("No. of genes", log[10](n)))
#   ) +
#   labs(
#     title = "KOs (ranked by mean correlation)",
#     y = "Pearson's correlation (in vivo vs ex vivo)",
#     x = "KOs"
#   ) +
#   optimized_theme_fig() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     legend.position = "right"
#   )
# 
# Fig2B
# 
# ggsave(
#   filename = basedir("Fig2B.pdf"),
#   plot = Fig2B,
#   width = 18,
#   height = 4.5,
#   units = "cm"
# )

#split 

##############
#Fig2C current

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


# Ifit1_Rcor1 <- run_and_extract(
#   goi = "Ifit1", ct = "Eo.Ba", kos = c("Rcor1"),
#   effect_labels = c("Rcor1" = "De-novo effect"),
#   geneset = "ISG", goi_exp = goi_exp, limmaRes = limmaRes
# )
# 
# Myc_GMP_Setdb1 <- run_and_extract(
#   goi = "Myc", ct = "GMP", kos = c("Setdb1"),
#   effect_labels = c("Setdb1" = "Effect not captured"),
#   geneset = "Growth regulator", goi_exp = goi_exp, limmaRes = limmaRes
# )
# 
# Atp7b <- run_and_extract(
#   goi = "Atp7b", ct = "Mono", kos = c("Cbx3","Brd9"),
#   effect_labels = c("Cbx3" = "Consistent trend"),  # FIXED label key
#   geneset = "Copper homeostasis", goi_exp = goi_exp, limmaRes = limmaRes
# )
# 
# Rcor1_Myc_GMP <- run_and_extract(
#   goi = "Myc", ct = "GMP", kos = c("Rcor1"),
#   effect_labels = c("Rcor1" = "No effect"),
#   geneset = "Growth regulator", goi_exp = goi_exp, limmaRes = limmaRes
# )
# Run for all genes of interest
genes_to_run <- list(
  list(goi="Ifit1", ct="Eo.Ba", kos=c("Brd9"), effect_labels=c("Brd9"="Opposite trend"), geneset="ISG"),
  list(goi="Ifit1", ct="Eo.Ba", kos=c("Rcor1"), effect_labels=c("Rcor1"="De-novo effect"), geneset="ISG"),
  list(goi="Myc", ct="GMP", kos=c("Setdb1"), effect_labels=c("Setdb1"="Effect not captured"), geneset="Growth regulator"),
  list(goi="Atp7b", ct="Mono", kos=c("Cbx3","Brd9"), effect_labels=c("Cbx3"="Consistent trend","Brd9"="Opposite trend"), geneset="Copper homeostasis"),
  list(goi="Myc", ct="GMP", kos=c("Rcor1"), effect_labels=c("Rcor1"="No effect"), geneset="Growth regulator"),
  list(goi="Klk1", ct="Mono", kos=c("Cbx3"), effect_labels=c("Cbx3"="Consistent trend"), geneset="Serine protease"),
  list(goi="Cebpa", ct="GMP", kos=c("Brd9"), effect_labels=c("Brd9"="Effect not captured"), geneset="Myeloid activator")
  
)

results_list <- lapply(genes_to_run, function(x) {
  analyze_kos(x$goi, x$ct, x$kos, x$effect_labels, goi_exp_limma, x$geneset)
})


# Example: Combine first KO plots into a multi-panel figure
Fig.2D <- results_list[[1]]$plots[[1]] + 
  results_list[[7]]$plots[[1]] + 
  results_list[[3]]$plots[[1]] + 
  results_list[[2]]$plots[[1]] + 
  plot_layout(ncol=4, guides="collect") &
  theme(legend.position="right")




Fig.2D <- Fig.2D + 
  plot_annotation(
    title = str_wrap("Gene expression representing consistent and inconsistent KO-effects between experimental models", 80)
  ) &
  optimized_theme_fig(
               )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  

  

Fig.2D
ggsave(
  filename = basedir(paste0("Fig2D",".pdf")),
  plot = Fig.2D,
  width = 18,
  height = 5 ,
  units = "cm"
)
Sup.Fig3C <- Ifit1_Rcor1$plots[[1]] + Rcor1_Myc_GMP$plots[[1]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "right"  # removes grid lines
  )
ggsave(
  filename = basedir(paste0("Sup.Fig3C",".pdf")),
  plot = Sup.Fig3C,
  width = 9,
  height = 5 ,
  units = "cm"
)


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

Fig2E <- ggplot(combined, aes(x = coef, y = log2.or.filtered)) +
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


Fig2E

ggsave(basedir("Fig2E.pdf"),plot=Fig2E,
       w=14.8,h=4, units = "cm")
#

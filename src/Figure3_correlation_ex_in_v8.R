#
#load libraries and functions--------------------------------------------------
#
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
source("src/Ag_ko_classification.R")

library(scales)
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
base <- "Figure3_correlation_ex_in_v8"
basedir <- dirout("Figure3_correlation_ex_in_v8")

ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%pull(value)
IFN_genes = union(enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`,
                  enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`) 

IFN_genes <- union(ISG_core, IFN_genes)
IFN_genes <-IFN_genes[IFN_genes != ""]
#Data and function

#Fig3A------------

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
correlation_deg_flagged %>%
  ungroup() %>%   # remove grouping
  summarise(
    greater_than_0.5 = sum(correlation > 0.5, na.rm = TRUE),
    less_than_0.5    = sum(correlation < 0.5, na.rm = TRUE)
  )


#alternate------------
cluster_colors <- c(
  "Mono" = "#E69F00",      # Orange
  "Eo.Ba" = "#56B4E9",     # Sky Blue
  "GMP" = "#009E73",       # Green
  "MEP.early" = "#F0E442", # Yellow
  "MkP" = "#CC79A7",       # Pink/Purple
  "Gran.P" = "#0072B2",   # Blue
  "Gran." = "#D55E00",     # Reddish Orange
  "HSC" = "#A020F0",       # Purple
  "GMP (early)" = "#999999",  # Light Gray (unchanged)
  "CLP" = "#D9D9D9",       # Light gray for CLP
  "unclear" = "#B0B0B0",    # Gray for unclear
  "Imm. B-cell" = "#8DA0CB", # Soft Blue
  "MEP" = "pink",        # Lighter gray for MEP
  "Ery" = "#A9A9A9",        # Slightly darker gray for Ery
  "Imm.B.cell" = "gray"     # Other
)
plot_data <- correlation_deg_flagged %>%
  filter(genotype %in% koi) %>%
  mutate(celltype_factor = factor(celltype, levels = names(cluster_colors)))
# Plot with exact cluster_colors mapping
Fig3A <- ggplot(plot_data, aes(x = genotype, y = correlation)) +
  
  geom_point(aes(
    size = pmin(3, log10(num_degs)),
    color = celltype_factor  # map color directly
  ),
  stroke = 0.7) +  # optional outline for visibility
  
  # reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", color = "darkred") +
  
  # apply your exact cluster colors
  scale_color_manual(values = cluster_colors, name = "Cell type") +
  
  # size scale for number of DEGs
  scale_size_continuous(
    range = c(0.1, 1.5),
    name = expression(log[10]("No. of DEGs"))
  ) +
  ylim(c(-1,1))+
  
  labs(
    x = "KOs",
    y = "Correlation",
    title = "Corralation reveals discordant KO effects in vivo vs ex vivo KO effect in hematopoietic cells"
  ) +
  
  optimized_theme_fig() +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

Fig3A
ggsave(basedir("Fig3A.pdf"))

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
  filename = basedir("Fig3A_example_wd.pdf"),
  plot = Fig3A_example,
  width = 4.5,
  height = 2.5,
  units = "cm"
)

n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$celltype))
ggsave(
  filename = basedir("Fig3A_perturb.pdf"),
  plot = Fig3A,
  width = n_col * 0.28,    # PDF becomes wider
  height = n_row * 0.6,
  units = "cm"
)

# Save the version with the legend
Fig3A_B_all <- Fig3A_example + Fig3A + plot_layout(widths = c(0.8, 5))
#paper--------------
ggsave(
  filename = basedir("Fig.3A_B.pdf"),
  plot = Fig3A_B_all,
  width = 12,
  height = 5.2,
  units = "cm"
)
############################
#KO target expression
limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
KOs <- limmaRes %>%
  pull(coef)%>%
  unique()
coef_sign <- limmaRes_NTC %>% filter(genes %in% KOs) %>%
  filter(group != "n.s")
coef_logFC <- limmaRes_NTC %>% filter(genes %in% KOs)%>%
  filter(genes %in% koi) %>%
  mutate(genes = factor(genes, levels = rev(ko_order)))

ggplot(coef_logFC, aes(
  y = celltype,
  x = genes,
  color = pmin(2, pmax(-2, logFC)),
  size = pmin(3, -log10(adj.P.Val))
)) +  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    limits = c(-2,2),
    name = TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0, 1.8),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    title = "Culture effect in target genes",
    x = "cell type",
    y = "Genes"
  ) +
  
  optimized_theme_fig() +
  theme(
    legend.position = "right",
    strip.text.x = element_text(angle = 45, hjust = 0, vjust = 0)
  )
ggsave(basedir("Fig3C.pdf"), w = 11, h = 4.5,  units = "cm")
############################
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
deg_plot_data <- deg_plot_data %>%
  #filter(genotype %in% c("IRF9KO","STAT1KO","STAT2KO","TYK2CMV"))%>%
  filter(cell_type == "T8") %>%
  mutate(genotype = ifelse(grepl("IRF9KO",genotype),"Irf9KO",ifelse(
    grepl("STAT1KO",genotype),"Stat1KO",ifelse(
      grepl("STAT2KO",genotype),"StatKO",ifelse(
        grepl("TYK2CMV",genotype),"Tyk2CMV",ifelse(
          grepl("STAT3VAV",genotype),"Stat3Vav",ifelse(
            grepl("STAT1BB",genotype),"Stat1BB","Tyk2KE"
          )))))))
# Compute mean correlation per genotype across cell types
ko_order <- deg_plot_data %>%
  group_by(genotype) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_corr)) %>%  # descending order
  pull(genotype)

# Update the factor levels for plotting
deg_plot_data <- deg_plot_data %>%
  mutate(genotype = factor(genotype, levels = ko_order))
unique(deg_plot_data$genotype)
# Plot with genotypes ordered by aggregated correlation

##################################################################
#alternate
Fig3D <- ggplot(deg_plot_data %>% filter(cell_type == "T8"), aes(x = genotype, y = correlation)) +
  
  geom_point(aes(
    size = pmin(3, log10(num_degs)),
    fill = "black"  # map color directly
  ),
  stroke = 0.7) +  # optional outline for visibility
  
  # reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", color = "darkred") +
  
  # apply your exact cluster colors
  #scale_color_manual(values = cluster_colors, name = "Cell type") +
  
  # size scale for number of DEGs
  scale_size_continuous(
    range = c(0.1, 1.5),
    name = expression(log[10]("No. of DEGs"))
  ) +
  ylim(c(-1,1))+
  
  labs(
    x = "KOs",
    y = "Correlation",
    title = "Corralation reveals discordant KO effects in vivo vs ex vivo KO effect in splenic T-cells"
  ) +
  
  optimized_theme_fig() +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )
Fig3D
# Wrap plot with extra space on the right for the legend


n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$cell_type))
# Save larger PDF; plot itself remains the same size
ggsave(
  filename = basedir("Fig3C_splenic_cells.pdf"),
  plot = Fig3D,
  width = 0.8 * n_col,   # bigger PDF width to include legend
  height = 5,   # keep plot height small
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
  summarise(correlation = cor(logFC_ex.vivo, 
                              logFC_in.vivo, use = "complete.obs")) %>%
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
deg_counts <- full_join(deg_ex_vivo, deg_in_vivo, 
                        by = c("genotype")) %>% na.omit()
deg_plot_data <- deg_counts %>%
  pivot_longer(cols = starts_with("degs"), names_to = "condition",
               values_to = "num_degs") %>%
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
Fig3E <- ggplot(deg_plot_data , aes(x = genotype, y = correlation)) +
  
  geom_point(aes(
    size = pmin(3, log10(num_degs)),
    fill = "black"  # map color directly
  ),
  stroke = 0.7) +  # optional outline for visibility
  
  # reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", color = "darkred") +
  
  # apply your exact cluster colors
  #scale_color_manual(values = cluster_colors, name = "Cell type") +
  
  # size scale for number of DEGs
  scale_size_continuous(
    range = c(0.1, 1.5),
    name = expression(log[10]("No. of DEGs"))
  ) +
  ylim(c(-1,1))+
  
  labs(
    x = "KOs",
    y = "Correlation",
    title = "Corralation reveals discordant KO effects in vivo vs ex vivo KO effect in glioblastoma model"
  ) +
  
  optimized_theme_fig() +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )
Fig3E
# Wrap plot with extra space on the right for the legend


n_col = length(unique(deg_plot_data$genotype))
n_row = length(unique(deg_plot_data$celltype))
# Save larger PDF; plot itself remains the same size
ggsave(
  filename = basedir("Fig3E_glioblastoma.pdf"),
  plot = Fig3D,
  width =16.5,   # bigger PDF width to include legend
  height = 5,   # keep plot height small
  units = "cm"
)
############################################################################

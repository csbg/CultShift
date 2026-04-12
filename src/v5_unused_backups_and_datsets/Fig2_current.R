
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

# splenic
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
#Glioblastoma
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

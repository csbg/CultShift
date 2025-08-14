#
#load libraries and functions--------------------------------------------------
#
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
source("src/Ag_ko_classification.R")
#source("src/Ag_enrichR_mouse_genes.R")
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

merged_logFC <- read_rds(InDir_cor("in.vivo_ex.vivo_logFC.rds"))
correlation_matrix_data <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(InDir_cor("DEGs_per_tissue.rds"))%>%
  select(-correlation)

# Prepare correlation data
correlation_matrix_no_na <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)

# Perform hierarchical clustering
hc_cols <- hclust(dist(t(correlation_matrix_no_na)), method = "ward.D2")
column_order <- colnames(correlation_matrix_data)[hc_cols$order]
write_rds(column_order,basedir("column_order.rds"))
# Reshape correlation data
correlation <- correlation_matrix_data %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(cols = colnames(correlation_matrix_data),
               names_to = "genotype",
               values_to = "correlation",
               values_drop_na = FALSE)


# Prepare DEG data
deg_plot_data <- deg_plot_data %>%
  na.omit() %>%
  group_by(genotype, celltype) %>%
  filter(num_degs == max(num_degs, na.rm = TRUE)) %>%
  dplyr::select(-condition) %>%
  unique()

# Merge DEG data with correlation data
correlation_deg <- inner_join(deg_plot_data, correlation,
                              by = c("celltype", "genotype"))
correlation_deg %>% write_rds(basedir("correlation_deg.rds"))
# Merge with KO flags to include valid KO status
correlation_deg_flagged <- correlation_deg %>%
  inner_join(ko_flags, by = c("genotype", "celltype")) %>%
  na.omit()%>%
  filter(valid_ko)
# Update the genotype factor levels for plotting
correlation_deg_flagged$genotype <- factor(correlation_deg_flagged$genotype,
                                           levels = column_order)

#fig-----
# Create the plot
Fig2A <- correlation_deg_flagged %>%
 #filter(genotype %in% koi) %>%
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
       title =  "Correlation of KO-effects (in vivo versus ex vivo)") +
  optimized_theme_fig()+
  theme(
    
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
    x = "logFC (Ex Vivo)",
    y = "logFC (In Vivo)",
    
    fill = "Gene Count"
  )+ optimized_theme_fig()
ggsave(
  filename = basedir("Fig2A_example.pdf"),
  plot = Fig2A_example,
  width = 4,
  height = 2.5,
  units = "cm"
)
# Save the version with the legend
Fig2A_all <- Fig2A_example + Fig2A + plot_layout(widths = c(1, 5))
#paper--------------
ggsave(
  filename = basedir("row1_for_fig.pdf"),
  plot = Fig2A_all,
  width = 18,
  height = 5,
  units = "cm"
)
ggsave(
  filename = basedir("row1_for_legend.pdf"),
  plot = Fig2A_all,
  width = 18,
  height = 6,
  units = "cm"
)

koi_t <- as.data.frame(koi,"koi")
write.table(koi_t,basedir("genotypes.tsv"), row.names = F)
################################################################################
InDir_cor <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/")
merged_logFC <- read_rds(InDir_cor("in.vivo_ex.vivo_logFC.rds"))
correlation_matrix_data <- read_rds(InDir_cor("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(InDir_cor("DEGs_per_tissue.rds"))%>%
  select(-correlation)

# Merge with KO flags to include valid KO status
correlation_deg_flagged <- correlation_deg %>%
  inner_join(ko_flags, by = c("genotype", "celltype")) %>%
  filter(valid_ko)%>%
  na.omit()%>%
  filter(valid_ko)
# Update the genotype factor levels for plotting
correlation_deg_flagged$genotype <- factor(correlation_deg_flagged$genotype,
                                           levels = column_order)

# Calculate mean correlation per KO for ordering
ko_order <- correlation_deg_flagged %>%
  group_by(genotype) %>%
  summarise(mean_corr = mean(correlation, na.rm = TRUE)) %>%
  arrange(desc(mean_corr)) %>%
  pull(genotype)

# Update factor levels of genotype based on mean correlation
correlation_deg_flagged$genotype <- factor(correlation_deg_flagged$genotype, levels = ko_order)


#
Fig2B <- correlation_deg_flagged %>%
  filter(valid_ko) %>%
  ggplot(aes(
    x = reorder(genotype, correlation, FUN = mean),
    y = correlation
  )) +
  
  geom_boxplot(outlier.shape = NA, fill = "white", color = "#213c75ff",
               width = 0.6,) +
  geom_jitter(aes(
    shape = celltype,
    size = pmin(3, log10(num_degs))
  ),
  width = 0.15,
  alpha = 0.7, color = "#808080ff"
  ) +
  scale_shape_manual(
    values = 0:6,
    name = expression("Cell type")
  )+
# Use distinct shapes (you can customize)
  scale_size_continuous(
    range = c(0, 1.8),
    breaks = c(1, 2, 3),
    name = expression(atop("No. of genes", log[10](n)))
  ) +
  labs(
    title = "KOs (ranked by mean correlation)",
    y = "Pearson's correlation (in vivo vs ex vivo)",
    x = "KOs"
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

Fig2B
ggsave(
  filename = basedir("Fig2B.pdf"),
  plot = Fig2B,
  width = 18,
  height = 4.5,
  units = "cm"
)

##############
#Fig2C current

InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")


limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
limmaRes_all <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))
#consistent
genes_consistent <- limmaRes_all %>%


# function


get_consistent_genes <- function(limma_df, coef, pval_thresh = 0.01, logfc_thresh = 1) {
  coefs <- c(paste0("in.vivo",coef), paste0("ex.vivo",coef))
  
  genes_consistent <- limma_df %>%
    filter(coef %in% coefs) %>%
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
  list_of_genes <- c("Oas2","Gbp3","Tnfaip6",
                     "Oas3","Irf7","Gvin1","Ifit1","Myc_GMP",
                     "Msmo1","Mthfd2","Idi1","Ccnd1","Myc",
                     "Dppa5a","Rbakdn","Pcbp4","Aqp1","Myo1b",
                     "Rgs13","Atp7b",
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
#function
stat_tests_all <- list()  # Initialize a list to store results

analyze_kos <- function(goi, ct, kos, effect_labels, goi_exp, limmaRes,geneset) {
  stat_tests_all <- list()  # Initialize storage for statistical results
  
  # Step 1: Perform statistical tests
  for (KO in kos) {
    filtered_data <- goi_exp %>%
      filter(gene == goi, celltype == ct, comparison == KO)
    
    filtered_data$genotype <- factor(filtered_data$genotype,
                                     levels = c("NTC", setdiff(filtered_data$genotype, "NTC")))
    filtered_data$E <- as.numeric(filtered_data$E)
    
    if (nrow(filtered_data) > 0) {
      for (tissue_type in unique(filtered_data$tissue)) {
        filtered_tissue <- filtered_data %>% filter(tissue == tissue_type)
        
        if (nrow(filtered_tissue) > 1 && length(unique(filtered_tissue$genotype)) > 1) {
          E_values <- filtered_tissue$E
          normality_p <- shapiro.test(E_values)$p.value
          test_method <- ifelse(normality_p > 0.05, "t.test", "wilcox.test")
          
          stat_test <- tryCatch({
            compare_means(E ~ genotype, data = filtered_tissue, method = test_method) %>%
              mutate(tissue = tissue_type, KO = KO)
          }, error = function(e) {
            message(paste("compare_means failed for", KO, "in", tissue_type, ":", e$message))
            NULL
          })
          
          if (!is.null(stat_test)) {
            stat_tests_all[[paste(KO, tissue_type, sep = "_")]] <- stat_test
          }
        }
      }
    } else {
      message(paste("No data for", goi, "in", ct, "KO:", KO))
    }
  }
  
  # Combine all stats
  stat_tests_combined <- bind_rows(stat_tests_all) %>%
    mutate(significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
  
  # Step 2: Generate plots
  plots <- lapply(kos, function(KO) {
    filtered_limma <- limmaRes %>%
      filter(ensg == goi, coef == KO, celltype == ct)
    
    effect_label <- effect_labels[KO]
    
    filtered_data <- goi_exp %>%
      filter(gene == goi, celltype == ct, comparison == KO)
    filtered_data$genotype <- factor(filtered_data$genotype,
                                     levels = c("NTC", setdiff(filtered_data$genotype, "NTC")))
    
    if (nrow(filtered_data) > 0) {
      stat_subset <- stat_tests_combined %>%
        filter(KO == !!KO) %>%
        select(tissue, significance)
      
      annotation_data <- filtered_data %>%
        group_by(tissue) %>%
        summarize(y_pos = max(E, na.rm = TRUE) * 0.8, .groups = "drop") %>%
        left_join(stat_subset, by = "tissue")
      
      p <- ggplot(filtered_data, aes(x = genotype, y = E, color = tissue)) + 
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
              panel.grid.minor = element_blank())+
        geom_text(data = annotation_data,
                  aes(x = 1.5, y = y_pos, label = significance),
                  inherit.aes = FALSE,
                  size = 4)
      
      return(p)
    } else {
      message(paste("No data available for", goi, "in", ct, "KO:", KO))
      return(NULL)
    }
  })
  names(plots) <- kos
  
  return(list(
    stat_tests = stat_tests_combined,
    plots = plots
  ))
}
run_and_extract <- function(goi, ct, kos, effect_labels, geneset, goi_exp, limmaRes) {
  result <- analyze_kos(
    goi = goi,
    ct = ct,
    kos = kos,
    effect_labels = effect_labels,
    goi_exp = goi_exp,
    limmaRes = limmaRes,
    geneset = geneset
  )
  
  list(
    stat = result$stat_tests,
    plots = lapply(names(result$plots), function(koname) {
      result$plots[[koname]] + theme(legend.position = "none")
    }) %>% setNames(names(result$plots))
  )
}
# Run all panels
Ifit1_Brd9 <- run_and_extract(
  goi = "Ifit1", ct = "Eo.Ba", kos = c("Brd9"),
  effect_labels = c("Brd9" = "Opposite effect"),
  geneset = "ISG", goi_exp = goi_exp, limmaRes = limmaRes
)
Ifit1_Rcor1 <- run_and_extract(
  goi = "Ifit1", ct = "Eo.Ba", kos = c("Rcor1"),
  effect_labels = c("Rcor1" = "No effect"),
  geneset = "ISG", goi_exp = goi_exp, limmaRes = limmaRes
)

Myc_GMP_Brd9 <- run_and_extract(
  goi = "Myc", ct = "GMP", kos = c("Brd9"),
  effect_labels = c("Brd9" = "De novo effect"),
  geneset = "Protease inhibitor", goi_exp = goi_exp, limmaRes = limmaRes
)

Atp7b <- run_and_extract(
  goi = "Atp7b", ct = "Mono", kos = c("Cbx3"),
  effect_labels = c("Cbx3" = "Consistent effect"),  # FIXED label key
  geneset = "Copper homeostasis", goi_exp = goi_exp, limmaRes = limmaRes
)

Pcbp4 <- run_and_extract(
  goi = "Pcbp4", ct = "HSC", kos = c("Chd4"),
  effect_labels = c("Chd4" = "Consistent effect"),  # FIXED label key
  geneset = "Copper homeostasis", goi_exp = goi_exp, limmaRes = limmaRes
)

Myc_GMP <- run_and_extract(
  goi = "Myc", ct = "GMP", kos = c("Brd9"),
  effect_labels = c("Brd9" = "De novo effect"),  # FIXED label key
  geneset = "growth/metabolism", goi_exp = goi_exp, limmaRes = limmaRes
)


goi_exp %>%
  filter(gene == "Myc", genotype == "Brd9")%>%
  pull(celltype)%>%
  unique()
# Stats
stat_results_Ifit1 <- Ifit1$stat
stat_results_Myc_Rcor1 <- Myc_Rcor1$stat
stat_results_Myc_GMP <- Myc_GMP$stat
stat_results_Atp7b <- Atp7b$stat
stat_results_Pcbp4 <- Pcbp4$stat
#stat_results_Myc <- Myc_Rcor1$stat
# Combine all into one data frame
all_stats <- bind_rows(
  stat_results_Ifit1,
  stat_results_Myc_GMP,
  stat_results_Myc_Rcor1,
  stat_results_Atp7b,
  stat_results_Pcbp4
)

# Write to a single CSV
write.csv(all_stats, basedir("all_stats.csv"), row.names = FALSE)
# Plots
Brd9_Ifit1 <- Ifit1_Brd9$plots[["Brd9"]]
Brd9_Myc_GMP <- Myc_GMP$plots[["Brd9"]]
Cbx3_Atp7b <- Atp7b$plots[["Cbx3"]]
Chd4_Pcbp4 <- Pcbp4$plots[["Chd4"]]
#

Fig.2C <- Brd9_Ifit1  + Brd9_Myc_GMP +
  Cbx3_Atp7b + Chd4_Pcbp4 +
  plot_layout(ncol = 4, guides = "collect") &
  theme(
    legend.position = "right"  # removes grid lines
  )
Fig.2C <- Fig.2C +
  plot_annotation(
    title = "Gene expression representing consistent and inconsistent KO-effects between experimental models",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 7,
        face = "bold",
        color = "black"
      )
    )
  )

Fig.2C

ggsave(
  filename = basedir(paste0("Fig.2C",".pdf")),
  plot = Fig.2C,
  width = 18,
  height = 5 ,
  units = "cm"
)

Sup.Fig3C <- Rcor1_Ifit1 + Rcor1_Myc_GMP +
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

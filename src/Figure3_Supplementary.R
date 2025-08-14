#load libraries and functions--------------------------------------------------
#
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
source("src/Ag_ko_classification.R")
#source("src/Ag_enrichR_mouse_genes.R")
library("scales")
library(tidyverse)

basedir <- dirout("Figure3_Supplementary")
InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")

limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
KOs <- limmaRes %>%
  pull(coef)%>%
  unique()
coef_sign <- limmaRes_NTC %>% filter(genes %in% KOs) %>%
  filter(group != "n.s")
coef_logFC <- limmaRes_NTC %>% filter(genes %in% KOs)
ggplot(coef_logFC, aes(
  x = celltype,
  y = genes,
  color = pmin(2, pmax(-2, logFC)),
  size = pmin(3, -log10(adj.P.Val))
)) +  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0, 2),
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
ggsave(basedir("culture_effect_in_Kogenes.pdf"), w = 6, h = 11,  units = "cm")
KO <- koi[1]
ct <- unique(meta$celltype)[1]


dataVoom_NTC_in_ex <- read_rds(InDir_NTC("dataVoom_perCTex.vivovsin.vivo.rds"))
NTC_meta_in_ex <- read_rds(InDir_NTC("NTC_meta.rds"))

example <- c("Cebpa","Brd9","Smc3")
example <- intersect(example, rownames(dataVoom_NTC_in_ex$E))
dat.list<-list()
for(gg in unique(example)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta_in_ex %>%
    mutate(E = dataVoom_NTC_in_ex$E[gg,]) %>%
    rownames_to_column("samples") %>%
    remove_rownames()
  
  dat.list[[gg]] <- gene_data
}
dat.list <- bind_rows(dat.list,.id="gene")
head(dat.list)
#function
unique(dat.list$celltype)
create_gene_plots_NTC <- function(data, geneset, remove_guides = FALSE) {
  
  # Perform Wilcoxon tests by gene & celltype
  stat_table <- data %>%
    group_by(gene, celltype) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(E ~ tissue)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      significance = case_when(
        is.na(p_value) ~ "n.s.",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "n.s."
      )
    ) %>%
    # Set label y-position based on max expression in each facet
    left_join(
      data %>%
        group_by(gene, celltype) %>%
        summarise(y_pos = max(E, na.rm = TRUE) * 1.05, .groups = "drop"),
      by = c("gene", "celltype")
    )
  
  # Create the plot
  plot <- ggplot(data, aes(x = celltype, y = E, color = tissue, group = tissue)) +
    geom_boxplot(
      fill = NA,
      outlier.shape = NA,
      position = position_dodge(width = 0.8),
      size = 0.2
    ) +
    # geom_jitter(
    #   position = position_jitterdodge(
    #     jitter.width = 0.3,
    #     dodge.width = 0.8
    #   ),
    #   alpha = 0.3,
    #   size = 0.5,
    #   show.legend = FALSE
    # ) +
    geom_text(
      data = stat_table,
      aes(x = celltype, y = y_pos, label = significance),
      inherit.aes = FALSE,
      size = 2.5
    ) +
    facet_grid(
      rows = vars(gene),
      cols = vars(celltype),
      space = "free_x",
      scales = "free",
      labeller = labeller(gene = label_wrap_gen(width = 18))
    ) +
    labs(
      title = "Representative gene expression patterns",
      y = "Scaled Gene Expression",
      x = NULL
    ) +
    scale_color_manual(
      values = c("ex.vivo" = "#6a3d9aff", "in.vivo" = "#d38d5fff"),
      name = "Experimental model",
      labels = c("ex.vivo" = "Ex vivo", "in.vivo" = "In vivo")
    ) +
    optimized_theme_fig() +
    theme(
      axis.text.x = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),    # Light grey major gridlines, subtle but visible
  panel.grid.minor = element_blank()
    )
  
  if (remove_guides) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}

# -------------------------------
# Prepare data for plotting
# -------------------------------

# Filter and combine genes
example_genes <- c("Cebpa", "Brd9", "Smc3")
combined_data <- purrr::map_dfr(example_genes, function(gg) {
  df <- dat.list %>%
    filter(gene == gg) %>%
    filter(celltype %in% c("Eo.Ba", "HSC", "MkP", "Mono", "GMP", "Gran.P", "Gran."))
  
  if (nrow(df) == 0) {
    warning(paste("No data available for gene:", gg))
    return(NULL)
  }
  
  df$gene <- gg
  return(df)
})

# Relabel if needed
combined_data$gene <- recode(combined_data$gene,
                             "Cebpa" = "Cebpa",
                             "Brd9" = "Brd9",
                             "Smc3" = "Smc3")

# Set factor levels for consistent ordering
combined_data$celltype <- factor(combined_data$celltype,
                                 levels = c("HSC", "MEP.early", "MkP", 
                                            "GMP", "Gran.P", "Gran.", 
                                            "Mono", "Eo.Ba"))

# Create the plot
Sup.Fig3B <- create_gene_plots_NTC(combined_data, "example_NTC_w.o_jitter", remove_guides = FALSE)

# Save the plot
ggsave(basedir("Sup.Fig3B_w.o_jitter.pdf"), Sup.Fig3B, width = 11, height = 7, units = "cm")

InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")


limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
limmaRes_all <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))
#
ko_expr <- data.table()  # or tibble() or data.frame()

for (coef in selected_KOs) {
  data <- limmaRes %>% 
    
    filter(ensg == coef)
  
  ko_expr <- bind_rows(ko_expr, data)
}

ko_expr %>%
  filter(group != "n.s") %>%
  ggplot(aes(x = ensg, y = -log10(adj.P.Val),
             fill = logFC 
  )) +  # Use alpha based on validity
  geom_col() +  # Use geom_point to create dots
  scale_fill_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(1,2),
    #limits = c(0,5),
    #breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(title = "Interaction effect of ISG core genes",
       x = "KOs",
       y = TeX("$-\\log_{10}(p_{adj})$"))+
  facet_grid(cols = vars(celltype)) +
  theme_bw() +
  optimized_theme_fig()+
  theme(

    legend.position = "right",
    strip.text.x = element_text(angle = 45, hjust = 0, vjust = 0)
  )




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

for (KO in c(selected_KOs)){
  list_of_genes <- selected_KOs
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


#function
stat_tests_all <- list()  # Initialize a list to store results

analyze_ko_across_cts <- function(goi, KO, celltypes, effect_labels, goi_exp, limmaRes, geneset) {
  stat_tests_all <- list()  # Initialize storage for statistical results
  
  # Step 1: Perform statistical tests for the same KO across all cell types
  for (ct in celltypes) {
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
              mutate(tissue = tissue_type, celltype = ct)
          }, error = function(e) {
            message(paste("compare_means failed for", ct, "in", tissue_type, ":", e$message))
            NULL
          })
          
          if (!is.null(stat_test)) {
            stat_tests_all[[paste(ct, tissue_type, sep = "_")]] <- stat_test
          }
        }
      }
    } else {
      message(paste("No data for", goi, "celltype:", ct, "KO:", KO))
    }
  }
  
  # Combine all stats safely
  if (length(stat_tests_all) > 0) {
    stat_tests_combined <- bind_rows(stat_tests_all) %>%
      mutate(significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "ns"
      ))
  } else {
    stat_tests_combined <- tibble()
  }
  
  # Step 2: Generate a single plot across cell types for the same KO
  filtered_data <- goi_exp %>%
    filter(gene == goi, celltype %in% celltypes, comparison == KO)
  
  filtered_data$genotype <- factor(filtered_data$genotype,
                                   levels = c("NTC", setdiff(filtered_data$genotype, "NTC")))
  
  effect_label <- if (!is.null(effect_labels)) effect_labels[KO] else ""
  
  if (nrow(filtered_data) > 0) {
    annotation_data <- filtered_data %>%
      group_by(tissue, celltype) %>%
      summarize(y_pos = max(E, na.rm = TRUE) * 0.8, .groups = "drop")
    
    if (nrow(stat_tests_combined) > 0) {
      annotation_data <- annotation_data %>%
        left_join(stat_tests_combined %>% select(tissue, celltype, significance),
                  by = c("tissue", "celltype"))
    } else {
      annotation_data$significance <- NA
    }
    
    p <- ggplot(filtered_data, aes(x = genotype, y = E, color = tissue)) +
      geom_boxplot(outlier.shape = NA,
                   position = position_dodge(width = 0.8),
                   size = 0.2) +
      geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                  alpha = 0.5) +
      facet_grid(
        cols = vars(celltype, tissue),
        scales = "free_y",
        labeller = labeller(
          celltype = label_value,
          tissue = label_value
        )
      )+
    
      scale_color_manual(
        values = c("ex.vivo" = "#6a3d9aff", "in.vivo" = "#d38d5fff"),
        name = expression("Culture model")
      ) +
      labs(
        title = bquote(atop(.(paste0(goi, ": ", geneset)), .(KO))),
        y = "Expression"
      ) +
      xlab(paste0(KO, " KO", if (effect_label != "") paste0(" (", effect_label, ")"))) +
      theme(legend.position = "none") +
      optimized_theme_fig() +
      geom_text(data = annotation_data,
                aes(x = 1.5, y = y_pos, label = significance),
                inherit.aes = FALSE,
                size = 4)
    
    return(list(
      stat_tests = stat_tests_combined,
      plot = p
    ))
  } else {
    message(paste("No data available for", goi, "KO:", KO))
    return(NULL)
  }
}

ko_list <- c(koi,"Cebpa","Smc3")
results <- map(ko_list, function(ko) {
  result <- analyze_ko_across_cts(
    goi = ko,
    KO = ko,
    celltypes = unique(limmaRes$celltype),
    effect_labels = NULL,
    goi_exp = goi_exp,
    limmaRes = limmaRes,
    geneset = "regulator expression"
  )
  
  # Save plot if it exists
  if (!is.null(result) && !is.null(result$plot)) {
    ggsave(
      filename = basedir(paste0(ko, "_int_expr.pdf")),
      plot = result$plot,
      width = 18,
      height = 6, units = "cm"
    )
  }
  
  return(result)
}) %>% set_names(ko_list)

result_brd9 <- analyze_ko_across_cts(
  goi = "Brd9",
  KO = "Brd9",
  celltypes = unique(limmaRes$celltype),
  effect_labels = NULL,
  goi_exp = goi_exp,
  limmaRes = limmaRes,
  geneset = "regulator expression"
)
ggsave(basedir("Brd9_int_expr.pdf"))

# To access the plot
result_brd9$plot

result_Cebpa <- analyze_ko_across_cts(
  goi = "Cebpa",
  KO = "Cebpa",
  celltypes = unique(limmaRes$celltype),
  effect_labels = NULL,
  goi_exp = goi_exp,
  limmaRes = limmaRes,
  geneset = "regulator expression"
)
ggsave(basedir("Cebpa_int_expr.pdf"))
result_Smc3 <- analyze_ko_across_cts(
  goi = "Smc3",
  KO = "Smc3",
  celltypes = unique(limmaRes$celltype),
  effect_labels = NULL,
  goi_exp = goi_exp,
  limmaRes = limmaRes,
  geneset = "regulator expression"
)

result_Prmt1 <- analyze_ko_across_cts(
  goi = "Prmt1",
  KO = "Prmt1",
  celltypes = unique(limmaRes$celltype),
  effect_labels = NULL,
  goi_exp = goi_exp,
  limmaRes = limmaRes,
  geneset = "regulator expression"
)

result_Setdb1 <- analyze_ko_across_cts(
  goi = "Setdb1",
  KO = "Setdb1",
  celltypes = unique(limmaRes$celltype),
  effect_labels = NULL,
  goi_exp = goi_exp,
  limmaRes = limmaRes,
  geneset = "regulator expression"
)

result_Smc3 <- analyze_ko_across_cts(
  goi = "Smc3",
  KO = "Smc3",
  celltypes = unique(limmaRes$celltype),
  effect_labels = NULL,
  goi_exp = goi_exp,
  limmaRes = limmaRes,
  geneset = "regulator expression"
)

#
create_gene_plots_NTC <- function(data, geneset, remove_guides = FALSE) {
  
  # Helper function to check normality and decide test for each gene-celltype
  data <- data %>%
    group_by(gene, celltype) %>%
    mutate(
      test_method = {
        # Split data by tissue group for this gene-celltype
        groups <- split(E, tissue)
        
        # Run Shapiro-Wilk test for normality on each group if >= 3 values
        normal1 <- if (length(groups[[1]]) >= 3) shapiro.test(groups[[1]])$p.value > 0.05 else FALSE
        normal2 <- if (length(groups[[2]]) >= 3) shapiro.test(groups[[2]])$p.value > 0.05 else FALSE
        
        # If both groups normal, use t-test, else wilcox.test
        if (normal1 & normal2) "t.test" else "wilcox.test"
      }
    ) %>%
    ungroup()
  
  # Extract unique test methods per gene-celltype for later use in plotting
  stat_methods <- data %>%
    distinct(gene, celltype, test_method)
  
  # Create the plot
  plot <- ggplot(data, aes(x = celltype, y = E, color = tissue, group = tissue)) + 
    geom_boxplot(
      fill = NA,  # Hollow boxplots
      outlier.shape = NA,
      position = position_dodge(width = 0.8),  # Match jitter
      size = 0.2
    ) + 
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.3,
        dodge.width = 0.8
      ), 
      alpha = 0.3,
      size = 0.5,
      show.legend = FALSE
    ) +
    facet_grid(rows = vars(gene), cols = vars(celltype), space = "free_x",
               scales = "free",
               labeller = labeller(gene = label_wrap_gen(width = 18))) +  
    labs(
      legend = "Experimental model",
      title = "Representative gene expression patterns"
    ) +
    xlab(NULL) +
    ylab("Scaled Gene Expression") +
    scale_color_manual(
      values = c("ex.vivo" = "#6a3d9aff", "in.vivo" = "#d38d5fff"),
      name = "Experimental model",
      labels = c("ex.vivo" = "Ex vivo", "in.vivo" = "In vivo")
    ) +
    optimized_theme_fig() +
    theme(
      axis.text.x = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank()
    )
  
  if (remove_guides) {
    plot <- plot + theme(legend.position = "none")
  }
  
  # Add significance using the right test method per facet
  for(i in seq_len(nrow(stat_methods))) {
    this_gene <- stat_methods$gene[i]
    this_celltype <- stat_methods$celltype[i]
    this_method <- stat_methods$test_method[i]
    
    plot <- plot + stat_compare_means(
      data = data %>% filter(gene == this_gene, celltype == this_celltype),
      aes(x = celltype, y = E, group = tissue),
      method = this_method,
      label = "p.signif",
      label.y = max(data$E[data$gene == this_gene & data$celltype == this_celltype], na.rm = TRUE) * 0.5,
      inherit.aes = FALSE
    )
  }
  
  return(plot)
}
# Example genes
example <- c("Cebpa", "Brd9", "Prmt5")

# Prepare combined_data filtered for your genes and celltypes
combined_data <- data.frame()
for (gg in example) {
  data <- dat.list %>% 
    filter(gene == gg) %>%
    filter(celltype %in% c("Eo.Ba", "HSC", "MkP", "Mono", "GMP", "Gran.P", "Gran."))
  
  if (nrow(data) == 0) {
    warning(paste("No data available for gene:", gg))
    next
  }
  
  data$gene <- gg
  combined_data <- rbind(combined_data, data)
}

# Apply custom gene labels if needed
custom_labels <- c(
  "Cebpa" = "Cebpa",
  "Brd9" = "Brd9",
  "Prmt5" = "Prmt5"
)
combined_data <- combined_data %>%
  mutate(gene = dplyr::recode(gene, !!!custom_labels))

# Set factor levels for consistent celltype ordering
combined_data$celltype <- factor(combined_data$celltype,
                                 levels = c("HSC", "MEP.early", "MkP", 
                                            "GMP", "Gran.P", "Gran.", 
                                            "Mono", "Eo.Ba"))

#Sup.Fig.3C--------------------------
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
              panel.grid.minor = element_blank())
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

Ifit1_Rcor1 <- run_and_extract(
  goi = "Ifit1", ct = "Eo.Ba", kos = c("Rcor1"),
  effect_labels = c("Rcor1" = "No effect"),
  geneset = "ISG", goi_exp = goi_exp, limmaRes = limmaRes
)

Myc_Rcor1 <- run_and_extract(
  goi = "Myc", ct = "GMP", kos = c("Rcor1"),
  effect_labels = c("Rcor1" = "No effect"),  # FIXED label key
  geneset = "growth/metabolism", goi_exp = goi_exp, limmaRes = limmaRes
)

Rcor1_Myc_GMP <- Myc_Rcor1$plots[["Rcor1"]]
Rcor1_Ifit1 <- Ifit1_Rcor1$plots[["Rcor1"]]
# Stats
stat_results_Myc_Rcor1 <- Myc_Rcor1$stat
stat_results_Myc_GMP <- Myc_GMP$stat
#stat_results_Myc <- Myc_Rcor1$stat
# Combine all into one data frame
all_stats <- bind_rows(

  stat_results_Myc_GMP,
  stat_results_Myc_Rcor1

)

# Write to a single CSV
write.csv(all_stats, basedir("all_stats.csv"), row.names = FALSE)
# Plots



#
Atp7b_Brd9 <- run_and_extract(
  goi = "Atp7b", ct = "Mono", kos = c("Brd9"),
  effect_labels = c("Brd9" = "No effect"),  # FIXED label key
  geneset = "Copper homeostasis", goi_exp = goi_exp, limmaRes = limmaRes
)
Brd9_Atp7b <- Atp7b_Brd9$plots[["Brd9"]]

Pcbp4_Brd9 <- run_and_extract(
  goi = "Pcbp4", ct = "HSC", kos = c("Brd9"),
  effect_labels = c("Brd9" = "No effect"),  # FIXED label key
  geneset = "Copper homeostasis", goi_exp = goi_exp, limmaRes = limmaRes
)
Brd9_Pcbp4 <- Pcbp4_Brd9$plots[["Brd9"]]
 
Sup.Fig3C <- Rcor1_Ifit1 + Rcor1_Myc_GMP +
  Brd9_Pcbp4 + Brd9_Atp7b +
  plot_layout(ncol = 4, guides = "collect") +
  theme(
    legend.position = "right"  # removes grid lines
  )
ggsave(
  filename = basedir(paste0("Sup.Fig3C_w.o.jitter",".pdf")),
  plot = Sup.Fig3C,
  width = 18,
  height = 5 ,
  units = "cm"
)


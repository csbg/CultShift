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
ggsave(basedir("Sup.Fig3A.pdf"), w = 6, h = 11,  units = "cm")
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
#func
create_gene_plots_NTC <- function(data, geneset, remove_guides = FALSE) {
  plot <- ggplot(data, aes(x = celltype, y = E, color = tissue, group = tissue)) + 
    geom_boxplot(
      fill = NA,
      outlier.shape = NA,
      position = position_dodge(width = 0.8),
      size = 0.2
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
      axis.ticks.x = NULL,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank()
    ) +
    # Add significance stars from precomputed limma results
    geom_text(
      data = data %>% distinct(gene, celltype, y_pos, sig_label), # one label per facet
      aes(x = 1, y = y_pos, label = sig_label),
      inherit.aes = FALSE,
      size = 2,               # thinner stars
      fontface = "plain",     # not bold
      #family = "Arial",       # optional: slimmer font
      alpha = 0.6
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




# Initialize empty list to store individual plots
all_gene_plots <- list()


# adjust as needed to match your limmaRes_NTC structure
sig_thresholds <- function(p) {
  ifelse(p < 0.0001, "****",  
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "n.s"))))
}

# Merge limma results into combined_data
combined_data <- combined_data %>%
  left_join(limmaRes_NTC %>% 
              dplyr::rename(gene=genes)%>%
              dplyr::select(gene, celltype, adj.P.Val), 
            by = c("gene", "celltype")) %>%
  mutate(sig_label = sig_thresholds(adj.P.Val),
         y_pos = max(E, na.rm = TRUE) * 1.1)   # position stars above boxplots

# Create the plot
Sup.Fig3B <- create_gene_plots_NTC(combined_data, "example_NTC_w.o_jitter", remove_guides = FALSE)

# Save the plot
ggsave(basedir("Sup.Fig3B_w.o_jitter.pdf"), Sup.Fig3B, width = 11, height = 7, units = "cm")




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
non_affected <- c("Chd4","Prmt5")
for (KO in c(selected_KOs,non_affected)){
  list_of_genes <- c("Oas2",
                     "Cebpa",
                     "Cpt1",
                     "Apoa1",
                     "Cd36",
                     "Rab44",
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
limmaRes_all <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))
#consistent

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
unique(goi_exp_limma$celltype)
Apoa1_Brd9 <- run_and_extract(
  goi = "Apoa1",
  ct = "MkP",
  kos = c("Brd9"),
  effect_labels = c("Brd9" = "Consistent trend"),
  geneset = "fatty acid metabolism",
  goi_exp_limma = goi_exp_limma
)
Cd36_Brd9 <- run_and_extract(
  goi = "Cd36",
  ct = "Mono",
  kos = c("Brd9"),
  effect_labels = c("Brd9" = "Consistent trend"),
  geneset = "fatty acid metabolism",
  goi_exp_limma = goi_exp_limma
)
Cpt1_Brd9 <- run_and_extract(
  goi = "Cpt1",
  ct = "MkP",
  kos = c("Brd9"),
  effect_labels = c("Brd9" = "Consistent trend"),
  geneset = "fatty acid metabolism",
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

Cebpa_Brd9 <- run_and_extract(
  goi = "Cebpa", ct = "GMP", kos = c("Brd9"),
  effect_labels = c("Brd9" = "Opposite trend"),
  geneset = "Myeloid differentiation",
  goi_exp_limma = goi_exp_limma
)

Rab44_Brd9 <- run_and_extract(
  goi = "Rab44", ct = "Eo.Ba", kos = c("Brd9"),
  effect_labels = c("Brd9" = "No effect"),
  geneset = "Rab GTPase",
  goi_exp_limma = goi_exp_limma
)




# Example: Combine first KO plots into a multi-panel figure
Sup.Fig.3C <- Cebpa_Brd9$plots[[1]]+
  Rab44_Brd9$plots[[1]]+
  plot_layout(ncol=, guides="collect") &
  theme(legend.position="right")





ggsave(
  filename = basedir(paste0("Sup.Fig3C_w.o.jitter",".pdf")),
  plot = Sup.Fig.3C,
  width = 7.5,
  height = 4 ,
  units = "cm"
)


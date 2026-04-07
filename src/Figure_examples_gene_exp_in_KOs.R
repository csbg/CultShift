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


basedir <- dirout("Figure_examples_gene_exp_in_KOs")
#Fig3D alternate

koi_t <- as.data.frame(koi,"koi")
write.table(koi_t,basedir("genotypes.tsv"), row.names = F)

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
Sup.Fig3D <- Ifit1_Rcor1$plots[[1]] + Rcor1_Myc_GMP$plots[[1]] +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "right"  # removes grid lines
  )
ggsave(
  filename = basedir(paste0("Sup.Fig3D",".pdf")),
  plot = Sup.Fig3D,
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


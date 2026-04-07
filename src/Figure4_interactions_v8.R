source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
library(tidyverse)
library(enrichR)
library(purrr)
library(scales)
library(purrr)
library(patchwork)
library(cowplot)
library(latex2exp)
library(readr)

#directories ------
base <-"Figure4_interactions_v8"
basedir <- dirout("Figure4_interactions_v8")

#InDir6 <- dirout("Ag_ScRNA_12_fgsea_overlap")


#Fig4B-D------------

#HSC----
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
# Step 1: Filter for significant genes (adj.P.Val < 0.05 and abs(logFC) > 1)
limmaRes_significant <- limmaRes %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)  # Only significantly altered genes


data <- summary_df %>%
  filter(coef %in% koi) %>%
  filter(Count >= 10) %>%
  inner_join(ko_flags, by = c("celltype", "coef")) %>%
  filter(valid_ko)

data %>% write_rds(basedir("kos_Fig4.rds"))
#for text
summary_total <- data %>%
  group_by(celltype, coef, genotype) %>%
  summarise(Total_Regulated = sum(Count, na.rm = TRUE), .groups = "drop")%>%
  filter(Total_Regulated >= 10)
  

median(summary_total$Total_Regulated)
summary_total %>%
  mutate(celltype_coef = paste0(celltype,coef))%>%
  pull(celltype_coef)%>%
  unique()
length(unique(summary_total$coef))
summary_total$data <-"HSC"

#Fig4C#Spleen-----------------------------

out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"

InDir1 <- dirout("Figure_JAK_STAT")

limmaRes <- read_rds(InDir1("combined_jakstat_diff_exp.rds"))

limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))
# Modify the 'coef' column for any KO
limmaRes <- limmaRes %>%
  mutate(coef = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(cell_type = gsub("M", "Macrophage", cell_type)) %>%
  mutate(cell_type = gsub("T8", "T-cells", cell_type)) 
limmaRes$coef <- gsub("treatmentex_vivo","WT",
                      limmaRes$coef)

DEG_spleen <- limmaRes %>%
  dplyr::filter(grepl("Interaction", coef)) %>%
  dplyr::mutate(genotype = sub("Interaction_", "", coef)) %>%
  filter(genotype %in% c("IRF9KO","STAT1KO","STAT2KO","TYK2CMV",
                         "TYK2KE","STAT3VAV", "STAT1BB"))%>%
  mutate(genotype = ifelse(grepl("IRF9KO",genotype),"Irf9KO",ifelse(
    grepl("STAT1KO",genotype),"Stat1KO",ifelse(
      grepl("STAT2KO",genotype),"StatKO",ifelse(
        grepl("TYK2CMV",genotype),"Tyk2CMV",ifelse(
          grepl("TYK2KE",genotype),"Tyk2KE",ifelse(
            grepl("STAT3VAV",genotype),"Stat3Vav","Stat1BB"
          )))))))%>%
  dplyr::filter(cell_type == "T-cells") %>%
  dplyr::mutate(celltype = cell_type) %>%
  dplyr::select(-cell_type) %>%
  dplyr::filter(abs(logFC) > 1, adj.P.Val < 0.01) %>%
  
  
  dplyr::group_by(coef, celltype, genotype) %>%
  dplyr::summarise(Total_Regulated = dplyr::n(), .groups = "drop") %>%  # 👈 avoid tally
  
  dplyr::filter(Total_Regulated >= 10) %>%
  dplyr::mutate(data = "Splenic T-cells") %>%
  dplyr::select(celltype, coef, genotype, Total_Regulated, data)
#summaries for text---
DEG_spleen %>%
  mutate(celltype_coef = paste0(celltype,coef))%>%
  pull(celltype_coef)%>%
  unique()



#############
#Fig4D---------------
#Glioblastoma---------------------
InDir2 <- dirout("Glioblastoma_limmaRes.3way.mod")

limmaRes <- read_rds(InDir2("limmaRes_threeway.rds"))
unique(limmaRes$coef)

limmaRes_NTC_noRT <- limmaRes%>%
  filter(coef == "tissueex.vivo")%>%
  dplyr::select(-coef)
limmaRes_int_glio <- limmaRes%>% 
  filter(grepl("tissueex.vivo.genotype[A-Za-z0-9]*$",coef))%>%
  mutate(coef = gsub("tissueex.vivo.genotype","", coef))

limmaRes_clean <- limmaRes %>%
  filter(
    coef %in% c(
      grep("tissueex.vivo.genotype[A-Za-z0-9]*$", coef, value = TRUE),
      "tissueex.vivo",
      "tissueex.vivo.RT_statusRT"
    )
  ) %>%
  mutate(
    coef = case_when(
      coef == "tissueex.vivo" ~ "culture_effect",
      coef == "tissueex.vivo.RT_statusRT" ~ "interaction_RT",
      str_detect(coef, "^tissueex\\.vivo\\.genotype") ~ str_replace(coef, "^tissueex\\.vivo\\.genotype", "Interaction_genotype_"),
      
      TRUE ~ coef
    )
  )


deg_counts_glio <- limmaRes_int_glio %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(coef) %>%
  tally(name = "Total_Regulated") %>%
  mutate(genotype = coef)%>%
  mutate(data = "Glioblastoma")%>%
  mutate(celltype = "GL261")%>%
  filter(Total_Regulated  >= 10)
# summaries for text-------
deg_counts_glio %>%
  mutate(celltype_coef = paste0(celltype,coef))%>%
  pull(celltype_coef)%>%
  unique()
#combined-----
DEG_table <- dplyr::bind_rows(
  summary_total,
  DEG_spleen,
  deg_counts_glio
)
DEG_table <- DEG_table %>%
  mutate(
    genotype = as.character(genotype),
    celltype = as.character(celltype),
    data = as.character(data)
  )
DEG_table <- DEG_table %>%
  mutate(
    facet_id = paste(celltype, data, sep = " | ")
  )

# Split by facet
DEG_list <- DEG_table %>% 
  mutate(data = factor(data))%>%
  group_split(facet_id)

# Count genotypes per facet (for patchwork widths)
facet_counts <- DEG_table %>% group_by(facet_id) %>% summarise(n = n())

# Make one plot per facet
plot_list <- lapply(DEG_list, function(df) {
  ggplot(df %>%
           mutate(data = factor(data,
                                levels = c("HSC","Glioblastoma","Splenic"))), 
         aes(x = reorder(genotype,Total_Regulated), y = log10(Total_Regulated))) +
    geom_col(color = "darkgrey", width = 0.35, fill = NA) +  # now fixed width
    labs(title = unique(df$facet_id), x = NULL, y = NULL) +
    optimized_theme_fig() +
    theme(axis.text.x = element_text(angle = 70, hjust = 0, vjust = 0))+
    optimized_theme_fig()+
    theme(
      plot.margin = margin(0, 0, 0, 0),       # remove outer margins
      panel.spacing = unit(0.0, "cm") ,       # reduce spacing between facets
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
        )
})

# Combine with variable widths proportional to number of genotypes
combined_plot <- wrap_plots(plot_list, nrow = 1, widths = facet_counts$n)
#plot-----------
combined_plot

ggsave(basedir("4B_DHSC_Glio_Spleen_N_DEGs_Interaction_all.pdf"), w=22,
       h = 3, units = "cm")

#############

#4Eexample genes----------------
InDir2 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")


limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
limmaRes_all <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_all_coef.rds"))
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
  
)

results_list <- lapply(genes_to_run, function(x) {
  analyze_kos(x$goi, x$ct, x$kos, x$effect_labels, goi_exp_limma, x$geneset)
})


# Example: Combine first KO plots into a multi-panel figure
Fig.4E <- results_list[[1]]$plots[[1]] + 
  results_list[[7]]$plots[[1]] + 
  results_list[[3]]$plots[[1]] + 
  results_list[[2]]$plots[[1]] + 
  plot_layout(ncol=4, guides="collect") &
  theme(legend.position="right")

Fig.4E <- Fig.4E + 
  plot_annotation(
    title = str_wrap("Gene expression representing consistent and inconsistent KO-effects between experimental models", 80)
  ) &
  optimized_theme_fig(
  )+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
Fig.4E
ggsave(
  filename = basedir(paste0("Fig4E",".pdf")),
  plot = Fig.4E,
  width = 18,
  height = 5 ,
  units = "cm"
)#

#4F---------
#Correlation of culture effect and KO effect
limmaRes_int <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))%>%
  mutate(genes = ensg)
#########
limmaRes_int %>%
  group_by(celltype,coef) %>%
  summar
#########
limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

correlation_results_perturb <- merged_data %>%
  inner_join(ko_flags, by = c("coef","celltype")) %>%
  filter(valid_ko) %>%
  group_by(coef, celltype) %>%
  filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )%>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

correlation_results_summarized <- correlation_results_perturb %>%
  filter(coef %in% koi) %>%                              # filter for your koi
  inner_join(summary_df, by = c("coef","celltype")) %>%  # join
  group_by(coef, celltype) %>%                           # group by the grouping columns
  summarise(num_degs = sum(Count, na.rm = TRUE),      # sum Count
            .groups = "drop")                            # ungroup after summarise


# Step 3: Add significance labels based on p-value thresholds
correlation_results_perturb <- correlation_results_perturb %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )

correlation_results_perturb$sign_level <- factor(
  correlation_results_perturb$significance,
  levels = c("***","**","*","")
)

Fig4F_heatmap <- ggplot(correlation_results_perturb,
                         aes(x = coef, y = celltype, fill = cor_abs,
                             size = pmin(-log10(p_adj),5))) +
  
  geom_point(shape = 21) +
  #geom_text(aes(label = significance), size = 1.5) +
  
  scale_fill_gradient2(
    low = "#9A9C39",    # muted blue-teal
    mid = "white",    # soft yellow
    high = "#E7298A",   # muted red-orange
    limits = c(-1, 1),
    name = str_wrap("Correlation", width = 20)
  )+
  
  scale_size_continuous(
    range = c(0,2),
    name = TeX("$\\log_{10}\\; (Padj)$")
  )+
  labs(title = str_wrap("Culture effects drive discrepancies of KO-effects in ex vivo models", 
                        width = 80),
       x = "KO",
       y = "Cell type") +
  
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
Fig4F_heatmap
ggsave(basedir("Fig4F_heatmap_perurb_correlation_cult_int.pdf"),plot=Fig4F_heatmap,
       w=9,h = 5, units = "cm")


######
#Fig4C--------density-----------

#Spleen-------------
out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"
#*
InDir1 <- dirout("Figure_JAK_STAT")

limmaRes <- read_rds(InDir1("combined_jakstat_diff_exp.rds"))

limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))

# Modify the 'coef' column for any KO
limmaRes <- limmaRes %>%
  mutate(coef = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(cell_type = gsub("M", "Macrophage", cell_type)) %>%
  mutate(cell_type = gsub("T8", "T-cells", cell_type)) 
limmaRes$coef <- gsub("treatmentex_vivo","WT",
                      limmaRes$coef)

limmaRes_int <- limmaRes %>%
  filter(grepl("Interaction", limmaRes$coef))%>%
  filter(cell_type == "T-cells")
deg_counts <- limmaRes_int %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.01) %>%
  group_by(coef, cell_type) %>%
  tally(name = "n_DEGs") %>%
  mutate(valid = n_DEGs >= 10)

limmaRes_WT <- limmaRes %>%
  filter(coef == "WT")%>%
  filter(cell_type == "T-cells")%>%
  dplyr::select(-coef)
merged_data <- limmaRes_int %>%
  inner_join(deg_counts %>% filter(valid), 
             by = c("coef", "cell_type"))%>%
  inner_join(limmaRes_WT, by = c("gene","cell_type"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)%>%
  mutate(celltype = cell_type)


correlation_results_spleen <- merged_data %>%
  group_by(coef, celltype) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results_spleen <- correlation_results_spleen %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Step 3: Add significance labels based on p-value thresholds
correlation_results_spleen <- correlation_results_spleen %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )

Sup.Fig.9_heatmap_spleen <- ggplot(correlation_results_spleen,
                                aes(x = gsub("Interaction_","",coef), 
                                    y = celltype,
                                    fill = cor_abs,
                                    size = pmin(-log10(p_adj),5))) +
  
  geom_point(shape = 21) +
  #geom_text(aes(label = significance), size = 1.5) +
  
  scale_fill_gradient2(
    low = "#9A9C39",    # muted blue-teal
    mid = "white",    # soft yellow
    high = "#E7298A",   # muted red-orange
    limits = c(-1, 1),
    name = str_wrap("Correlation", width = 20)
  )+
  
  scale_size_continuous(
    range = c(0,2),
    limits = c(0,5),
    name = TeX("$\\log_{10}\\; (Padj)$")
  )+
  labs(title = str_wrap("Culture effects drive discrepancies of KO-effects in ex vivo models", 
                        width = 50),
       x = "KO",
       y = "Cell type") +
  
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(basedir("Sup.Fig.9_NTC_inter_cor.pdf"), w = 3, h = 4,units = "cm")
#
#Glioblastoma---------------------
InDir2 <- dirout("Glioblastoma_limmaRes.3way.mod")

limmaRes <- read_rds(InDir2("limmaRes_threeway.rds"))

limmaRes_NTC_noRT <- limmaRes%>%
  filter(coef == "tissueex.vivo")%>%
  dplyr::select(-coef)
limmaRes_int_glio <- limmaRes%>% 
  filter(grepl("tissueex.vivo.genotype[A-Za-z0-9]*$",coef))%>%
  mutate(coef = gsub("tissueex.vivo.genotype","", coef))


# Extract only in.vivo_noRT and ex.vivo_noRT genotype coefficients
genotype_noRT <- limmaRes %>%
  filter(str_detect(coef, "_noRT$")) %>%                  # keep only noRT rows
  filter(str_detect(coef, "in\\.vivo|ex\\.vivo")) %>%    # keep in.vivo or ex.vivo
  filter(!str_detect(coef, "tissueex\\.vivo|RT_status")) # remove tissueex or other unwanted coeffs

deg_counts <- limmaRes_int_glio %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.01) %>%
  group_by(coef) %>%
  tally(name = "n_DEGs") %>%
  mutate(valid = n_DEGs >= 10) 
merged_data <- limmaRes_int_glio %>%
  inner_join(deg_counts %>% filter(valid), 
             by = c("coef"))%>%
  inner_join(limmaRes_NTC_noRT, by = c("ensg"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)%>%
  mutate(celltype = "glio")

correlation_results_glio <- merged_data %>%
  group_by(coef, celltype) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results_glio <- correlation_results_glio %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Step 3: Add significance labels based on p-value thresholds
correlation_results_glio <- correlation_results_glio %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )
Sup.Fig.9_heatmap_glio <- ggplot(correlation_results_glio,
                              aes(x = gsub("Interaction","",coef), 
                                  y = celltype,
                                  fill = cor_abs,
                                  size = pmin(-log10(p_adj),5))) +
  
  geom_point(shape = 21) +
  #geom_text(aes(label = significance), size = 1.5) +
  
  scale_fill_gradient2(
    low = "#9A9C39",    # muted blue-teal
    mid = "white",    # soft yellow
    high = "#E7298A",   # muted red-orange
    limits = c(-1, 1),
    name = str_wrap("Correlation", width = 20)
  )+
  scale_size_continuous(
    range = c(0,2),
    limits = c(0,5),
    name = TeX("$\\log_{10}\\; (Padj)$")
  )+
  labs(title = str_wrap("Culture effects drive discrepancies of KO-effects in ex vivo models", 
                        width = 70),
       x = "KO",
       y = "Cell type") +
  
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
ggsave(basedir("Sup.Fig.9_heatmap_glio.pdf"),plot=Sup.Fig.9_heatmap_glio,
       w=9.5,h=3.5, units = "cm")
#Combined plot --------------------
correlation_results_glio$dataset <- "Glioblastoma"

correlation_results_perturb$dataset <- str_wrap("Perturb-seq(Hematiopoietic)", width = 11)
correlation_results_perturb$sign_level <- NULL
correlation_all <- rbind(correlation_results_perturb,
                         correlation_results_spleen, 
                         correlation_results_glio)

condition_colors <- c("#4E92D0" ,
                      "#478E5C",
                      "#AD1F69"
)
ggplot(correlation_all, aes(x = cor_abs, color = dataset)) +
  geom_density(alpha = 0.4, size = 0.2) +        # semi-transparent fill
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype ="dotted", color = "red", size =0.2)+
  #scale_fill_manual(values = condition_colors)+
  scale_color_manual(values = condition_colors)+
  labs(
    title = str_wrap("Correlation between culture effects and interaction effects",
                     width = 40),
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
                                panel.grid.minor = element_blank(),
                                axis.text.x = element_text(angle = 0))
ggsave(basedir("4G_Correlation_density_interaction.NTC.only_Tcells.pdf"), w = 5.9, h = 3.5, units = "cm")
##############
#Fig4H---------correlation Observed vs shuffled-----------


InDir7  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_per_pathway_fgsea_in.vivo")


limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))%>%
  mutate(genes = ensg)
limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
merged_data <- limmaRes %>%
  inner_join(limmaRes_NTC, by = c("genes","celltype"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)


# Step 1: Calculate observed correlations with p-values
observed_correlations <- merged_data %>%
  filter(coef %in% koi) %>%
  merge(ko_flags, by = c("coef", "celltype")) %>%
  filter(valid_ko)%>%
  group_by(coef, celltype) %>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,
    .groups = 'drop'
  ) %>%
  mutate(type = "Observed")  # Label as observed correlations

# Step 2: Shuffle values and calculate correlations again
set.seed(42)  # Set seed for reproducibility
shuffled_correlations <- merged_data %>%
  merge(ko_flags, by = c("coef", "celltype")) %>%
  filter(valid_ko)%>%
  filter(coef %in% koi) %>%
  group_by(coef, celltype) %>%
  mutate(shuffled_logFC_KO = sample(logFC_KO)) %>%  # Shuffle logFC_KO within each group
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(shuffled_logFC_KO), method = "pearson"),
    .groups = 'drop'
  ) %>%
  mutate(type = "Shuffled")  # Label as shuffled correlations

# Step 3: Combine observed and shuffled correlations
combined_correlations <- bind_rows(
  observed_correlations %>% dplyr::select(coef, celltype, cor_abs, type),
  shuffled_correlations %>% dplyr::select(coef, celltype, cor_abs, type)
)

Fig4H <- ggplot(combined_correlations, aes(x = cor_abs, y = celltype, fill = type)) +
  geom_density_ridges(alpha = 0.4, scale = 1) +
  labs(
    title = "Correlation of culture effects to observed vs 
shuffled interaction effects across KOs",
    x = "Correlation",
    y = "Cell type",
    fill = "Correlation Type"
  ) +
  scale_x_continuous(breaks = c(0,0.5, 0.8, 1)) + 
  # xlim(-1,2)+
  theme_minimal() +
  scale_fill_manual(values = c("Observed" = "blue", "Shuffled" = "gray"))+
  optimized_theme_fig()


ggsave(basedir("Fig4H_Correlation.pdf"),plot = Fig4H,
       w=7,h=5, units = "cm")



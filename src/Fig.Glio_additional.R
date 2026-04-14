##############################
# Modular GEO RNA-seq pipeline
##############################
#Aarathy
###############
source("src/00_init.R")
#source("src/Ag_enrichR_mouse_genes.R")
source("src/Ag_ko_classification.R")
source("src/Ag_Optimized_theme_fig.R")
library(edgeR)
library(limma)
library(biomaRt)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
library(GEOquery)
library(R.utils)
library(msigdbr)
library(fgsea)
library(readxl)
library(latex2exp)
library(Biobase)
library(Matrix)
library(circlize)
library(grid)

InDir1 <- dirout("Ag_top_pathway_genes")
InDir <- dirout("Glioblastoma_limmaRes.3way.mod")
out <- dirout("Fig.Glio")
limmaRes <- read_rds(InDir("limmaRes_threeway.rds"))

threeway_terms <- unique(grep("^tissue.*.genotype.*.RT_status", limmaRes$coef, value = TRUE))

triple_interaction_summary <- limmaRes %>%
  filter(coef %in% threeway_terms) %>%       # keep only three-way interactions
  group_by(coef) %>%                         # group by interaction term
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")  # count up/down

unique(limmaRes$coef)

limmaRes_NTC <- limmaRes %>%
  filter(coef == "tissueex.vivo")

###########
#enrichment to ntc----KO interaction

limmaRes_invivo <- limmaRes %>% filter(coef %in% grep("[A-Za-z0-9]*_in.vivo_noRT",limmaRes$coef, value = T))

limmaRes_exvivo <- limmaRes %>% filter(coef %in% grep("[A-Za-z0-9]*_ex.vivo_noRT",limmaRes$coef, value = T))
#num_degs------

# Apply to each dataset
deg_invivo <- count_deg(limmaRes_invivo)
deg_exvivo <- count_deg(limmaRes_exvivo)

deg_combined <- full_join(deg_exvivo, deg_invivo, by = "coef")

#correlation in vivo vs ex vivo


limmaRes_invivo <- limmaRes_invivo %>%
  mutate(KO = str_replace(coef, "_in\\.vivo_noRT", ""))

limmaRes_exvivo <- limmaRes_exvivo %>%
  mutate(KO = str_replace(coef, "_ex\\.vivo_noRT", ""))


count_deg <- function(df, padj_cutoff = 0.05, lfc_cutoff = 1) {
  df %>%
    group_by(KO) %>%   # <<-- changed
    summarize(
      n_sig = sum(adj.P.Val < padj_cutoff & abs(logFC) > lfc_cutoff, na.rm = TRUE),
      n_up  = sum(adj.P.Val < padj_cutoff & logFC >  lfc_cutoff, na.rm = TRUE),
      n_down = sum(adj.P.Val < padj_cutoff & logFC < -lfc_cutoff, na.rm = TRUE)
    ) %>%
    arrange(KO)
}

deg_exvivo <- count_deg(limmaRes_exvivo) %>% rename_with(~paste0(.x, "_exvivo"), -KO)
deg_invivo <- count_deg(limmaRes_invivo) %>% rename_with(~paste0(.x, "_invivo"), -KO)

deg_combined <- full_join(deg_exvivo, deg_invivo, by = "KO") %>%
  arrange(KO)

# set rownames
deg_combined_df <- as.data.frame(deg_combined)
rownames(deg_combined_df) <- deg_combined_df$KO
deg_combined_df$KO <- NULL



merged <- limmaRes_invivo %>%
  dplyr::select(KO, ensg, logFC_in = logFC) %>%
  inner_join(
    limmaRes_exvivo %>% dplyr::select(KO, ensg, logFC_ex = logFC),
    by = c("KO", "ensg")
  )

cor_results <- merged %>%
  group_by(KO) %>%
  summarize(
    cor = cor(logFC_in, logFC_ex, method = "pearson"),
    n = n()
  )
KO_order <- merged %>%
  group_by(KO) %>%
  summarize(
    cor = cor(logFC_in, logFC_ex, method = "pearson"),
    n = n()
  ) %>%
  arrange(cor) %>%
  pull(KO)

combined <- merge(deg_combined[,c("KO","n_sig_exvivo","n_sig_invivo")],cor_results, by = "KO")
combined <- combined %>%
  mutate(n_sig_max = pmax(n_sig_exvivo, n_sig_invivo))
combined$KO <- factor(combined$KO, levels = KO_order)

combined <- combined %>%
  arrange(KO)


# Make KO a factor with your desired order if not already
combined$KO <- factor(combined$KO, levels = KO_order)



# Reference gene sets: significant genes per cell type
ref_gene_sets <- limmaRes_NTC %>%
  filter(group != "n.s") %>%
  summarise(ref_genes = list(unique(ensg)), .groups = "drop")

input_gene_sets <- limmaRes_invivo %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  
  group_by(coef) %>%
  summarise(sig_genes = list(unique(ensg)), .groups = "drop")


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
enrichment_ntc_in.vivo <- enrichment_results%>%
  mutate(KO = gsub("_in.vivo_noRT","",coef))

combined <- combined %>%
  distinct()%>%
  left_join(enrichment_ntc_in.vivo, by = c("KO"))

# Step 1: Combine the correlation summary and enrichment summary
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
  mutate(
    log2.or.filtered = ifelse(overlap >= 5, pmin(log2.odds.ratio, 4), NA)  # keep value only if overlap>5
  )

combined <- combined %>%
  mutate(log2.or.filtered.plot = ifelse(is.na(log2.or.filtered), 0.01, log2.or.filtered))
combined <- combined %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))  # FDR correction


# Make KO a factor in the order of signif_count or custom order
combined$KO <- factor(combined$KO, levels = KO_order)
#
combined %>% write_rds(out("overlap_enr_corr_inex.rds"))
# Create bubble plot
Fig_overlap_bubble <- ggplot(combined%>%
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
heatmap <- ggplot(combined%>%
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

ggsave(out("Overlap_enrichment_mean_cor_Glio.pdf"),
       plot = final_plot,
       width = 10, height = 4.5, units = "cm")




################
Fig_overlap <- ggplot(combined, aes(x =  KO, y = log2.or.filtered.plot)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  # significance stars (for overlap > 5)
  geom_text(
    data = combined %>% filter(overlap >= 5),
    aes(label = significance_en, y = 5), # place text just above bar
    color = "black",
    size = 1.5
  ) +
  
  # add "NA" label for overlap <= 5
  geom_text(
    data = combined %>% filter(overlap < 5),
    aes(label = "NA", y = 3.5),
    color = "black",
    size = 1.5
  ) +
  
  labs(
    x = "KOs",
    y = expression(atop(log[2]*"(Odds ratio)", 
                        "in vivo KO effects and culture effects")),
    title = "Overlap of genesets with culture effect and in vivo KO effect"
  ) +
  optimized_theme_fig() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

Fig_overlap
ggsave(out("overlap_NTC_to_invivo.pdf"),
       w= 12, h =4, units = "cm")
Fig_interaction / Fig_overlap 
#####################
limmaRes_int <-  limmaRes %>% filter(coef %in% grep("tissueex.vivo.genotype[A-Za-z0-9]*$",
                                                    limmaRes$coef, value = T))

gsea.res <- data.table() 

run_gsea <- function(limmaRes, enr.terms, celltypes = NULL, coefs = NULL) {
  
  # Initialize the result table
  gsea_res <- data.table() 
  
  # Determine cell types to process
  if (is.null(celltypes)) {
    celltypes <- unique(limmaRes$celltype)
  }
  
  # Determine coefficients to process
  if (is.null(coefs)) {
    coefs <- unique(limmaRes$coef)
  }
  
  # Loop through each cell type
  for (ct in celltypes) {
    
    # Loop through each coefficient
    for (de_grp in coefs) {
      
      # Loop through each database in the enrichment terms
      for (dbx in names(enr.terms)) {
        
        # Subset the limma results based on the current cell type and coefficient
        subset_limmaRes <- limmaRes[limmaRes$celltype == ct & limmaRes$coef == de_grp, ]
        
        # Extract statistics (logFC) and assign gene names as names
        stats <- with(subset_limmaRes, setNames(logFC, nm = ensg))
        
        # Skip this iteration if there are missing values in stats
        if (any(is.na(stats))) {
          next
        }
        
        # Perform fgsea analysis
        fgsea_output <- fgsea(
          pathways = enr.terms[[dbx]],
          stats = stats
          #minSize = 15,   # Example additional arguments, adjust as necessary
          #maxSize = 500,  # Example additional arguments, adjust as necessary
          #nperm = 1000    # Example additional arguments, adjust as necessary
        )
        
        # Check if fgsea output is not empty and append the results to gsea_res
        if (length(fgsea_output) > 0) {
          gsea_res <- rbind(gsea_res, data.table(fgsea_output,
                                                 coef = de_grp,
                                                 celltype = ct,
                                                 db = dbx))
        }
      }
    }
  }
  
  # Return the combined GSEA results
  return(gsea_res)
}
limmaRes_int$celltype <- "glio"
source("src/Ag_enrichR_mouse_genes.R")
gsea.res_glio_KO_int <- run_gsea(limmaRes_int , enr.terms, celltypes = unique(limmaRes_int$celltype),
                     coefs =unique(limmaRes_int$coef))
gsea.res_glio_KO_int$coef <- gsub("tissueex.vivo.genotype","",gsea.res_glio_KO_int$coef)
write_rds(gsea.res_glio_KO_int,
          out("enrichment_interactionKO_at_noRT.rds"))

gsea.res_glio_KO_int <- read_rds(out("enrichment_interactionKO_at_noRT.rds"))

db = "MSigDB_Hallmark_2020"
summary_interaction_KO <- limmaRes_int %>%
  # group by interaction term
  mutate(KO = gsub("tissueex.vivo.genotype","", coef))%>%
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")


create_gsea_plot <- function(db) {
  
  # Filter the data for the given database
  pDT <- gsea.res_glio_KO_int %>%
    #filter(coef %in% koi) %>%
    filter(db == !!db)  # Correctly reference the current database
    #left_join(ko_flags, by = c("coef", "celltype")) %>%
    #left_join(summary_interaction_KO, by = c("coef", "celltype")) %>%
    #filter(Count > 5) %>%
    
   # dplyr::select(-c( "group"))
  # Merge with KO flags
  
  # Step 3 continued: Keep only valid KOs for the specific cell type
  #pDT <- pDT %>% filter(valid_ko == TRUE, padj < 0.05)
  
  # Select the pathways for plotting (both positive and negative)
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n = 7), by = c("coef", "celltype", "pathway")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n = 7), by = c("coef", "celltype", "pathway")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(union(pw.display.pos, pw.display.neg))
  # Remove duplicate rows
  pDT <- pDT %>% distinct()
  
  # Filter pDT to include only selected pathways
  # **Convert list-columns to character format**
  dat <- pDT %>%
    mutate(across(where(is.list), ~sapply(., toString)))  # Converts lists to comma-separated strings
  
  # Save the filtered data table
  write.table(dat, file = basedir(paste0("fgsea_", db, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # Create the plot for the current database
  fig <- ggplot(pDT, aes(x = coef, y = pathway, color = NES, size = pmin(3, -log10(padj)))) +
    geom_point() + 
    scale_color_gradient2(low = "#4C889C", mid = "white", high = "#D0154E", name = TeX("NES")) +
    geom_point(data = pDT[padj < 0.05], shape = 1) +
    scale_size_continuous(
      range = c(0, 1.5),
      name = TeX("$-\\log_{10}(p_{adj})$")
    ) +
    theme_bw() +
    xRot() +
    labs(
      x = "KOs",
      title = "Enriched pathways (Interaction effect)"
    ) +
    facet_grid(cols = vars(celltype), scales = "free", space = "free") +
    optimized_theme_fig() +
    theme(
      strip.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom",
      legend.box = "vertical" ,
      legend.text = element_text(angle = 0, hjust = 1)# Stack legends vertically
    ) +
    guides(
      color = guide_colorbar(title.position = "top"),
      size = guide_legend(title.position = "top")
    )
  #strip.text.y = element_text(angle = 0)) +
  
  # Save the plot for the current database
  ggsave(out("Sup.Fig.interactionKO_noRT", db, "_3rep.pdf"), fig,
         w = 18, h = 15, units = "cm")
  
  
  
}

# Example usage: 
# Create the plot for "MSigDB_Hallmark_2020"
create_gsea_plot("MSigDB_Hallmark_2020")


summary_invivo_KO <- limmaRes_invivo %>%
  # group by interaction term
  mutate(KO = gsub("_in.vivo_noRT","", coef))%>%
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")

Fig_interaction <-  ggplot(summary_interaction_KO,aes(
  x = reorder(KO,signif_count),
  y = log10(signif_count)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
# Display the plot
Fig_interaction
ggsave(out("interaction_KO_number.pdf"),
       w= 12, h = 4, units = "cm")

#overlap with ntc genes--------------------

#############################---------------------------------
#correlation NTC to KO interaction
merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("ensg"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

# Step 1: Calculate correlation with p-value for each KO and celltype
correlation_results <- merged_data %>%
  mutate(KO = gsub("tissueex.vivo.genotype","",coef.x))%>%
  inner_join(summary_interaction_KO, by = c("KO")) %>%
  #filter(signif_count > 10) %>%
  group_by(KO)%>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results <- correlation_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Step 3: Add significance labels based on p-value thresholds
correlation_results <- correlation_results %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )
#order
correlation_results <- correlation_results %>%
  mutate(KO = factor(KO, levels = KO_order))
Correlation <- ggplot(correlation_results, aes(x = KO, y = cor_abs)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  geom_text(aes(label = significance), 
            y =  0.5, 
            color = "black",
            size = 1.5) + 
  labs(
    title = "Correlation of interaction effect to culture effect",
    x = NULL)+
  # fill = "Cell Type")+
  ylab(expression(atop("Correlation of interaction effect", "to culture effects")))+
  
  
  optimized_theme_fig()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

Correlation
ggsave(out("correlation_KO_int_to_ntc.pdf"),
       w= 12, h = 4, units = "cm")
############
#######

#######
#3way interaction summary
limmaRes_3way <-  limmaRes %>%
  filter(coef %in% threeway_terms)

limmaRes_NTC <- limmaRes %>% filter(coef == "tissueex.vivo")

summary_interaction_threeway <- limmaRes_3way %>%
  # group by interaction term
  mutate(KO = gsub("tissueex.vivo.genotype|.RT_statusRT","", coef))%>%
  
  group_by(KO) %>%
  summarise(signif_count = sum(group %in% c("up", "down")), .groups = "drop")

#order
KO_order <- KO_order <- summary_interaction_KO %>%
  arrange(signif_count) %>%
  pull(KO)
# Make KO a factor with levels in your desired order
summary_interaction_threeway <- summary_interaction_threeway %>%
  mutate(KO = factor(KO, levels = KO_order))
#plot
Fig_interaction_threeway <-  ggplot(summary_interaction_threeway,aes(
  x = KO,
  y = log10(signif_count)
)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  labs(
    title = "No. of genes with interaction effects per KO",
    x = NULL,
    y = expression(atop("Number of genes with", 
                        paste("interaction effects ", log[10](n))))
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    
  )
# Display the plot
Fig_interaction_threeway
ggsave(out("No.of_genes_3way_interaction.pdf"))
#############
#correlation 3way------------
unique(limmaRes_3way$coef)
merged_data <- limmaRes_3way %>%
  inner_join(limmaRes_NTC, by = c("ensg"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y,
         adj.P.Val_KO = adj.P.Val.x,
         adj.P.Val_NTC = adj.P.Val.y)

# Step 1: Calculate correlation with p-value for each KO and celltype
correlation_results <- merged_data %>%
  mutate(KO = gsub("tissueex.vivo.genotype|.RT_statusRT","",coef.x))%>%
  inner_join(summary_interaction_KO, by = c("KO")) %>%
  #filter(signif_count > 10) %>%
  group_by(KO)%>%
  summarize(
    cor_abs = cor(abs(logFC_NTC), abs(logFC_KO), method = "pearson"),
    p_value = cor.test(abs(logFC_NTC), abs(logFC_KO), method = "pearson")$p.value,  # Get p-value
    .groups = 'drop'
    
  )

correlation_results <- correlation_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# Step 3: Add significance labels based on p-value thresholds
correlation_results <- correlation_results %>%
  mutate(
    significance = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE             ~ ""
    )
  )
#order
correlation_results <- correlation_results %>%
  mutate(KO = factor(KO, levels = KO_order))
Correlation_3way_to_NTC <- ggplot(correlation_results, aes(x = KO, y = cor_abs)) +
  geom_col(color = "darkgrey", fill = NA, width = 0.6) +
  
  geom_text(aes(label = significance), 
            y =  0.25, 
            color = "black",
            size = 1.5) + 
  labs(
    title = "Correlation of 3-way interaction effect to culture effect",
    x = NULL)+
  # fill = "Cell Type")+
  ylab(expression(atop("Correlation of interaction effect", "to culture effects")))+
  
  
  optimized_theme_fig()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing = unit(0.02, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
Correlation_3way_to_NTC
ggsave(out("correlation_3way_to_ntc.pdf"),
       w= 8, h =4)
#############################

#########
# ISGs with interaction effects-----------
InDir <- dirout("Ag_top_pathway_genes")

genesets <- fread(InDir("goi_logFC_perturb_seq.tsv"))
genesets <- genesets %>%
  dplyr::select(c("genes","logFC","adj.P.Val","pathway","celltype","group"))

# function for selecting genes
get_genes_for_pathway <- function(pathway_pattern, genesets) {
  unique(genesets[grepl(pathway_pattern, genesets$pathway, ignore.case = TRUE), ]$genes)
}

####
source("src/Ag_enrichR_mouse_genes.R")
pathways <- list(
  ISG_core_all = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
    filter(L1=="ISG_Core")%>%pull(value),
  Cholesterol_all = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  mTORC1_all = enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`,
  ISGs_all = union(enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`,
                   enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`),
  ROC_all = enr.terms$MSigDB_Hallmark_2020$`Reactive Oxygen Species Pathway`,
  
  #Cholesterol = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  Hypoxia_all = enr.terms$MSigDB_Hallmark_2020$`Hypoxia`,
  Glycolysis_all = enr.terms$MSigDB_Hallmark_2020$`Glycolysis`,
  Protein_loc_all = Reduce(union, list(
    enr.terms$MSigDB_Hallmark_2020$`Myc Targets V1`,
    enr.terms$GO_Biological_Process_2023$`Protein Import (GO:0017038)`,
    enr.terms$GO_Biological_Process_2023$`Protein Insertion Into ER Membrane (GO:0045048)`,
    enr.terms$GO_Biological_Process_2023$`Protein Insertion Into Membrane (GO:0051205)`,
    enr.terms$GO_Biological_Process_2021$`RNA biosynthetic process (GO:0032774)`
  ))
)
annotate_pathway <- function(df, pathways_list) {
  df %>%
    mutate(
      pathway = case_when(
        ensg %in% pathways_list$ISGs_all ~ "ISGs_all",
        ensg %in% pathways_list$ISG_core_all ~ "ISG_core_all",
        ensg %in% pathways_list$mTORC1_all ~ "mTORC1_all",
        ensg %in% pathways_list$Cholesterol_all ~ "Cholesterol_all",
        ensg %in% pathways_list$ROC_all ~ "ROS_all",
        ensg %in% pathways_list$Hypoxia_all ~ "Hypoxia_all",
        ensg %in% pathways_list$Glycolysis_all ~ "Glycolysis_all",
        ensg %in% pathways_list$Protein_loc_all ~ "Protein_localization_all",
        TRUE ~ "Other"
      )
    )
}

limmaRes_NTC$KO <- "NTC"
limmaRes_NTC$coef <- "NTC"

head(limmaRes_int)

top_NTC <- limmaRes_NTC %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(desc(abs(logFC)))%>%
  slice_head(n = 50)
 
top_genes <-  limmaRes_NTC %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(desc(abs(logFC)))%>%
  slice_head(n = 50)%>%
  pull(ensg)%>%
  unique() 

colnames(limmaRes_NTC) %in% colnames(limmaRes_int)
# Original gene sets
# ISGs <- get_genes_for_pathway("ISGs", genesets)
# ISG_core <- get_genes_for_pathway("ISG_core", genesets)
# mTORC1 <- get_genes_for_pathway("mTORC1", genesets)
# Cholesterol <- get_genes_for_pathway("Cholesterol", genesets)
# Housekeeping <- c("Tbp", "Ubc")
# NTC_ISG <- limmaRes_NTC %>%
#   filter(adj.P.Val < 0.05)%>%
#   arrange(abs(logFC))%>%
#   slice_head(n = 20)
#   
limmaRes_RT_int <- limmaRes %>%
  filter(coef == "tissueex.vivo.RT_statusRT")
limmaRes_RT_int$KO <- "Radiotherapy"
limmaRes_int$KO <- gsub("tissueex.vivo.genotype","",limmaRes_int$coef)
top_NTC_int <- top_NTC %>%
  rbind(limmaRes_int %>% filter(ensg %in% top_NTC$ensg))%>%
  rbind(limmaRes_RT_int %>% filter(ensg %in% top_NTC$ensg))

top_NTC_int <- annotate_pathway(top_NTC_int, pathways)
top_NTC_int$KO <- factor(top_NTC_int$KO, 
                         levels = c("NTC",
                                    setdiff(top_NTC_int$KO, "NTC")))

# fig---------
plot <- ggplot(top_NTC_int , aes(x = KO, y = ensg,
                                           color = pmin(2, pmax(-2, logFC)),
                                 size = pmin(3, -log10(adj.P.Val))))+  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0,2),
    #limits = c(0,5),
    #breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(
  title = str_wrap("Interaction effect (KO/RT) - culture model) of top culture effect genes", width = 80),
       x = "NTC/Interaction",
       y = "Genes")+
  #facet_grid(cols = vars(pathway), scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()+theme(
    legend.position = "right",
    strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
    axis.text.x = element_text(angle = 45)
    )

plot
ggsave(out("top_NTCs_interction.pdf"), w = 13, h = 12.5, units = "cm")
#axis.text.y = element_blank()
#interaction RT-------------------

limmaRes_RT_int$KO <- "NTC"
top_NTC$coef <- "culture effect"
top_NTC_RT_int <- top_NTC %>%
  rbind(limmaRes_RT_int %>% filter(ensg %in% top_NTC$ensg))
plot <- ggplot(top_NTC_RT_int , aes(x = coef, y = genes,
                                 color = pmin(2, pmax(-2, logFC)),
                                 size = pmin(3, -log10(adj.P.Val))))+  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0,1.5),
    #limits = c(0,5),
    #breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(
    title = str_wrap("Interaction effect (RT-culture model) of top culture effect genes", width = 40),
    x = "KOs",
    y = "Genes"
  )+
  #facet_grid(cols = vars(pathway), scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()+theme(
    legend.position = "right",
    strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
    axis.text.x = element_text(angle = 45)
  )

plot
ggsave(out("top_NTCs_RT_interction.pdf"), w = 4, h = 12, units = "cm")

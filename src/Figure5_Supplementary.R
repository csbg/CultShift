source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification_Mye.R")


basedir <- dirout("Figure5_Supplementary_Mye")
#Supplementary-------------------

InDir2 <- dirout("Figure1_Mye")
InDir3 <- dirout("Figure1_Supplementary_Mye")
InDir7  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_per_pathway_fgsea_in.vivo_Mye")

gsea.res <- read_rds(InDir7("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction_invivo.rds"))
gsea.res <- gsea.res %>%
  filter(coef %in% grep("interaction",gsea.res$coef, value = T))
gsea.res$coef <- gsub("interaction","",gsea.res$coef )
# Step 2: Summarize to find KOs with at least one valid cell type
valid_ko_summary <- ko_flags %>%
  group_by(genotype) %>%
  summarize(has_valid_celltype = any(valid_ko), .groups = 'drop')

# Step 3: Filter GSEA results based on valid KOs

db = "MSigDB_Hallmark_2020"
create_gsea_plot <- function(db) {
  
  # Filter the data for the given database
  pDT <- gsea.res %>%
    filter(coef %in% koi) %>%
    filter(db == !!db) %>%  # Correctly reference the current database
    left_join(ko_flags, by = c("coef", "celltype")) %>%
    left_join(summary_df, by = c("coef", "celltype")) %>%
    filter(Count > 5) %>%
    filter(valid_ko)%>%
    dplyr::select(-c("Count", "Regulation"))
  # Merge with KO flags
  
  # Step 3 continued: Keep only valid KOs for the specific cell type
  pDT <- pDT %>% filter(valid_ko == TRUE, padj < 0.05)
  
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
      legend.text = element_text(angle = 55, hjust = 1)# Stack legends vertically
    ) +
    guides(
      color = guide_colorbar(title.position = "top"),
      size = guide_legend(title.position = "top")
    )
  #strip.text.y = element_text(angle = 0)) +
  
  # Save the plot for the current database
  ggsave(basedir("Sup.Fig.5A_fgsea", db, "_3rep.pdf"), fig,
         w = 18, h = 15, units = "cm")
  
  
  
}

# Example usage: 
# Create the plot for "MSigDB_Hallmark_2020"
create_gsea_plot("MSigDB_Hallmark_2020")
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
# You can also loop through the databases to generate plots for all of them
lapply(databases, create_gsea_plot)

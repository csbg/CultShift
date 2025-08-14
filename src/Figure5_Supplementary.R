source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
source("src/Ag_ko_classification.R")
base <- "Figure5_Supplementary"
basedir <- dirout("Figure5_Supplementary")
#Supplementary-------------------
InDir1  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide")
InDir2 <- dirout("Figure1")
InDir3 <- dirout("Figure1_Supplementary")
gsea.res <- read_rds(InDir1("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
gsea.res$coef <- gsub("interaction","",gsea.res$coef )


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
      strip.text.x = element_text(angle = 90),
      legend.position = "bottom",
      legend.box = "vertical" ,
      legend.text = element_text(angle = 45, hjust = 1)# Stack legends vertically
    ) +
    guides(
      color = guide_colorbar(title.position = "top"),
      size = guide_legend(title.position = "top")
    )
  #strip.text.y = element_text(angle = 0)) +
  
  # Save the plot for the current database
  ggsave(basedir("Sup.Fig.5A_fgsea", db, "_3rep.pdf"), fig,
         w = 9, h = 12, units = "cm")
  
  
  
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
#Sup.Fig5B-----------------
limmaRes <- read_rds(InDir_int("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
#exclude fig1 genes
genes_fig1 <- read_rds(InDir2("genes_fig1.rds"))
# dplyr::select significant genes fron interaction, present in results

supp_fig1_genes <- read_rds(InDir3("supp_fig1_genes.rds"))

genes_ntc <- genes_fig1 %>%
  mutate(geneset = pathway)%>%
  dplyr::select(genes,geneset) %>%
  filter(!(geneset %in% "mTORC1_or_Cholesterol"))%>%
  rbind(supp_fig1_genes)
top_int <- limmaRes %>%
  filter(ensg %in% genes_ntc$genes)%>%
  mutate(genes = ensg) %>%
  inner_join(genes_ntc, by = "genes")


unique(top_int$geneset)

top_int <- limmaRes %>%
  filter(ensg %in% genes_ntc$genes)%>%
  mutate(genes = ensg) %>%
  inner_join(genes_ntc, by = "genes")%>%
  inner_join(ko_flags, by = c("coef","celltype"))%>%
  filter(valid_ko)%>%
  #filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","celltype")) %>%
  filter(Count > 10) %>% 
  filter(geneset != "ISG core")
# Wrap function for facet labels
wrapped_labeller <- labeller(
  geneset = label_wrap_gen(width = 10) # wrap text after ~8 characters
)

Sup.Fig.5B <- ggplot(top_int, 
                     aes(x = coef, y = ensg,
                         color = pmin(1.5, pmax(-1.5, logFC)),
                         size = pmin(3, -log10(adj.P.Val))
                     )) +
  geom_point() + 
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(0, 1.5),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    title = "Differentially expressed genesets",
    x = "KOs",
    y = "Genes"
  ) +
  facet_grid(
    rows = vars(geneset), 
    cols = vars(celltype),
    scales = "free", 
    space = "free", 
    labeller = wrapped_labeller
  ) +
  theme_bw() +
  optimized_theme_fig() + 
  theme(
    axis.text = element_text(size = 5),
    strip.text.y = element_text(angle = 90, hjust = 0),
    strip.text.x = element_text(angle = 90, hjust = 0),
    legend.position = "bottom",
    legend.text = element_text(angle = 45, hjust = 1)  # Rotates legend labels
  )

Sup.Fig.5B



ggsave(basedir("Sup.Fig.5B.pdf"), w=9, h= 16, units = "cm")

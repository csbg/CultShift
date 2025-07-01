###############
source("src/00_init.R")
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
require(fgsea)
library(latex2exp)

source("src/Ag_Optimized_theme_fig.R")
#####################################################################
inDir  <-  dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
InDir1 <- dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide")
InDir5 <- dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
base  <-  "Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_per_pathway_fgsea/"
basedir  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide_per_pathway_fgsea")
########################################################################
limmaRes  <-  read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
limmaRes  <-  limmaRes%>% filter(celltype != "MEP")
################################
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
################################################################################
adj_p_cutoff <- 0.05
logfc_cutoff <- 1
gsea.res <- read_rds(InDir1("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))

meta <- fread(InDir5("meta_cleaned.tsv")) # Read data
meta <- as.data.frame(meta)               # Convert to dataframe (optional)
rownames(meta) <- meta[[1]]   

meta <- meta[, -1, drop = FALSE] 
colnames(meta) <- gsub("rowname","sample1", colnames(meta))
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop")%>%
  mutate(coef = genotype)


selected_KOs <- meta %>%
  group_by(genotype, tissue, celltype) %>%                  # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample1), .groups = 'drop') %>% # Count distinct samples for each group
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>% # Spread tissue to separate columns (in.vivo and ex.vivo)
  group_by(genotype) %>%                                    # Regroup by genotype
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) %>%              # Keep genotypes that have at least one celltype with 3+ samples in both tissues
  pull(genotype) %>% unique()


summary_df <- limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")


count_threshold = 10
coefficients  <-  summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)%>%
  pull(coef)%>%
  unique()
correlation_deg <- read_rds(InDir1("correlation_deg.rds"))
KO_list <- correlation_deg %>% filter(correlation < 0.5, num_degs >= 10) %>%
  pull(genotype)%>%
  unique()
koi <- Reduce(intersect, list(selected_KOs,  KO_list))#, coefficients)) #only valid_k
################################################################################
#fgsea--------------------------------------------------------------------------
################################################################################

enr.terms <- enrichrGetGenesets(databases)
# # save(enr.terms, file=out("Genesets_Human.RData"))
# Convert to mouse --------------------------------------------------------
hm.map <- fread(PATHS$RESOURCES$HM.MAP, check.names = T)
hm <- unique(hm.map[Human.gene.name != "",c("Gene.name", "Human.gene.name")])
names(hm) <- c("Mouse", "Human")
enr.terms <- lapply(enr.terms, function(dbl){
  dbl <- lapply(dbl, function(gs){
    unique(hm[Human %in% gs]$Mouse)
  })
  dbl[sapply(dbl, length) > 0]
})



process_db <- function(dbx, gsea_res, base, size_n = 0, coefs) {
  print(paste("Processing database:", dbx))
  print(paste("Number of coefficients:", length(coefs)))
  
  # Filter data based on the conditions
  pDT <- gsea_res[gsea_res$db == dbx & size > size_n & padj < 0.05 & coef %in% coefs,]
  
  # Check if pDT has any data after filtering
  if (nrow(pDT) == 0) {
    print(paste("No data found for database:", dbx))
  } else {
    print(paste("Found", nrow(pDT), "entries for", dbx))
  }
  
  ## Splitting the task to handle both ends of the NES spectrum - positive and negative
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5), by=c("celltype", "coef")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("celltype", "coef")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  pDT <- pDT[pathway %in% pw.display]
  
  # Hierarchical ordering
  pDT <- hierarch.ordering(pDT, "celltype", "pathway", "NES", TRUE)
  
  # Plotting
  ggplot(pDT, aes(x=gsub("ex.vivo:", "", coef), y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
    scale_color_gradient2(low="blue", mid="white", high="red") +
    geom_point(data=pDT[padj < 0.05]) +
    scale_size_continuous(range=c(0, 5), limits = c(0, 5)) +
    theme_bw(12) +
    xRot() +
    facet_wrap(vars(celltype)) +
    labs(x = "interaction-KO") +
    theme(axis.text = element_text(size = 10)) +
    optimized_theme_fig()+
    theme(strip.text.x  = element_text(angle = 90))
  
  # Save the plot
  dat <- dirout(paste0(base,"/", dbx))
  ggsave(dat("GSEA_plot_", dbx,"_",size_n, ".pdf"), width = 30,
         height = length(unique(pDT$pathway)) * 0.3 + 3, limitsize = FALSE)
}

# Apply the function to each unique `db` using purrr::walk
unique_dbs <- unique(gsea.res$db)
walk(unique_dbs, process_db, gsea_res = gsea.res, base = base,size_n =5,coefs=coefficients)


################################

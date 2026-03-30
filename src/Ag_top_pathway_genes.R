source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_ko_classification.R")
source("src/Ag_enrichR_mouse_genes.R")
library(tidyverse)
library(enrichR)
library(purrr)
basedir <- dirout("Ag_top_pathway_genes")

limmaRes_NTC <- read_rds(InDir_NTC("limma_perCTex.vivovsin.vivo.rds"))
pathways <- list(
  ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
    filter(L1=="ISG_Core")%>%pull(value),
  Cholesterol = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  mTORC1 = enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`,
  ISGs = union(enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`,
               enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`),
  ROS = enr.terms$MSigDB_Hallmark_2020$`Reactive Oxygen Species Pathway`,
  
  #Cholesterol = enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`,
  Hypoxia = enr.terms$MSigDB_Hallmark_2020$`Hypoxia`,
  Glycolysis = enr.terms$MSigDB_Hallmark_2020$`Glycolysis`,
  Protein_loc = Reduce(union, list(
    enr.terms$MSigDB_Hallmark_2020$`Myc Targets V1`,
    enr.terms$GO_Biological_Process_2023$`Protein Import (GO:0017038)`,
    enr.terms$GO_Biological_Process_2023$`Protein Insertion Into ER Membrane (GO:0045048)`,
    enr.terms$GO_Biological_Process_2023$`Protein Insertion Into Membrane (GO:0051205)`,
    enr.terms$GO_Biological_Process_2021$`RNA biosynthetic process (GO:0032774)`
  ))
)
#IFN genes
get_top_genes_up <- function(pathway_name,
                             pathway_genes,
                             limma_results, 
                             logFC_threshold = 1,
                             pval_threshold = 0.05, 
                             top_n = Inf,
                             min_celltypes = 3) {
  
  
  
  # Filter for significant upregulated genes in the pathway
  filtered <- limma_results %>%
    filter(group == "up") %>%
    filter(adj.P.Val < pval_threshold) %>%
    filter(toupper(genes) %in% toupper(pathway_genes))
  
  # Count in how many cell types each gene is upregulated
  gene_counts <- filtered %>%
    group_by(genes) %>%
    summarise(n_celltypes = n_distinct(celltype)) %>%
    filter(n_celltypes >= min_celltypes)
  
  # Keep only those genes
  filtered_genes <- filtered %>%
    filter(toupper(genes) %in% toupper(gene_counts$genes)) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = top_n) %>%
    pull(genes) %>%
    unique()
  
  # Data for plotting
  pathway_plot <- limma_results %>%
    filter(toupper(genes) %in% toupper(filtered_genes))
  
  return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
}

get_top_genes_down <- function(pathway_name,
                               pathway_genes,
                               limma_results, 
                               logFC_threshold = 1,
                               pval_threshold = 0.05, 
                               top_n = Inf,
                               min_celltypes = 3) {
  
  
  
  # Filter for significant downregulated genes in the pathway
  filtered <- limma_results %>%
    filter(group == "down") %>%
    filter(adj.P.Val < pval_threshold) %>%
    filter(toupper(genes) %in% toupper(pathway_genes))
  
  # Count how many cell types each gene is down in
  gene_counts <- filtered %>%
    group_by(genes) %>%
    summarise(n_celltypes = n_distinct(celltype)) %>%
    filter(n_celltypes >= min_celltypes)
  
  # Keep only those genes
  filtered_genes <- filtered %>%
    filter(toupper(genes) %in% toupper(gene_counts$genes)) %>%
    arrange(logFC) %>%
    slice_head(n = top_n) %>%
    pull(genes) %>%
    unique()
  
  # Data for plotting
  pathway_plot <- limma_results %>%
    filter(toupper(genes) %in% toupper(filtered_genes))
  
  return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
}



# Initialize an empty list to store results

results_up <- map(names(pathways), function(pathway) {
  get_top_genes_up(
    pathway_name = pathway,
    pathway_genes = pathways[[pathway]],
    limma_results = limmaRes_NTC
  )
})
results_down <- map(names(pathways), function(pathway) {
  get_top_genes_down(
    pathway_name = pathway,
    pathway_genes = pathways[[pathway]],
    limma_results = limmaRes_NTC
  )
})
# Optionally, name the elements of the list using pathway names
names(results_up) <- names(pathways)
names(results_down) <- names(pathways)


# Access the results for each pathway
ISGs <- results_down$ISGs$pathway_plot
ISG_core <- results_down$ISG_core$pathway_plot
mTORC1 <- results_up$mTORC1$pathway_plot
Cholesterol <- results_up$Cholesterol$pathway_plot
ROS <- results_up$ROS$pathway_plot


#Combine results into a single data frame and add 'pathway' column
genes_of_interest <- bind_rows(
  list(
    ISGs       = results_down$ISGs$pathway_plot,
    ISG_core   = results_down$ISG_core$pathway_plot,
    Cholesterol = results_up$Cholesterol$pathway_plot,
    mTORC1     = results_up$mTORC1$pathway_plot,
    ROS = results_up$ROS$pathway_plot
  ),
  .id = "pathway"  # This creates a new column 'pathway' with the list names
)

# Write to file
write.table(
  genes_of_interest,
  file = basedir("goi_logFC_perturb_seq.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#write.table(genes_of_interest,basedir("goi_logFC_perturb_seq.tsv"))

# # Access the results for each pathway
# ISGs <- results_down$ISGs$top_genes
# ISG_core <- results_down$ISG_core$top_genes
# mTORC1 <- results_up$mTORC1$top_genes
# Cholesterol <- results_up$Cholesterol$top_genes






# # library(dplyr)
# # library(purrr)
# # library(tibble)
# # 
# # # Suppose 'results' is your list of lists (top_genes + pathway_plot)
# df_results <- map2_df(
#   names(pathways),    # pathway names
#   results_down,
#   ~ {
#     # .x = pathway name, .y = list(top_genes, pathway_plot)
#     pathway <- .x
#     top_genes <- .y$top_genes
#     pathway_plot <- .y$pathway_plot
# 
#     # Add a column for the pathway
#     pathway_plot %>%
#       filter(toupper(genes) %in% toupper(top_genes)) %>%
#       mutate(pathway = pathway)
#   }
# )
# df_results_up <- map2_df(
#   names(pathways),    # pathway names
#   results_up,
#   ~ {
#     # .x = pathway name, .y = list(top_genes, pathway_plot)
#     pathway <- .x
#     top_genes <- .y$top_genes
#     pathway_plot <- .y$pathway_plot
#     
#     # Add a column for the pathway
#     pathway_plot %>%
#       filter(toupper(genes) %in% toupper(top_genes)) %>%
#       mutate(pathway = pathway)
#   }
# )
# 
# # # Optional: reorder columns
# df_results_down <- df_results %>% dplyr::select(pathway, everything())
# df_results_up <- df_results %>% dplyr::select(pathway, everything())
# # 
# # # Check
# # head(df_results)
# ISGs <- df_results_down %>%
#   filter(pathway == "ISGs")
# ISG_core <- df_results_down %>%
#   filter(pathway == "ISG_core")
# mTORC1 <- df_results_up %>%
#   filter(pathway == "mTORC1")
# Cholesterol <- df_results_up %>%
#   filter(pathway == "Cholesterol")
# 
# dev.off()
# # Diverging bar plot
# ggplot(df_results %>%
#          filter(pathway == "ISGs"), aes(x = logFC, y = genes, fill = group)) +
#   geom_bar(stat = "identity", width = 0.7) +
#   scale_fill_manual(values = c("up" = "red", "down" = "blue", "n.s" = "grey")) +
#   labs(title = "ISG Gene Expression Changes",
#        x = "log2 Fold Change",
#        y = "Gene") +
#   facet_grid(cols = vars(celltype))+
#   theme(axis.text.y = element_text(size = 10))
# 
# ###########
# ISG_core <- df_results %>%
#   filter(pathway == "ISG_core")
# 
# # Diverging bar plot
# ggplot(ISG_core, aes(x = logFC, y = genes, fill = group)) +
#   geom_bar(stat = "identity", width = 0.7) +
#   scale_fill_manual(values = c("up" = "red", "down" = "blue", "n.s" = "grey")) +
#   labs(title = "ISG Gene Expression Changes",
#        x = "log2 Fold Change",
#        y = "Gene") +
#   facet_grid(cols = vars(celltype))+
#   theme(axis.text.y = element_text(size = 10))
# ########################################
# unique(df_results$pathway)
# Cholesterol <- df_results %>%
#   filter(pathway == "Cholesterol")
# 
# # Diverging bar plot
# ggplot(Cholesterol, aes(x = logFC, y = genes, fill = group)) +
#   geom_bar(stat = "identity", width = 0.7) +
#   scale_fill_manual(values = c("up" = "red", "down" = "blue", "n.s" = "grey")) +
#   labs(title = "Cholesterol Gene Expression Changes",
#        x = "log2 Fold Change",
#        y = "Gene") +
#   facet_grid(cols = vars(celltype))+
#   theme(axis.text.y = element_text(size = 10))
# ######################
# mTORC1 <- df_results %>%
#   filter(pathway == "mTORC1")
# 
# # Diverging bar plot
# ggplot(mTORC1, aes(x = logFC, y = genes, fill = group)) +
#   geom_bar(stat = "identity", width = 0.7) +
#   scale_fill_manual(values = c("up" = "red", "down" = "blue", "n.s" = "grey")) +
#   labs(title = "mTORC1 Gene Expression Changes",
#        x = "log2 Fold Change",
#        y = "Gene") +
#   facet_grid(cols = vars(celltype))+
#   theme(axis.text.y = element_text(size = 10))
# 
# genes_of_interest <- rbind(ISG_core,ISGs, Cholesterol, mTORC1)
# write.table(genes_of_interest,basedir("goi_logFC_perturb_seq.tsv"))

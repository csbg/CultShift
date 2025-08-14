#Figure2
#Invivo vs exvivo 
#Aarathy
#24-09-2024
###############################################################################
#load libraries and functions--------------------------------------------------
###############################################################################
source("src/00_init.R")
source("src/Ag_Optimized_theme_fig.R")
source("src/Ag_top_genes_per_pathway.R")
library(tidyverse)
library(enrichR)
library(purrr)
library(pheatmap)
library(ComplexHeatmap)
library("scales")
#
#directories ------
#
base<-"Figure2"
basedir<-dirout("Figure2")
#
#enrichr database-----
#
databases = c("KEGG_2019_Mouse",
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse",
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019",
              "Reactome_2022",
              "GO_Molecular_Function_2023",
              "GO_Biological_Process_2023",
              "CellMarker_2024")
enr.terms <- enrichrGetGenesets(databases)

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
#


#fig2.2----------------
InDir5 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
meta <- read_rds(InDir5("meta.rds"))
meta$sample1 <- rownames(meta)


# Check if there are at least 3 distinct samples per tissue for each genotype and celltype
ko_flags <- meta %>%
  group_by(genotype, celltype, tissue) %>%
  summarize(num_samples = n_distinct(sample1), .groups = 'drop') %>%
  pivot_wider(names_from = tissue, values_from = num_samples, values_fill = 0) %>%
  mutate(valid_ko = (in.vivo >= 3 & ex.vivo >= 3)) %>%
  group_by(genotype, celltype) %>%
  summarize(valid_ko = any(valid_ko), .groups = "drop") 



selected_KOs <- meta %>%
  group_by(genotype, tissue, celltype) %>%                  # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample1), .groups = 'drop') %>% # Count distinct samples for each group
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>% # Spread tissue to separate columns (in.vivo and ex.vivo)
  group_by(genotype) %>%                                    # Regroup by genotype
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) %>%              # Keep genotypes that have at least one celltype with 3+ samples in both tissues
  pull(genotype) %>% unique()

# Read in the correlation matrix and DEG data
Indir1 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/")
correlation_matrix_data <- read_rds(Indir1("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(Indir1("DEGs_per_tissue.rds"))

# Prepare correlation data
correlation_matrix_no_na <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)

# Perform hierarchical clustering
hc_cols <- hclust(dist(t(correlation_matrix_no_na)), method = "ward.D2")
column_order <- colnames(correlation_matrix_data)[hc_cols$order]

# Reshape correlation data
correlation <- correlation_matrix_data %>%
  as_tibble(rownames = "celltype") %>%
  pivot_longer(cols = colnames(correlation_matrix_data),
               names_to = "genotype",
               values_to = "correlation",
               values_drop_na = FALSE)

# Prepare DEG data
deg_plot_data <- deg_plot_data %>%
  group_by(genotype, celltype) %>%
  filter(num_degs == max(num_degs, na.rm = TRUE))

# Merge DEG data with correlation data
correlation_deg <- inner_join(deg_plot_data, correlation,
                              by = c("celltype", "genotype"))

# Merge with KO flags to include valid KO status
correlation_deg_flagged <- correlation_deg %>%
  left_join(ko_flags, by = c("genotype", "celltype")) %>%
  filter(genotype %in% selected_KOs)%>%
  na.omit()%>% 
  filter(valid_ko)# Keep only selected KOs

# Update the genotype factor levels for plotting
correlation_deg_flagged$genotype <- factor(correlation_deg_flagged$genotype, levels = column_order)
#fig-----
# Create the plot
Fig2.2 <- correlation_deg_flagged %>%
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
    high = "#D0154E"
                       ) +
  scale_size_continuous(
    limits = c(1,4),
    breaks = c(1,2,3),
    name =TeX("$\\log_{10}\\; (\\No.\\; of \\;DEGs)$")
                        )+
  labs(x = "KOs",
       y = "Celltype") +
  optimized_theme_fig()+
  theme(
    legend.position = "right",
       )
Fig2.2



#
#Fig2.3 DEG interaction logFC ---------------------------------------------
#

InDir2 <- dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide/")
InDir3 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")

gsea.res <- read_rds(InDir2("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))

limmaRes<- read_rds(InDir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))

adj_p_cutoff <- 0.05
logfc_cutoff <- 1

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
KO_list <- correlation_deg %>% filter(correlation < 0.5, num_degs >= 20) %>%
  pull(genotype)%>%
  unique()

#kos of interest (koi)
#KO selection
#atleast 3 rep in both tissues, 

koi <- Reduce(intersect, list(selected_KOs, coefficients, KO_list))
#fig ------------
Fig2.4 <- summary_df %>%
  filter(Count > 10, coef %in% koi) %>%
  # Plot with log10 scale on y-axis
  ggplot(aes(x = coef, y = Regulation, color = Regulation)) +
  geom_point(aes(size=log10(Count))) +    
  scale_size_continuous(
    limits = c(1,4),
    breaks = c(1,2,3),
    name =TeX("$-\\; \\log_{10}\\; (\\No.\\; of \\;DEGs)$")
    ) + 
  scale_color_manual(
    values = c("Downregulated" =  "#4C889C", "Upregulated" = "#D0154E")
    ) +
  facet_grid(rows = vars(celltype), space = "free") +
  labs(title="No. of genes with interaction effect of culture condition in KOs",
       x="KOs")+
 
  #scale_y_log10() +  # Apply log10 scale to the y-axis
  optimized_theme_fig()+theme(legend.position = "none")
Fig2.4
ggsave(basedir("Fig2.4.pdf"),plot=Fig2.4_legend,w=5,h=4, units = "cm")

#from 12_enrichr
###############################################################################
#Fig2.3.1_supplementary
InDir4 <- dirout("Figure1")
# grouping gene sets specifically for ISG and cholesterol
genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
colnames(genes_fig_1) <- c("pathway","ensg")

filtered_genes <- limmaRes %>%
  #filter(group != "n.s")%>%
  filter(ensg %in% genes_fig_1$ensg, coef %in% koi) %>%
  group_by(celltype, coef) %>%
  merge(genes_fig_1, by = "ensg") %>%
  left_join(ko_flags, by = c("coef" = "genotype", "celltype")) %>%  # Merge with KO flags per cell type
  filter(valid_ko == TRUE)  # Keep only valid KOs for the specific cell type

# Recode pathways for better labeling
filtered_genes$pathway <- recode(filtered_genes$pathways,
                                 "ISG_core" = "ISG core",
                                 "mTORC1_or_Cholesterol" = "mTORC1 or Cholesterol")

# Create the plot
Fig2.3.1_supplementary <- ggplot(filtered_genes, aes(x = coef, y = ensg,
                                     color = pmin(2, pmax(-2, logFC)) ,
                                     size = pmin(5, -log10(adj.P.Val))
)) +  # Use alpha based on validity
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(1,3),
    limits = c(0,5),
    breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(title = "Differentially expressed genesets",
       x = "KOs",
       y = "Genes")+
  facet_grid(rows = vars(pathway), cols = vars(celltype), scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()

Fig2.3.1_supplementary
#paper
ggsave(basedir(paste0("Fig2.3.1_supplementary_with_legend.pdf")),plot = Fig2.3.1_supplementary, 
       width = 18.3,
       height = 12, units = "cm")
ggsave(basedir(paste0("Fig2.3.1_supplementary_without_legend.pdf")),plot = Fig2.3.1_supplementary+theme(
  legend.position = "none"
), 
       width = 18.3,
       height = 12, units = "cm")
# #ppt
# ggsave(basedir(paste0("Fig2.3.1_supplementary_ppt.pdf")),plot = Fig2.3.1_supplementary+optimized_theme(), 
#        width =12,
#        height = 6)

#################################################################################
# Pathways of interest (you can adjust these as per your needs)
get_top_genes <- function(pathway_name,
                          pathway_genes,
                          limma_results,
                          logFC_threshold = 1,
                          pval_threshold = 0.05,
                          top_n = 1) {
  # Filter limma results for the pathway genes
  
  filtered_genes <- limma_results %>%
    filter(group != "n.s") %>%  # Skip non-significant genes
    filter(adj.P.Val < pval_threshold, abs(logFC) > logFC_threshold) %>%
    filter(toupper(ensg) %in% toupper(pathway_genes)) %>%  # Match pathway genes
    group_by(celltype, coef) %>%
    arrange(adj.P.Val) %>%  # Sort by adjusted p-value
    slice_head(n = top_n) %>%  # Get the top 5 genes per group (celltype)
    pull(ensg) %>% 
    unique()
  
  # Filter the main table to return data for plotting
  pathway_plot <- limma_results %>%
    filter(toupper(ensg) %in% toupper(filtered_genes), coef %in% koi)%>%
    mutate(pathway = pathway_name)
  
  return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
}

pathway_list <- list(
  "EMT" = enr.terms$MSigDB_Hallmark_2020$`Epithelial Mesenchymal Transition`,
 # "Myc Targets" = enr.terms$MSigDB_Hallmark_2020$`Myc Targets V2`,
  "ROS" = enr.terms$MSigDB_Hallmark_2020$`Reactive Oxygen Species Pathway`,
  "E2F targets" = enr.terms$MSigDB_Hallmark_2020$`E2F Targets`,
  "TNF" = enr.terms$MSigDB_Hallmark_2020$`TNF-alpha Signaling via NF-kB`)
 # "KRAS Dn" = enr.terms$MSigDB_Hallmark_2020$`KRAS Signaling Dn`


# Loop through each pathway to extract and plot the top genes
# by default, top 5 genes per celltype per KO are obtained 
pathway_results <- lapply(names(pathway_list), function(pathway_name) {
  pathway_genes <- pathway_list[[pathway_name]]
  get_top_genes(pathway_name, pathway_genes, limmaRes)
})


# Combine the results from all pathways for plotting

combined_results <- map_dfr(pathway_results, "pathway_plot")
combined_results <- combined_results %>%
  left_join(ko_flags, by = c("coef" = "genotype", "celltype")) %>%  # Merge with KO flags per cell type
  filter(valid_ko == TRUE)
pathway_gene <- rbind(combined_results,filtered_genes)
pathway_gene %>% write_rds(basedir("pathway_gene.rds"))
# Plot the results
Fig2.3.2_supplementary <- ggplot(combined_results, aes(x = coef, y = ensg,
                                       color = pmin(2, pmax(-2, logFC)),  # Clamp logFC between -2 and 2
                                       size = pmin(10, -log10(adj.P.Val))  # Set size by adjusted p-value
)) +
  geom_point() +  # Create points
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name =TeX("$\\log_{2}\\; (FC)$")
  ) +
  scale_size_continuous(
    range = c(1,3),
    limits = c(0,5),
    breaks = c(1,3,5),
    name =TeX("$-\\log_{10}(p_{adj})$")
  )+
  labs(title = "Differentially expressed genesets",
       x = "KOs",
       y = "Genes")+
  facet_grid(rows = vars(pathway), cols = vars(celltype),
             scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()+theme(
    axis.text = element_text(size=5)
  )
Fig2.3.2_supplementary

# Save the plot
ggsave(basedir(paste0("Fig2.3.2_supplementary.pdf")), plot = Fig2.3.2_supplementary,
       width = 18.3,
       height = 17, units = "cm" )

################################################################################
#only top 10 genes per pathway
pathway_results <- lapply(names(pathway_list), function(pathway_name) {
  pathway_genes <- pathway_list[[pathway_name]]
  get_top_genes_pathway(pathway_name, pathway_genes, limmaRes)
})


# Combine the results from all pathways for plotting

combined_results <- map_dfr(pathway_results, "pathway_plot")
combined_results <- combined_results %>% filter(coef %in% koi)%>%
  left_join(ko_flags, by = c("coef" = "genotype", "celltype")) %>%  # Merge with KO flags per cell type
  filter(valid_ko == TRUE)

# Plot the results
# Filter combined_results to remove rows with very low p-value significance
# and missing or unplotable values
filtered_results <- combined_results  # Set size threshold

# Plot with adjusted filtered data
Fig2.3.2_supplementary <- ggplot(filtered_results, aes(
  x = coef, 
  y = ensg,
  color = pmin(2, pmax(-2, logFC)) ,
  size = pmin(5, -log10(adj.P.Val))  # Set size by adjusted p-value
)) +
  geom_point() +  # Create points
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +  # Color scale
  scale_size_continuous(
    range = c(2, 5),  # Set the actual size range from 2 to 30
    breaks = c(1, 2, 4, 5)  # Set specific breaks to create distinct point sizes
  ) +
  labs(
    title = "Differentially Expressed Genesets",
    x = "KOs",
    y = "Genes",
    color = "logFC",
    size = "-log10(adj.P.Val)"
  ) +
  facet_grid(rows = vars(pathway), cols = vars(celltype),
             scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()
  
Fig2.3.2_supplementary

# Save the plot
ggsave(basedir(paste0("Fig2.3.2_supplementary.pdf")), plot = Fig2.3.2_supplementary,
       width = 13 * 0.8 *2,
       height = length(unique(filtered_results$ensg)) * 0.23 )


#
#Fig2.4----------------
#
InDir7  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide")
gsea.res <- read_rds(InDir7("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
gsea.res$coef <- gsub("interaction","",gsea.res$coef )
# Step 2: Summarize to find KOs with at least one valid cell type
valid_ko_summary <- ko_flags %>%
  group_by(genotype) %>%
  summarize(has_valid_celltype = any(valid_ko), .groups = 'drop')


# Step 3: Filter GSEA results based on valid KOs

db = "MSigDB_Hallmark_2020"
pDT <- gsea.res %>%
  filter(coef %in%koi) %>%
  filter(db == "MSigDB_Hallmark_2020") %>%
  left_join(ko_flags, by = c("coef" , "celltype"))%>%
  left_join(summary_df ,c("coef", "celltype"))%>%
  filter(Count > 10) %>%
  filter(valid_ko)# Merge with KO flags

# Step 3 continued: Keep only valid KOs for the specific cell type
pDT <- pDT %>% filter(valid_ko == TRUE, padj < 0.05)
# Step 4: Set alpha value based on validity
#pDT <- pDT %>%
 # mutate(alpha_value = if_else(valid_ko, 1, 0))  # Set alpha based on validity
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5),by=c("coef", "celltype","pathway")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("coef", "celltype","pathway")]$pathway)
# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(c(pw.display.pos, pw.display.neg))

pDT <- pDT[pathway %in% pw.display] 
# Create the plot
#fig-----------
fig2.4 <- ggplot(pDT, aes(x = coef, y = pathway,
                          color = NES,
                          size = pmin(5, -log10(padj)))) +
  geom_point() + 
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name=TeX("log_{2}(FC)")) +
  geom_point(data = pDT[padj < 0.05], shape = 1) +
  scale_size_continuous(
    range = c(0, 3),
    limits = c(0, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  theme_bw() +
  xRot() +
  labs(x="KOs")+
  facet_grid(cols= vars(celltype),scales = "free",space = "free") +  # Create facets for each cell type
  theme(strip.text.y = element_text(angle = 0)) +
  optimized_theme_fig() 
fig2.4
ggsave(basedir("Fig2.4_fgsea.pdf"), w=18,
       h=17, units = "cm")
fig2.4 <- ggplot(pDT, aes(x = coef, y = pathway,
                          color = NES,
                          size = pmin(5, -log10(padj)))) +
  geom_point() + 
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name=TeX("log_{2}(FC)")) +
  geom_point(data = pDT[padj < 0.05], shape = 1) +
  scale_size_continuous(
    range = c(0, 4),
    limits = c(1, 5),
    name=TeX("$-\\log_{10}(p_{adj})$"))+
  theme_bw() +
  xRot() +
  labs(x="KOs")+
  facet_grid(cols= vars(celltype),scales = "free",space = "free") +  # Create facets for each cell type
  theme(strip.text.y = element_text(angle = 0)) +
  optimized_theme_fig() +theme(legend.position = "right")
ggsave(basedir("Fig2.4_legend",db,"selected_pathways",".pdf"), w=16,
       h=9)

#combined ----------------
# Define Row 1: Spacer + Fig2.2
row1 <- (plot_spacer() + Fig2.2) + plot_layout(widths = c(1, 1))

# Define Row 2 Column 1: Fig2.3 + Spacer
row2 <- Fig2.3 + plot_spacer() + plot_layout(widths = c(0.8, 2))

# Define Row 2: Combine row2_col1 with Fig2.4
row3 <- (fig2.4)

# Combine Row 1 and Row 2 in final layout
final_plot <- row1 / row2 / row3 + plot_layout(heights = c(0.8, 0.7,1.4))
final_plot
# Display the final plot
ggsave(basedir("Fig2_06_11.pdf"),plot =final_plot, w=18, h=17)
################################################################################
# Step 1: Ensure the summary of differential genes is per celltype and KO



#################################################################################

#
#fig 2.5----------------
#

InDir6 <- dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")

#
#load data and clean metadata
#

meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))


###########################
dataVoom_Eo.Ba<-read_rds(InDir5("Eo.Ba_dataVoom.rds"))
dataVoom_Mono<-read_rds(InDir5("Mono_dataVoom.rds"))
dataVoom_MkP<-read_rds(InDir5("MkP_dataVoom.rds"))
dataVoom_GMP<-read_rds(InDir5("GMP_dataVoom.rds"))
dataVoom_HSC<-read_rds(InDir5("HSC_dataVoom.rds"))
dataVoom_MEP<-read_rds(InDir5("MEP_dataVoom.rds"))
dataVoom_Gran.<-read_rds(InDir5("Gran._dataVoom.rds"))

dat.list <-list()
for (KO in koi){
  list_of_genes <- c("Oas2","Gbp3","Gvin1","Msmo1","Mthfd2","Idi1","Ccnd1")
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
          gene_data <- meta[meta$celltype == ct,] %>%
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
goi_exp %>% write_rds(basedir("Norm_exp_goi.rds"))
goi_exp <- read_rds(basedir("Norm_exp_goi.rds"))


# Summarize the data to calculate mean scaled expression for each combination of guide, gene, and condition
create_gene_plots <- function(data, gene, KO, ct, limmaRes, significance_threshold = 0.05) {
  # Filter limmaRes for the given gene, KO, and celltype
  filtered_limma <- limmaRes %>%
    filter(gene == gene, coef == KO, celltype == ct)
  
  # Check if the gene is significant in this KO and celltype
  is_significant <- any(filtered_limma$adj.P.Val < significance_threshold)
  
  if (!is_significant) {
    message(paste("Gene", gene, "is not significant in", KO, "KO for celltype:", ct))
    return(NULL)  # Skip the plot if the gene is not significant
  }
  
  # Filter goi_exp data for the gene, KO, and celltype
  filtered_data <- data[data$gene == gene & data$celltype == ct & data$comparison == KO,]
  # Check if filtered_data is empty after filtering
  if (nrow(filtered_data) == 0) {
    message(paste("No data available for gene:", gene, "in celltype:", ct, "for KO:", KO))
    return(NULL)  # Skip plotting if no data is available
  }
  
  ## Check if both tissues have at least 3 replicates
  tissue_counts <- filtered_data %>%
    group_by(tissue) %>%
    summarise(unique_samples = n_distinct(sample1), .groups = "drop")  # Count unique samples by tissue
  
  # Verify the replicate counts for each tissue
  if (any(tissue_counts$unique_samples[tissue_counts$tissue == "ex.vivo"] < 3) ||
      any(tissue_counts$unique_samples[tissue_counts$tissue == "in.vivo"] < 3)) {
    message(paste("Not enough replicates for gene:", gene, "in celltype:", ct, "for both tissues (ex.vivo and in.vivo)"))
    return(NULL)  # Skip plot if either tissue has fewer than 3 replicates
  }
  # Proceed to generate the plot if both tissues have enough replicates and the gene is significant
  p <- ggplot(filtered_data, aes(x = genotype, y = scaled_E, fill = tissue)) + 
    geom_boxplot(
      outlier.shape = NA, 
      position = position_dodge(width = 0.8),  # Adjust width to create space between tissues
      color = "black", 
      size = 0.5,
      coef = Inf
    ) + 
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.2, 
        dodge.width = 0.8  # Adjust dodge width to match the boxplot
      ), 
      alpha = 0.5
    ) +  # Jittered points
    facet_grid(cols = vars(tissue), scales = "free") +
    scale_fill_manual(values = c("#C1A0AC", "#87B1D6"),
                      name = "Experimental model") +
    labs(title = paste0(gene, " (", ct, ")")) +
    xlab(paste0(KO, " KO")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "black"), # Centered, larger plot title
      axis.title = element_text(size = 18, face = "bold", color = "black"),              # Bold, black axis titles
      axis.text = element_text(size = 18, color = "black"),                              # Clear axis text with larger size
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),  # Angled X-axis labels
      legend.title = element_text(size = 16, face = "bold"),  # Bold legend title
      legend.text = element_text(size = 16),  # Clear legend text
      legend.position = "right",  # Legend positioned on the right
      legend.key = element_blank(),
      strip.text  = element_text(size = 18)
    ) 
  
  
  # Get pathway information (if available in the limmaRes table)
  
  dir <-dirout(paste0(base,"/","Interested_genes"))
  # Save the plot
  ggsave(dir(paste0("Fig2.5_", ct, "_", KO, "_", gene, "_", pathway, ".pdf")), 
         plot = p, device = "pdf", width = 10)
}

# Loop over KOs and cell types, getting limmaRes info for significance and plotting
for (comp in koi) {
  for (ct in unique(goi_exp$celltype)) {  
    data <- goi_exp %>%
      filter(comparison == comp, celltype == ct)
    
    # Create gene plots for each gene in the list of genes
    gene_plots <- lapply(unique(list_of_genes), function(gene) {
      create_gene_plots(data, gene, comp, ct, limmaRes)
    })
  }
}
################################################################################
#additional genes

pathway_gene <-  pathway_gene
# Get top genes with their pathway information
top_genes_per_ko_ct <- pathway_gene %>%
  group_by(celltype, coef) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  dplyr::select(ensg, pathway) %>%
  unique()

# Assuming KO_list, cell_types, and limmaRes are already defined
genes <- top_genes_per_ko_ct$ensg %>% unique()


process_ko_celltype <- function(KO, ct, genes) {
  # Get significant genes and their pathways for the current KO and cell type
  significant_genes_info <- pathway_gene %>%
    filter(coef == KO, celltype == ct, ensg %in% genes, adj.P.Val < 0.05) %>%
    dplyr::select(ensg, pathway) %>%
    distinct()
  
  # If no significant genes are found, return NULL
  if (nrow(significant_genes_info) == 0) {
    message(paste("No significant genes found for", KO, "in cell type:", ct))
    return(NULL)
  }
  
  # Extract the genes and pathways
  significant_genes <- significant_genes_info$ensg
  pathways <- significant_genes_info$pathway
  
  # Get the dataVoom object for the current cell type
  dataVoom_ct <- get(paste0("dataVoom_", ct))
  
  # Process each significant gene using map_dfr
  gene_data_list <- map_dfr(significant_genes, function(gene) {
    # Subset metadata and expression values for the specified gene
    if (gene %in% rownames(dataVoom_ct$E)) {
      gene_data <- meta %>%
        filter(celltype == ct) %>%
        mutate(E = dataVoom_ct$E[gene, ]) %>%
        rownames_to_column("samples") %>%
        filter(genotype %in% c(KO, "NTC")) %>%
        mutate(scaled_E = scale(E),
               genes = gene,
               celltype = ct,
               comparison = KO)
      
      # Check for tissue replicates for both tissues
      tissue_counts <- gene_data %>%
        group_by(tissue) %>%
        summarise(unique_samples = n_distinct(samples), .groups = "drop")
      
      # Extract replicate counts for ex.vivo and in.vivo
      ex_vivo_replicates <- tissue_counts %>%
        filter(tissue == "ex.vivo") %>%
        pull(unique_samples)
      
      in_vivo_replicates <- tissue_counts %>%
        filter(tissue == "in.vivo") %>%
        pull(unique_samples)
      
      # Default to 0 if no replicates are found
      if (length(ex_vivo_replicates) == 0) ex_vivo_replicates <- 0
      if (length(in_vivo_replicates) == 0) in_vivo_replicates <- 0
      
      # Check if both tissues have at least 3 replicates
      if (ex_vivo_replicates < 3 || in_vivo_replicates < 3) {
        message(paste("Not enough replicates for gene:", gene, "in celltype:", ct, "for both tissues (ex.vivo and in.vivo)"))
        return(NULL)  # Skip to the next gene if not enough replicates
      }
      
      # Include pathway information in the processed gene data
      pathway_info <- significant_genes_info %>%
        filter(ensg == gene) %>%
        dplyr::select(pathway) %>%
        mutate(genes = gene)  # Ensure we have the gene in the pathway data
      
      # Combine gene data with pathway information
      return(gene_data %>% left_join(pathway_info, by = "genes"))
    } else {
      return(NULL)  # Skip if the gene is not in dataVoom
    }
  })
  
  # Return the combined gene data for the KO and cell type
  return(gene_data_list)
}

# Initialize an empty list to store results
all_results_list <- list()

# Iterate through each KO
for (KO in koi) {
  # Initialize an empty data frame to store results for the current KO
  current_ko_results <- data.frame()
  
  # Iterate through each cell type
  for (ct in unique(limmaRes$celltype)) {
    # Call the processing function and get the gene data
    gene_data <- process_ko_celltype(KO, ct, genes)
    
    # If gene_data is not NULL, combine it with current results
    if (!is.null(gene_data)) {
      current_ko_results <- rbind(current_ko_results, gene_data)
    }
  }
  
  # If current KO results are not empty, add it to the results list
  if (nrow(current_ko_results) > 0) {
    current_ko_results$KO <- KO  # Add the KO information to the data
    all_results_list[[KO]] <- current_ko_results
  }
}

# Combine all results into one data frame
all_results <- do.call(rbind, all_results_list)
# Inspect the combined results
print(head(all_results))
str(all_results)


####################################################################
create_gene_plots <- function(data, gene, KO, ct, limmaRes, significance_threshold = 0.05, pathway) {
  # Filter limmaRes for the given gene, KO, and celltype
  filtered_limma <- limmaRes %>%
    filter(gene == gene, coef == KO, celltype == ct)
  
  # Check if the gene is significant in this KO and celltype
  is_significant <- any(filtered_limma$adj.P.Val < significance_threshold)
  
  if (!is_significant) {
    message(paste("Gene", gene, "is not significant in", KO, "KO for celltype:", ct))
    return(NULL)  # Skip the plot if the gene is not significant
  }
  
  # Filter goi_exp data for the gene, KO, and celltype
  # Filter data for the gene, KO, and celltype
  filtered_data <- data[data$gene == gene & data$celltype == ct & data$comparison == KO,]
  
  # Check tissue replicates
  tissue_counts <- filtered_data %>%
    group_by(tissue) %>%
    summarise(unique_samples = n_distinct(sample1), .groups = "drop")
  
  # Verify replicate counts
  if (any(tissue_counts$unique_samples[tissue_counts$tissue == "ex.vivo"] < 3) ||
      any(tissue_counts$unique_samples[tissue_counts$tissue == "in.vivo"] < 3)) {
    message(paste("Not enough replicates for gene:", gene, "in celltype:", ct, "for both tissues"))
    return(NULL)
  }
  
  # Skip if no data or valid values for faceting
  if (nrow(filtered_data) == 0 || all(is.na(filtered_data$tissue))) {
    message("No valid data available for plotting.")
    return(NULL)
  }
  # Proceed to generate the plot if both tissues have enough replicates and the gene is significant
  p <- ggplot(filtered_data, aes(x = genotype, y = scaled_E, fill = tissue)) + 
    geom_boxplot(
      outlier.shape = NA, 
      position = position_dodge(width = 0.8),  # Adjust width to create space between tissues
      color = "black", 
      size = 0.5,
      coef = Inf
    ) + 
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.2, 
        dodge.width = 0.8  # Adjust dodge width to match the boxplot
      ), 
      alpha = 0.5
    ) +  # Jittered points
    facet_grid(cols = vars(tissue), scales = "free") +
    scale_fill_manual(values = c("#C1A0AC", "#87B1D6"),
                      name = "Experimental model") +
    labs(title = paste0(gene, " (", ct, ")")) +
    xlab(paste0(KO, " KO")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "black"), # Centered, larger plot title
      axis.title = element_text(size = 18, face = "bold", color = "black"),              # Bold, black axis titles
      axis.text = element_text(size = 18, color = "black"),                              # Clear axis text with larger size
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),  # Angled X-axis labels
      legend.title = element_text(size = 16, face = "bold"),  # Bold legend title
      legend.text = element_text(size = 16),  # Clear legend text
      legend.position = "right",  # Legend positioned on the right
      legend.key = element_blank(),
      strip.text  = element_text(size = 18)
    ) 
  
  # Get pathway information (if available in the limmaRes table)
  pathway <- all_results %>% filter(genes == gene &
                                      celltype == ct &
                                      comparison == KO) %>% pull(pathway)%>%
    unique()
  pathway <- pathway[1]
  dir <-dirout(paste0(base,"/",pathway))
  # Save the plot
  ggsave(dir(paste0("Fig2.5_", ct, "_", KO, "_", gene, ".pdf")), 
         plot = p, device = "pdf", width = 10)
}

# Loop over KOs and cell types, getting limmaRes info for significance and plotting
# Loop over KOs and cell types, plotting
for (KO in koi) {
  for (ct in unique(all_results$celltype)) {  
    for (pathway in unique(all_results$pathway)){
      # Filter all_results for the current KO and celltype
      data <- all_results %>% filter(comparison == KO,
                                     celltype == ct,
                                     pathway == pathway)
      
      # Create gene plots for each unique gene in all_results for the specific KO and celltype
      unique_genes <- unique(data$genes)  # Extract unique genes for the current cell type
      gene_plots <- lapply(unique_genes, function(gene) {
        create_gene_plots(data, gene, KO, ct, limmaRes, pathway)
      })
    }
  }
}
#####################################################################

library(dplyr)
library(ggplot2)
library(hexbin)
library(tidyr)
# Define genes to highlight
limmaRes_NTC <- read_rds(dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")("limma_perCTex.vivovsin.vivo.rds"))%>%
  mutate(ensg = genes)%>%
 dplyr::select(-genes)
head(limmaRes_NTC)

# Thresholds for defining DEGs
adj_p_cutoff <- 0.05
logfc_cutoff <- 2

Brd9_genes <- limmaRes %>%
  filter(coef == "Brd9")%>%
  mutate(gr=group)%>%
  dplyr::select(celltype, ensg, gr)
# Load necessary libraries

# Example list of genes for each cell type (replace with actual data)



# Convert named list to a data frame to facilitate joining

# Merge with the main data
volcano_data <- limmaRes_NTC %>%
  left_join(Brd9_genes, by = c("ensg","celltype")) %>%
  mutate(
    highlight = ifelse(gr != "n.s", "yes", "no")
  )

# Create the volcano plot with cell type highlighting
# Load necessary libraries
library(ggplot2)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create the volcano plot with hex bins, highlighted genes, and a vertical line at logFC = 1
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  # Hex bin layer with grey-to-blue color gradient for density
  geom_hex(bins = 30, aes(fill = ..count..)) +
  scale_fill_gradient(low = "grey", high = "blue") +
  
  # Overlay red points for highlighted genes
  geom_point(data = filter(volcano_data, highlight == "yes"), color = "red", alpha = 0.8, size = 1.5) +
  
  # Labels and theme adjustments
  labs(
    title = "Volcano Plot of logFC vs adj.P.Val",
    x = "Log Fold Change (logFC)",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme_minimal() +
  facet_wrap(vars(celltype), scales = "free") +
  theme(legend.position = "right") +
  
  # Optional significance threshold line
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  
  # Vertical line at logFC = 1
  geom_vline(xintercept = c(-1,1), linetype = "dotted", color = "black")

# Save the plot to a PDF file
ggsave(basedir("test.pdf"))
#
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create the volcano plot with two hex bin layers for density
ggplot(volcano_data %>%, aes(x = logFC, y = -log10(adj.P.Val))) +
  # Hex bin layer for all genes, with grey-to-blue color gradient for density
 # geom_hex(bins = 10, aes(fill = ..count..), color = NA) +
  #scale_fill_gradient(low = "grey", high = "blue", name = "All Genes Density") +
  
  # Hex bin layer for highlighted genes only, with grey-to-red color gradient
  geom_hex(data = filter(volcano_data, highlight == "yes"), bins = 50, aes(fill = ..count..), color = NA) +
  scale_fill_gradient(low = "grey", high = "red", name = "Highlighted Genes Density") +
  
  # Labels and theme adjustments
  labs(
    title = "Volcano Plot of logFC vs adj.P.Val",
    x = "Log Fold Change (logFC)",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme_minimal() +
  facet_wrap(vars(celltype), scales = "free") +
  theme(legend.position = "right") +
  
  # Optional significance threshold line
  #geom_hline(yintercept = -log10(0.05), 
     #        linetype = "dashed", 
      #       color = "blue") +
  
  # Vertical lines at logFC = -1 and 1
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dotted",
             color = "black")

# Save the plot to a PDF file
ggsave(basedir("test1.pdf"))
##########################
#interferon genes in volcanoplot
ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
  filter(L1=="ISG_Core")%>%
  pull(value)
# Modify NTC_ISG to include a 'highlighted' column based on the conditions
NTC_ISG <- limmaRes_NTC %>%
  group_by(celltype) %>%
  mutate(highlighted = ifelse(ensg %in% ISG_core & group != "n.s", "yes", "no")) %>%
  dplyr::select(ensg,celltype,highlighted)%>%
  ungroup()

volcano_data <- limmaRes %>%
  filter(coef %in% koi)%>%
  inner_join(ko_flags, by = c("coef", "celltype"))%>%
  filter(valid_ko == 'TRUE')%>%
  inner_join(summary_df, by = c("coef", "celltype"))%>%
  filter(Count >10)
# Join NTC_ISG with volcano_data to get celltype-specific genes to highlight
volcano_data_highlight <- volcano_data %>%
  inner_join(NTC_ISG, by = c("celltype","ensg")) %>%
  filter(highlighted == "yes")  # Filter only genes in the highlight list

# Volcano plot
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  
  # Hex bin layer for all genes, with grey-to-blue color gradient for density
  geom_hex(bins = 100, aes(fill = ..count..), color = NA) +
  scale_fill_gradient(low = "grey", high = "blue", name = "Counts") +
  
  # Highlight points from NTC_ISG in red for each celltype
 # geom_point(data = volcano_data_highlight,
   #         aes(x = logFC, y = -log10(adj.P.Val)), color = NA)+
  
  geom_point(data = volcano_data_highlight,
             aes(x = logFC, y = -log10(adj.P.Val)),
             color = "red", size = 0.5) +
  
  # Labels and theme adjustments
  labs(
    title = "Volcano Plot of logFC vs adj.P.Val",
    x = "Log Fold Change (logFC)",
    y = "-log10(Adjusted P-Value)"
  ) +
  theme_minimal() +
  facet_wrap(coef ~celltype, scales = "free") +
  theme(legend.position = "right") +
  
  # Optional significance threshold line
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  
  # Vertical lines at logFC = -1 and 1
  geom_vline(xintercept = c(-1, 1),
             linetype = "dotted",
             color = "black")

# Save the plot to a PDF file
ggsave(basedir("ISG_core_in_KO.pdf"))
#waste----------------------
library(ggplot2)
library(dplyr)
library(hexbin)

# Thresholds for defining DEGs
adj_p_cutoff <- 0.05
logfc_cutoff <- 2

# Step 1: Get the list of DEGs for each KO-celltype combination
highlighted_genes_list <- limmaRes %>%
  filter(adj.P.Val < adj_p_cutoff, abs(logFC) > logfc_cutoff) %>%
  group_by(celltype, coef) %>%
  summarize(genes_to_highlight = list(unique(ensg)), .groups = 'drop')

# Step 2: Join limmaRes_NTC with the highlighted genes list and add a highlight flag
highlighted_genes <- limmaRes_NTC %>%
  left_join(highlighted_genes_list, by = c("celltype")) %>%
  mutate(highlight = ensg %in% unlist(genes_to_highlight)) # Flag genes for highlighting

# Step 3: Create a hexbin volcano plot with separate layers for NTC and KO-celltype specific DEGs
library(ggplot2)
library(dplyr)
library(hexbin)

library(ggplot2)
library(dplyr)
library(hexbin)

volcano_plot_per_KO_celltype <- ggplot() +
  # Layer 1: Base hexbin for NTC genes, with count-based color gradient (grey to blue)
  
  
  # Layer 2: Overlay hexbin plot for KO-celltype specific DEGs with a semi-transparent red
  stat_binhex(
    data = highlighted_genes %>% filter(highlight),
    aes(x = logFC, y = -log10(P.Value), fill = ..count..),
    bins = 50,
    color = NA,
    alpha = 0.4  # Semi-transparent red to distinguish this layer
  ) +
  
  # Single color scale for both layers, going from grey to blue
  scale_fill_gradient(low = "grey", high = "blue", name = "Count") +
  
  # Facet by each KO and celltype combination
  facet_wrap(celltype ~ coef, scales = "free") +
  theme_minimal() +
  labs(
    title = "Faceted Volcano Plot of NTC and KO-Specific DEGs by Cell Type",
    x = "log2 Fold Change",
    y = "-log10 P-value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    strip.text = element_text(size = 10)  # Adjust facet label size if needed
  )

# Display the plot
print(volcano_plot_per_KO_celltype)

# Optionally save the plot
ggsave(basedir("Volcano_Plot_Per_KO_Celltype_Hexbin.pdf"), plot = volcano_plot_per_KO_celltype, width = 20, h=30)

       
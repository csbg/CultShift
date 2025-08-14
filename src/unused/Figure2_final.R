#Figure2
#Invivo vs exvivo 
#Aarathy
#24-09-2024
###############################################################################
#load libraries and functions--------------------------------------------------
###############################################################################
source("src/00_init.R")
source("src/Ag_Optimized_theme.R")
library(tidyverse)
library(enrichR)
#library(purrr)
#library(gridExtra)
library(pheatmap)
library(ComplexHeatmap)
library("scales")
#library(circlize)  # For color palettes
########################
#directories ------
########################
base<-"Figure2"
basedir<-dirout("Figure2")
################################################################################
InDir5 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")


################################################

meta <- read_rds(InDir5("meta.rds"))
colnames(meta)
selected_KOs <- meta %>%
  group_by(genotype, tissue, celltype) %>%                  # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample), .groups = 'drop') %>% # Count distinct guides for each group
  filter(num_sample >= 3) %>%                               # Keep groups with at least 3 guides
  pull(genotype) %>% unique()


selected_KOs1 <- meta %>%
  group_by(genotype, tissue, celltype) %>%                  # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample), .groups = 'drop') %>% # Count distinct samples for each group
  pivot_wider(names_from = tissue, values_from = num_sample, values_fill = 0) %>% # Spread tissue to separate columns (in.vivo and ex.vivo)
  filter(in.vivo >= 3 & ex.vivo >= 3) %>%                   # Keep groups with at least 3 samples in both tissues
  group_by(genotype) %>%                                    # Regroup by genotype
  filter(any(in.vivo >= 3 & ex.vivo >= 3)) %>%              # Keep genotypes that have at least one celltype with 3+ samples in both tissues
  pull(genotype) %>% unique()


###############################
#enrichr database-----
###############################

enr.terms <- enrichrGetGenesets(ENRICHR.DBS)
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

###############################################################################
#Fig2.1 Correlation between logFC ---------------------------------------------
###############################################################################
Indir1<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_correlation/")
correlation_matrix_data <- read_rds(Indir1("correlation_ex.vivo_vs_in.vivo.rds"))
deg_plot_data <- read_rds(Indir1("DEGs_per_tissue.rds"))
correlation_matrix_no_na <- correlation_matrix_data %>%
  replace(is.na(.), 0) %>%
  replace(is.nan(.), 0) %>%
  replace(is.infinite(.), 0)
hc_cols <- hclust(dist(t(correlation_matrix_no_na)), method = "ward.D2")
column_order <- colnames(correlation_matrix_data)[hc_cols$order]

correlation <- correlation_matrix_data%>%as_tibble(rownames = "celltype")%>%
  pivot_longer(cols = colnames(correlation_matrix_data),
               names_to ="genotype",
               values_to ="correlation",
               values_drop_na =F)

deg_plot_data <- deg_plot_data %>%
  group_by(genotype, celltype) %>%   # Group by genotype and celltype
  filter(num_degs == max(num_degs, na.rm = TRUE)) %>%   # Keep the row with max num_degs in each group
  select(-condition)
correlation_deg <- inner_join(deg_plot_data,correlation,by = c("celltype","genotype"))

correlation_deg$genotype <- factor(correlation_deg$genotype, levels = column_order)
genotypes <- unique(correlation_deg$genotype)
celltypes <- unique(correlation_deg$celltype)
Fig2.1 <- correlation_deg %>% 
  filter(genotype %in% selected_KOs) %>%
  ggplot() +
  geom_point(aes(
    x = genotype,
    y = celltype,
    size = log10(num_degs),
    fill = correlation  # Use 'fill' for inside color
  ),
  shape = 21,           # Use shape 21 to enable fill and color
  color = "black",       # Black outline
  stroke = 0.5          # Adjust the width of the outline
  ) +
  scale_fill_gradient2(low = "#4C889C", 
                       mid = "white", 
                       high = "#D0154E") +  # Use 'fill' for inside color
  labs(x = "KOs",
       y = "Celltype") +
  optimized_theme()

Fig2.1

ggsave(basedir("Fig2.1.pdf"),
       w=length(genotypes)*0.3+3,
       h=length(celltypes)*0.3+3)
#############################
#Fig2.1.1

InDir2 <- dirout("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/")
InDir3 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")

gsea.res <- read_rds(InDir2("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))
limmaRes<- read_rds(InDir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))%>%
  mutate(coef = gsub("interaction","",coef))

KO_list <- correlation_deg %>% filter(correlation < 0.5, num_degs >= 10) %>%
  pull(genotype)%>%
  unique()

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
#kos of interest (koi)
coefficients <- intersect(coefficients, selected_KOs)
koi <- KO_list[KO_list %in% coefficients ]

Fig2.1.2 <- summary_df %>%
  filter(Count > 10, coef %in% koi) %>%
  # Plot with log10 scale on y-axis
  ggplot(aes(x = coef, y = Regulation, color = Regulation)) +
  geom_point(aes(size=log10(Count))) +    
  #scale_size_continuous(limits = c(10,2000),breaks = c(10,50,100,500,1000,2000)) + 
  scale_color_manual(values = c("Downregulated" =  "#4C889C", "Upregulated" = "#D0154E")) +
  facet_grid(rows = vars(celltype)) +
  labs(title="No: of genes with interaction effect of culture condition in KOs",
       x="KOs")+
  guides(
    
    size = guide_legend(title = "log10(Count)")
  )+
  #scale_y_log10() +  # Apply log10 scale to the y-axis
  optimized_theme()

ggsave(basedir("Fig2.1.2.pdf"),plot=Fig2.1.2, h=7)
#from 12_enrichr
###############################################################################
#Fig2.3
InDir4 <- dirout("Figure1")
genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
colnames(genes_fig_1) <- c("pathways","ensg")
mTORC1_or_Cholesterol <- c("Idi1", "Cyp51","Stard4", "Scd2","Mthfd2","Sqle","Fads2","Dhcr24",
                           "Hmgcs1","Ldlr","Plscr1","Acat1",
                           "Acat2","Msmo1")
cholesterol_df <- data.frame(
  # Cholesterol genes
  pathways = rep("mTORC1_or_Cholesterol", length(mTORC1_or_Cholesterol)),
  ensg = mTORC1_or_Cholesterol # Cholesterol pathway
)
genes_fig_1 <- rbind(genes_fig_1, cholesterol_df) %>% distinct()

filtered_genes = limmaRes %>%
  filter(ensg %in% genes_fig_1$ensg, coef %in% koi) %>%
  group_by(celltype, coef)%>%
  filter(coef %in% coefficients)%>%
  merge(genes_fig_1, by = "ensg" )# Group by celltype and coefficient


filtered_genes$pathways <- recode(filtered_genes$pathways,
                                  "ISG_core"= "ISG core",
                                  "mTORC1_or_Cholesterol" = "mTORC1/Cholesterol")
fig2.3 <- ggplot(filtered_genes, aes(x = coef, y = ensg,
                                       color = logFC,
                                       size = pmin(30,-log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(
    range = c(2, 6),  # Set the actual size range from 2 to 30
    breaks = c(1,2,4, 5)  # Set specific breaks to create distinct point sizes
  ) +
  labs(title = "Differentially expressed genesets",
       x = "KOs",
       y = "Genes",
       color = "logFC",
       size = "-log10(adj.P.Val)") +
  facet_grid(rows = vars(pathways), cols = vars(celltype),scales = "free_y", space = "free") +
  theme_bw() +
  optimized_theme()
fig2.3 
ggsave(basedir(paste0("fig2.3.pdf")),plot = fig2.3, 
       width = 13*0.3*6,
       height = length(unique(filtered_genes$ensg))*0.25)
################################################################################
#fig 2.4
####################################################################
#***Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01***#
####################################################################
#limma with all 
#data
#
#Aarathy
###############
###############

#####################################################################

InDir6<-dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")

#####################################################################
#load data and clean metadata
################################################

meta <- read_rds(InDir5("meta.rds"))
colnames(meta)
selected_KOs <- meta %>%
  group_by(genotype, tissue, celltype) %>%                  # Group by genotype, tissue, and celltype
  summarize(num_sample = n_distinct(sample), .groups = 'drop') %>% # Count distinct guides for each group
    filter(num_sample >= 3) %>%                               # Keep groups with at least 3 guides
  distinct(genotype) 


################################################
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))


###########################
dataVoom_Eo.Ba<-read_rds(InDir5("Eo.Ba_dataVoom.rds"))
dataVoom_Mono<-read_rds(InDir5("Mono_dataVoom.rds"))
dataVoom_MkP<-read_rds(InDir5("MkP_dataVoom.rds"))
dataVoom_GMP<-read_rds(InDir5("GMP_dataVoom.rds"))
dataVoom_HSC<-read_rds(InDir5("HSC_dataVoom.rds"))
dataVoom_MEP<-read_rds(InDir5("MEP_dataVoom.rds"))

dat.list <-list()
for (KO in coefficients){
  list_of_genes <- c("Tbp", "Oas2","Gbp3","Gvin1","Msmo1","Hmgcs1","Stard4","Mthfd2","Idi1")
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
            rownames_to_column("sample1") %>%
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
# Function to create plots for each gene
create_gene_plots <- function(data, gene, KO, ct) {
  filtered_data <- data[data$gene == gene & data$celltype == ct,]
  
  # Check if there's any data to plot
  if (nrow(filtered_data) == 0) {
    message(paste("No data for gene:", gene, "and celltype:", ct))
    return(NULL)  # Skip the plot if no data
  }
  
  # Generate the plot if data is present
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
    theme_bw()+
    theme(
      # Text Elements
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"), # Centered, larger plot title
      axis.title = element_text(size = 16, face = "bold", color = "black"),              # Bold, black axis titles
      axis.text = element_text(size = 14, color = "black"),                              # Clear axis text with larger size
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"),  # Angled X-axis labels for better readability
      
      # Legend
      legend.title = element_text(size = 14, face = "bold"),            # Bold legend title
      legend.text = element_text(size = 12),                            # Clear legend text
      legend.position = "right",                                        # Legend positioned on the right
      legend.key = element_blank())+
    optimized_theme()
  
  # Get pathway information
  pathway <- genes_fig_1 %>% filter(ensg == gene) %>% pull(pathways)
  
  # Save the plot
  ggsave(basedir(paste0("fig2.4_", ct, "_", KO, "_", gene, "_", pathway, ".pdf")), 
         plot = p, device = "pdf")
}

# Loop over KOs and cell types
for (comp in c("Brd9", "Cebpa", "Bcl11a", "Smc3","Rcor1", "Wdr82", "Spi1", "Phf10","Kmt2d")) {
  for (ct in unique(goi_exp$celltype)) {  
    data <- goi_exp %>%
      filter(comparison == comp, celltype == ct)
    
    # Create gene plots for each gene in the list of genes
    gene_plots <- lapply(unique(list_of_genes), function(gene) {
      create_gene_plots(data, gene, comp, ct)
    })
  }
}

############################################
InDir7  <-  dirout("Ag_ScRNA_12_Pseudobulk_FGSEA_per_celltype_guide.R")
gsea.res <- read_rds(InDir7("fgsea_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
gsea.res$coef <- gsub("interaction","",gsea.res$coef )
# pathways <- list(
#   "Cholesterol Homeostasis",
#   "mTORC1 Signaling",
#   "Interferon Alpha Response",
#   "Interferon Gamma Response",
#   "Inflammatory Response",
#   "p53 Pathway",
#   "Myc Targets V1",
#   "IL−6/JAK/STAT3 Signaling",
#   "Oxidative Phosphorylation",
#   "TGF−beta Signaling",
#   "Unfolded Protein Response",
#   "E2F Targets",
#   "Protein Secretion"
# )
#pDT  <-  gsea.res[db == dbx & pathway %in% pathways]
#pathways <- names(enr.terms$MSigDB_Hallmark_2020)
dbx <- "MSigDB_Hallmark_2020"

pDT <- gsea.res 
# Step 3: Filter GSEA results based on valid KOs
pDT <- gsea.res %>%
  filter(coef %in% koi) %>%
  filter(db == "MSigDB_Hallmark_2020") %>%
  left_join(ko_flags, by = c("coef" = "genotype", "celltype"))  # Merge with KO flags

# Step 3 continued: Keep only valid KOs for the specific cell type
pDT <- pDT %>% filter(valid_ko == TRUE)
# Step 4: Set alpha value based on validity
pDT <- pDT %>%
  mutate(alpha_value = if_else(valid_ko, 1, 0.2))  # Set alpha based on validity

pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5),by=c("coef", "celltype","pathway")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("coef", "celltype","pathway")]$pathway)
# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(c(pw.display.pos, pw.display.neg))
pDT <- pDT[pathway %in% pw.display]
# Create the plot
# Step 4: Set alpha value based on validity
pDT <- pDT %>%
  mutate(alpha_value = if_else(valid_ko, 1, 0.2))  # Set alpha based on validity

# Create the plot
ggplot(pDT, aes(x = coef, y = pathway, color = NES, size = pmin(5, -log10(padj)), alpha = alpha_value)) +
  geom_point() + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  geom_point(data = pDT[padj < 0.05], shape = 1) +
  scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +
  theme_bw(12) +
  xRot() +
  facet_wrap(vars(celltype)) +  # Create facets for each cell type
  theme(strip.text.y = element_text(angle = 0)) +
  optimized_theme() +
  labs(alpha = "Validity of KO")


ggsave(basedir("Fig2.01_GSEA_plot_","selected_pathways",dbx,".png"), w=16,
       h=length(unique(pDT$pathway)) * 0.2*2 + 3, limitsize = FALSE)


################################################################################
# Step 1: Ensure the summary of differential genes is per celltype and KO


# Step 1: Summarize differential genes for each celltype and KO
summary_df <- limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count") %>%
  group_by(celltype, coef) %>%
  summarise(total_differential = sum(Count), .groups = 'drop')

# Step 2: Filter KOs per celltype where total differential genes >= 10
valid_KOs_per_celltype <- summary_df %>%
  filter(total_differential >= 10) %>%
  select(celltype, coef)  # Keep only relevant columns for joining

# Step 3: Filter GSEA results to include only valid KOs per celltype
pDT_filtered <- gsea.res %>%
  inner_join(valid_KOs_per_celltype, by = c("celltype", "coef")) %>%  # Ensure celltype-specific KO filtering
  filter(padj < 0.05) %>%                                             # Significant GSEA pathways
  filter(db == "MSigDB_Hallmark_2020")                                # Only keep MSigDB Hallmark pathways

# Step 4: Select top 5 upregulated and downregulated pathways per KO, celltype, and pathway
pw.display.pos <- unique(pDT_filtered[order(-NES)][, head(.SD, n = 5), by = c("coef", "celltype", "pathway")]$pathway)
pw.display.neg <- unique(pDT_filtered[order(NES)][, head(.SD, n = 5), by = c("coef", "celltype", "pathway")]$pathway)

# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(c(pw.display.pos, pw.display.neg))

# Step 5: Filter the GSEA result for the selected pathways
pDT_final <- pDT_filtered %>%
  filter(pathway %in% pw.display)

# Step 6: Plot the GSEA results for the selected pathways and KOs
ggplot(pDT_final, aes(x = coef, y = pathway, color = NES, size = pmin(5, -log10(padj)))) +
  geom_point() + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +  # NES color gradient
  geom_point(data = pDT_final[padj < 0.05], shape = 1) +              # Highlight significant points
  scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +          # Limit the size of points
  theme_bw(base_size = 12) + 
  xRot() + 
  labs(x= "KOs")+
# Rotate x-axis labels
  facet_wrap(vars(celltype)) +                                         # Facet by celltype
  theme(strip.text.y = element_text(angle = 0)) +                      # Customize facet text angle
  optimized_theme()

# Step 7: Save the plot
ggsave(basedir("Fig2.01_GSEA_plot_", "selected_pathways", dbx, ".pdf"), 
       width = length(valid_KOs_per_celltype$coef)*.3*2, 
       height = length(unique(pDT_final$pathway)) * 0.2 * 2 + 3, 
       limitsize = FALSE)
#################################################################################
#plot per celltype
# Step 1: Summarize differential genes for each celltype and KO
summary_df <- limmaRes %>%
  group_by(celltype, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count") %>%
  group_by(celltype, coef) %>%
  summarise(total_differential = sum(Count), .groups = 'drop')

# Step 2: Filter KOs per celltype where total differential genes >= 10
valid_KOs_per_celltype <- summary_df %>%
  filter(total_differential >= 10) %>%
  select(celltype, coef)  # Keep only relevant columns for joining

# Step 3: Loop over each celltype to create a plot per celltype
celltypes <- unique(valid_KOs_per_celltype$celltype)

for (ct in celltypes) {
  
  # Step 4: Filter GSEA results for the current celltype
  pDT_filtered <- gsea.res %>%
    filter(celltype == ct) %>%  # Filter for the current celltype
    inner_join(valid_KOs_per_celltype, by = c("celltype", "coef")) %>%  # Keep only valid KOs
    filter(padj < 0.05) %>%                                               # Significant GSEA pathways
    filter(db == "MSigDB_Hallmark_2020")                                  # Only keep MSigDB Hallmark pathways
  
  # Step 5: Select top 5 upregulated and downregulated pathways per KO
  pw.display.pos <- unique(pDT_filtered[order(-NES)][, head(.SD, n = 5), by = c("coef", "pathway")]$pathway)
  pw.display.neg <- unique(pDT_filtered[order(NES)][, head(.SD, n = 5), by = c("coef", "pathway")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  
  # Step 6: Filter the GSEA result for the selected pathways
  pDT_final <- pDT_filtered %>%
    filter(pathway %in% pw.display)
  
  # Step 7: Plot the GSEA results for the current celltype
  plot <- ggplot(pDT_final, aes(x = coef, y = pathway, color = NES, size = pmin(5, -log10(padj)))) +
    geom_point() + 
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +  # NES color gradient
    geom_point(data = pDT_final[padj < 0.05], shape = 1) +              # Highlight significant points
    scale_size_continuous(range = c(0, 5), limits = c(0, 5)) +          # Limit the size of points
    theme_bw(base_size = 12) + 
    xRot() +                                                            # Rotate x-axis labels
    labs(title = paste("GSEA Plot for Celltype:", ct)) +           # Add title with celltype name
    theme(strip.text.y = element_text(angle = 0)) +                      # Customize facet text angle
    optimized_theme()
  
  # Step 8: Save the plot for the current celltype
  ggsave(basedir(paste0("Fig2.01_GSEA_plot_", ct, "_", dbx, ".pdf")), 
         plot = plot, width = 10, 
         height = length(unique(pDT_final$pathway)) * 0.2 * 2 , 
         limitsize = FALSE)
}


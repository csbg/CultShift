#Figure2_initial
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
base<-"Figure2_initial"
basedir<-dirout("Figure2_initial")
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
fig2.1 <- correlation_deg %>% 
  ggplot()+
  geom_point(aes(x=genotype,
                 y= celltype,
                 size =num_degs,
                 color=correlation))+
  scale_color_gradient2(low = "#4C889C",#muted("blue"),
                        mid = "white",
                        high = "#D0154E")+#muted("red"))+
  scale_size_continuous(breaks = c(10,100,500,1000,5000))+
  optimized_theme()
fig2.1
ggsave(basedir("fig2.1.pdf"),
       w=length(genotypes)*0.3+3,
       h=length(celltypes)*0.3+3)
unique(deg_plot_data$deg_category)
KO_list <- correlation_deg %>% filter(correlation < 0.5, num_degs >= 100) %>%
  pull(genotype)%>%
  unique()

  
 

##################
#what pathways
#fgsea
InDir2 <- dirout("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype_guide/")
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
#kos of interest (koi)
koi <- KO_list[KO_list %in% coefficients]
length(koi)
fig2.1.2<-summary_df %>%
  filter(Count > 10, coef %in% koi) %>%
  # Plot with log10 scale on y-axis
  ggplot(aes(x = gsub("interaction", "", coef), y = Count, fill = Regulation)) +
  geom_col(position = position_dodge(preserve = "single")) +    # Separate bars for Upregulated and Downregulated
  scale_fill_manual(values = c("Downregulated" =  "#4C889C", "Upregulated" = "#D0154E")) +
  facet_grid(rows = vars(celltype)) +
  labs(title="No: of genes with interaction effect of culture condition in KOs",
       x="KOs")+
  scale_y_log10() +  # Apply log10 scale to the y-axis
  optimized_theme()

summary_df %>%
  filter(Count > 10, coef %in% koi) %>%
  # Plot with log10 scale on y-axis
  ggplot(aes(x = coef, y = Regulation, color = Regulation)) +
  geom_point(aes(size=pmin(2000,Count))) +    
  scale_size_continuous(limits = c(10,2000),breaks = c(10,50,100,500,1000,2000)) + 
  scale_color_manual(values = c("Downregulated" =  "#4C889C", "Upregulated" = "#D0154E")) +
  facet_grid(rows = vars(celltype)) +
  labs(title="No: of genes with interaction effect of culture condition in KOs",
       x="KOs")+
  #scale_y_log10() +  # Apply log10 scale to the y-axis
  optimized_theme()
ggsave(basedir("fig2.1.2.pdf"),plot=fig2.1.2)
#from 12_enrichr
################################################################################
#Fig2.3 ------------------------------------------------------------------------
################################################################################
dbx <- "MSigDB_Hallmark_2020"
pathways <- list(
  "Cholesterol Homeostasis",
  "mTORC1 Signaling",
  "Interferon Alpha Response",
  "Interferon Gamma Response",
  "Inflammatory Response")
pDT <- gsea.res[gsea.res$db == dbx & size > 10 & #padj < 0.05 &
                  coef %in% koi & pathway %in% pathways,]

## Splitting the task to handle both ends of the NES spectrum - positive and negative
pw.display.pos <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("celltype", "coef")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5), by=c("celltype", "coef")]$pathway)

# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(c(pw.display.pos, pw.display.neg))
pDT <- pDT[pathway %in% pw.display]

if (nrow(pDT) > 0) {
  # Proceed with your operation
  pDT <- hierarch.ordering(pDT, "celltype", "pathway", "NES", TRUE)
  
  #pDT <- pDT %>% filter(padj < 0.05)
  # Plotting
fig2.3 <- ggplot(pDT, aes(x=coef, y=celltype, color=NES, size=pmin(5, -log10(padj)))) +
    scale_color_gradient2(low="blue", mid="white", high="red") +
    geom_point(data=pDT)+#[padj < 0.05]) +
    scale_size_continuous(range=c(0, 5), limits = c(0, 5)) +
    theme_bw(12) +
    xRot() +
    facet_wrap(vars(pathway)) +
    labs(y = "Interaction-KO", x="Celltypes") +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))+
  optimized_theme()
  
  # Save the plot
  #dat <- dirout(paste0(base,"/", dbx))
  ggsave(basedir("Fig2.3.pdf"), plot = fig2.3,
         height = length(unique(pDT$celltype)) * 0.25 + 4,
         width = length(koi)* 0.35 + 14)
} else {
  warning("Data table is empty. Skipping operation.")
} 
################################################################################



################################################################################
#genes from pathways -----------------------------------------------------------
################################################################################
inDir3 <- dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype_guide/")
limmaRes<- read_rds(inDir3("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
# Clean and split the genes from the leadingEdge column
pathway_genes_list <- pDT %>% filter(padj < 0.05)%>%
  # Group by the pathway column
  group_by(pathway) %>%
  # Summarize to extract and clean the leadingEdge genes
  summarise(
    unique_genes = list(
      leadingEdge %>%
        # Replace extraneous characters like parentheses, quotes, and "c"
        str_replace_all("[\\(\\)\"c]", "") %>%
        # Split by commas (with optional space after)
        str_split(",\\s*") %>%
        # Flatten the nested lists into one vector
        unlist() %>%
        # Trim whitespace from each gene
        str_trim() %>%
        # Keep only unique gene names
        unique()
    )
  ) %>%
  # Convert the result to a named list with pathways as keys and gene lists as values
  deframe()

# Print the cleaned list of genes per pathway
pathway_genes_list



# Apply the function to each unique `db` using purrr::walk
IFN_response <- union(pathway_genes_list$`Interferon Alpha Response`,
                      pathway_genes_list$`Interferon Gamma Response`)
#Option1
pathways <- list(
  "Cholesterol Homeostasis" = pathway_genes_list$`Cholesterol Homeostasis`,
  "mTORC1 Signaling" = pathway_genes_list$`mTORC1 Signaling`,
  "IFN_response" = union(pathway_genes_list$`Interferon Alpha Response`,
                         pathway_genes_list$`Interferon Gamma Response`),
  "Inflammatory Response" = pathway_genes_list$`Inflammatory Response`)

#" TNF-alpha Signaling via NF-kB",
#"p53 Pathway",
#"Myc Targets V1",
#"IL−6/JAK/STAT3 Signaling",
# "Oxidative Phosphorylation",
# "TGF−beta Signaling",
#"Unfolded Protein Response",
#"E2F Targets",
#"Protein Secretion"
#Option2
# pathways <- list(
#   ISG_core = read.delim(paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Mostafavi_Cell2016.tsv"))%>%
#     filter(L1=="ISG_Core")%>%pull(value),
#   mTORC1_or_Cholesterol = union(enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`, 
#                                 enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`),
#   IFN_response = union(enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`, 
#                        enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`))

# Define the function to get top genes for a given pathway
get_top_genes <- function(pathway_name,
                          pathway_genes,
                          limma_results, 
                          logFC_threshold = 1,
                          pval_threshold = 0.05, 
                          top_n = 3) {
  # Filter limma results for the pathway genes
  filtered_genes <- limma_results %>%
    filter(group != "n.s") %>%
    group_by(celltype, coef) %>%   # Group by celltype and coefficient
    filter(adj.P.Val < pval_threshold) %>%
    filter(toupper(ensg) %in% toupper(pathway_genes)) %>%
    arrange(desc(abs(logFC))) %>%  # Sort by absolute logFC in descending order
    slice_head(n = top_n) %>%      # Get the top n genes per group (celltype)
    pull(ensg) %>%                 # Extract the ensg (gene names)
    unique()    
  
  # Extract and return data for the filtered genes
  pathway_plot <- limma_results %>%
    filter(toupper(ensg) %in% toupper(filtered_genes),
           coef %in% coefficients)
  #%>%
  #filter(group != "n.s") 
  
  return(list(top_genes = filtered_genes, pathway_plot = pathway_plot))
}

# Initialize an empty list to store results
results <- list()

# Apply the get_top_genes function to each pathway using map()
results <- map(
  names(pathways),  # Apply the function to each pathway
  ~ get_top_genes(
    pathway_name = .x,  # .x is the current pathway name in map
    pathway_genes = pathways[[.x]],  # Access the corresponding genes for that pathway
    limma_results = limmaRes
  )
)

# Set the names of the results list to match the pathway names
names(results) <- names(pathways)
# Extract the `pathway_plot` from each result and combine them into one data frame
combined <- bind_rows(
  map(results, "pathway_plot"),  # Extract pathway_plot from each element
  .id = "Geneset"                # Create a column to indicate the pathway
)
################################################################################
# Define a function to create a plot for each gene set
create_plot <- function(data, geneset) {
  
  filt <-  data %>% filter(Geneset == geneset)%>%
    group_by(celltype,coef)%>%
    arrange(adj.P.Val)%>%
    slice_head(n=5)%>%
    pull(ensg)%>%
    unique()
  data <- data %>% filter(Geneset == geneset,ensg %in% filt)
  len <-length(unique(data$ensg))
  ggplot(data,
         aes(x = gsub("interaction","",coef), y = ensg,
             color = logFC,
             size = pmin(30, -log10(adj.P.Val)))) +  # pmin caps the size
    geom_point() +  # Scatter plot
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +  # Color scale for logFC
    scale_size_continuous(range = c(2, 6), breaks = c(2, 5, 10, 20, 30)) +  # Size for adj.P.Val
    labs(
      title = paste("LogFC and adj.P.Val for Gene Set:", geneset),
      x = "Coefficient",
      y = "Genes",
      color = "logFC",
      size = "-log10(adj.P.Val)"
    ) +
    facet_grid(cols = vars(celltype))+
    theme_bw() +  # Use a clean theme
    optimized_theme()+
    theme(
      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5))
  ggsave(basedir(paste0(geneset, ".pdf")), width = 18, height = len*0.15+2)
}

# Generate plots for each unique gene set
unique_gene_sets <- unique(combined$Geneset)
plots <- map(unique_gene_sets, ~ create_plot(combined, .x))
#walk2(plots, unique_gene_sets, ~ ggsave(basedir(paste0(.y, ".pdf")), plot = .x, width = 10,))
#################################################################################
#genes of interest
Cholesterol <-c("Idi1", "Cyp51","Stard4", "Scd2","Mthfd2","Sqle","Fads2","Dhcr24",
                "Hmgcs1","Ldlr","Plscr1","Acat1",
                "Acat2")


# Plotting the dot plot
f1_3<-ggplot(limmaRes%>% filter(coef %in% coefficients,ensg %in% Cholesterol),
             aes(x = gsub("interaction","",coef), y = ensg,
                                                                   
                  color = logFC,
                                                                    
                 size = pmin(30,-log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(
    range = c(2, 6),  # Set the actual size range from 2 to 30
    breaks = c(2, 5, 10, 20, 30)  # Set specific breaks to create distinct point sizes
  ) +
  labs(title = "Differentially expressed genesets:Cholesterol biosynthesis",
       x = "Cell Type",
       y = "ensg",
       color = "logFC",
       size = "-log10(adj.P.Val)") +
  facet_grid(cols = vars(celltype),scales = "free_y", space = "free") +
  theme_bw() + optimized_theme()+
  theme(
    axis.text.x = element_text(size = 12,angle = 90,vjust = 0.5))
ggsave(basedir("Cholesterol_selected.pdf"),h=13*0.15+2,w=18)

f1_3

#############################
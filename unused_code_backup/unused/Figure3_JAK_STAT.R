source("src/00_init.R")
require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)
require(enrichR)
library(dplyr)
# renv::snapshot(lockfile = "renv_NF.lock")

source("~/code/resources/RFunctions/Basics.R")
source("src/Ag_Optimized_theme.R")

out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/"
base <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Enrichr"

base_fgsea<-"/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/FGSEA"

basedir <- dirout("Figure3")
################################################################################
dirout_jak <- function(out, ext="", init=TRUE){
  out.dir <- paste0("/media/AGFORTELNY/PROJECTS/TfCf_AG/", "/JAKSTAT/", out, "/")
  if(init){
    dir.create(out.dir,showWarnings=FALSE); 
    message("Setting output directory: ", out.dir)
  }
  function(...){
    paste0(out.dir, paste0(...), ext)
  }
}

#function
extract_unique_genes <- function(enrichr_data,
                                 terms,
                                 db) {
  unique_genes <- enrichr_data %>%
    filter(db == db & Term %in% terms) %>%
    pull(Genes) %>%
    unique() %>%
    strsplit(";") %>%
    unlist() %>%
    toupper() %>%
    unique()
  return(unique_genes)
}
################################################################################
res <- read.delim(paste0(out,"/DEG.tsv"))
res$probe<-res$rn
res$rn<-NULL
res <- as.data.frame(res)
gmap <- as.data.frame(read_rds(file = file.path(out, "DEG_GMP.RDS")))
res <- merge(res,gmap[,c("probe","gene")],by="probe")
################################################################################
# Visualize results ---------------------------------------------------------
limmaRes <- res %>% filter(grepl("treatmentex_vivo",res$genotype))
limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))



# Modify the 'genotype' column dynamically for any KO
limmaRes <- limmaRes %>%
  mutate(genotype = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", genotype))

################################################################################


#for all KOs
# Calculate the number of up and downregulated genes for each coefficient and cell type
adj_p_cutoff <- 0.05
logfc_cutoff <- 1
summary_df <- limmaRes %>%
  group_by(cell_type, genotype) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")
summary_df$genotype <- gsub("treatmentex_vivo","WT",
                            summary_df$genotype)
summary_df$genotype <- factor(summary_df$genotype , levels = c("WT",
                                                               setdiff(unique(summary_df$genotype),"WILDTYPE")))

# Filter out rows with count equal to 0 and based on count threshold
filtered_data <-summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= 10) %>%
  group_by(cell_type) %>% 
  filter(genotype %in% unique(genotype)) %>%  # Keep only relevant KOs for each cell type
  ungroup()
filtered_data$genotype <- factor(filtered_data$genotype ,
                                 levels = c("WT",
                                            setdiff(unique(filtered_data$genotype),"WILDTYPE")))
  
fig3.1 <- ggplot(filtered_data, 
              aes(y = fct_relevel(gsub("Interaction_", "", genotype), "WT"),
                  x = ifelse(Regulation == "Downregulated",
                             -log10(Count), log10(Count)), 
                  fill = Regulation)) + geom_col() +
    scale_fill_manual(values = c("Upregulated" = "#D0154E", "Downregulated" = "#4C889C")) +
    labs(y = "Interaction-KO",
         x = "Log10(No: of Genes)",
         title = paste("Upregulated and Downregulated Genes")) +
    theme_bw(12) +
    theme(axis.text = element_text(size = 15)) +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    facet_grid(cols = vars(cell_type), scales = "free",space = "free") +
    theme(strip.text = element_text(size = 15)) +
  coord_flip()+
    optimized_theme()
fig3.1 
# Save plot
ggsave(basedir(paste0("fig3.1", ".pdf")), plot = fig3.1,w=7,
         height = 3.5)


################################################################################
InDir4 <- dirout("Figure1")
genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
colnames(genes_fig_1) <- c("pathways","gene")
# mTORC1_or_Cholesterol <- c("Idi1", "Cyp51","Stard4", "Scd2","Mthfd2","Sqle","Fads2","Dhcr24",
#                            "Hmgcs1","Ldlr","Plscr1","Acat1",
#                            "Acat2","Msmo1")
# cholesterol_df <- data.frame(
#   # Cholesterol genes
#   pathways = rep("mTORC1_or_Cholesterol", length(mTORC1_or_Cholesterol)),
#   gene = mTORC1_or_Cholesterol # Cholesterol pathway
# )
#genes_fig_1 <- rbind(genes_fig_1, cholesterol_df) %>% distinct()


#IFN_genes <- genes_fig_1%>%filter(pathways =="ISG_core")%>% pull(ensg)
# Subset limmaRes for specific genes
limma_subset <- limmaRes %>%
  filter(genotype == "treatmentex_vivo")
 
limma_subset <- merge(limma_subset,genes_fig_1, by ="gene")
#limma_subset$pathways <- gsub("mTORC1_or_Cholesterol","mTORC1/Cholesterol",limma_subset$pathways)
#limma_subset$pathways <- gsub("ISG_core","ISG core",limma_subset$pathways)
limma_subset$pathways <- recode(limma_subset$pathways,
       "ISG_core"= "ISG core",
       "mTORC1_or_Cholesterol" = "mTORC1/Cholesterol")
fig3.2 <- limma_subset %>%
  ggplot(aes(x = cell_type, y = gene, color = logFC,
             size = pmin(10,-log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "blue",#muted("blue"),
                        mid = "white",
                        high = "red") +
  scale_size_continuous(
    range = c(2, 6),  # Set the actual size range from 2 to 30
    breaks = c(1,2,4, 5, 10)  # Set specific breaks to create distinct point sizes
  ) +
  labs(title = "Differentially expressed genesets",
       x = "Celltype",
       y = "Genes",
       color = "logFC",
       size = "-log10(adj.P.Val)") +
 facet_grid(rows = vars(pathways), scales = "free_y", space = "free") +
  theme_bw() + NULL+optimized_theme()
fig3.2
ggsave(basedir("Fig3.2.pdf"),fig3.2, width = 7, height = 7.5)
#################################################################################

###################################################
#fig3.3 

# Calculate the number of up and downregulated genes for each coefficient and cell type
summary_df <- limmaRes %>%
  group_by(cell_type, genotype) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")

# Define your cutoffs and filter the dataframe
adj_p_cutoff <- 0.05
logfc_cutoff <- 1


#filtered based on KOs with atleast 10 differentially regulated genes between tissue types
count_threshold = 10
coefficients <- summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= count_threshold)%>%
  pull(genotype)%>%
  unique()
##########
limma_KO <- limmaRes %>%
  filter(gene %in% genes_fig_1$gene)%>%
  filter(genotype != "treatmentex_vivo", genotype %in% coefficients)%>%
  merge(genes_fig_1, by = "gene" )
limma_KO$pathways <- recode(limma_KO$pathways,
                                  "ISG_core"= "ISG core",
                                  "mTORC1_or_Cholesterol" = "mTORC1/Cholesterol")
fig3.3 <- ggplot(limma_KO, aes(x = gsub("Interaction_","",genotype), y = gene,
                                     
                               color = pmin(2, pmax(-2, logFC)) ,
                               size = pmin(5, -log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  scale_size_continuous(
    range = c(2, 6),  # Set the actual size range from 2 to 30
    breaks = c(1,2,4, 5)  # Set specific breaks to create distinct point sizes
  ) +
  labs(title = "Differentially expressed genesets",
       x = "Cell Type",
       y = "Genes",
       color = "logFC",
       size = "-log10(adj.P.Val)") +
  facet_grid(rows = vars(pathways), cols = vars(cell_type),scales = "free_y", space = "free") +
  theme_bw() +
  optimized_theme()
fig3.3 
ggsave(basedir(paste0("fig3.3.pdf")),plot = fig3.3, 
       width = 10,
       height = 7.5 )
##############################################

# Arrange the plots
combined_plot <- (
  (fig3.1 + plot_spacer() + fig3.3) + 
    plot_layout(widths = c(7.5, 0.5, 10))
) / 
  (fig3.2 +  plot_spacer()+plot_layout(widths = c(7,10.5)))

# Adjusting the heights to match your requirement
combined_plot <- combined_plot + plot_layout(heights = c(3.5, 3.5))

# Save the arranged figure
ggsave(basedir("combined_plot.pdf"), combined_plot, width = 18, height = 17)

###############################################
#Normalized_reads
dataVoom_M <- read.csv("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Normalized_reads/Normalized_readsM.csv",
                   header = T,
                   row.names = 1)
dataVoom_T8 <- read.csv("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Normalized_reads/Normalized_readsT8.csv",
                   header = T,
                   row.names = 1)


# Step 1: Ensure gmap is a data frame if it's a data.table
gmap <- as.data.frame(gmap)

# Step 2: Replace rownames of dataVoom_M
dataVoom_M <- dataVoom_M %>%
  rownames_to_column(var = "ensG") %>%     # Convert rownames to a column for merging
  left_join(gmap[, c("ensG", "gene")], by = "ensG") %>%  # Merge with gmap to get the gene names
  column_to_rownames(var = "gene")         # Set the 'gene' column as rownames

# Step 3: Remove the 'ensG' column after merging
dataVoom_M <- dataVoom_M %>%
  select(-ensG)

# Step 4: Repeat for dataVoom_T8
dataVoom_T8 <- dataVoom_T8 %>%
  rownames_to_column(var = "ensG") %>%
  left_join(gmap[, c("ensG", "gene")], by = "ensG") %>%
  column_to_rownames(var = "gene") %>%
  select(-ensG)

# Check the first few rows of dataVoom_M to confirm
head(dataVoom_M)

meta <- read_rds("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/DEG_Annotation.RDS")
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample_name
colnames(meta)
colnames(meta) <- c("sample_name",
                    "experiment_id",
                    "genotype",
                    "cell_type",
                    "tissue" 
                    )

unique(meta$cell_type)
list_of_genes <- c("Gbp3", "Oas3","Oas1a")
# Initialize an empty list to store data
dat.list <- list()
colnames(meta)

# Loop over each KO genotype and cell type
for (KO in gsub("Interaction_", "", coefficients)) {
  list_of_genes <- c("Gbp3", "Oas3","Oas1a")
 # , "Msmo1", "Hmgcs1", "Stard4", "Mthfd2"
  for (ct in c("M", "T8")) {
    # Get the dataVoom object corresponding to the current cell type
    dataVoom_ct <- get(paste0("dataVoom_", ct))
    
    # Check if any of the genes of interest (goi) exists in the row names of dataVoom_ct
    if (any(rownames(dataVoom_ct) %in% list_of_genes)) {
      for (goi in list_of_genes) {
        # Proceed only if the current gene exists in the row names of dataVoom_ct
        if (goi %in% rownames(dataVoom_ct)) {
          # Filter meta_ct to only include samples that have expression data in dataVoom_ct
          meta_ct <- meta[meta$cell_type == ct, ]
          meta_ct <- meta_ct[rownames(meta_ct) %in% colnames(dataVoom_ct), ]
          meta_ct[!(rownames(meta_ct) %in% colnames(dataVoom_ct)),]
          # Extract the expression values for each sample in the filtered meta_ct
          sample_ids <- rownames(meta_ct)
          expression_values <- dataVoom_ct[goi, sample_ids, drop = FALSE] 
          
          # Prepare gene data by adding the expression values for the current gene
          gene_data <- meta_ct %>%
            mutate(E = as.numeric(expression_values)) %>% # Add expression values as column E
            rownames_to_column("sample1") %>%
            filter(genotype %in% c(KO, "WT")) %>%
            mutate(scaled_E = scale(E)) %>% # Scale the expression data
            mutate(gene = goi,              # Add gene name
                   celltype = ct,           # Add cell type
                   comparison = KO)         # Add KO comparison label
          
          # Store the gene data in the list with a unique key
          dat.list[[paste0(ct, "_", goi, "_", KO)]] <- gene_data
        }
      }
    }
  }
}

# Combine all gene data into a single data frame
goi_exp <- bind_rows(dat.list, .id = "celltype_gene_genotype")
# View the final combined data frame

goi_exp %>% write_rds(basedir("Norm_exp_goi.rds"))
#goi_exp <- read_rds(basedir("Norm_exp_goi.rds"))


# Summarize the data to calculate mean scaled expression for each combination of guide, gene, and condition
# Function to create plots for each gene for each KO

create_gene_plots <- function(data, gene, KO, ct) {
  filtered_data <- data[data$gene == gene & data$celltype == ct,]
  
  filtered_data$genotype <- factor(filtered_data$genotype,levels = c("WT",setdiff(filtered_data$genotype,"WT")))
  # Check if there's any data to plot
  if (nrow(filtered_data) == 0) {
    message(paste("No data for gene:", gene, "and celltype:", ct))
    return(NULL)  # Skip the plot if no data
  }
  # Check for missing values
  if (any(is.na(filtered_data$scaled_E)) || any(is.na(filtered_data$genotype))) {
    message(paste("Missing values found for gene:", gene, "and celltype:", ct))
    return(NULL)  # Skip the plot if missing values
  }
  # Generate the plot if data is present
  p <- ggplot(filtered_data, aes(x = genotype, y = scaled_E, fill = tissue)) + 
      # Boxplot with tissue color fill
    geom_boxplot(
      outlier.colour = NA, 
      position = position_dodge(width = 0.8),  # Adjust width to create space between tissues
      color = "black", 
      size = 0.5
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
    theme_bw()+optimized_theme()+
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
      legend.key = element_blank())                                    # Remove background behind legend items
 
  
  # Get pathway information
  pathway <- genes_fig_1 %>% filter(gene == gene) %>% pull(pathways)
  
  # Save the plot
  ggsave(basedir(paste0("fig3.4_", ct, "_", KO, "_", gene, "_", pathway, ".pdf")), 
         plot = p, device = "pdf")
}
coefficients
# Loop over KOs and cell types
for (comp in gsub("Interaction_","",coefficients)) {
  for (ct in unique(goi_exp$celltype)) {  
    data <- goi_exp %>%
      filter(comparison == comp, celltype == ct)
    
    # Create gene plots for each gene in the list of genes
    gene_plots <- lapply(unique(list_of_genes), function(gene) {
      create_gene_plots(data, gene, comp, ct)
    })
  }
}
#pergene but for all KO


create_gene_plots_all <- function(data, gene,  ct) {
  filtered_data <- data[data$gene == gene & data$celltype == ct,]
  filtered_data$genotype <-factor(filtered_data$genotype,
                                  levels = c("WT",
                                             setdiff(filtered_data$genotype,"WT")))
  # Check if there's any data to plot
  if (nrow(filtered_data) == 0) {
    message(paste("No data for gene:", gene, "and celltype:", ct))
    return(NULL)  # Skip the plot if no data
  }
  # Check for missing values
  if (any(is.na(filtered_data$scaled_E)) || any(is.na(filtered_data$genotype))) {
    message(paste("Missing values found for gene:", gene, "and celltype:", ct))
    return(NULL)  # Skip the plot if missing values
  }
  # Generate the plot if data is present
  boxplot_jitter <- ggplot(filtered_data,
                           aes(x = genotype, y = scaled_E, fill = tissue)) + 
    # Boxplot with tissue color fill
    geom_boxplot(
      outlier.colour = NA,
      position = position_dodge(width = 0.8),  # Adjust width to create space between tissues
      color = "black", 
      size = 0.5
    ) + 
    geom_jitter(
      position = position_jitterdodge(
        jitter.width = 0.2,
        dodge.width = 0.8  # Adjust dodge width to match the boxplot
      ),
      alpha = 0.5)+
     # Jittered points
    facet_grid(cols = vars(tissue), scales = "free") +
    scale_fill_manual(values = c("#C1A0AC", "#87B1D6"),
                      name = "Experimental model") +
    labs(title = paste0(gene, " (", ct, ")")) +
    xlab(paste0(KO, " KO")) +
    theme_bw()+optimized_theme()+
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
      legend.key = element_blank())                                    # Remove background behind legend items
  
  
  # Get pathway information
  pathway <- genes_fig_1 %>% filter(gene == gene) %>% pull(pathways)
  summary_data <- filtered_data %>%
    group_by(genotype, tissue) %>%
    summarize(median_scaled_E = median(scaled_E, na.rm = TRUE)) %>%
    ungroup()
  
  diff_to_wt <- summary_data %>%
    group_by(tissue) %>%
    mutate(diff_to_WT = median_scaled_E - median_scaled_E[genotype == "WT"]) %>%
    ungroup() %>%
    filter(genotype != "WT")
  
  # Set up the difference-to-WT bar plot
  diff_to_wt_plot <- ggplot(diff_to_wt, aes(x = genotype, y = diff_to_WT, fill = tissue)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    labs(title = paste( gene, "(", ct, ")"),
         x = "Genotype (KO)", y = "KO-WT") +
    scale_fill_manual(values = c("ex_vivo" = "#C1A0AC", "in_vivo" = "#87B1D6"),
                      name = "Experimental Model",guide = "none") +
    theme_minimal() +
    optimized_theme()+guides()
   
  
  # Combine the two plots side by side and collect guides
  combined_plot <- (boxplot_jitter / diff_to_wt_plot) + 
    plot_layout(guides = "collect",
                widths = c(2,0.5),
                heights = c(1.5,1))
  combined_plot# Collect legends into a single guide
      # Add an overall title if needed
  ggsave(basedir(paste0("fig3.4_KO-WT", ct, "_", gene, "_", pathway, ".pdf")), 
         plot = combined_plot, device = "pdf",)
  # Display the combined plot
}  
# Loop over KOs and cell types
for (ct in unique(goi_exp$celltype)) {  
    data <- goi_exp %>%
      filter(celltype == ct)
    
    # Create gene plots for each gene in the list of genes
    gene_plots <- lapply(unique(list_of_genes), function(gene) {
      create_gene_plots_all(data, gene, ct)
    })
}



source("src/00_init.R")
require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)
require(enrichR)
library(dplyr)
require(latex2exp)
library(patchwork)
# renv::snapshot(lockfile = "renv_NF.lock")
source("~/code/resources/RFunctions/Basics.R")
source("src/Ag_Optimized_theme_fig.R")

out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/"
base <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Enrichr"

base_fgsea<-"/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/FGSEA"

basedir <- dirout("Figure3_paper")

#*******************#
res <- read.delim(paste0(out,"/DEG.tsv"))
res$probe <- res$rn
res$rn<-NULL
res <- as.data.frame(res)
gmap <- as.data.frame(read_rds(file = file.path(out, "DEG_GMP.RDS")))
res <- merge(res,gmap[,c("probe","gene")],by="probe")
limmaRes <- res %>% filter(grepl("treatmentex_vivo",res$genotype))
limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))

# Modify the 'genotype' column for any KO
limmaRes <- limmaRes %>%
  mutate(genotype = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", genotype))%>%
  mutate(cell_type = gsub("M", "Macrophage", cell_type)) %>%
  mutate(cell_type = gsub("T8", "T-cells", cell_type)) 
limmaRes$genotype <- gsub("treatmentex_vivo","WT",
                          limmaRes$genotype)


unique(limmaRes$genotype)

#
#fig3.1-------------------------------------------------------------------------
#

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

# Filter out rows with count equal to 0 and based on count threshold
filtered_data <-summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= 10) %>%
  group_by(cell_type) %>% 
  filter(genotype %in% unique(genotype)) %>%  # Keep only relevant KOs for each cell type
  ungroup()

filtered_data$genotype <- factor(filtered_data$genotype ,
                                 levels = c("WT",
                                            setdiff(unique(filtered_data$genotype),"WT")))

#fig------------- 
fig3.1 <- ggplot(filtered_data, 
              aes(y = fct_relevel(gsub("Interaction_", "", genotype), "WT"),
                  x = ifelse(Regulation == "Downregulated",
                             -log10(Count), log10(Count)), 
                  fill = Regulation)) + geom_col() +
  scale_fill_manual(values = c("Upregulated" = "#D0154E", "Downregulated" = "#4C889C")) +
  labs(y = "Interaction-KO",
         x = TeX("log_{10}(No. of genes)"),
         title = paste("Upregulated and Downregulated Genes")) +
  facet_grid(cols = vars(cell_type), scales = "free",space = "free")  +
  coord_flip()+
  optimized_theme_fig()+
  theme(legend.position = "bottom") 
fig3.1 
# Save plot
ggsave(basedir(paste0("fig3.1", ".pdf")), plot = fig3.1,w=7,
         height = 3.5)

#
#fig3.2--------
#

InDir4 <- dirout("Figure1")
genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
colnames(genes_fig_1) <- c("pathways","gene")

limma_subset <- limmaRes %>%
  filter(genotype == "WT")
 
limma_subset <- merge(limma_subset,genes_fig_1, by ="gene")

limma_subset$pathways <- recode(limma_subset$pathways,
       "ISG_core"= "ISG core",
       "mTORC1_or_Cholesterol" = "mTORC1/Cholesterol")
#fig----
fig3.2 <- limma_subset %>%
  ggplot(aes(x = cell_type, y = gene, 
             color = pmin(2, pmax(-2, logFC)) ,
             size = pmin(5, -log10(adj.P.Val))))+
  geom_point() + 
  
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name=TeX("log_{2}(FC)")) +
  scale_size_continuous(
    range = c(1, 3),  # Set the actual size range from 2 to 30
    breaks = c(1,2,4, 5, 10),# Set specific breaks 
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(title = "DEGs",
       x = "Celltype",
       y = "Genes") +
 facet_grid(rows = vars(pathways), scales = "free_y", space = "free") +
  theme_bw() + optimized_theme_fig()
fig3.2
ggsave(basedir("Fig3.2.pdf"),fig3.2, width = 7, height = 7.5)

#
#fig3.3--------
#

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

head(limmaRes)
limma_KO <- limmaRes %>%
  filter(gene %in% genes_fig_1$gene)%>%
  filter(genotype != "WT", genotype %in% coefficients)%>%
  inner_join(summary_df, by = c("genotype","cell_type"))%>%
  filter(Count >10)%>%
  merge(genes_fig_1, by = "gene" )
limma_KO$pathways <- recode(limma_KO$pathways,
                                  "ISG_core"= "ISG core",
                                  "mTORC1_or_Cholesterol" = "mTORC1/Cholesterol")
#fig----
fig3.3 <- ggplot(limma_KO, aes(x = gsub("Interaction_","",genotype), y = gene,
                               color = pmin(2, pmax(-2, logFC)) ,
                               size = pmin(5, -log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name = TeX("$log_{2}(FC)$"))+
  scale_size_continuous(
    range = c(1, 3),  # Set the actual size range from 2 to 30
    breaks = c(1,2,4, 5, 10),
    name = TeX("$-\\log_{10}(p_{adj})$")# Set specific breaks to create distinct point sizes
  ) +
  labs(title = "DEGs",
       x = "Cell Type",
       y = "Genes"
       ) +
  facet_grid(rows = vars(pathways), cols = vars(cell_type),scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()
fig3.3 
ggsave(basedir(paste0("fig3.3.pdf")),plot = fig3.3, 
       width = 10,
       height = 7.5 ,units = "cm")

#
#fig3.4-----
#

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
  column_to_rownames(var = "gene") %>% select(-ensG)

# Step 4: Repeat for dataVoom_T8
dataVoom_T8 <- dataVoom_T8 %>%
  rownames_to_column(var = "ensG") %>%
  left_join(gmap[, c("ensG", "gene")], by = "ensG") %>%
  column_to_rownames(var = "gene") %>%
  select(-ensG)

meta <- read_rds("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/DEG_Annotation.RDS")
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample_name
colnames(meta) <- c("sample_name",
                    "experiment_id",
                    "genotype",
                    "cell_type",
                    "tissue" 
                    )
#make it custom
list_of_genes <- c("Gbp3", "Oas3","Oas1a")
# Initialize an empty list to store data
dat.list <- list()
colnames(meta)

# Loop over each KO genotype and cell type
for (KO in gsub("Interaction_", "", coefficients)) {
  list_of_genes <- c("Oas1a")
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

#Macrophage fig----
summary <- summary_df %>%
  mutate( celltype = cell_type)%>%
  select(-cell_type)

summary$genotype <- gsub("Interaction_","",summary$genotype)
summary <- summary %>%
  mutate(comparison = genotype)%>%
  select(-genotype)
summary$celltype <- gsub("Macrophage","M", summary$celltype)
summary$celltype <- gsub("T-cells","T8", summary$celltype)
ct <- "M"
data <- goi_exp
gene <-"Oas1a"
filtered_data <- data[data$gene == gene & data$celltype == ct,]
filtered_data$genotype <-factor(filtered_data$genotype,
                                levels = c("WT",
                                           setdiff(filtered_data$genotype,"WT")))
filtered_data <- filtered_data%>%
  inner_join(summary, by = c("comparison", "celltype"))%>%
  filter(Count > 10)
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
boxplot_jitter_M <- ggplot(filtered_data,
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
  theme_bw()+
  optimized_theme_fig()#+

# Get pathway information
pathway <- genes_fig_1 %>% filter(gene == gene) %>% pull(pathways)
summary_data <- filtered_data %>%
  #inner_join(summary, by = c("genotype","celltype"))%>%
 # filter(Count > 10)%>%
  group_by(genotype, tissue) %>%
  summarize(median_scaled_E = median(scaled_E, na.rm = TRUE)) %>%
  ungroup()

diff_to_wt <- summary_data %>%
  group_by(tissue) %>%
  mutate(diff_to_WT = median_scaled_E - median_scaled_E[genotype == "WT"]) %>%
  ungroup() %>%
  filter(genotype != "WT")

# Set up the difference-to-WT bar plot
diff_to_wt_plot_M <- ggplot(diff_to_wt, aes(x = genotype, y = diff_to_WT, fill = tissue)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = paste( gene, "(", ct, ")"),
       x = "Genotype (KO)", y = "KO-WT") +
  scale_fill_manual(values = c("ex_vivo" = "#C1A0AC", "in_vivo" = "#87B1D6"),
                    name = "Experimental Model") +
  theme_minimal() +
  optimized_theme_fig()
fig3.4.1 <- diff_to_wt_plot_M
fig3.4.1.sup <- boxplot_jitter_M
# Combine the two plots side by side and collect guides

  # Display the combined plot
# T8
ct <- "T8"
data <- goi_exp
gene <-"Oas1a"
filtered_data <- data[data$gene == gene & data$celltype == ct,]
filtered_data$genotype <-factor(filtered_data$genotype,
                                levels = c("WT",
                                           setdiff(filtered_data$genotype,"WT")))
filtered_data <- filtered_data%>%
  inner_join(summary, by = c("comparison", "celltype"))%>%
  filter(Count > 10)
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
#fig----
# Generate the plot if data is present
boxplot_jitter_T8 <- ggplot(filtered_data,
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
  theme_bw()+optimized_theme_fig()#+

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
diff_to_wt_plot_T8 <- ggplot(diff_to_wt, aes(x = genotype, y = diff_to_WT, fill = tissue)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = paste( gene, "(", ct, ")"),
       x = "Genotype (KO)", y = "KO-WT") +
  scale_fill_manual(values = c("ex_vivo" = "#C1A0AC", "in_vivo" = "#87B1D6"),
                    name = "Experimental Model") +
  theme_minimal() +
  optimized_theme_fig()
fig3.4.2 <- diff_to_wt_plot_T8 
fig3.4.2.sup <- boxplot_jitter_T8
# Combine the two plots side by side and collect guides
combined_plot_T8 <- (boxplot_jitter_T8 / diff_to_wt_plot_T8) + 
  plot_layout(guides = "collect",
              widths = c(2,0.5),
              heights = c(2,1))
combined_plot_T8
ggsave(basedir(paste0("fig3.4_KO-WT", ct, "_", gene, "_", pathway, ".pdf")), 
       plot = combined_plot_T8, device = "pdf")  

##combined_plots-------------


# Assuming fig3.1, fig3.2, fig3.3, combined_plot_T8, and combined_plot_M are defined

fig3_sup_2 <- (fig3.4.1.sup + fig3.4.2.sup) + 
  plot_layout(widths = c(1, 0.7),guides = "collect")

ggsave(basedir("Fig3_supplementary.pdf"), 
       plot = fig3_sup_2, 
       width = 18, 
       height = 17, 
       units = "cm")

# Create the left column and combine with fig3.2

second_column <-  (diff_to_wt_plot_M | diff_to_wt_plot_T8) + 
  plot_layout(widths = c(1, 0.7), guides = "collect") &
  theme(legend.position = "bottom") 
# Combine the rows according to your specification
# fig3 <- (top_column + second_column) +  # '/' arranges them in rows
#   plot_layout(heights = c(5, 1), nrow = 2) 

# fig3 <- (top_column/second_column )+ 
#   plot_layout(heights = c(5, 1),widths = c(1,1))

# Set the overall dimensions

# Set the output size (in cm)
ggsave(basedir("fig3.4.pdf"), 
       plot = second_column, 
       width = 18, 
       height = 17, 
       units = "cm")

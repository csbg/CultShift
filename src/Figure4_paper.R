source("src/00_init.R")
require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)
require(enrichR)
library(dplyr)
require(latex2exp)
# renv::snapshot(lockfile = "renv_NF.lock")
source("~/code/resources/RFunctions/Basics.R")
source("src/Ag_Optimized_theme_fig.R")

out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT_Ar/"
basedir <- dirout("Figure4_paper")

#*******************#
res <- rbind(read.delim(paste0(out,"/DEG_ResultsT8.tsv")),
             read.delim(paste0(out,"/DEG_ResultsM.tsv"))
        )
res$probe <- res$rn
res$rn<-NULL
res <- as.data.frame(res)
gmap <- as.data.frame(read_rds(file = file.path(out, "DEG_GMP.RDS")))
res <- merge(res,gmap[,c("probe","gene")],by="probe")
limmaRes <- res %>% filter(grepl("treatmentex_vivo",res$coef))


limmaRes$group <- ifelse(limmaRes$logFC >= 1 & 
                           limmaRes$adj.P.Val <= 0.05, "up", 
                         ifelse(limmaRes$logFC <= -1 & 
                                  limmaRes$adj.P.Val <= 0.05, "down", "n.s"))

# Modify the 'coef' column for any KO
limmaRes <- limmaRes %>%
  mutate(coef = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(cell_type = gsub("M", "Macrophage", cell_type)) %>%
  mutate(cell_type = gsub("T8", "T-cells", cell_type)) 
limmaRes$coef <- gsub("treatmentex_vivo","WT",
                          limmaRes$coef)


unique(limmaRes$coef)

#
#Fig4.2-------------------------------------------------------------------------
#

# Calculate the number of up and downregulated genes for each coefficient and cell type
adj_p_cutoff <- 0.05
logfc_cutoff <- 1
summary_df <- limmaRes %>%
  group_by(cell_type, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")

ko_flag <- limmaRes %>%
  group_by(cell_type, coef) %>%
  summarise(
    Count = sum(adj.P.Val < adj_p_cutoff & abs(logFC) > logfc_cutoff)
  ) %>% group_by(coef, cell_type)%>%
  mutate(valid_ko = ifelse(Count >10, TRUE, FALSE))
# Filter out rows with count equal to 0 and based on count threshold
filtered_data <-summary_df %>% 
  filter(Count != 0) %>% 
  filter(Count >= 10) %>%
  group_by(cell_type) %>% 
  filter(coef %in% unique(coef)) %>%  # Keep only relevant KOs for each cell type
  ungroup()

  
  
filtered_data$coef <- factor(filtered_data$coef ,
                                 levels = c("WT",
                                            setdiff(unique(filtered_data$coef),"WT")))
WT <- filtered_data %>%
  filter(coef == "WT")
#fig------------- 
Fig4.2_1 <- ggplot(WT, 
              aes(y = cell_type,
                  x = ifelse(Regulation == "Downregulated",
                             -log10(Count), log10(Count)), 
                  fill = Regulation)) + geom_col() +
  scale_fill_manual(values = c("Upregulated" = "#D0154E", "Downregulated" = "#4C889C")) +
  labs(y = "Interaction-KO",
         x = TeX("log_{10}(No. of genes)"),
         title = paste("Upregulated and Downregulated Genes")) +
 # facet_grid(cols = vars(cell_type), scales = "free",space = "free")  +
  coord_flip()+
  optimized_theme_fig()+
    theme(legend.position = "right") 
Fig4.2_1
# Save plot
ggsave(basedir(paste0("Fig4.2_1", ".pdf")), plot = Fig4.2_1,w=5,
         height = 5, units = "cm")

Interaction <- filtered_data %>%
  filter(coef != "WT")
Interaction <-Interaction %>%
  group_by(cell_type, coef) %>%
  summarize(Total_Regulated = sum(Count), .groups = "drop")
Fig4.2_2 <- ggplot(Interaction,aes(
  x = gsub("Interaction_","",coef),
  y = log10(Total_Regulated)
)) +
  geom_col() +
  # Custom colors for upregulated and downregulated bars
  scale_fill_manual(values = c("Upregulated" = "#D0154E", "Downregulated" = "#4C889C")) +
  # Facet by celltype with free space for flexibility in cell widths
  facet_grid(cols = vars(cell_type), space = "free", scales = "free") +
  labs(
    #title = "No. of Genes with interaction effect of culture condition in KOs",
    x = NULL,
    y = TeX("$\\log_{10}\\; (No.of\\; genes\\; with\\; interaction\\; effects)$")
  ) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    strip.text.x = element_text(angle = 90,
                                hjust = 0))


Fig4.2_2
ggsave(basedir(paste0("Fig4.2_2", ".pdf")), plot = Fig4.2_2,w=7,
       height = 5, units = "cm")
#Fig4.1--------
#

InDir4 <- dirout("Figure1")
genes_fig_1 <- read_rds(InDir4("genes_fig1.rds"))
colnames(genes_fig_1) <- c("gene","pathways")

limma_subset <- limmaRes %>%
  filter(coef == "WT")
#only want ISG
limma_subset <- merge(limma_subset,genes_fig_1, by ="gene")%>%
  filter(pathways %in% c("ISG core","mTORC1/Cholesterol","mTORC1_or_Cholesterol"))

limma_subset$pathways <- dplyr::recode(limma_subset$pathways,
       "ISG core"= "ISG core",
       "mTORC1_or_Cholesterol" = "mTORC1/Cholesterol")
#fig----
Fig4.1 <- limma_subset %>%
  ggplot(aes(x = cell_type, y = gene, 
             color = pmin(2, pmax(-2, logFC)) ,
             size = pmin(5, -log10(adj.P.Val))))+
  geom_point() + 
  
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name=TeX("log_{2}(FC)")) +
  scale_size_continuous(
    range = c(1, 2.5),  # Set the actual size range from 2 to 30
    #breaks = c(1,2,4, 5, 10),# Set specific breaks 
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(title = "DEGs",
       x = "cell_type",
       y = "Genes") +
 facet_grid(rows = vars(pathways), scales = "free_y", space = "free") +
  theme_bw() + optimized_theme_fig()
#+theme(axis.text.x = element_text(angle = 0))
Fig4.1
ggsave(basedir("Fig4.1.pdf"),Fig4.1, width = 5, height = 13, units = "cm")

#
#Fig4.3--------
#

# Calculate the number of up and downregulated genes for each coefficient and cell type
summary_df <- limmaRes %>%
  group_by(cell_type, coef) %>%
  summarise(
    Upregulated = sum(adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff),
    Downregulated = sum(adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff)
  ) %>%
  pivot_longer(cols = c(Upregulated, Downregulated),
               names_to = "Regulation", values_to = "Count")%>%
  mutate(comparison = gsub("Interaction_","",coef))%>%
  select(-coef)%>%
  mutate(cell_type =cell_type)

# Define your cutoffs and filter the dataframe
adj_p_cutoff <- 0.05
logfc_cutoff <- 1


#filtered based on KOs with atleast 10 differentially regulated genes between tissue types
count_threshold = 10
coefficients <- ko_flag%>%
  filter(Count >= count_threshold)%>%
  pull(coef)%>%
  unique()

limma_KO <- limmaRes %>%
  filter(gene %in% genes_fig_1$gene)%>%
  filter(coef != "WT", coef %in% coefficients)%>%
  merge(genes_fig_1, by = "gene" )%>%
  filter(pathways == "ISG core")
  
limma_KO$pathways <- recode(limma_KO$pathways,
                                  "ISG_core"= "ISG core",
                                  "mTORC1_or_Cholesterol" = "mTORC1/Cholesterol")
limma_KO <- limma_KO%>%
  inner_join(ko_flag, by = c("coef","cell_type"))%>%
  filter(Count > 10)
#fig----
Fig4.3 <- ggplot(limma_KO, aes(x = gsub("Interaction_","",coef), y = gene,
                               color = pmin(2, pmax(-2, logFC)) ,
                               size = pmin(5, -log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name = TeX("$log_{2}(FC)$"))+
  scale_size_continuous(
    range = c(1, 2.5),  # Set the actual size range from 2 to 30
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
Fig4.3 
ggsave(basedir(paste0("Fig4.3.pdf")),plot = Fig4.3, 
       width = 8,
       height = 7 ,units = "cm")

#
#Fig4.4-----
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
list_of_genes <- c("Gbp3", "Oas3","Oas1a","Irf7")
# Initialize an empty list to store data
dat.list <- list()
colnames(meta)

# Loop over each KO genotype and cell type
for (KO in gsub("Interaction_", "", coefficients)) {
 # list_of_genes <- c("Oas3")
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
           # mutate(E = E))%>% # Scale the expression data
            mutate(gene = goi,              # Add gene name
                   cell_type = ifelse(ct == "M",gsub("M","Macrophage",ct), "T-cells"),          # Add cell type
                   comparison = KO)         # Add KO comparison label
          
          # Store the gene data in the list with a unique key
          dat.list[[paste0(ct, "_", goi, "_", KO)]] <- gene_data
        }
      }
    }
  }
}

# Combine all gene data into a single data frame
goi_exp <- bind_rows(dat.list, .id = "cell_type_gene_genotype")

#goi_exp <- read_rds(basedir("Norm_exp_goi.rds"))
goi <- "Oas1a"
# fig4D----and supplementary 4B
koi <-  ko_flag %>% # Remove grouping from previous steps
  mutate(comparison = gsub("Interaction_", "", coef)) %>%
  mutate(cell_type = cell_type) %>%
  select(cell_type, coef, valid_ko, Count,comparison)

filtered_data <- goi_exp %>%
  
  select(-cell_type)%>%
  filter(gene == goi)%>%
  inner_join(koi, by =c("comparison"))%>%
  filter(valid_ko == TRUE)
#%>%group_by(cell_type) %>%
 
   
filtered_data$coef <-factor(filtered_data$coef,
                                levels = c("WT",
                                           setdiff(filtered_data$coef,"WT")))
# Check if there's any data to plot
if (nrow(filtered_data) == 0) {
  message(paste("No data for gene:", gene, "and cell_type:", ct))
  return(NULL)  # Skip the plot if no data
}
# Check for missing values
if (any(is.na(filtered_data$E)) || any(is.na(filtered_data$genotype))) {
  message(paste("Missing values found for gene:", gene, "and cell_type:", ct))
  return(NULL)  # Skip the plot if missing values
}
#fig supplementary 4B----

# Generate the plot if data is present
boxplot_jitter <- ggplot(filtered_data,
                            aes(x = genotype, y = E, fill = tissue)) + 
  # Boxplot with tissue color fill
  geom_boxplot(
    outlier.colour = NA,
    position = position_dodge(width = 0.8),  # Adjust width to create space between tissues
    color = "black", 
    size = 0.3
  ) + 
  geom_jitter(data = filtered_data%>%filter(genotype != "WT"),
              
              position = position_jitterdodge(
                jitter.width = 0.2,
                dodge.width = 0.8  # Adjust dodge width to match the boxplot
              ),
              alpha = 0.5)+
  
  facet_grid(rows = vars(tissue),cols = vars(cell_type), scales = "free") +
  scale_fill_manual(values = c("#6F987B", "#764BAB"),
                    name = "Experimental model") +
  labs(title = paste0(goi)) +
  xlab(paste0(KO, " KO")) +
  theme_bw()+
    optimized_theme_fig()#+
#supplementaryfig4B----
ggsave(basedir("SupplementaryFig6B.pdf"), 
       plot = boxplot_jitter, 
       width = 12, 
       height = 6, 
       units = "cm")
# Get pathway information
pathway <- genes_fig_1 %>% filter(gene == gene) %>% pull(pathways)
summary_data <- filtered_data %>%
  group_by(genotype, tissue, cell_type) %>%
  summarize(median_E = median(E, na.rm = TRUE)) %>%
  ungroup()

diff_to_wt <- summary_data %>%
  group_by(tissue, cell_type) %>%
  mutate(diff_to_WT = median_E - median_E[genotype == "WT"]) %>%
  ungroup() %>%
  filter(genotype != "WT")

# fig 4D
diff_to_wt_plot <- ggplot(diff_to_wt, 
                          aes(
                            x = genotype,
                              y = diff_to_WT,
                              fill = tissue
                              )
                          ) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    color = "black") +
  geom_hline(
    yintercept = 0,
             linetype = "dashed",
    color = "gray") +
  labs(
    title = paste( goi),
       x = "Genotype (KO)", y = "KO-WT") +
  scale_fill_manual(
    values = c("ex_vivo" = "#6F987B",
               "in_vivo" = "#764BAB"),
                    name = "Experimental Model"
    ) +
  facet_grid(cols = vars(cell_type), 
             scale= "free",
             space = "free"
             )+
  theme_minimal() +
  optimized_theme_fig()


# Set the output size (in cm)
ggsave(basedir("Fig4D.pdf"), 
       plot = diff_to_wt_plot, 
       width = 12, 
       height = 5, 
       units = "cm")

################
diff_to_wt_plot2<- ggplot(diff_to_wt, aes(x = genotype, y = diff_to_WT, group = genotype, color = tissue)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_line(position = position_dodge(width = 0.3), aes(group = genotype), color = "gray70") +
  facet_wrap(~ cell_type, scales = "free_x") +
  
  scale_color_manual(values = c("ex_vivo" = "#6F987B", "in_vivo" = "#764BAB")) +
  labs(x = "Genotype (KO)", y = "KO - WT", color = "Model") +
  optimized_theme_fig() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Set the output size (in cm)
ggsave(basedir("Fig4D_option2.pdf"), 
       plot = diff_to_wt_plot2, 
       width = 12, 
       height = 5, 
       units = "cm")
##############
#option 3
# Filter to just the gene of interest
oas1a_data <- filtered_data %>%
  filter(gene == "Oas1a") %>%
  select(sample_name, genotype, cell_type, tissue, E) %>%
  mutate(rep_id = row_number())  # gives each row a unique ID

paired_data <- oas1a_data %>%
  pivot_wider(names_from = tissue, values_from = E, id_cols = c(rep_id, sample_name, genotype, cell_type)) %>%
  filter(!is.na(in_vivo) & !is.na(ex_vivo))

ggplot(paired_data, aes(x = ex_vivo, y = in_vivo)) +
  geom_point(size = 3, aes(color = genotype)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ genotype + cell_type) +
  labs(x = "Ex vivo Oas1a expression",
       y = "In vivo Oas1a expression",
       title = "Oas1a Expression: In Vivo vs Ex Vivo") +
  theme_minimal(base_size = 13)

#correlation-------------------------
res$group <- ifelse(res$logFC >= 1 & 
                      res$adj.P.Val <= 0.05, "up", 
                    ifelse(res$logFC <= -1 & 
                             res$adj.P.Val <= 0.05, "down", "n.s"))



# summary_df <-in.vivo_degs %>%
#   group_by(cell_type, coef) %>%
#   summarise(
#     Upregulated = sum(adj.P.Val < 0.05 & logFC > 1),
#     Downregulated = sum(adj.P.Val < 0.05 & logFC < -1)
#   ) %>%
#   pivot_longer(cols = c(Upregulated, Downregulated),
#                names_to = "Regulation", values_to = "Count")

# Prepare ex.vivo data
ex_vivo_data <- res %>%
  filter(coef != "treatmentex_vivo") %>%
  filter(str_detect(coef, "^ex.vivo")) %>%
  mutate(coef_clean = str_replace(coef, "ex.vivo\\.", "")) %>%
  mutate(genotype = gsub("ex.vivo_genotype", "", coef)) %>%
  dplyr::select("gene", "logFC", "cell_type", "genotype", "adj.P.Val")

# Prepare in.vivo data
in_vivo_data <- res %>%
  filter(!(coef %in% grep("treatmentex_vivo", res$coef, value = T))) %>%
  filter(str_detect(coef, "^genotype")) %>%
  mutate(coef_clean = str_replace(coef, "genotype\\.", "")) %>%
  mutate(genotype = gsub("genotype", "", coef)) %>%
  dplyr::select("gene", "logFC", "cell_type", "genotype", "adj.P.Val")
unique(ex_vivo_data$cell_type)
# Combine datasets
merged_data <- ex_vivo_data %>%
  inner_join(in_vivo_data, by = c("gene", "cell_type", "genotype"), 
             suffix = c("_ex.vivo", "_in.vivo"))
saveRDS(merged_data,basedir("in.vivo_ex.vivo_logFC_JAK_stat.rds"))

merged_data %>%
  filter(genotype == "IRF9KO")%>%
  filter(cell_type=="T8")%>%
  ggplot(aes(x = logFC_ex.vivo, y = logFC_in.vivo)) +
  geom_hex(bins = 40, aes(fill = after_stat(count))) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") +  # Blue gradient
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c")+  # Linear regression line
  #scale_fill_viridis_c() +
  #scale_fill_viridis_c() +  # Optional: nice continuous color scale
  
  labs(
    title = "Hexbin Plot: logFC ex vivo vs in vivo",
    x = "logFC (Ex Vivo)",
    y = "logFC (In Vivo)",
    fill = "Gene Count"
  )+ optimized_theme_fig()


# Calculate the correlation for each genotype within each cell type, excluding NA
correlation_results <- merged_data %>%
  group_by(cell_type, genotype) %>%
  summarise(correlation = cor(logFC_ex.vivo, logFC_in.vivo, use = "complete.obs")) %>%
  ungroup()

#deg
deg_ex_vivo <- ex_vivo_data %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype, cell_type) %>%
  summarise(degs_ex_vivo = n())

deg_in_vivo <- in_vivo_data %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05) %>%
  group_by(genotype, cell_type) %>%
  summarise(degs_in_vivo = n())

# Merge the DEG counts
deg_counts <- full_join(deg_ex_vivo, deg_in_vivo, by = c("genotype", "cell_type")) %>% na.omit()
deg_plot_data <- deg_counts %>%
  pivot_longer(cols = starts_with("degs"), names_to = "condition", values_to = "num_degs") %>%
  mutate(condition = ifelse(condition == "degs_ex_vivo", "Ex Vivo", "In Vivo")) %>%
  left_join(correlation_results, by = c("genotype", "cell_type")) %>%
  group_by(cell_type, genotype) %>%
  slice_max(order_by = num_degs, n = 1, with_ties = FALSE) %>%  # Retain only the row with max num_degs
  ungroup()


write_rds(deg_plot_data,basedir("DEGs_per_tissue.rds"))


# Create the plot with adjusted correlation color scale
Fig4.5 <- deg_plot_data %>%
  ggplot() +
  geom_point(aes(
    x = genotype,
    y = cell_type,
    size = pmin(3, log10(num_degs)),
    fill = correlation  # Set transparency based on KO validity
  ),
  shape = 21,           # Use shape 21 to enable fill and color
  color = "black",       # Black outline
  stroke = 0.5          # Adjust the width of the outline
  ) +
  scale_fill_gradient2(
    low = "#4C889C",    # Color for low correlation
    mid = "white",      # Midpoint color
    high = "#D0154E",   # Color for high correlation
    midpoint = 0,       # Set the midpoint of the color scale to 0 for more balanced color mapping
    limits = c(-1, 1)   # Set limits to cover the full range of correlation values (assuming correlation ranges from -1 to 1)
  ) +
  scale_size_continuous(
    range = c(0,3.5),
    limits = c(0,3),
    breaks = c(1,2,3),
    name = TeX("$\\log_{10}\\; (\\No.\\; of \\;DEGs)$")
  ) +
  labs(x = "KOs",
       y = "cell_type") +
  optimized_theme_fig() +
  theme(
    legend.justification = "right"
  )
example <- merged_data %>%
  filter(genotype == "IRF9KO")%>%
  filter(cell_type=="T8")%>%
  ggplot(aes(x = logFC_ex.vivo, y = logFC_in.vivo)) +
  geom_hex(bins = 40, aes(fill = after_stat(count))) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") +  # Blue gradient
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c")+  # Linear regression line
  #scale_fill_viridis_c() +
  #scale_fill_viridis_c() +  # Optional: nice continuous color scale
  
  labs(
    title = "Hexbin Plot: logFC ex vivo vs in vivo",
    x = "logFC (Ex Vivo)",
    y = "logFC (In Vivo)",
    fill = "Gene Count"
  )+ optimized_theme_fig()

Fig4.5
Fig4.5 <- example+ Fig4.5 + plot_layout(widths = c(1,3))
ggsave(basedir("correlation_homeo_ex.pdf"), w=12, h= 5, units = "cm")
##########################
# correlation


oas1a_data <- merged_data %>%
  filter(gene == "Oas1a")

ggplot(oas1a_data %>% filter(cell_type== "M"), aes(x = logFC_ex.vivo, y = logFC_in.vivo)) +
  geom_point(size = 1, color = "#4C889C") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  facet_grid(cols = vars(genotype), scales = "free") +
  labs(
    title = "Oas1a: logFC in Ex Vivo vs In Vivo",
    x = "logFC (Ex vivo)",
    y = "logFC (In vivo)"
  ) +
  optimized_theme_fig()+
  theme(
    strip.text.x = element_text(angle = 45)
  )

ggsave(basedir("Fig.4E.logFC_scatter_M.pdf"), w = 6.42 , h =2, units = "cm" )


ggplot(oas1a_data, aes(x = logFC_ex.vivo, y = logFC_in.vivo)) +
  geom_point(size = 1, color = "#4C889C") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  facet_wrap(genotype ~ cell_type) +
  labs(
    title = "Oas1a: logFC in Ex Vivo vs In Vivo",
    x = "logFC (Ex vivo)",
    y = "logFC (In vivo)"
  ) +
  optimized_theme_fig()
 # theme(
 #   strip.text.x = element_text(angle = 45)
#)

ggsave(basedir("Fig.4E.logFC_scatter.pdf"), w = 9 , h = 8, units = "cm" )

limmaRes_int <- limmaRes %>%
  #filter(group != "n.s") %>%
  filter(coef %in% grep("Interaction", coef, value = TRUE)) %>%
  mutate(
    coef = gsub("interaction_", "", coef, ignore.case = TRUE),
    cell_type = if_else(cell_type == "T-cells", "T8", "M")
  )


limmaRes_NTC <- res %>% filter(coef == "treatmentex_vivo")%>%
  #filter(group != "n.s") %>%
  dplyr::select(gene,logFC,cell_type,group)


merged_data <- limmaRes_int %>%
  inner_join(limmaRes_NTC, by = c("gene","cell_type"))%>%  # Adjust "gene" to your actual column name for joining
  mutate(logFC_KO = logFC.x,
         logFC_NTC = logFC.y)
         #adj.P.Val_KO = adj.P.Val.x,
         #adj.P.Val_NTC = adj.P.Val.y)
summary_df$coef <- summary_df$comparison
summary_df <- summary_df %>%
  mutate(cell_type = if_else(cell_type == "T-cells", "T8", "M")
  )
# Step 1: Calculate correlation with p-value for each KO and celltype
correlation_results <- merged_data %>%
  #inner_join(ko_flags, by = c("coef","celltype")) %>%
  #filter(valid_ko) %>%
  group_by(coef, cell_type) %>%
  
 # filter(coef %in% koi) %>%
  inner_join(summary_df, by = c("coef","cell_type")) %>%
  filter(Count > 10) %>%
  distinct(coef, cell_type, gene, .keep_all = TRUE)%>%

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
Fig4.5 <- ggplot(correlation_results, aes(x = coef, y = cor_abs)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = significance), 
            vjust = -0.4, 
            color = "black",
            size = 2.5) + 
  # Add enrichment dots just below the bars
  # geom_point(
  #   aes(x = coef, y = -0.2, size =  pmin(10,-log10(fisher_p_adj)), color = perc.overlap),
  #   shape = 16  # Dot shape
  # ) +
  # Customize the color and size scales
  # scale_color_gradient2(
  #   low = "#4C889C",
  #   mid = "white",
  #   high = "#D0154E")+
  # scale_size(range = c(0, 3), name = "adj.p.value") +  # Size by number of overlapping genes
  labs(
    #title = TeX("$Correlation of $\\log_{2}$(\\FC)$(ex-vivo vs in-vivo in NTC) vs KO-tissue interaction"),
    
    x = "KO",
    y = "Correlation (abs logFC)",
    fill = "Cell Type"
  ) +
  facet_grid(cols = vars(cell_type), scales = "free_x", space = "free_x") +
  optimized_theme_fig()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        strip.text.x = element_blank()) 
Fig4.5
###################################################################

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

################################################################################
# fgsea -------------------------------------------------------------------
# Initialize the result table
################################################################################

# Initialize the result table
gsea.res <- data.table() 

run_gsea <- function(limmaRes, enr.terms, celltypes = NULL, coefs = NULL) {
  
  # Initialize the result table
  gsea_res <- data.table() 
  
  # Determine cell types to process
  if (is.null(celltypes)) {
    celltypes <- unique(limmaRes$cell_type)
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
        subset_limmaRes <- limmaRes[limmaRes$cell_type == ct & limmaRes$coef == de_grp, ]
        
        # Extract statistics (logFC) and assign gene names as names
        stats <- with(subset_limmaRes, setNames(logFC, nm = gene))
        
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
gsea.res <- run_gsea(res, enr.terms, celltypes = unique(res$cell_type),
                     coefs =unique(res$coef))
gsea.res %>% write_rds(basedir("fgsea_hom_vs_ex.vivo_per_CT.rds"))

#############################################
gsea.res <- read_rds(basedir("fgsea_hom_vs_ex.vivo_per_CT.rds"))
unique(res$coef)
unique(gsea.res$coef)
FGSEA <- gsea.res
pDT <- FGSEA %>%
  #filter(coef %in%)
  mutate(genotype = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(celltype = gsub("M", "Macrophage", celltype)) %>%
  mutate(celltype = gsub("T8", "T-cells", celltype))%>%
  filter(coef != "(Intercept)")
unique(pDT$genotype)  
#pDT$genotype <- gsub("Interaction_","",pDT$genotype)
pDT$genotype <- gsub("treatmentex_vivo","WT",pDT$genotype )
# Step 2: Summarize to find KOs with at least one valid cell type
pDT <- pDT %>%
  filter(genotype %in% c(grep("Interaction",pDT$genotype, value = T),"WT"))
db = "MSigDB_Hallmark_2020"
pDT <- pDT %>%
  filter(db == "MSigDB_Hallmark_2020") 

# Step 3 continued: Keep only valid KOs for the specific cell type
pDT <- pDT %>% filter( padj < 0.05)

# mutate(alpha_value = if_else(valid_ko, 1, 0))  # Set alpha based on validity
pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5),by=c("coef", "celltype","pathway")]$pathway)
pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("coef", "celltype","pathway")]$pathway)
# Combine and remove duplicates across both positive and negative selections
pw.display <- unique(union(pw.display.pos, pw.display.neg))
pDT <- pDT[pathway %in% pw.display]
# Remove duplicate rows
pDT <- pDT %>% distinct()
pDT_agg <- pDT %>%
  group_by(pathway) %>%
  summarize(average_NES = mean(NES, na.rm = TRUE)) %>%
  arrange(desc(average_NES))  # Ordering pathways by the average NES, highest first

# Step 2: Create a factor for pathway that reflects the aggregated NES order
pDT$pathway <- factor(pDT$pathway, levels = pDT_agg$pathway)




fig4.1_supplementary <- ggplot(pDT, aes(x = gsub("Interaction_","",genotype), y = pathway,
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
  optimized_theme_fig() +
  theme(
    strip.text.x = element_text(angle = 90)
  )
fig4.1_supplementary
ggsave(basedir("Fig4.1_supplementary.pdf"), 
       plot = fig4.1_supplementary, 
       width = 8, 
       height = 13, 
       units = "cm")
#################################
# Assuming limmaRes is already loaded and contains a column named 'logFC'


unique(res$cell_type)
limmaRes <- res %>%
  mutate(genotype = gsub("genotype(.*):treatmentex_vivo", "Interaction_\\1", coef))%>%
  mutate(celltype = gsub("M", "Macrophage", cell_type)) %>%
  mutate(celltype = gsub("T8", "T-cells", cell_type))%>%
  filter(coef != "(Intercept)")
unique(limmaRes$genotype)  
test_irf9_num <- limmaRes %>%
  filter(celltype == "T-cells") %>%
  filter(genotype == "Interaction_IRF9KO")%>%
  filter(adj.P.Val < 0.05)%>%
  filter(abs(logFC) > 1)

# Basic histogram
hist(test_irf9$logFC,
     breaks = 50,
     col = "skyblue",
     main = "Distribution of logFC in Brd9 KO",
     xlab = "logFC",
     ylab = "Frequency")

# Add a vertical line at 0
abline(v = 0, col = "red", lwd = 2)

# Optional: Add thresholds (e.g. |logFC| > 1)
abline(v = c(-1, 1), col = "darkgreen", lwd = 2, lty = 2)

ggplot(limmaRes, aes(x = logFC)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", size = 1) +
  geom_vline(xintercept = c(-1, 1), color = "darkgreen", linetype = "dashed") +
  labs(title = "Distribution of logFC in Brd9 KO",
       x = "logFC",
       y = "Count") +
  theme_minimal()

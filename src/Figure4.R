source("src/00_init.R")
source("~/code/resources/RFunctions/Basics.R")
source("src/Ag_Optimized_theme_fig.R")
require(tidyverse)
require(data.table)
require(edgeR)
require(variancePartition)
require(pheatmap)
require(enrichR)
library(dplyr)
require(latex2exp)
library(patchwork)
library(fgsea)
# renv::snapshot(lockfile = "renv_NF.lock")


out <- "/media/AGFORTELNY/PROJECTS/TfCf_AG/Ag_ScRNA_22_JAKSTAT_Ar/"
basedir <- dirout("Figure4")

#*******************#
selected_coef <- c("Interaction_STAT1KO","Interaction_STAT2KO",
                   "Interaction_TYK2CMV", "Interaction_IRF9KO")
selected_KO <- gsub("Interaction_","",selected_coef)
#*
#*
res <- rbind(read.delim(paste0(out,"/DEG_ResultsT8.tsv")),
             read.delim(paste0(out,"/DEG_ResultsM.tsv"))
        )
res$probe <- res$rn
res$rn<-NULL
res <- as.data.frame(res)
gmap <- as.data.frame(read_rds(file = file.path(out, "DEG_GMP.RDS")))
res <- merge(res,gmap[,c("probe","gene")],by="probe")
res %>% write_rds(basedir("combined_jakstat_diff_exp.rds"))
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



#
#Fig4.Aa prep--------------------------------------------------------------
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
  mutate(valid_ko = ifelse(Count > 10, TRUE, FALSE))
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
WT$cell_type_label <- recode(WT$cell_type,
                             "Macrophage" = "M",
                             "T-cells" = "T")
Fig4Aa <- ggplot(WT, 
                  aes(y = cell_type,
                      x = ifelse(Regulation == "Downregulated",
                                 -log10(Count), log10(Count)), 
                      fill = Regulation)) + 
  geom_col() +
  scale_fill_manual(
    values = c("Upregulated" = "#D0154E", "Downregulated" = "#4C889C"),
    labels = c("Upregulated" = "Higher expr. ex vivo", 
               "Downregulated" = "Higher expr. in vivo")
  ) +
  labs(
    y = "Cell type",
    x = expression(atop("No. of genes", paste(log[10](n + 1)))),
    title = "No. of DEGs\n(culture effect)",
    fill = NULL
  ) +
  coord_flip() +
  optimized_theme_fig() +
  theme(
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0,          # Left-align
      size = 5, color = "black"    # Or "Times", "Courier", etc. (must be installed)
    )
  )



ggsave(basedir(paste0("Fig4Aa", ".pdf")), plot = Fig4Aa,w=6,
         height = 6, units = "cm")


#Fig4Ab--------
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

# Ensure pathway recoding is done
limma_subset$pathways <- dplyr::recode(limma_subset$pathways,
                                       "ISG core"= "ISG core",
                                       "mTORC1_or_Cholesterol" = "mTORC1/Cholesterol")

# Set as factor to control gene order if needed
limma_subset$pathways <- factor(limma_subset$pathways, 
                                levels = c("ISG core", "mTORC1/Cholesterol"))
# limma_subset$cell_type <- recode(limma_subset$cell_type,"Macrophage" = "M",
#                                  "T-cells"  = "T")
Fig4Ab <- limma_subset %>%
  ggplot(aes(
    x = cell_type, 
    y = gene, 
    color = pmin(2, pmax(-2, logFC)),
    size = pmin(5, -log10(adj.P.Val))
  )) +
  geom_point() + 
  scale_color_gradient2(
    low = "#4C889C",
    mid = "white",
    high = "#D0154E",
    name = TeX("log_{2}(FC)")
  ) +
  scale_size_continuous(
    range = c(0, 1.8),
    name = TeX("$-\\log_{10}(p_{adj})$")
  ) +
  labs(
    title = "Differentially expressed ISGs and\ngrowth/metabolic genes",
    x = "Cell type",
    y = "Genes"
    #caption = "M: Macrophages    T: T-cells"
  ) +
  facet_grid(
    rows = vars(pathways),       # << put facets side-by-side
    scales = "free_y", 
    space = "free"
  ) +
  optimized_theme_fig()+
  theme(
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0,          # Left-align
      size = 5, color = "black"    # Or "Times", "Courier", etc. (must be installed)
    ))
Fig4Ab
ggsave(basedir(paste0("Fig4Ab", ".pdf")), plot = Fig4Ab,w=4.5,
       height = 11.5, units = "cm")

#Fig4B--------
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
coefficients <- ko_flag %>%
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
selected_coef <- c("Interaction_STAT1KO","Interaction_STAT2KO",
                   "Interaction_TYK2CMV", "Interaction_IRF9KO")

limma_KO <- limma_KO %>%
  inner_join(ko_flag, by = c("coef","cell_type"))%>%
  filter(Count > 10) %>%
  filter(coef %in% c("Interaction_STAT1KO","Interaction_STAT2KO",
                     "Interaction_TYK2CMV", "Interaction_IRF9KO"))

Fig4B <- ggplot(limma_KO, aes(x = gsub("Interaction_","",coef), y = gene,
                              color = pmin(2, pmax(-2, logFC)) ,
                              size = pmin(5, -log10(adj.P.Val))))+
  geom_point() +  # Use geom_point to create dots
  scale_color_gradient2(low = "#4C889C",
                        mid = "white",
                        high = "#D0154E",
                        name = TeX("$log_{2}(FC)$"))+
  scale_size_continuous(
    range = c(0, 1.8),  # Set the actual size range from 2 to 30
    #breaks = c(1,2,4, 5, 10),
    name = TeX("$-\\log_{10}(p_{adj})$")# Set specific breaks to create distinct point sizes
  ) +
  labs(title = "ISG core (Interaction effect)",
       x = "KOs",
       y = "Genes"
  ) +
  facet_grid(rows = vars(pathways), cols = vars(cell_type),scales = "free", space = "free") +
  theme_bw() +
  optimized_theme_fig()
Fig4B 
ggsave(basedir(paste0("Fig4D.pdf")),plot = Fig4D, 
       width = 6,
       height =6 ,units = "cm")



#Fig4C----------
Interaction <- filtered_data %>%
  filter(coef != "WT") %>%
  filter(coef %in% selected_coef) %>%
  group_by(cell_type, coef) %>%
  summarize(Total_Regulated = sum(Count), .groups = "drop")
Fig4C <- ggplot(Interaction,aes(
  x = gsub("Interaction_","",coef),
  y = log10(Total_Regulated)
)) +
  geom_col() +
  
  # Facet by celltype with free space for flexibility in cell widths
  facet_grid(cols = vars(cell_type), space = "free", scales = "free") +
  labs(
    title = "No. of DEGs (Interaction effect)",
    x = "KOs",
    y = expression(atop("No. of genes", 
                        paste(log[10](n)))
  )) +
  # Custom theme with no legend if not needed
  optimized_theme_fig() + 
  theme(
    strip.text.x = element_text(angle = 0,
                                hjust = 0))


ggsave(basedir(paste0("Fig4C", ".pdf")), plot = Fig4B,w=5,
       height = 5, units = "cm")


#Fig4D-------

res$group <- ifelse(res$logFC >= 1 & 
                      res$adj.P.Val <= 0.05, "up", 
                    ifelse(res$logFC <= -1 & 
                             res$adj.P.Val <= 0.05, "down", "n.s"))

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

# Combine datasets
merged_data <- ex_vivo_data %>%
  inner_join(in_vivo_data, by = c("gene", "cell_type", "genotype"), 
             suffix = c("_ex.vivo", "_in.vivo"))
saveRDS(merged_data,basedir("in.vivo_ex.vivo_logFC_JAK_stat.rds"))

merged_data %>%
  filter(genotype == "IRF9KO")%>%
  filter(cell_type =="T8")%>%
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

#Fig.4D----------
Fig.4D <- merged_data %>%
  filter(genotype == "IRF9KO")%>%
  filter(cell_type=="T8")%>%
  ggplot(aes(x = logFC_ex.vivo, y = logFC_in.vivo)) +
  geom_hex(bins = 40, aes(fill = after_stat(count))) +
  scale_fill_gradient(low = "#d0e1f2", high = "#08306b") +  # Blue gradient
  geom_smooth(method = "lm", linetype = "dashed", se = FALSE, size = 0.8, color ="#e41a1c")+  # Linear regression line
  #scale_fill_viridis_c() +
  #scale_fill_viridis_c() +  # Optional: nice continuous color scale
  
  labs(
    title = "IRF9KO (T-cells)\nlogFC ex vivo versus in vivo",
    x = "logFC (Ex Vivo)",
    y = "logFC (In Vivo)",
    fill = "Gene Count"
  )+ 
  optimized_theme_fig() +
  theme(
    plot.title = element_text(hjust = 1, vjust = 1)
  )
ggsave(basedir("Fig.4D.pdf"),plot = Fig.4D, w=4.2, h= 3.5, units = "cm")

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
data <- deg_plot_data %>%
  filter(genotype %in% gsub("Interaction_","", selected_coef)) %>%
  filter(!(cell_type == "T8" & genotype == "STAT1KO"))
#Fig.4E-------
Fig.4E <- ggplot(data, aes(
  x = genotype,
  y = correlation)) +
  geom_col() +
  facet_grid(
    cols = vars(cell_type),
    labeller = labeller(cell_type = c("M" = "Macrophages", "T8" = "T cells")),
  
    scales = "free"
  ) +
  labs(
    title = "KO effects (in vivo versus ex vivo)",
    x = "KOs",
    y = "Pearson's correlation"
  ) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  optimized_theme_fig() +
  theme(
    legend.justification = "right"
  )
ggsave(basedir("Fig.4E.pdf"),plot = Fig.4E, width = 5,
       height =4,units = "cm")

#######Fig.4D+$4Cb---------------------


#Fig4F-----
#

#Normalized_reads
dataVoom_M <- read.csv("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Normalized_reads/Normalized_readsM.csv",
                   header = T,
                   row.names = 1)
dataVoom_T8 <- read.csv("/media/AGFORTELNY/PROJECTS/TfCf_AG/JAKSTAT/Normalized_reads/Normalized_readsT8.csv",
                   header = T,
                   row.names = 1)


gmap <- as.data.frame(gmap)

dataVoom_M <- dataVoom_M %>%
  rownames_to_column(var = "ensG") %>%     # Convert rownames to a column for merging
  left_join(gmap[, c("ensG", "gene")], by = "ensG") %>%  # Merge with gmap to get the gene names
  column_to_rownames(var = "gene") %>% 
  dplyr::select(-ensG)


dataVoom_T8 <- dataVoom_T8 %>%
  rownames_to_column(var = "ensG") %>%
  left_join(gmap[, c("ensG", "gene")], by = "ensG") %>%
  column_to_rownames(var = "gene") %>%
  dplyr::select(-ensG)

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
# koi <-  ko_flag %>% # Remove grouping from previous steps
#   mutate(comparison = gsub("Interaction_", "", coef)) %>%
#   mutate(cell_type = cell_type) %>%
#   dplyr::select(cell_type, coef, valid_ko, Count,comparison)

filtered_data <- goi_exp %>%
  filter(comparison %in% gsub("Interaction_","",selected_coef))%>%
  filter(gene == goi)

filtered_data$coef <- factor(filtered_data$coef,
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
write_rds(filtered_data, basedir("filtered_data.rds"))

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

oas1a_data <- filtered_data %>%
  filter(gene == "Oas1a") %>%
  dplyr::select(sample_name, genotype, cell_type, tissue, E)
oas1a_data %>%
  dplyr::count(sample_name, genotype, cell_type, tissue) %>%
  filter(n > 1)

oas1a_summary <- oas1a_data %>%
  group_by(genotype, cell_type, tissue) %>%
  summarise(E = mean(E, na.rm = TRUE), .groups = "drop")

paired_data <- oas1a_summary %>%
  pivot_wider(names_from = tissue, values_from = E) %>%
  filter(!is.na(in_vivo) & !is.na(ex_vivo))

paired_data <- paired_data %>%
  mutate(celltype = paste(cell_type,genotype,sep=":"))
my_color <- c("#004949", "#009292", "#490092", "#006DDB", "#B66DFF")

    
Fig4F <- ggplot(paired_data, aes(
  x = ex_vivo, 
  y = in_vivo, 
  color = genotype,
  shape = cell_type
)) +
  geom_point(size = 1, position = position_jitter(width = 0.02, height = 0.02)) +  # add small jitter if you want
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = my_color) +
  labs(
    x = "Ex vivo expression",
    y = "In vivo expression",
    title = "Oas1a Expression\nIn Vivo versus Ex Vivo (Mean)",
    color = "Genotype",
    shape = "Cell type"
  ) +
  optimized_theme_fig() +
  theme(legend.position = "right",
        strip.text = element_text(size = 7))

ggsave(basedir("Fig4F_combined.pdf"), plot = Fig4F, width = 5.5,
       height = 3.5 ,units = "cm")


row2 <- (Fig4D | Fig4E) +
  plot_layout(widths = c(1,1.2))
ggsave(basedir("row2.pdf"), plot = row2, width = 15, height = 7, units = "cm")


##########################


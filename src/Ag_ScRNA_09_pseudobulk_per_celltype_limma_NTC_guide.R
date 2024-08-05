#limma_NTC
#Comparison:ex.vivo-in.vivo

###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(ggrepel)

##################################################################################
inDir<-dirout("Ag_ScRNA_08_Pseudobulk_limma_guide")

base<-"Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/"
basedir<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
########################
##################################################################################
#load data
########################
#metadata
meta<-read.delim(inDir("metadata_guide.tsv"),row.names=1)
counts <- read.delim(inDir("combined_in_ex_counts_guide.tsv"), row.names = 1)

#meta data
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","T-cell","Gran.","MEP")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]
#only select the genotypes present in both tissue conditions
meta<-meta[meta$genotype %in% meta[meta$tissue=="ex.vivo",]$genotype,]

# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)-]", ".", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta<-meta%>%filter(!grepl("NA",rownames(meta)))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))
#selecting only NTC
NTC_counts<-counts[,grep("NTC",colnames(counts),value = T)]
NTC_meta<-meta[grep("NTC",rownames(meta),value = T),]
inDir1<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide/")
NTC_meta%>%write_rds(inDir1("NTC_meta.rds"))
counts<-counts[,rownames(NTC_meta)]
stopifnot(all(colnames(counts)==rownames(NTC_meta)))
meta<-NULL

##################################
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group <- interaction(NTC_meta$celltype, NTC_meta$tissue)
unique(group)

mm <- model.matrix(~0 + group)
rownames(mm) <- rownames(NTC_meta)

# Normalization
dataVoom <- voom(d, mm)
dataVoom%>% write_rds(basedir("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaFit <- lmFit(dataVoom, mm)
limmaFit <- eBayes(limmaFit)

# Generate comparisons based on tissue types
tissue_type1<-"ex.vivo"
tissue_type2<-"in.vivo"
comparisons <- lapply(unique(NTC_meta$celltype), function(celltype) {
  comp1 <- paste0("group", celltype, ".",tissue_type1)
  comp2 <- paste0("group", celltype, ".", tissue_type2)
  list(comp1, comp2, celltype)
})

# Initialize list to store top table results
top_table <- list()

# Perform differential expression analysis for each comparison
for (comp_list in comparisons) {
  group1 <- comp_list[[1]]
  group2 <- comp_list[[2]]
  id <- comp_list[[3]]
  
  # Create contrast matrix
  contrast <- paste0(group1, "-", group2)
  
  # Debugging: Print contrast matrix
  print(contrast)
  
  # Ensure that contrast inputs are accepted
  if (!(all(c(group1, group2) %in% colnames(limmaFit)))) {
    stop("Invalid column names for comparison:", group1, " and ", group2)
  }
  
  # Fit model and perform differential expression analysis
  tmp <- contrasts.fit(limmaFit, makeContrasts(contrast, levels = colnames(coef(limmaFit))))
  tmp <- eBayes(tmp)
  top_table[[id]] <- topTable(tmp, sort.by = "P", n = Inf)
}

# Bind results and add NTC_metadata
ex_in_NTC_per_ct <- bind_rows(top_table, .id = 'celltype') %>%
  mutate(genes = gsub("\\...\\d+", "", rownames(.)))

# Create a column 'group' to group the genes as up, down, or ns
ex_in_NTC_per_ct$group <- ifelse(ex_in_NTC_per_ct$logFC >= 1 & 
                                ex_in_NTC_per_ct$adj.P.Val <= 0.05, "up", 
                              ifelse(ex_in_NTC_per_ct$logFC <= -1 & 
                                       ex_in_NTC_per_ct$adj.P.Val <= 0.05, "down", "n.s"))

ex_in_NTC_per_ct %>% write_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
################################################################################
# Plotting and saving-----------------------------------------------------------
################################################################################
ex_in_NTC_per_ct <- read_rds(basedir("limma_perCTex.vivovsin.vivo.rds"))
ex_in_NTC_per_ct %>%
  ggplot(aes(x = logFC, y = pmin(10, -log10(adj.P.Val)), col = group)) +
  geom_point() +
  scale_color_manual(values = c("#5782A7", "#B9B8B6", "#8A264A")) +
  ggtitle(paste0(tissue_type1, " ", tissue_type2, " NTC D.E")) +
  facet_wrap(vars(celltype), scales = "free_x") +
  theme_bw()+
  theme(strip.text = element_text(size = 18))


ggsave(basedir(paste0(tissue_type1, "_", tissue_type2, "_NTC_D.E_percelltype_trial.pdf")))

# Assuming ex_in_NTC_per_ct is already loaded and 'group' column is created
# Assuming ex_in_NTC_per_ct is already loaded and 'group' column is created
# Filter out 'n.s' group
filtered_data <- ex_in_NTC_per_ct %>%
  filter(group != "n.s")

# Count the number of up and down genes for each cell type
gene_counts <- filtered_data %>%
  group_by(celltype, group) %>%
  summarise(count = n()) %>%
  ungroup()

# Log10-transform the counts, adding 1 to avoid log(0)
gene_counts <- gene_counts %>%
  mutate(log10_count = log10(count + 1))

# Plotting
ggplot(gene_counts, aes(x = celltype, y = log10_count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("down" = "#5782A7", "up" = "#8A264A")) +
  labs(title = "NTC D.E_genes",
       x = "Cell Type",
       y = "Log10(Number of Genes + 1)",
       fill = "Regulation") +
  theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15),
        title = element_text(size=18))+
  theme(axis.text.x = element_blank()) +
  facet_wrap(~ celltype, scales = "free_x")+
  theme(strip.text = element_text(size = 18))
ggsave(basedir(paste0("NTC_D.E_genes_percelltype.pdf")))
##########################################

#fig:1
# Identify top genes for labeling (you can adjust criteria here)
top_genes <- ex_in_NTC_per_ct %>%
  filter(group != "n.s") %>%
  group_by(celltype, group) %>%
  top_n(3, wt = -adj.P.Val) # Select top 3 by adjusted p-value for each group

# Plot
ggplot() +
  # Hexbin plot for the "others" group
  stat_bin_hex(data = filter(ex_in_NTC_per_ct, group == "n.s"), 
               aes(x = logFC, y = -log10(adj.P.Val), fill = ..count..), 
               bins = 20, color = NA, alpha = 0.7) +
  scale_fill_gradient(low = "lightgrey", high = "black",
                      limits = c(1, 5000), name = "Gene Count") +
  
  # Overlay points for the main groups
  geom_point(data = filter(ex_in_NTC_per_ct, group != "n.s"), 
             aes(x = logFC, y = -log10(adj.P.Val), color = group), 
             alpha = 0.9, size = 2.5) +
  
  # # Add text labels for top genes with ggrepel
  # geom_text_repel(data = top_genes, 
  #                 aes(x = logFC, y = pmin(10, -log10(adj.P.Val)), label = genes), 
  #                 size = 3, color = "black", 
  #                 max.overlaps = 10, # Adjust the number of overlaps allowed
  #                 box.padding = 0.5, 
  #                 point.padding = 0.5) + # Adjust padding for label positions
  # 
  # Manually setting colors for groups
  scale_color_manual(values = c(
    "up" = "#D0154E", 
    "down" = "#4C889C"
  ), name = "Group") +
  
  labs(title = "Ex-vivo vs In-vivo",
       x = "logFC",
       y = "-log10(adj.P)") +
  facet_wrap(vars(celltype), scales = "free") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
  
ggsave(basedir(paste0("NTC_D.E_genes_percelltype.pdf")))

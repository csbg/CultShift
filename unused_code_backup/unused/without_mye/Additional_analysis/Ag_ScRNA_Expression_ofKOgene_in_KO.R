####################################################################
#***Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01***#
####################################################################
#limma with all 
#data
#
#Aarathy
###############
###############
source("src/00_init.R")
source("src/Ag_Optimized_theme.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
library(ComplexHeatmap)
#####################################################################
InDir<-dirout("/Ag_ScRNA_08_Pseudobulk_limma_guide")
base<-"Ag_ScRNA_Expression_in_KO/"
basedir<-dirout("Ag_ScRNA_Expression_in_KO/")
inDir1<-dirout_load("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC_guide")
#####################################################################
#load data and clean metadata
################################################
genes_to_exclude <- c("B2m","S100a11","Actg1","Sri","Ly6e","Vamp8","Mt1","Hba-a1",
                      "Hba-a2","Pim1","Fabp5","Fdps","Cd9")
#metadata
meta<-read.delim(InDir("metadata_guide.tsv"),row.names=1)

#exclude given celltypes
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear","Gran.","MEP")
meta <- meta[!(meta$celltype %in% celltypes_to_exclude), ]

#only select the genotypes present in both tissue conditions
#meta<-meta[meta$genotype %in% meta[meta$tissue=="ex.vivo",]$genotype,]
#meta %>%filter(tissue =="ex.vivo")%>% pull(genotype)%>%unique()
#meta<-meta[meta$genotype %in% meta[meta$tissue=="in.vivo",]$genotype,]

# Replace space (\\s), left parenthesis (\\(), right parenthesis (\\)), or hyphen (-)
rownames(meta) <- gsub("[\\ \\(\\)-]", ".", rownames(meta))
rownames(meta) <- gsub("[\\ \\(\\)-]", ".", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in row names
rownames(meta) <- gsub("Eo/Ba", "Eo.Ba", rownames(meta))

# Replace "Eo/Ba" with "Eo.Ba" in all relevant columns
meta[] <- lapply(meta, gsub, pattern = "Eo/Ba", replacement = "Eo.Ba")
meta<-meta%>%filter(!grepl("NA",rownames(meta)))


meta%>%write_rds(basedir("meta.rds"))
#counts
counts <- read.delim(InDir("combined_in_ex_counts_guide.tsv"), row.names = 1)
counts<-counts[!(rownames(counts) %in% genes_to_exclude),rownames(meta)]

stopifnot(all(colnames(counts)==rownames(meta)))

stopifnot(setequal(unique(meta[meta$tissue == "ex.vivo", ]$celltype), 
                   unique(meta[meta$tissue == "in.vivo", ]$celltype)))

################################################
#factors and levels
################################################
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))

unique_celltypes <- unique(meta$celltype)
for (ct in unique_celltypes) {
  cat("Performing DE analysis for cell type:", ct, "\n")
  
  # Subset metadata and counts data for the current cell type
  meta_subset <- meta %>% filter(celltype == ct)
  counts_subset <- counts[, rownames(meta_subset)]
  stopifnot(all(colnames(counts_subset)==rownames(meta_subset)))
# Prepare the data for differential expression analysis
  d0 <- DGEList(counts_subset)
  d0 <- calcNormFactors(d0,method = "TMM")
#threshold <- 30
#drop <- which(apply(cpm(d0), 1, max) < threshold)
#d <- d0[-drop, ]
#########################

###############################################################################
#model matrix
  meta_subset <- as.data.frame(meta_subset)
  modelMatrix <- model.matrix(~ tissue * genotype, data = meta_subset)
  dataVoom <- voom(d0, modelMatrix)
  saveRDS(dataVoom, basedir(paste0(ct, "_dataVoom.rds")))
}

###########################
dataVoom_Eo.Ba<-read_rds(basedir("Eo.Ba_dataVoom.rds"))
dataVoom_Mono<-read_rds(basedir("Mono_dataVoom.rds"))
dataVoom_MkP<-read_rds(basedir("MkP_dataVoom.rds"))
dataVoom_GMP<-read_rds(basedir("GMP_dataVoom.rds"))
dataVoom_HSC<-read_rds(basedir("HSC_dataVoom.rds"))
dataVoom_MEP<-read_rds(basedir("MEP_dataVoom.rds"))


meta<-read_rds(basedir("meta.rds"))

#Theses are the KO with significant upregulated interaction terms
list_guide <- unique(meta$guide)

dat.list <- list()
KO  <- list_guide[2]
for (KO in list_guide){
  
  for (ct in unique(meta$celltype)) {
    # Get the dataVoom object corresponding to the current cell type
    dataVoom_ct <- get(paste0("dataVoom_", ct))
    #CAPITALIZE
    goi <- gsub("_.*","",KO)
    # Check if goi exists in the row names of dataVoom_ct$E
    if (goi %in% rownames(dataVoom_ct$E)){
    # Subset the metadata and E values for the current gene and cell type
      gene_data <- meta[meta$celltype == ct,] %>%
      mutate(E = dataVoom_ct$E[goi,]) %>%
      rownames_to_column("sample1") %>%
      filter(guide %in% c(KO, "NTC")) %>%
      mutate(scaled_E = scale(E)) %>%
      mutate(gene = goi)%>%
      mutate(celltype=ct)%>%
      mutate(comparison=KO)
          
      # Store the gene data in the list
    dat.list[[paste0(ct, "_", goi,KO)]] <- gene_data
    }
  }
}

KO_list<-bind_rows(dat.list,.id = "celltype_gene_guide")
KO_list %>% write_rds(basedir("Norm_exp.rds"))
KO_list <- read_rds(basedir("Norm_exp.rds"))

list_guide <- list_gene %>% na.omit()
list_gene <- gsub("_.*","",list_guide) %>% na.omit()


# Summarize the data to calculate mean scaled expression for each combination of guide, gene, and condition
KO_summary <- KO_list %>%
  group_by(tissue,gene, guide) %>%
  summarize(mean_scaled_E = mean(scaled_E, na.rm = TRUE)) %>%
  ungroup()

# Plot the data
for(tis in unique(KO_summary$tissue)){
    data <- KO_summary %>% filter(tissue == tis)
    ggplot(data, aes(x = guide, y = gene)) +
      geom_point(aes(color = mean_scaled_E), size = 4) +  # Use color for points instead of fill
      scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                            name = "Mean Scaled Expression") +  # Color scale for mean expression
      theme_minimal() +
      labs(
        title = paste("Mean Scaled Expression ","(", tissue, ")"),
        x = "Guide",
        y = "Gene"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      )+
      optimized_theme()
    ggsave(basedir(paste0(tis,"summary_KO_plot.pdf")), w=30,h=20)
    
}


  


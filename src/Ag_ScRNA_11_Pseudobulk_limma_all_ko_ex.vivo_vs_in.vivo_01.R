####################################################################
#***Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01***#
####################################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#14-03-24
#Aarathy
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
#####################################################################
inDir<-dirout("/Ag_ScRNA_08_Pseudobulk")
base<-"Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01/"
basedir<-dirout("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_01/")
#####################################################################
#source DE function
#####################################################################
source("src/Ag_ScRNA_11_invivo_exvivo_KO_limma_function.R")
################################################
#load data and clean metadata
################################################
#metadata
meta<-read.csv(inDir("metadata.csv"),row.names=1)

#exclude given celltypes
celltypes_to_exclude <- c("B-cell", "CLP", "Ery", "EBMP", "unclear")
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
#counts
counts <- read.csv(inDir("ex_counts.csv"), row.names = 1)
counts<-cbind(counts,read.csv(inDir("in_counts.csv"), row.names = 1))
counts<-counts[,rownames(meta)]
stopifnot(all(colnames(counts)==rownames(meta)))
################################################
#factors and levels
################################################
meta$genotype <- factor(meta$genotype, levels=c("NTC", unique(setdiff(meta$genotype,"NTC"))))
meta$tissue <- factor(meta$tissue, levels=c("in.vivo", "ex.vivo"))
################################################
#perform DE independent of celltype
model_formula <- "~tissue*genotype"
limmaRes <- performDE(meta, counts,model_formula)

limmaRes%>%write_rds(basedir("limma_ex.vivo_vs_in.vivo_all_CT.rds"))
pmap(unique(limmaRes$coef),~{
  coefx<-.x
  data<-limmaRes[limmaRes$coef==coefx,]
  ggplot(data,aes(x=logFC,y=-log10(adj.P.Val),col=group))+
    geom_point()+
    scale_color_manual(values =  c("#5782A7", "#B9B8B6", "#8A264A"))+
    ggtitle(paste0(gsub("ex.vivo:","interaction:",coefx),"D.E"))
  # facet_wrap(vars(data$celltype),scale="free")
  out <- dirout(paste0(base, coefx))
  ggsave(out("DE.pdf"))
})

#############################################################
#ex.vivo condition
#############################################################
#volcano plot
coefx<-"ex.vivo"
data<-limmaRes[limmaRes$coef=="ex.vivo",]
ggplot(data,aes(x=logFC,y=-log10(adj.P.Val),col=group))+
  geom_point()+
  scale_color_manual(values = c("#5782A7", "#B9B8B6", "#8A264A"))+
  ggtitle(paste0("ex.vivo","_D.E"))
# facet_wrap(vars(data$celltype),scale="free")

out <- dirout(paste0(base, coefx))
ggsave(out("DE.pdf"))

#enrichment

source("src/Ag_enrichment.R")

perform_enrichment_analysis(data=data,
                            logFC_cutoff_up= 1,
                            logFC_cutoff_down= -3,
                            databases=c("KEGG_2019_Mouse", 
                                        "MSigDB_Hallmark_2020", 
                                        "WikiPathways_2019_Mouse", 
                                        "GO_Biological_Process_2021"),
                            output_file_prefix="all_CT")


####################################################

#interaction
data<-limmaRes[limmaRes$coef %in% grep("ex.vivo:",limmaRes$coef,value=T),]

# Filter data based on the condition
filtered_data <- data %>%
  filter(adj.P.Val < 0.01) %>%
  group_by(coef) %>%
  summarize(genes = list(unique(ifelse(abs(logFC) > threshold, genes, NA)))) %>%
  unnest(genes)

# Count occurrences of each gene
gene_counts <- unnested_data %>%
  count(genes)

# Filter genes that occur in multiple coefficients
genes_overlap <- gene_counts %>%
  filter(n > 1) %>%
  pull(genes)

# Resulting list of genes that overlap between coefficients
genes_overlap

goi.all<-data %>%
  filter(adj.P.Val<0.01) %>%
  group_by(coef) %>%
  top_n(20,abs(logFC)) %>%
  pull(genes)%>%
  unique()
length(goi.all)

(p.coef <- data %>%
    filter(genes %in% goi.all) %>%
    ggplot(aes(y=genes, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 5)) 
dat.list <- list()
for(gg in goi.all){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
head(dataVoom$E)

(p.vals <- bind_rows(dat.list, .id="genes") %>%
    
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ cell, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
ggsave(out(""))
################################################
#perform DE per celltype
################################################

# Use purrr functions to iterate over unique KO and celltypes
ct<-unique(meta$celltype)[1]

map(unique(meta$celltype), ~ {
  ct <- .x
  meta_ct <- meta[meta$celltype==ct,] 
  
  if(nrow(meta_ct[meta_ct$tissue=="in.vivo",])>1 &
     nrow(meta_ct[meta_ct$tissue=="ex.vivo",])>1){
    meta_ct$genotype <- factor(meta_ct$genotype,
                               levels=c("NTC",unique(setdiff(meta_ct$genotype,"NTC"))))
    
    meta_ct$tissue <- factor(meta_ct$tissue, levels=c("in.vivo", "ex.vivo"))
    
    # Call the function for differential expression analysis
    limmaRes_ct <- performDE(meta_ct, counts,model_formula)
    
    # Add celltype column to the combined results
    limmaRes_ct$celltype <- ct
    
    # Combine with the existing limmaRes data frame
    limmaRes <<- rbind(limmaRes, limmaRes_ct)
    limmaRes$group <- ifelse(limmaRes$logFC >= 1 &
                               limmaRes$adj.P.Val <= 0.05, "up",
                             ifelse(limmaRes$logFC <= -1 &
                                      limmaRes$adj.P.Val <= 0.05, "down", "n.s"))
  }else{
    print("null")
  }
  limmaRes
})

# Assuming your data frame is called 'data'

# Filter data based on the condition abs(logFC) > 1
gene_counts <- data %>%
  filter(abs(logFC) > 1)%>%
  group_by(genes) %>%
  summarize(unique_coefs = n_distinct(coef))

# Filter genes where the number of unique coefficients is greater than 1
overlap_genes <- gene_counts %>%
  filter(unique_coefs > 20) %>%
  pull(genes)

# Resulting list of genes that overlap across multiple coefficients with abs(logFC) > 1
overlap_genes
head(data)

pheatmap(data[genes %in% overlap_genes,])

###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
require(fgsea)
library(msigdbr)

#load data
inDir<-dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype/")
inDir1<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC/")
inDir2<-dirout_load("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment/ENRICHR")
base<-"Ag_ScRNA_16_per_celltype_gene_plots/"
basedir<-dirout("Ag_ScRNA_16_per_celltype_gene_plots")


IFN_inter<-dirout(paste0(base,"IFN_inter"))
IFN_base<-"Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/"
################################################################################
#load limma results from NTC
dataVoom_NTC<-read_rds(inDir1("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaRes_NTC<-read_rds(inDir1("limma_perCTex.vivovsin.vivo.rds"))
#load enrichment_results from NTC
enr.res.all_up<-read_rds(inDir2("enr.res.all_NTC_up.rds"))
enr.res.all_down<-read_rds(inDir2("enr.res.all_NTC_down.rds"))

#load_limma across KO
limmaRes<-read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
dataVoom_Eo.Ba<-read_rds(inDir("Eo.Ba_dataVoom.rds"))
dataVoom_Mono<-read_rds(inDir("Mono_dataVoom.rds"))
dataVoom_MkP<-read_rds(inDir("MkP_dataVoom.rds"))
dataVoom_GMP<-read_rds(inDir("GMP_dataVoom.rds"))
dataVoom_HSC<-read_rds(inDir("HSC_dataVoom.rds"))
dataVoom_MEP<-read_rds(inDir("MEP_dataVoom.rds"))
NTC_meta<-read_rds(inDir1("NTC_meta.rds"))
meta<-read_rds(inDir1("meta.rds"))
################################################################################
#Function to extract genes from terms
extract_unique_genes <- function(dataframe, terms) {
  # Filter the dataframe for rows where the Term is in the specified list of terms
  filtered_table <- subset(dataframe, Term %in% terms)
  
  # Extract the genes column from the filtered dataframe
  genes <- filtered_table$Genes
  
  # Split the genes column by ";" and flatten the resulting list
  all_genes <- unlist(strsplit(genes, ";"))
  
  # Select only unique genes
  unique_genes <- unique(all_genes)
  
  return(unique_genes)
}
################################################################################
Cholesterol_related<-
  enr.res.all_up%>%
  filter(Term %in% grep("cholesterol",enr.res.all_up$Term,value = T))%>%
  pull(Genes)%>%
  gsub(";",",")%>%unlist%>%
  unique()
################################################################################
#IFN_genes
IFN_genes <- extract_unique_genes(enr.res.all_down, c("Interferon Alpha Response",
                                                                               "Interferon Gamma Response"),
                                                           "MSigDB_Hallmark_2020")
#Cholesterol homeostasis
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, c("Cholesterol Homeostasis"))
#
################################################################################
#Interaction upregulated------------------------------------------------------
 ################################################################################
 # logFC-----------
 ################################################################################
list_ko <- limmaRes%>%
  filter(logFC>1 & adj.P.Val<0.05)%>%
  pull(coef)%>%
  str_replace("ex.vivo:","")%>%
  unique()

for (KO in list_ko){
  IFN_genes_KO<-limmaRes %>%
    filter(coef == paste0("ex.vivo:",KO)) %>%
    #filter(adj.P.Val < 0.05) %>%
    filter(logFC > 1) %>%
    filter(toupper(ensg) %in% IFN_genes)
  
  
  ggplot(IFN_genes_KO,aes(x=celltype,y=ensg,size=pmin(-log10(adj.P.Val),5),
                            color=logFC))+
    
    geom_point()+scale_color_gradient2(high="red", low="blue")+
    theme(axis.text = element_text(size = 15)) +
    #facet_wrap(vars(celltype), scales = "free_x")+
    theme_bw(12)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
    scale_size_continuous(range=c(0,5), limits = c(0,5))+
    ggtitle("Brd9KO IFN genes-Interaction(ex.vivo-in.vivo)")+
    ylab("IFN response genes")
  ko_dir<-dirout(paste0(base,"/IFN_inter/",KO))
  ggsave(ko_dir("Brd9KO_interaction_upregulated_IFN_genes.pdf"),h=nrow(Brd9_IFN_genes)*0.04+1)
  
}

###################
#Expression Upregulated
##############
#NOTE:CHECK IF CAPITALIZED
gene_list<-IFN_genes
#Theses are the KO with significant upregulated interaction terms
list_ko <- limmaRes%>%
  filter(logFC>1 & adj.P.Val<0.05)%>%
  pull(coef)%>%
  str_replace("ex.vivo:","")%>%
  unique()

dat.list <- list()
for (KO in list_ko){
  list_of_genes<-limmaRes%>%
    #FILTERING ONLY SIGNIFICANT-UP IFN genes for that particular KO
    filter(coef==paste0("ex.vivo:",KO))%>%
             # ONLY UPREGULATED GENES
             filter(logFC>1 & adj.P.Val<0.05)%>%
    filter(toupper(ensg) %in% gene_list)%>%
    pull(ensg)
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
#Here we have sign.upregulated genes for each KO in all celltypes(not only the one with sig.up) of that KO.
limma_KO_up_IFN_genes<-bind_rows(dat.list,.id = "celltype_gene_KO")
limma_KO_up_IFN_genes%>%write_rds(basedir("IFN_genes_in_each_ko_upregulated_for all_celltypes"))
head(limma_KO_up_IFN_genes)
######################################################
# IFN genes-------------------------------------------
######################################################
IFN_UP<-dirout(paste0(base,"IFN_inter/UP/"))
dirout_load(IFN_UP)
# Function to create plots for each gene
create_gene_plots <- function(data, gene,KO) {
  ggplot(data[data$gene == gene,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype~tissue, scales = "free") +
    labs(title = gene)+
    xlab(paste0(KO,"KO"))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/UP/",KO))
  ggsave(ko_dir(paste0(gene,".pdf")))
}
#For each KO where there were significantly up genes, make the plot
for (comp in unique(limma_KO_up_IFN_genes$comparison)){
  data<-limma_KO_up_IFN_genes%>%
    filter(comparison==comp)
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots(data, gene,comp)
    
  })
}
# Generate plots for each gene


#####################################################################
#Barplopt
# Group data by celltype, tissue, and genotype, calculate the mean scaled_E
# Define a custom function to calculate median and return a data frame
calculate_median <- function(x) {
  return(data.frame(y = median(x)))
}
for (comp in unique(limma_KO_up_IFN_genes$comparison)){
  data<-limma_KO_up_IFN_genes%>%
    filter(comparison==comp)
  #calculate mean  
  mean_data <- data %>%
    group_by(celltype, tissue, genotype) %>%
    summarise(mean_scaled_E = mean(scaled_E))
    
  # Calculate median scaled_E for each celltype, tissue, and genotype
  median_data <- data %>%
      group_by(celltype, tissue, genotype) %>%
      summarise(median_scaled_E = median(scaled_E))
  # Create the violin plot
  ggplot(data, aes(x = tissue, y = scaled_E, fill = genotype)) +
    #geom_violin(trim = FALSE) +
    geom_boxplot(aes(middle = mean(scaled_E)))+
    #geom_jitter()+
    facet_wrap(~ celltype, scales = "free") +
    labs(x = "Tissue", y = "Scaled E", fill = "Genotype", title = paste0("IFN genes in upregulated_",comp,"_interaction")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/UP/",comp,"Summary_IFN/"))
  ggsave(ko_dir("Boxplot_IFN_genes_in_Brd9_up.pdf"),h=min(20,nrow(data)*0.04+1))
    
}   

# Plot mean scaled_E for each celltype in each tissue
ggplot(mean_data, aes(x = celltype, y = mean_scaled_E, fill = genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ tissue, scales = "free") +
  labs(x = "Cell Type", y = "Mean Scaled E", fill = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(IFN_inter(paste0("mean_of_scaled_E_across_IFN_gene",".pdf")))
###############################################################
#####################################################################
#
#
#####################################################################
#library(ggplot2)

# Define a custom function to calculate median and return a data frame
calculate_median <- function(x) {
  return(data.frame(y = median(x)))
}

# Calculate median scaled_E for each celltype, tissue, and genotype
median_data <- data_Brd9 %>%
  group_by(celltype, tissue, genotype) %>%
  summarise(median_scaled_E = median(scaled_E))


# Create the violin plot
ggplot(data_Brd9, aes(x = tissue, y = scaled_E, fill = genotype)) +
  geom_violin(trim = FALSE) +
  #geom_hline(data = median_data, aes(yintercept = median_scaled_E),
  #color = "black", linetype = "solid", size = 1) +
  #geom_point(data = median_data, aes(y = median_scaled_E), 
  # color = "black", size = 3, shape = 95) +
  #geom_point(data = median_data, aes(x = celltype, y = median_scaled_E),
  #color = "black", size = 3, shape = 95) +
  facet_wrap(~ celltype, scales = "free") +
  labs(x = "Cell Type", y = "Scaled E", fill = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(IFN_inter(""))

#####
# Boxplot
# Calculate median scaled_E for each genotype, tissue, and celltype
median_data <- data_Brd9 %>%
  group_by(genotype, tissue, celltype) %>%
  summarise(median_scaled_E = median(scaled_E))

# Create the violin plot
ggplot(data_Brd9, aes(x = tissue, y = scaled_E, fill = genotype)) +
  #geom_violin(trim = FALSE) +
  geom_boxplot(aes(middle = mean(scaled_E)))+
  #geom_jitter()+
  facet_wrap(~ celltype, scales = "free") +
  labs(x = "Tissue", y = "Scaled E", fill = "Genotype", title = "IFN genes in upregulated Brd9KO interaction.pdf") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(IFN_inter("Boxplot_IFN_genes_in_Brd9_up.pdf"))
################################################################################
#
Cholesterol_inter<-dirout(paste0(base,"Cholesterol_inter"))
#Cholesterol homeostasis
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, c("Cholesterol Homeostasis"))
################################################################################
#Interaction upregulated------------------------------------------------------
################################################################################
# logFC-----------

###################
#Expression Upregulated
##############
#NOTE:CHECK IF CAPITALIZED
gene_list<-Cholesterol_genes
#Theses are the KO with significant upregulated interaction terms
list_ko <- limmaRes%>%
  filter(logFC< -1 & adj.P.Val<0.05)%>%
  pull(coef)%>%
  str_replace("ex.vivo:","")%>%
  unique()

dat.list <- list()
for (KO in list_ko){
  list_of_genes<-limmaRes%>%
    #FILTERING ONLY SIGNIFICANT-UP IFN genes for that particular KO
    filter(coef==paste0("ex.vivo:",KO))%>%
    # ONLY UPREGULATED GENES
    filter(logFC< -1 & adj.P.Val<0.05)%>%
    filter(toupper(ensg) %in% gene_list)%>%
    pull(ensg)
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
#Here we have sign.upregulated genes for each KO in all celltypes(not only the one with sig.up) of that KO.
limma_KO_down_Cholesterol_genes<-bind_rows(dat.list,.id = "celltype_gene_KO")
limma_KO_down_Cholesterol_genes%>%write_rds(basedir("Cholesterol_genes_in_each_ko_downregulated_for all_celltypes"))

######################################################
# Cholesterol_genes-------------------------------------------
######################################################
Cholesterol_down<-dirout(paste0(base,"Cholesterol_Inter/DOWN/"))

# Function to create plots for each gene
create_gene_plots <- function(data, gene,KO) {
  ggplot(data[data$gene == gene,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype~tissue, scales = "free") +
    labs(title = gene)+
    xlab(paste0(KO,"KO"))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/Cholesterol_Inter/DOWN/",KO))
  ggsave(ko_dir(paste0(gene,".pdf")))
}
#For each KO where there were significantly up genes, make the plot
for (comp in unique(limma_KO_down_Cholesterol_genes$comparison)){
  data<-limma_KO_down_Cholesterol_genes%>%
    filter(comparison==comp)
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots(data, gene,comp)
    
  })
}
# Generate plots for each gene


#####################################################################
#Barplopt
# Group data by celltype, tissue, and genotype, calculate the mean scaled_E
# Define a custom function to calculate median and return a data frame
calculate_median <- function(x) {
  return(data.frame(y = median(x)))
}
for (comp in unique(limma_KO_down_Cholesterol_genes$comparison)){
  data<-limma_KO_down_Cholesterol_genes%>%
    filter(comparison==comp)
  #calculate mean  
  mean_data <- data %>%
    group_by(celltype, tissue, genotype) %>%
    summarise(mean_scaled_E = mean(scaled_E))
  
  # Calculate median scaled_E for each celltype, tissue, and genotype
  median_data <- data %>%
    group_by(celltype, tissue, genotype) %>%
    summarise(median_scaled_E = median(scaled_E))
  # Create the violin plot
  ggplot(data, aes(x = tissue, y = scaled_E, fill = genotype)) +
    #geom_violin(trim = FALSE) +
    geom_boxplot(aes(middle = mean(scaled_E)))+
    #geom_jitter()+
    facet_wrap(~ celltype, scales = "free") +
    labs(x = "Tissue", y = "Scaled E", fill = "Genotype", title = paste0("Cholesterol genes in downregulated_",comp,"_interaction")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/IFN_inter/UP/",comp,"Summary_IFN/"))
  ggsave(ko_dir("Boxplot_Cholesterol_in",comp,"down.pdf"),h=min(20,nrow(data)*0.04+1))
  
}   
##################################################################
##################################################################
basedir<-dirout_load("Ag_ScRNA_12_Pseudobulk_enrichr_per_celltype/FGSEA")
gsea.res<-read_rds(basedir("FGSEA_Interaction_across_KO.RDS"))
head(gsea.res)
gsea.res%>%filter(pathway %in% grep("Cholesterol Homeostasis",gsea.res$pathway,value = T))%>%
  pull(leadingEdge)%>% unlist()%>%unique()

gsea.res%>%filter(pathway %in% grep("Cholesterol Homeostasis",gsea.res$pathway,value = T))%>%
  pull(leadingEdge)%>% unlist()%>%unique()
##################################################################
limmaRes%>%filter(coef=="ex.vivo:Brd9")%>%
  filter(ensg %in% c("Actb", "Gapdh", "Rpl13a", "Hprt", "Tbp"))
list_of_genes<-c("Actb", "Gapdh", "Rpl13a", "Hprt", "Tbp")
dat.list <- list()
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
          filter(genotype %in% c("Brd9", "NTC")) %>%
          mutate(scaled_E = scale(E)) %>%
          mutate(gene = goi)%>%
          mutate(celltype=ct)%>%
          mutate(comparison="Brd9")
        
        # Store the gene data in the list
        dat.list[[paste0(ct, "_", goi,"Brd9")]] <- gene_data
      }
    }
  }}
control_genes<-bind_rows(dat.list,.id = "celltype_gene_KO")
create_gene_plots <- function(data, gene,KO) {
  ggplot(data[data$gene == gene,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype~tissue, scales = "free") +
    labs(title = gene)+
    xlab(paste0(KO,"KO"))
  
  ko_dir<-dirout(paste0("/Ag_ScRNA_16_per_celltype_gene_plots/Brd9_Control_genes",KO))
  ggsave(ko_dir(paste0(gene,".pdf")))
}
head(control_genes)
#For each KO where there were significantly up genes, make the plot
for (comp in unique(control_genes$comparison)){
  data<-control_genes%>%
    filter(comparison==comp)
  gene_plots <- lapply(unique(data$gene), function(gene) {
    create_gene_plots(data, gene,comp)
    
  })
}
# Generate plots for each gene

IFN_genes_logFC<-limmaRes%>%filter(toupper(ensg)%in% IFN_genes)%>%
  filter(abs(logFC)>1)

ggplot(IFN_genes_logFC,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                        color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5))+
  ggtitle("IFN genes-Interaction(ex.vivo-in.vivo)")+
  ylab("IFN response genes")
ggsave(basedir("IFN_genes_in_ko_interaction.pdf"),w=25,h=10)

#############################################################################################

Cholesterol_Hom_logFC<-limmaRes%>%filter(toupper(ensg)%in% Cholesterol_Hom_genes )%>%
  filter(abs(logFC)>1)

ggplot(Cholesterol_Hom_logFC,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                           color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5))+
  ggtitle("Cholesterol_Hom-Interaction(ex.vivo-in.vivo)")+
  ylab("Cholesterol_Hom")
ggsave(basedir("Cholesterol_Hom_in_ko_interaction.pdf"),w=25,h=10)

############################################################################################

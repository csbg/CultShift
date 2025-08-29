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
#####################################################################
inDir<-dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype/")
inDir1<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC/")
inDir2<-dirout_load("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment/ENRICHR")
base<-"Ag_ScRNA_14_per_celltype_gene_plots/"
basedir<-dirout("Ag_ScRNA_14_per_celltype_gene_plots")
########################################################################
#load limma results from NTC
dataVoom_NTC<-read_rds(inDir1("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaRes_NTC<-read_rds(inDir1("limma_perCTex.vivovsin.vivo.rds"))
#load enrichment_results from NTC
enr.res.all_up<-read_rds(inDir2("enr.res.all_NTC_up.rds"))
enr.res.all_down<-read_rds(inDir2("enr.res.all_NTC_down.rds"))

#######################################################################
#load_limma across KO
limmaRes<-read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
dataVoom_Eo.Ba<-read_rds(inDir("Eo.Ba_dataVoom.rds"))
dataVoom_Mono<-read_rds(inDir("Mono_dataVoom.rds"))
dataVoom_MkP<-read_rds(inDir("MkP_dataVoom.rds"))
dataVoom_GMP<-read_rds(inDir("GMP_dataVoom.rds"))
dataVoom_HSC<-read_rds(inDir("HSC_dataVoom.rds"))
dataVoom_MEP<-read_rds(inDir("MEP_dataVoom.rds"))
###########################################################################
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

#Enrichr from NTC--------------------------------------------
############################################
#Cholesterol Homeostasis----------------
############################################
chol<-dirout(paste0(base,"NTC_Cholesterol"))
terms_of_interest <- c("Cholesterol Homeostasis")
#filtered
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, terms_of_interest)
Cholesterol_filtered<-limmaRes_NTC%>%filter(toupper(genes) %in% Cholesterol_genes)%>%
  filter(adj.P.Val<0.05)
ggplot(Cholesterol_filtered,aes(x=celltype,y=genes,size=pmin(-log10(adj.P.Val),5),
                                color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  #facet_wrap(vars(celltype), scales = "free_x")+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5)) +
  #theme(strip.text = element_text(size = 15))+
  #scale_size(range = c(0.5,7))+
  ggtitle(paste0(terms_of_interest,"_NTC_logFC"))
ggsave(chol(paste0("Cholesterol_filtered.pdf")), 
       w=8,h=10)
#unfiltered
Cholesterol_unfiltered<-limmaRes_NTC%>%filter(toupper(genes) %in% Cholesterol_genes)
ggplot(Cholesterol_unfiltered,aes(x=celltype,y=genes,size=pmin(-log10(adj.P.Val),5),
                                  color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  #facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5))+
  ggtitle(paste0(terms_of_interest,"_NTC_logFC"))
ggsave(chol(paste0("Cholesterol_unfiltered.pdf")), 
       w=8,h=10)
############################################
#Cholesterol Homeostasis-INTERACTION-----------------------
Cholesterol<-dirout(paste0(base,"Cholesterol"))
############################################
terms_of_interest <- c("Cholesterol Homeostasis")
#filtered
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, terms_of_interest)
Cholesterol_filtered<-limmaRes%>%filter(toupper(ensg) %in% Cholesterol_genes)%>%
  filter(adj.P.Val<0.05)
ggplot(Cholesterol_filtered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(Cholesterol(paste0("Cholesterol_filtered.pdf")), 
       w=10,h=length(Cholesterol_genes) * 0.2 + 3,limitsize = FALSE)
#unfiltered
Cholesterol_unfiltered<-limmaRes%>%filter(toupper(ensg) %in% Cholesterol_genes)
ggplot(Cholesterol_unfiltered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                  color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5)) 
ggsave(Cholesterol(paste0("Cholesterol_unfiltered.pdf")), 
       w=17,h=length(Cholesterol_genes) * 0.25 + 5,limitsize = FALSE)

library(dplyr)
###########################################################################

#########
#Plot gene expression
#########
NTC_Expr<-dirout(paste0(base,"NTC_Expr"))
Inter_Expr<-dirout(paste0(base,"Inter_Expr"))
#
# Initialize an empty list to store the results
dat.list <- list()
############################
#NTC
############################
# Iterate over each unique gene
for(gg in unique(Cholesterol_unfiltered$ensg)) {
  # Subset the metadata and E values for the current gene
  gene_data <- NTC_meta %>%
    mutate(E = scale(dataVoom_NTC$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
  
  # Group the data by tissue and celltype, and calculate the average E for each group
  avg_gene_data <- gene_data %>%
    group_by(tissue, celltype) %>%
    summarise(avg_E = mean(E))
  
  # Store the average gene data in the list
  dat.list[[gg]] <- avg_gene_data
}
dat.list<-bind_rows(dat.list,.id="gene")
head(dat.list)
ggplot(dat.list,aes(x=tissue, y=gene, fill=avg_E)) + 
  geom_tile() +
  facet_grid(~ celltype, space ="free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))
ggsave(NTC_Expr(paste0(terms_of_interest,"_NTC_averagedsample.png")),w=10,h=8)
####################################################
meta<-read_rds(inDir1("meta.rds"))
NTC_meta<-read_rds(inDir1("NTC_meta.rds"))
dat.list <- list()
for(gg in unique(Cholesterol_unfiltered$ensg)){
  dat.list[[gg]] <- NTC_meta %>%
    mutate(E=scale(dataVoom_NTC$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}

p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(celltype))
head(p.vals)    
(ggplot(p.vals,aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(~ cell, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                 hjust = 1))
ggsave(NTC_Expr(paste0(terms_of_interest,"_NTC_persample.png")),w=18,h=15)  
###############################################################################
#Cholesterol homeostasis genes in INTERACTION----------------------------------
###############################################################################
# # Initialize an empty list to store the results
# dat.list <- list()
# head(meta)
# # Iterate over each unique gene
# for (gg in unique(Cholesterol_unfiltered$ensg)) {
#   # Iterate over each unique cell type
#   for (ct in unique(meta$celltype)) {
#     # Subset the metadata and E values for the current gene and cell type
#     gene_data <- meta %>%
#       filter(celltype == ct) %>%
#       mutate(E = scale(get(paste0("dataVoom_", ct))$E[gg,])) %>%
#       rownames_to_column("sample1") %>%
#       remove_rownames()
#     
#     # Group the data by tissue and coefficient, and calculate the average E for each group
#     avg_gene_data <- gene_data %>%
#       group_by(tissue, genotype) %>%
#       summarise(avg_E = mean(E)) %>%
#       mutate(gene = gg) %>%
#       bind_rows()
#     
#     
#     # Store the average gene data in the list
#     dat.list[[paste0(ct)]] <- avg_gene_data
#   }
# }
# 
#   
# dat.list<-bind_rows(dat.list,.id="celltype")
# head(dat.list)
# ggplot(dat.list,aes(x=genotype, y=gene, fill=avg_E)) + 
#   geom_tile() +
#   facet_grid(celltype~ tissue, space ="free", scales = "free") +
#   scale_fill_gradient2(low="blue", high="red")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
#                                    hjust = 1))
# ggsave(NTC_Expr(paste0(terms_of_interest,"_NTC_averagedsample.png")),w=10,h=8)
selected_genes<-limmaRes%>%filter(ensg%in%c("Sc5d","Idi1","Cd63","Cd9","Fads2","Plscr1"))
ggplot(selected_genes,aes(x=gsub("ex.vivo:","",coef),y=ensg,size=pmin(-log10(adj.P.Val),5),
  color=pmax(logFC,-4)))+
  geom_point()+scale_color_gradient2(high="red", low="blue",limits = c(-6,6))+
  ggtitle("Selected_genes_Choles.KO-Interaction")+
  facet_wrap(vars(celltype), scales = "free_x")+
  theme(axis.text = element_text(size = 12)) +
  theme(strip.text = element_text(size = 15))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                  hjust = 1))
  
ggsave(Cholesterol("selected_genes.png"),w=15) 
####################
SC5D<-limmaRes%>%filter(ensg=="Sc5d" & coef=="ex.vivo:Brd9")
ggplot(SC5D,aes(x=celltype,y=ensg,size=pmin(-log10(adj.P.Val),5),
                color=logFC))+
  geom_point()+scale_color_gradient2(high="red", low="blue")
ggtitle("Scd5_in Brd9KO-Interaction")
#facet_wrap(vars(celltype), scales = "free_x")+
#theme(axis.text = element_text(size = 15)) +
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
#                                hjust = 1))+
ggsave() 
  #scale_size_continuous(range=c(0,5), limits = c(0,5)) 

###############################################
# IFN response-------------
###############################################
terms_of_interest <- c("Interferon Alpha Response", "Interferon Gamma Response")

# Call the function with your dataframe and terms
IFN_genes <- extract_unique_genes(enr.res.all_down, terms_of_interest)

#Interferon_genes_across_interaction_terms in KO
#filtered
IFN_filtered_p_value<-limmaRes%>%
  filter(toupper(ensg) %in% IFN_genes)%>%filter(adj.P.Val<0.05)

ggplot(IFN_filtered_p_value,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
             color=logFC))+
geom_point()+scale_color_gradient2(high="red", low="blue")+
facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(basedir(paste0("Interferon_genes_filtered.pdf")), w=20,h=length(IFN_genes) * 0.05 + 3, limitsize = FALSE)
#unfiltered 
IFN_unfiltered_p_value<-limmaRes%>%
  filter(toupper(ensg) %in% IFN_genes)

ggplot(IFN_unfiltered_p_value,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
             color=logFC))+
  
geom_point()+scale_color_gradient2(high="red", low="blue")+
    
facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
    hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
    scale_size(range = c(0.2,7))
    
ggsave(basedir(paste0("Interferon_genes_unfiltered.pdf")), w=20,h=length(IFN_genes) * 0.3 + 3, limitsize = FALSE)
############################################
#Cholesterol Homeostasis----------------
############################################
terms_of_interest <- c("Cholesterol Homeostasis")
#filtered
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, terms_of_interest)
Cholesterol_filtered<-limmaRes%>%filter(toupper(ensg) %in% Cholesterol_genes)%>%
  filter(adj.P.Val<0.05)
ggplot(Cholesterol_filtered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                  color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(basedir(paste0("Cholesterol_filtered.pdf")), 
       w=20,h=length(Cholesterol) * 0.8 + 3,limitsize = FALSE)
#unfiltered
Cholesterol_unfiltered<-limmaRes%>%filter(toupper(ensg) %in% Cholesterol_genes)
ggplot(Cholesterol_unfiltered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                  color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5)) 
ggsave(basedir(paste0("Cholesterol_unfiltered.pdf")), 
       w=20,h=length(Cholesterol_genes) * 0.8 + 3,limitsize = FALSE)
#######################################################################
#Fatty Acid Metabolism-------------------------------
terms_of_interest <- c("Fatty Acid Metabolism")
#filtered
Fatty_Acid_Metabolism_genes <- extract_unique_genes(enr.res.all_down, terms_of_interest)
Fatty_Acid_Metabolism_filtered<-limmaRes%>%filter(toupper(ensg) %in% Fatty_Acid_Metabolism_genes)%>%
  filter(adj.P.Val<0.05)
ggplot(Fatty_Acid_Metabolism_filtered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                          color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5)) 
ggsave(basedir(paste0("Fatty_Acid_Metabolism_genes_filtered.pdf")), 
       w=20,h=length(Cholesterol) * 0.8 + 3,limitsize = FALSE)
#unfiltered
Fatty_Acid_Metabolism_unfiltered<-limmaRes%>%filter(toupper(ensg) %in% Fatty_Acid_Metabolism_genes)
ggplot(Fatty_Acid_Metabolism_unfiltered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                            color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(basedir(paste0("Fatty_Acid_Metabolism_unfiltered.pdf")), 
       w=20,h=length(Cholesterol) * 0.8 + 3,limitsize = FALSE)
list_up<-list()
#mTORC1 Signaling--------------
terms_of_interest <- c("mTORC1 Signaling")
#filtered
mTORC1_Signaling_genes <- extract_unique_genes(enr.res.all_down, terms_of_interest)
mTORC1_Signaling_filtered<-limmaRes%>%filter(toupper(ensg) %in% mTORC1_Signaling_genes)%>%
  filter(adj.P.Val<0.05)
ggplot(mTORC1_Signaling_filtered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                     color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(basedir(paste0("mTORC1_Signaling_genes_filtered.pdf")), 
       w=20,h=length(mTORC1_Signaling_genes) * 0.8 + 3,limitsize = FALSE)
#unfiltered
mTORC1_Signaling_unfiltered<-limmaRes%>%filter(toupper(ensg) %in% mTORC1_Signaling_genes)
ggplot(mTORC1_Signaling_unfiltered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                       color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(basedir(paste0("mTORC1_Signaling_unfiltered.pdf")), 
       w=20,h=length(Cholesterol) * 0.8 + 3,limitsize = FALSE)
#Idl1---------------
Idi1<-limmaRes%>%filter(toupper(ensg) == "IDI1")
ggplot(Idi1,aes(x=coef,y=celltype,size=pmin(-log10(adj.P.Val),5),
                                  color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  #facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  ggtitle("IDI1")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(basedir(paste0("IDI1.pdf")))

########################################################
#
########################################################
head(dataVoom_Mono)
dat.list <- list()
for(gg in unique(Cholesterol_filtered$ensg)){
  dat.list[[gg]] <- meta %>%filter(celltype=="Mono")%>%filter(genotype=="Brd9")%>%
    mutate(E=scale(dataVoom_Mono$E[gg,grep("Brd9",colnames(dataVoom_Mono$E))])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}

p.vals <- bind_rows(dat.list, .id="genes") %>%
  mutate(cell = as.character(celltype))
head(p.vals)    
(ggplot(p.vals,aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    #facet_grid(~ , space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))
ggsave(Inter_Expr(paste0(terms_of_interest,"_NTC_persample_Mono.png")),w=10,h=12)  
###############################################################
for(gg in unique(Cholesterol_filtered$ensg)){
  dat.list[[gg]] <- meta %>%filter(celltype=="Mono")%>%filter(genotype=="Smarcb1")%>%
    mutate(E=scale(dataVoom_Mono$E[gg,grep("Smarcb1",colnames(dataVoom_Mono$E))])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}

p.vals <- bind_rows(dat.list, .id="genes") %>%
  mutate(cell = as.character(celltype))
head(p.vals)    
(ggplot(p.vals,aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    #facet_grid(~ , space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  ggtitle("Smarcb1")
ggsave(Inter_Expr(paste0(terms_of_interest,"_NTC_persample_Mono_Smarcb1.png")),w=10,h=12) 
#
dat.list <- list()
for(gg in unique(Cholesterol_filtered$ensg)){
  dat.list[[gg]] <- meta %>%filter(celltype=="Eo.Ba")%>%filter(genotype=="Brd9")%>%
    mutate(E=scale(dataVoom_Eo.Ba$E[gg,grep("Brd9",colnames(dataVoom_Eo.Ba$E))])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}

p.vals <- bind_rows(dat.list, .id="genes") %>%
  mutate(cell = as.character(celltype))
head(p.vals)    
(ggplot(p.vals,aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    #facet_grid(~ , space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  ggtitle("Brd9")
ggsave(Inter_Expr(paste0(terms_of_interest,"_Inter_persample_Eo.Ba_Brd9.png")),w=10,h=12) 
##################################################
#Brd9-IFN
##################################################
#Brd9---------------------





Brd9_50 <- limmaRes %>%
  filter(coef == "ex.vivo:Brd9") %>%
  filter(adj.P.Val < 0.05) %>%
  filter(logFC > 1) %>%
  #filter(toupper(ensg) %in% IFN_genes)%>%
  group_by(celltype) %>%
  top_n(50, wt = logFC)%>%pull(ensg)%>%unique()
Brd9<-enrichr(Brd9_50,"MSigDB_Hallmark_2020")%>%bind_rows()

Brd9_IFN_genes<-limmaRes %>%
  filter(coef == "ex.vivo:Brd9") %>%
  #filter(adj.P.Val < 0.05) %>%
  filter(logFC > 1) %>%
  filter(toupper(ensg) %in% IFN_genes)


ggplot(Brd9_IFN_genes,aes(x=celltype,y=ensg,size=pmin(-log10(adj.P.Val),5),
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
ggsave(IFN_inter("Brd9KO_interaction_upregulated_IFN_genes.pdf"),h=nrow(Brd9_IFN_genes)*0.04+1)

###################
#Expression
##############
meta<-read_rds(inDir1("meta.rds"))
NTC_meta<-read_rds(inDir1("NTC_meta.rds"))

KO <- "Brd9"
dat.list <- list()

for (ct in unique(meta$celltype)) {
  # Get the dataVoom object corresponding to the current cell type
  dataVoom_ct <- get(paste0("dataVoom_", ct))
  
  # Check if goi exists in the row names of dataVoom_ct$E
  if (any(rownames(dataVoom_ct$E) %in% unique(Brd9_IFN_genes$ensg))) {
    for (goi in unique(Brd9_IFN$ensg)) {
      # Proceed only if goi exists in the row names of dataVoom_ct$E
      if (goi %in% rownames(dataVoom_ct$E)) {
        # Subset the metadata and E values for the current gene and cell type
        gene_data <- meta[meta$celltype == ct,] %>%
          mutate(E = dataVoom_ct$E[goi,]) %>%
          rownames_to_column("sample1") %>%
          filter(genotype %in% c(KO, "NTC")) %>%
          mutate(scaled_E = scale(E)) %>%
          mutate(gene = goi)%>%
          mutate(celltype=ct)
        
        # Store the gene data in the list
        dat.list[[paste0(ct, "_", goi)]] <- gene_data
      }
    }
  }
}

data_Brd9<-bind_rows(dat.list,.id = "celltype_gene")

##########################################################
#create dir
IFN_inter<-dirout(paste0(base,"IFN_inter"))
# Function to create plots for each gene
create_gene_plots <- function(data, gene,KO) {
  ggplot(data[data$gene == gene,], aes(x = genotype, y = scaled_E)) + 
    geom_boxplot() +
    geom_jitter() +
    facet_grid(celltype~tissue, scales = "free") +
    labs(title = gene)+
    xlab(paste0(KO,"KO"))
  ggsave(IFN_inter(paste0(gene,".pdf")))
}

# Generate plots for each gene
gene_plots <- lapply(unique(data_Brd9$gene), function(gene) {
  create_gene_plots(data_Brd9, gene,"Brd9")
  
})

#####################################################################
#Barplopt
# Group data by celltype, tissue, and genotype, calculate the mean scaled_E
mean_data <- data_Brd9 %>%
  group_by(celltype, tissue, genotype) %>%
  summarise(mean_scaled_E = mean(scaled_E))

# Plot mean scaled_E for each celltype in each tissue
ggplot(mean_data, aes(x = celltype, y = mean_scaled_E, fill = genotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ tissue, scales = "free") +
  labs(x = "Cell Type", y = "Mean Scaled E", fill = "Genotype") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(IFN_inter(paste0("mean_of_scaled_E_across_IFN_gene",".pdf")))

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
#####################################################################
Ifit3_data <- meta[meta$celltype=="Eo.Ba",]%>%
  mutate(E = dataVoom_Eo.Ba$E["Ifit3",])%>%
  rownames_to_column("sample1")%>%
  filter(genotype %in% c("Brd9","NTC"))%>%
  mutate(scaled_E=scale(E))%>%
  mutate(gene="Ifit3")


##############################################################
#
ggplot(Ifit3_data,aes(x=genotype, y=scaled_E)) + 
  geom_boxplot()+
  geom_jitter() +
  facet_grid(~tissue, space ="free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))
#################################################################
Rplp0_data <- meta[meta$celltype=="Eo.Ba",]%>%
  mutate(E = dataVoom_Eo.Ba$E["Rplp0",])%>%
  rownames_to_column("sample1")%>%
  filter(genotype %in% c("Brd9","NTC"))%>%
  mutate(scaled_E=scale(E))%>%
  mutate(gene="Rplp0")
ggplot(Rplp0_data,aes(x=genotype, y=scaled_E)) + 
  geom_boxplot()+
  geom_jitter() +
  facet_grid(~tissue, space ="free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))

#mutate(gene=gg)%>%
  #remove_rownames()
  # 
  #   
   dat.list<-bind_rows(dat.list,.id="celltype")
   head(dat.list)
   ggplot(dat.list,aes(x=genotype, y=gene, fill=avg_E)) + 
     geom_tile() +
     facet_grid(celltype~ tissue, space ="free", scales = "free") +
     scale_fill_gradient2(low="blue", high="red")+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                      hjust = 1))
#########################################################
terms_of_interest <- c("Interferon Alpha Response", "Interferon Gamma Response")

# Call the function with your dataframe and terms
IFN_genes <- extract_unique_genes(enr.res.all_down, terms_of_interest)

#Interferon_genes_across_interaction_terms in KO
#filtered
IFN_filtered_p_value<-limmaRes%>%
  filter(toupper(ensg) %in% IFN_genes)%>%filter(adj.P.Val<0.05)

ggplot(IFN_filtered_p_value,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                color=logFC))+
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(basedir(paste0("Interferon_genes_filtered.pdf")), w=20,h=length(IFN_genes) * 0.05 + 3, limitsize = FALSE)
#unfiltered 
IFN_unfiltered_p_value<-limmaRes%>%
  filter(toupper(ensg) %in% IFN_genes)

ggplot(IFN_unfiltered_p_value,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                  color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))

ggsave(basedir(paste0("Interferon_genes_unfiltered.pdf")), w=20,h=length(IFN_genes) * 0.3 + 3, limitsize = FALSE)
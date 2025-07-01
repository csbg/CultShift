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
inDir2<-dirout_load("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment/ENRICHR")
base<-"Ag_ScRNA_15_per_celltype_gene_plots_specific/"
basedir<-dirout("Ag_ScRNA_15_per_celltype_gene_plots_specific")
#####################################################################
#loading enrichr terms
source("src/00_init.R")
out <- dirout("EXT_02_EnrichR_Genesets/")
# # # Download gene sets ------------------------------------------------------
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
#save(enr.terms, file=out("Genesets_Mouse.RData"))
#function
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
#

########################################################################
#load data
limmaRes<-read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))
enr.res.all_up<-read_rds(inDir2("enr.res.all_NTC_up.rds"))
enr.res.all_down<-read_rds(inDir2("enr.res.all_NTC_down.rds"))
#########################################################################
ct<-"Mono"
dir<-dirout(paste0(base,ct))

Mono<-limmaRes%>%
  filter(adj.P.Val<adj_p_cutoff &celltype ==ct )
  ggplot(Mono,aes(x =gsub("ex.vivo:","",coef),
             y = ensg, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),3))) +
  scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5)) +
  geom_point() +
  scale_size_continuous(range=c(0,8), limits = c(0,3)) +
  theme_bw(12) +
  #xRot() +
  #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
  labs(x = "interaction-KO")+
  
  theme(axis.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
ggsave(dir("Mono_KO-interaction.png"))
#############################################################################

terms_of_interest <- c("Cholesterol Homeostasis")
#filtered
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, terms_of_interest)
Cholesterol_filtered<-limmaRes%>%filter(toupper(ensg) %in% Cholesterol_genes)%>%
  filter(adj.P.Val<0.05 &celltype==ct)
ggplot(Cholesterol_filtered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                color=logFC))+
  
  geom_point()+scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(dir(paste0("NTC_overl.Cholesterol_filtered.png")), 
       w=16,h=length(Cholesterol_genes) * 0.3 + 3,limitsize = FALSE)
#unfiltered
Cholesterol_unfiltered<-limmaRes%>%filter(toupper(ensg) %in% Cholesterol_genes & celltype==ct)
ggplot(Cholesterol_unfiltered,aes(x=coef,y=ensg,size=-log10(adj.P.Val),
                                  color=pmax(logFC,-5)))+
  
  geom_point()+
  scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5)) +
  theme(strip.text = element_text(size = 15))+
  ggtitle(paste0(terms_of_interest))
ggsave(dir(paste0("NTC_overl.Cholesterol_unfiltered.png")), 
       w=15,h=length(Cholesterol_genes) * 0.2 + 3,limitsize = FALSE)
##############################################################################
#Any cholesterol homeostasis
cholesterol_homestasis<-enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`
cholesterol_homeostasis<-unlist(strsplit(cholesterol_homestasis, split=","))
cholesterol_homeostasis <- trimws(cholesterol_homeostasis)
cholesterol_homeostasis <- cholesterol_homeostasis[cholesterol_homeostasis != ""] 


Cholesterol_filtered<-limmaRes%>%filter(ensg %in% cholesterol_homeostasis)%>%
  filter(adj.P.Val<0.05 &celltype==ct)
ggplot(Cholesterol_filtered,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                color=logFC))+
  
  geom_point()+scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  scale_size_continuous(range=c(1,5), limits = c(1,5)) +
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(dir(paste0("All_Cholesterol_homeostasis_filtered.png")), 
       w=16,h=length(cholesterol_homeostasis) * 0.3 + 3,limitsize = FALSE)
#unfiltered
Cholesterol_unfiltered<-limmaRes%>%filter(ensg %in% cholesterol_homeostasis & celltype==ct)
ggplot(Cholesterol_unfiltered,aes(x=coef,y=ensg,size=-log10(adj.P.Val),
                                  color=pmax(logFC,-5)))+
  
  geom_point()+
  scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5)) +
  theme(strip.text = element_text(size = 15))+
  ggtitle(paste0(terms_of_interest))
ggsave(dir(paste0("All_Cholesterol_homeostasis_unfiltered.png")), 
       w=10,h=length(cholesterol_homeostasis) * 0.1 + 2,limitsize = FALSE)
################################################################################
#############
#Eo.Ba-------
#############
ct<-"Eo.Ba"
dir<-dirout(paste0(base,ct))

Eo.Ba<-limmaRes%>%
  filter(adj.P.Val<adj_p_cutoff &celltype ==ct )
ggplot(Eo.Ba,aes(x =gsub("ex.vivo:","",coef),
                 y = ensg, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),3))) +
  scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5)) +
  geom_point() +
  scale_size_continuous(range=c(0,8), limits = c(0,3)) +
  theme_bw(12) +
  #xRot() +
  #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
  labs(x = "interaction-KO")+
  
  theme(axis.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
ggsave(dir("Eo.Ba_KO-interaction.png"))
#############################################################################
terms_of_interest <- c("Cholesterol Homeostasis")
#filtered
Cholesterol_genes <- extract_unique_genes(enr.res.all_up, terms_of_interest)
Cholesterol_filtered<-limmaRes%>%filter(toupper(ensg) %in% Cholesterol_genes)%>%
  filter(adj.P.Val<0.05 &celltype==ct)
ggplot(Cholesterol_filtered,aes(x =gsub("ex.vivo:","",coef),y=ensg,size=pmin(-log10(adj.P.Val),5),
                                color=logFC))+
  
  geom_point()+scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(dir(paste0("NTC_overl.Cholesterol_filtered.png")))
#unfiltered
Cholesterol_unfiltered<-limmaRes%>%filter(toupper(ensg) %in% Cholesterol_genes & celltype==ct)
ggplot(Cholesterol_unfiltered,aes(x =gsub("ex.vivo:","",coef),y=ensg,size=-log10(adj.P.Val),
                                  color=pmax(logFC,-5)))+
  
  geom_point()+
  scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,8), limits = c(0,4)) +
  theme(strip.text = element_text(size = 15))+
  #scale_size(range = c(0.5,7))+
  ggtitle(paste0(terms_of_interest))
ggsave(dir(paste0("NTC_overl.Cholesterol_unfiltered.png")), 
       w=12,h=length(Cholesterol_genes) * 0.2 + 3,limitsize = FALSE)
##############################################################################
#Any cholesterol homeostasis
cholesterol_homestasis<-enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`
cholesterol_homeostasis<-unlist(strsplit(cholesterol_homestasis, split=","))
cholesterol_homeostasis <- trimws(cholesterol_homeostasis)
cholesterol_homeostasis <- cholesterol_homeostasis[cholesterol_homeostasis != ""] 
cholesterol_homeostasis

Cholesterol_filtered<-limmaRes%>%filter(ensg %in% cholesterol_homeostasis)%>%
  filter(adj.P.Val<0.05 &celltype==ct)
ggplot(Cholesterol_filtered,aes(x =gsub("ex.vivo:","",coef),y=ensg,size=pmin(-log10(adj.P.Val),5),
                                color=logFC))+
  
  geom_point()+scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  scale_size_continuous(range=c(1,5), limits = c(1,5)) +
  
  facet_wrap(vars(celltype), scales = "free_x")+theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 0.5))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size(range = c(0.2,7))
ggsave(dir(paste0("All_Cholesterol_homeostasis_filtered.png")), 
       w=16,h=length(cholesterol_homeostasis) * 0.3 + 3,limitsize = FALSE)
#unfiltered
Cholesterol_unfiltered<-limmaRes%>%filter(ensg %in% cholesterol_homeostasis & celltype==ct)
ggplot(Cholesterol_unfiltered,aes(x =gsub("ex.vivo:","",coef),y=ensg,size=-log10(adj.P.Val),
                                  color=pmax(logFC,-5)))+
  
  geom_point()+
  scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,8), limits = c(0,5)) +
  theme(strip.text = element_text(size = 15))+
  ggtitle(paste0(terms_of_interest))
ggsave(dir(paste0("All_Cholesterol_homeostasis_unfiltered.png")), 
       w=16,h=length(cholesterol_homeostasis) * 0.2 + 2,limitsize = FALSE)
################################################################################
exosome_genes<-limmaRes%>%filter(ensg %in% c("Pdcd6ip", "Tsg101", "Cd63", "Rab27a", "Rab27b", "Smpd3", "Plp2") & celltype==ct)
ggplot(exosome_genes,aes(x=coef,y=ensg,size=-log10(adj.P.Val),
                                  color=pmax(logFC,-5)))+
  
  geom_point()+
  scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5)) +
  theme(strip.text = element_text(size = 15))+
  ggtitle(paste0(terms_of_interest))
ggsave(dir(paste0("Exosome_related_genes.png")))
################################################################################

##################################
selected_chol<-c("Gstm7","Gpx8","Gnai1","Gldc")

#unfiltered
selected_chol_unfiltered<-limmaRes%>%filter(ensg %in% selected_chol )
ggplot(selected_chol_unfiltered,aes(x =gsub("ex.vivo:","",coef),
                                    y=ensg,size=-log10(adj.P.Val),
                                   color=pmax(logFC,-5)))+
  
  geom_point()+
  scale_color_gradient2(low="blue", mid="white", high="red", limits = c(-5, 5))+
  
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+
  geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5)) +
  theme(strip.text = element_text(size = 15))
 # ggtitle(paste0(terms_of_interest))
ggsave(basedir(paste0("Selected_Cholesterol_unfiltered_across_celltypes.png")), 
       w=18,h=length(selected_chol_unfiltered) * 0.1 + 3,limitsize = FALSE)

##################################
inflammatory_response_genes<-enr.terms$MSigDB_Hallmark_2020$`Inflammatory Response`
interferon_alpha_genes<-enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`
interferon_gamma_genes<-enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`

oxphos<-enr.terms$MSigDB_Hallmark_2020$`Oxidative Phosphorylation`
glycolysis<-enr.terms$MSigDB_Hallmark_2020$Glycolysis
random<-enr.terms$WikiPathways_2019_Mouse$`Dysregulated miRNA Targeting in Insulin/PI3K-AKT Signaling WP3855`
mtorc<-enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`
cholesterol<-enr.terms$WikiPathways_2019_Mouse$`Cholesterol Biosynthesis WP103`
adipogenesis<-enr.terms$WikiPathways_2019_Mouse$`Adipogenesis genes WP447`
fatty_acid<-enr.terms$MSigDB_Hallmark_2020$`Fatty Acid Metabolism`
Tnf<-enr.terms$MSigDB_Hallmark_2020$`TNF-alpha Signaling via NF-kB`

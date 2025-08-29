####################################################################
#******#
####################################################################
#Create pseudobulk object and corresponding metada from corresponding single cell 
#data
#
#Aarathy
###############
source("src/00_init.R")
library(edgeR)
library(limma)
library(tidyverse)
library(enrichR)
library(purrr)
library(gridExtra)
#####################################################################
inDir<-dirout_load("Ag_ScRNA_11_Pseudobulk_limma_all_ko_ex.vivo_vs_in.vivo_per_celltype/")
base<-"Ag_ScRNA_13_Pseudobulk_limma_all_ko_perCT_plots/"
basedir<-dirout("Ag_ScRNA_13_Pseudobulk_limma_all_ko_perCT_data_plots/")
###########################################################
limmaRes<-read_rds(inDir("limma_ex.vivo_vs_in.vivo_per_CT_interaction.rds"))

data<-limmaRes[limmaRes$coef %in% grep("ex.vivo:",limmaRes$coef,value=T),]
######################
#upregulated
#######################
#top_pos_list<-list()
cutoff<-10
for (ct in unique(data$celltype)) {
  for (coefx in unique(data$coef)) {
    top_pos <- data %>%
      filter(adj.P.Val < 0.05 & coef == coefx & celltype == ct & logFC > 0) %>%
      arrange(desc(logFC)) %>%
      head(cutoff) %>%
      pull(ensg) %>%
      unique()
    top_pos_list[[paste0(coefx,"_",ct)]]<-top_pos
  }
}
top_pos_list_table <- data.frame(
  list_name = names(top_pos_list),
  genes = sapply(top_pos_list, function(x) toString(x))
)

top_pos_list_table_filtered <- top_pos_list_table %>%
  filter(grepl("\\S", genes))

write.csv(top_pos_list_table_filtered,
          basedir("top_positive_interaction_genes_coef_table.csv"))
goi.all<-unique(unlist(top_pos_list))
genes_to_exclude <- c("Hbb-bs", "Hba-a2", "Hba-a1")  # Genes to exclude from goi.all

goi.all <- goi.all[!goi.all %in% genes_to_exclude]

length(goi.all)
#unfiltered
(top_10_pos <- data %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    facet_wrap(vars(celltype))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 8)) +
  ggtitle(paste0("Top_",cutoff,"_up_INTERACTION"))
ggsave(basedir(paste0("top_pos_logFC_",cutoff,"_genes_across_KO_unfiltered.pdf")),plot=top_10_pos,
       w = 20, h = 30, limitsize = FALSE)
#filtered
data_sign<-data%>%
  filter(adj.P.Val<0.05)%>%
  filter(ensg%in%goi.all)
(top_10_pos <- data_sign %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    facet_grid(cols=vars(celltype))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 10)) +
  ggtitle(paste0("Top_",cutoff,"_up_INTERACTION"))
ggsave(basedir(paste0("top_pos_logFC_",cutoff,"_genes_across_KO_filtered_out_adj.p.pdf")),plot=top_10_pos,
       w = 20, h = 20, limitsize = FALSE)

#################################################
#enrichr
enrich_top_genes<-enrichr(goi.all,"MSigDB_Hallmark_2020")%>%bind_rows()
#################################################
#interferon
#IFN PROCR;IFITM2;B2M;GBP4
IFN<-c("Procr","Ifitm2","B2m","Gbp4")
#unfiltered
data_ifn<-data%>%
  filter(ensg%in%IFN)
ggplot(data=data_ifn,aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  facet_wrap(vars(celltype))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  theme(axis.text = element_text(size = 10)) +
  coord_flip()+
  ggtitle(paste0("Interferon_up_interaction"))
ggsave(basedir(paste0("top_up_interaction_IFN_genes_unfiltered.pdf")),
      w = 20, h = 20, limitsize = FALSE)
#filtered
data_ifn<-data%>%
  filter(ensg%in%IFN)%>%
  filter(adj.P.Val<0.05)
ggplot(data=data_ifn,aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
  geom_point() +
  scale_color_gradient2(high="red", low="blue") +
  theme_bw()+
  facet_wrap(vars(celltype))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  theme(axis.text = element_text(size = 10)) +
  coord_flip()+
  ggtitle(paste0("Interferon_up_interaction"))
ggsave(basedir(paste0("top_up_interaction_IFN_genes_filtered.pdf")),
       w = 20, h = 20, limitsize = FALSE)

######################################
#downregulated
######################################
top_neg_list<-list()
cutoff<-10

for (ct in unique(data$celltype)) {
  for (coefx in unique(data$coef)) {
    top_neg <- data %>%
      filter(adj.P.Val < 0.05 & coef == coefx & celltype == ct) %>%
      arrange(logFC) %>%
      head(cutoff) %>%
      pull(ensg) %>%
      unique()
    top_neg_list[[paste0(coefx,"_",ct)]]<-top_neg
  }
}
top_neg_list_table <- data.frame(
  list_name = names(top_neg_list),
  genes = sapply(top_neg_list, function(x) toString(x))
)

top_neg_list_table_filtered <- top_neg_list_table %>%
  filter(grepl("\\S", genes))

write.csv(top_neg_list_table_filtered,
          basedir("top_negative_interaction_genes_coef_table.csv"))
goi.all<-NULL
goi.all<-unique(unlist(top_neg_list))
genes_to_exclude <- c("Hbb-bs", "Hba-a2", "Hba-a1")  # Genes to exclude from goi.all

goi.all <- goi.all[!goi.all %in% genes_to_exclude]

length(goi.all)
#unfiltered
(top_10_neg <- data %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    facet_wrap(vars(celltype))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 8)) +
  ggtitle(paste0("Top_",cutoff,"_down_INTERACTION"))
ggsave(basedir(paste0("top_neg_logFC_",cutoff,"_genes_across_KO_unfiltered_edgeR_limma.pdf")),plot=top_10_neg,
       w = 20, h = 30, limitsize = FALSE)
#filtered
data_sign<-data%>%
  filter(adj.P.Val<0.05)%>%
  filter(ensg%in%goi.all)
(top_10_neg <- data_sign %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    facet_grid(cols=vars(celltype))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 10)) +
  ggtitle(paste0("Top_",cutoff,"_down_INTERACTION"))
ggsave(basedir(paste0("top_neg_logFC_",cutoff,"_genes_across_KO_filtered_out_adj.pedgeR_limma.pdf")),plot=top_10_neg,
       w = 35, h = 30, limitsize = FALSE)
enriched_neg<-enrichr(goi.all,"MSigDB_Hallmark_2020")%>%bind_rows()
########################
#Y linked
#unfiltered
(top_10_neg <- data %>%
   filter(ensg %in% y_linked) %>%
   ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
   geom_point() +
   scale_color_gradient2(high="red", low="blue") +
   theme_bw()+
   facet_wrap(vars(celltype))+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 8)) +
  ggtitle(paste0("Top_",cutoff,"_down_INTERACTION"))
ggsave(basedir(paste0("ylinked_",cutoff,"_genes_across_KO_unfiltered.pdf")),plot=top_10_neg,
       w = 10, h = 10, limitsize = FALSE)
#filtered
y_linked<-c("Kdm5d","Ddx3y","Eif2s3y")
data_sign<-data%>%
  filter(adj.P.Val<0.05)%>%
  filter(ensg%in%y_linked)
(top_10_neg <- data_sign %>%
    filter(ensg %in% goi.all) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    facet_grid(cols=vars(celltype))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 10)) +
  ggtitle(paste0("Top_",cutoff,"_down_INTERACTION"))
ggsave(basedir(paste0("ylinked_",cutoff,"_genes_across_KO_filtered_out_adj.p.pdf")),plot=top_10_neg,
       w = 10, h = 10, limitsize = FALSE)
#Cholesterol homeostasis genes
cholesterol<-c("Idl1","Sc5d","Cd9","Decr1","S100a10")
(top_10_neg <- data %>%
    filter(ensg %in% cholesterol) %>%
    ggplot(aes(y=ensg, x=coef, color=pmax(-5,logFC), size=pmin(-log10(adj.P.Val),5))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw()+
    facet_wrap(vars(celltype))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)))+
  theme(axis.text = element_text(size = 8)) +
  ggtitle(paste0("Top_",cutoff,"_down_INTERACTION"))
ggsave(basedir(paste0("cholesterol_",cutoff,"_genes_across_KO_unfiltered.pdf")),plot=top_10_neg,
       w = 10, h = 10, limitsize = FALSE)
####################################
load(PATHS$RESOURCES$Enrichr.mouse)
dbx<-unique(names(enr.terms))[1]


#pathway genes
inflammatory_response_genes<-enr.terms$MSigDB_Hallmark_2020$`Inflammatory Response`
interferon_alpha_genes<-enr.terms$MSigDB_Hallmark_2020$`Interferon Alpha Response`
interferon_gamma_genes<-enr.terms$MSigDB_Hallmark_2020$`Interferon Gamma Response`
cholesterol_homestasis<-enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`
oxphos<-enr.terms$MSigDB_Hallmark_2020$`Oxidative Phosphorylation`
glycolysis<-enr.terms$MSigDB_Hallmark_2020$Glycolysis
random<-enr.terms$WikiPathways_2019_Mouse$`Dysregulated miRNA Targeting in Insulin/PI3K-AKT Signaling WP3855`
mtorc<-enr.terms$MSigDB_Hallmark_2020$`mTORC1 Signaling`
cholesterol<-enr.terms$WikiPathways_2019_Mouse$`Cholesterol Biosynthesis WP103`
adipogenesis<-enr.terms$WikiPathways_2019_Mouse$`Adipogenesis genes WP447`
fatty_acid<-enr.terms$MSigDB_Hallmark_2020$`Fatty Acid Metabolism`
Tnf<-enr.terms$MSigDB_Hallmark_2020$`TNF-alpha Signaling via NF-kB`
#############################################
dat.list <- list()


for(gg in goi.all){
  dat.list[[gg]] <- meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
p.vals <- bind_rows(dat.list, .id="genes")

(p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ celltype, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(out(paste0(selection,"_normalized_expression_across_sample.pdf")),
       w = 30, h = 15)
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(rows = vars(celltype), space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text = element_text(size = 5)) 
ggsave(out(paste0(selection,"_normalized_expression_across_sample_cell_wise.pdf")),
       w = 40, h = 15)
unique(p.vals$sample1)
unique(p.vals$sample)
unique(p.vals$cell)
p.vals <- bind_rows(dat.list, .id="genes") %>%
  mutate(cell = as.character(cell))
dat<-ggplot(limmaRes=p.vals,aes(x=tissue, y=genes, fill=E)) + 
  geom_tile() +
  facet_grid(. ~ celltype, space ="free", scales = "free") +
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dat1<-dat$limmaRes
ggsave(out(paste0(selection,"_normalized_expression_across_sample.pdf")),
       w = 30, h = 15)
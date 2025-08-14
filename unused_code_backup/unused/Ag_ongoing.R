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
library(readr)
################################################################################
inDir<-dirout("Ag_ScRNA_09_pseudobulk_per_celltype_limma_NTC/")
basedir<-dirout("Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment/")
base<-"Ag_ScRNA_10_Pseudobulk_ex_in_NTC_Enrichment/"
FGSEA<-dirout(paste0(base,"FGSEA"))
ENRICHR<-dirout(paste0(base,"ENRICHR"))
################################################################################
#load limma results
dataVoom_NTC<-read_rds(inDir("dataVoom_perCTex.vivovsin.vivo.rds"))
limmaRes_NTC<-read_rds(inDir("limma_perCTex.vivovsin.vivo.rds"))
################################################################################
db = c("KEGG_2019_Mouse" ,
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse", 
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019")
################################################################################
# enrichment up-regulated genesin NTC ex.vivo
enr.res.all_up<-read_rds(ENRICHR("enr.res.all_NTC_up.rds"))
# enrichmentdown-regulated genesin NTC ex.vivo
enr.res.all_down<-read_rds(ENRICHR("enr.res.all_NTC_down.rds"))
################################################################################
# Create an empty list to store the data frames
dbx<-"KEGG_2019_Mouse"
# Loop through unique values of dbx
for(dbx in db) {
  dat <- dirout(paste0(base, "FGSEA/", dbx))
  # Read the TSV file and store it in the list
  fgsea_NTC_list[[dbx]] <- read_tsv(dat(paste0("GSEA_significant_", dbx, ".tsv")))
}

################################################################################
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
#############################################
#genes_from_pathways
############################################
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

################################################################################
#function
################################################################################
extract_unique_genes <- function(dataframe, terms,dbx) {
  # Filter the dataframe for rows where the Term is in the specified list of terms
  filtered_table <- subset(dataframe, Term %in% terms)%>%filter(db==dbx)
  
  # Extract the genes column from the filtered dataframe
  genes <- filtered_table$Genes
  
  # Split the genes column by ";" and flatten the resulting list
  all_genes <- unlist(strsplit(genes, ";"))
  
  # Select only unique genes
  unique_genes <- unique(all_genes)
  
  return(unique_genes)
}
################################################################################
#Enrichr results----------------------------------------------------------------
################################################################################
#plot
################################################################################

dfs<-list(enr.res.all_up,enr.res.all_down)
names(dfs)<-c("up","down")
databases = c("KEGG_2019_Mouse" ,
              "MSigDB_Hallmark_2020",
              "WikiPathways_2019_Mouse", 
              "GO_Biological_Process_2021",
              "TRRUST_Transcription_Factors_2019")
df_filtered <- enr.res.all_down[enr.res.all_down$db == "MSigDB_Hallmark_2020", ]

pmap(list(names(dfs), dfs), function(name, df) {
  walk(databases, ~{
    db <- .x
    df_filtered <- df[df$db == db, ]
    
    # Check if subset is not empty
    if (nrow(df_filtered) > 0) {
      # Remove the prefix from column names
      colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
      
      # Perform filtering and mutation
      plotting_enr <- df_filtered %>%
        filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
        mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
      
      ggplot(plotting_enr, aes(x = celltype, y = Term, size = log2(Odds.Ratio),
                               color = neg.log10.Adjusted.P.value)) +
        geom_point() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
        theme_bw(12) 
      #facet_wrap(vars(db))
      
      ggsave(ENRICHR(paste0(name, "_per.celltype_", db, ".pdf")), 
             w = 20, h = length(unique(plotting_enr$Term)) * 0.2 + 3, limitsize = FALSE)
      
      # Do something with 'result'
    } else {
      warning("Subset is empty for iteration ")
    }
  })
})
##################
#without filter
##################
# pmap(list(names(dfs), dfs), function(name, df) {
#   walk(databases, ~{
#     db <- .x
#     df_filtered <- df[df$db == db, ]
#     df_filtered <- df_filtered %>% mutate(neg.log10.Adjusted.P.value = -log10(Adjusted.P.value))
#     
#     # Check if subset is not empty
#     if (nrow(df_filtered) > 0) {
#       # Remove the prefix from column names
#       colnames(df_filtered) <- gsub(paste0("^", name, "_"), "", colnames(df_filtered))
#       
#       # Perform filtering and mutation
#       plotting_enr <- df_filtered %>%
#         filter(Odds.Ratio > 5 & Adjusted.P.value < 0.01) %>%
#         pull(Term)
#       enr<-df_filtered %>%filter(Term %in% plotting_enr)
#       ggplot(enr, aes(x = celltype, y = Term, size = log2(Odds.Ratio),
#                       color = neg.log10.Adjusted.P.value)) +
#         geom_point() +
#         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
#         theme_bw(12) 
#       #facet_wrap(vars(db))
#       
#       ggsave(ENRICHR(paste0(name, "_per.celltype_", db, "without_filter.pdf")), 
#              w = 20, h = length(unique(enr$Term)) * 0.2 + 3, limitsize = FALSE)
#       
#       # Do something with 'result'
#     } else {
#       warning("Subset is empty for iteration ")
#     }
#   })
# })




################################################################################
#IFN_genes
IFN_genes_overlapping <- extract_unique_genes(enr.res.all_down, c("Interferon Alpha Response",
                                                      "Interferon Gamma Response"),
                                  "MSigDB_Hallmark_2020")
#All IFN genes
IFN_genes_all_MSig<-union(interferon_alpha_genes,interferon_gamma_genes)
#Cholesterol homeostasis
Cholesterol_genes_overlapping <- extract_unique_genes(enr.res.all_up, c("Cholesterol Homeostasis","Cholesterol Biosynthesis WP103")
                                          ,c("MSigDB_Hallmark_2020","WikiPathways_2019_Mouse"))
#All cholesterol genes
cholesterol_homestasis<-enr.terms$MSigDB_Hallmark_2020$`Cholesterol Homeostasis`
cholesterol<-enr.terms$WikiPathways_2019_Mouse$`Cholesterol Biosynthesis WP103`

Cholesterol_genes_all<-union(cholesterol_homestasis,cholesterol)

################################################################################
#Plots
selected_genes_plot<-function(limma_data,input_genes,data_name){
  if (all(grepl("^[A-Z]+$", input_genes))) {
    selected <- limma_data %>%
      filter(toupper(ensg) %in% input_genes) %>%
      filter(abs(logFC) > 1)
  } else {
    selected <- limma_data %>%
      filter(ensg %in% input_genes) %>%
      filter(abs(logFC) > 1)
  }
  ggplot(selected,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                                         color=logFC))+
    
    geom_point()+scale_color_gradient2(high="red", low="blue")+
    theme(axis.text = element_text(size = 15)) +
    facet_wrap(vars(celltype), scales = "free_x")+
    theme_bw(12)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                     hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
    scale_size_continuous(range=c(0,5), limits = c(0,5))+
    ggtitle(paste0(names(input_genes),data_name))+
    ylab(paste0(names(input_genes)))
  ggsave(basedir(paste0(names(input_genes),data_name,".pdf")),w=25,h=10)
         }



IFN_genes_overlapping_logFC<-limmaRes%>%filter(toupper(ensg)%in% IFN_genes_overlapping)%>%
  filter(abs(logFC)>1)

ggplot(IFN_genes_overlapping_logFC,aes(x=coef,y=ensg,size=pmin(-log10(adj.P.Val),5),
                           color=logFC))+
  
  geom_point()+scale_color_gradient2(high="red", low="blue")+
  theme(axis.text = element_text(size = 15)) +
  facet_wrap(vars(celltype), scales = "free_x")+
  theme_bw(12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                   hjust = 1))+geom_point(alpha = 0.6) +  # Adjust transparency to see overlapping points better
  scale_size_continuous(range=c(0,5), limits = c(0,5))+
  ggtitle("IFN genes-Interaction(ex.vivo-in.vivo)")+
  ylab("IFN response genes overlapping in NTC")
ggsave(basedir("IFN_genes_in_NTC_logFC_in_ko_interaction.pdf"),w=25,h=10)

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

################################################################################
# fgsea ------------------------------------------------------------------------
################################################################################
gsea.res <- data.table()  # Initialize the result table
for (ct in unique(limmaRes_NTC$celltype)) {
  #for (de.grp in unique(limmaRes_NTC[limmaRes_NTC$celltype == ct, ]$coef)) {
  for (dbx in names(enr.terms)) {
    # limmaRes_NTC_subset <- limmaRes_NTC[limmaRes_NTC$celltype == ct, ]
    # Correct subsetting with both conditions
    subset_limmaRes_NTC <- limmaRes_NTC[limmaRes_NTC$celltype == ct, ]
    # Now using `with` correctly on the subsetted data
    stats <- with(subset_limmaRes_NTC, setNames(logFC, nm = genes))
    
    # Check for missing values in stats
    if (any(is.na(stats))) {
      next  # Skip this iteration if there are missing values in stats
    }
    
    # Perform fgsea analysis
    fgsea_output <- fgsea(
      pathways = enr.terms[[dbx]],
      stats = stats
      #minSize = 15,   # Example additional arguments, adjust as necessary
      #maxSize = 500,  # Example additional arguments, adjust as necessary
      #nperm = 1000    # Example additional arguments, adjust as necessary
    )
    
    # Check if fgsea output is not empty
    if (length(fgsea_output) > 0) {
      gsea.res <- rbind(gsea.res, data.table(fgsea_output,
                                             #coef = de.grp,
                                             celltype = ct,
                                             db = dbx))
    }
  }
}

################################################################################
#plot fgsea---------------------------------------------------------------------
################################################################################
gsea.res[is.nan(NES), NES := 0]
gsea.res.export <- gsea.res[padj < 0.05][,-c("log2err", "NES", "size", "pval"),with=F]
gsea.res.export$leadingEdge <- sapply(gsea.res.export$leadingEdge, function(vec) paste(vec[1:10], collapse = ","))

for(dbx in unique(gsea.res$db)){
  dat<-dirout(paste0(base,"FGSEA/",dbx))
  write.tsv(gsea.res.export[db == dbx], dat("GSEA_significant_",dbx,".tsv"))
}

dbx<-"MSigDB_Hallmark_2020"
# Prepare for plotting
for(dbx in unique(gsea.res$db)){
  
  pDT <- gsea.res[db == dbx]
  ## Splitting the task to handle both ends of the NES spectrum-positive and negative
  pw.display.pos <- unique(pDT[padj < 0.05][order(-NES)][, head(.SD, n=5), by=c("celltype")]$pathway)
  pw.display.neg <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=5), by=c("celltype")]$pathway)
  
  # Combine and remove duplicates across both positive and negative selections
  pw.display <- unique(c(pw.display.pos, pw.display.neg))
  pDT <- pDT[pathway %in% pw.display]
  #pDT <- hierarch.ordering(pDT, "pathway", "celltype", "NES", TRUE)
  pDT <- hierarch.ordering(pDT, "celltype", "pathway", "NES", TRUE)
  ggplot(pDT, aes(x=celltype, y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
    
    scale_color_gradient2(low="blue", mid="white", high="red") +
    geom_point() +
    scale_size_continuous(range=c(0,5), limits = c(0,5)) +
    theme_bw(12) +
    xRot() +
    #facet_wrap(vars(celltype))+#,space="free", scales="free") +)+
    labs(x = "celltype")+
    
    theme(axis.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
  #theme(strip.text.y=element_text(angle=0))
  dat<-dirout(paste0(base,"FGSEA/",dbx))
  ggsave(dat("GSEA_plot_",dbx,".png"), w=20,h=length(unique(pDT$pathway)) * 0.2 + 3, limitsize = FALSE)
}

################################


################################
enr.res.all <- bind_rows(enr.res.list, .id="coefx")
enr.res.all%>% filter(Odds.Ratio > 10 & 
                        Adjusted.P.value<0.01 & 
                        db=="MSigDB_Hallmark_2020")%>%
  ggplot(aes(x=coefx,y=Term,size=log2(Odds.Ratio),
             color=-log10(Adjusted.P.value)))+
  geom_point()+scale_color_gradient2(high="red", low="blue") 
ggsave(out("NTC_DEG_enrichment.pdf"))
############################
goi.enr <- enr.res.all %>%
  filter(Adjusted.P.value < 0.01 & Odds.Ratio > 10) %>%
  pull(Genes) %>%
  strsplit(";") %>%
  unlist()%>%
  unique()
enriched<-top_table_res %>%
  #mutate(gene = gmap[ensg,]$external_gene_name) %>%
  filter(toupper(genes) %in% goi.enr)

set1<-enriched %>% 
  filter(adj.P.Val<0.01)

ggplot(set1[200:450,],aes(x=id,y=genes,size=-log10(adj.P.Val),
                          color=logFC))+
  geom_point()+
  scale_color_gradient2(high="red", low="blue")+
  theme(text = element_text(size=10))
###################################################
goi.all<-top_table_res %>%
  filter(adj.P.Val<0.01) %>%
  group_by(id) %>%
  top_n(10,logFC) %>%
  pull(genes)

(p.coef <- top_table_res %>%
    filter(genes %in% goi.all) %>%
    ggplot(aes(y=genes, x=id, color=logFC, size=-log10(adj.P.Val))) + 
    geom_point() +
    scale_color_gradient2(high="red", low="blue") +
    theme_bw())
#############################################
#
dat.list <- list()
for(gg in goi.all){
  dat.list[[gg]] <- NTC_meta %>%
    mutate(E=scale(dataVoom$E[gg,])) %>%
    rownames_to_column("sample1") %>%
    remove_rownames()
}
(p.vals <- bind_rows(dat.list, .id="genes") %>%
    
    mutate(cell = as.character(cell)) %>%
    ggplot(aes(x=sample, y=genes, fill=E)) + 
    geom_tile() +
    facet_grid(. ~ cell, space ="free", scales = "free") +
    scale_fill_gradient2(low="blue", high="red"))
############################################
t<-top_table_res %>% filter(adj.P.Val<0.01)%>%
  top_n(1,abs(logFC))%>%
  pull(genes)

head(t)
test2<-NTC_meta%>%
  mutate(exp=dataVoom$E["Gbp4",])
ggplot(test2,aes(x=tissue,y=exp))+
  geom_point()+facet_grid(.~celltype)

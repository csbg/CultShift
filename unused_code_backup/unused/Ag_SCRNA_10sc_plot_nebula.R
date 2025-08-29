source("src/00_init.R")

require(nebula)
require(doMC)
source("src/FUNC_Monocle_PLUS.R")

basedir <- "Ag_SCRNA_10_single_Cluster_expression/"
outB <- dirout(basedir)


# Sample annotation -------------------------------------------------------
SANN <- fread(PATHS$SCRNA$ANN)
TISSUES <- PATHS$SCRNA$MONOCLE.NAMES

# load datasets -----------------------------------------------------------
mobjs <- list()
for(tissuex in PATHS$SCRNA$MONOCLE.NAMES){
  (load(PATHS$SCRNA$MONOCLE.DIR(tissuex)))
  mobjs[[tissuex]] <- monocle.obj
}


# Which analysis to take Clusters from? -----------------------------------
ANALYSIS <- "monocle.singleR"

# # Analysis of most/all Clusters in vivo ---------------------------------------------------------------
Clusters <- list(
  in.vivo = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_in.vivo_noMixscape.tsv"))$Clusters,
  ex.vivo = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_ex.vivo_noMixscape.tsv"))$Clusters,
  leukemia = fread(dirout_load("SCRNA_21_02_ClusterEnrichments_simple")("Guides_Fisher_Mixscape_basic_leukemia_noMixscape.tsv"))$Clusters
)

ct.use <- "everythingMonocle"
(tissue.name <- "in.vivo")
combined_list<-list()
for(tissue.name in TISSUES[1:2]){

  (timex <- mobjs[[tissue.name]]$timepoint[1])
  for(timex in unique(mobjs[[tissue.name]]$timepoint)){

    # Monocle object
    monocle.obj <- mobjs[[tissue.name]]
    monocle.obj <- monocle.obj[, monocle.obj$timepoint == timex]
    monocle.obj <- monocle.obj[, monocle.obj$sample != "WT-LSK_OP0_NM_7d_1"]
    #monocle.obj <- monocle.obj[, colnames(monocle.obj) %in% Clusters[[tissue.name]]]
    head(monocle.obj)
    # select Mixscape effects
    monocle.obj@colData$mixscape_class.global
    monocle.obj[,monocle.obj@colData$mixscape_class.global %in% c("KO", "NTC")]
    #monocle.obj <- monocle.obj[,monocle.obj@colData$mixscape_class.global %in% c("KO", "NTC")]
    head(monocle_obj)
    if(ncol(monocle.obj) == 0) next

    
    ann <- fread(dirout_load(paste0("SCRNA_20_Summary/",tissue.name, "_", ANALYSIS))("Annotation.tsv"))
    ann <- ann[rn %in% colnames(monocle.obj)]
    sort(unique(ann$Clusters))
    for(x in c("Ery", "B.Cluster", "CLP", "MEP(pert.)", "T.Cluster")){
      ann <- ann[!grepl(x, Clusters)]}
    #   sort(unique(ann$Clusters))
    
    # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
    monocle.obj <- monocle.obj[, ann$rn]
    monocle.obj$clusterDE <- ann$Clusters
    # Keep only KO and NTCs
    monocle.obj <- monocle.obj[,monocle.obj$mixscape_class.global %in% c("KO", "NTC")]
    monocle.obj <- monocle.obj[,monocle.obj$clusterDE %in% names(which(table(monocle.obj$clusterDE) > 30))]
    
    # Remove lowly expressed genes
    monocle.obj <- monocle.obj[Matrix::rowSums(counts(monocle.obj)) > 20,]
    
    head(monocle_obj)   
    # Assign clusters to use as covariate (to avoid seeing shifts in populations but changes within populations)
    
    combined_list[[tissue.name]] <-monocle.obj
    
    
    
    #############################################
    
  }
}
monocle_ex.in<-combine_cds(
  combined_list,
  # #keep_all_genes = TRUE,
  # #Cluster_names_unique = FALSE,
  sample_col_name = "samples"
  
)
#############################
#Filter Clusters of interest
#############################
guides<-unique(monocle_ex.in[,monocle_ex.in$tissue=="ex.vivo"]$guide)

monocle_ex.in<-monocle_ex.in[,monocle_ex.in$guide %in% guides]
monocle_ex.in$guide<-gsub("_.+$","",monocle_ex.in$guide)

source("src/Ag_ScRNA_09_Limma_on_Pseudobulk_NTC_enrichment.R")
head(enr.res.down.all$db)
t1<-enr.res.down.all[enr.res.down.all$db=="MSigDB_Hallmark_2020" &
                       enr.res.down.all$Term=="Interferon Alpha Response",]
gene_list<-arrange(t1,desc(Odds.Ratio))
gene_list<-str_split(gene_list[1,]$Genes,";")
#gene_data<-gene_data[gene_data$Cluster %in% c("GMP","Mono","MEP","Eo_Ba","MkP","EBMP","HSC","Gran."),]

rownames(exprs(monocle_ex.in))[toupper(rownames(exprs(monocle_ex.in)))%in% gene_list]

gene_of_interest<-genes[1]

for(gene_of_interest in genes){
  names(exprs(monocle_ex.in)[gene_of_interest, ])
  Tissue = monocle_ex.in[gene_of_interest,]$tissue
  names(Tissue)<-names(exprs(monocle_ex.in)[gene_of_interest, ])
  Cluster=monocle_ex.in[gene_of_interest, ]$clusterDE
  names(Cluster)<-names(exprs(monocle_ex.in)[gene_of_interest, ])
  Guide=monocle_ex.in[gene_of_interest, ]$guide
  names(Guide)<-names(exprs(monocle_ex.in)[gene_of_interest, ])
  gene_data <- data.table(GeneExpression = exprs(monocle_ex.in)[toupper(gene_of_interest), ],
                          Tissue,Cluster,Guide)

  gene_data<-gene_data[gene_data$Cluster!="B-Cluster" & 
                       gene_data$Cluster!="CLP" & gene_data$Cluster!="Ery" & 
                       gene_data$Cluster!="Gran." & gene_data$Cluster!="unclear",]

#change names of celltype with sp.character iexpression data

  gene_data$Cluster<-gsub("Eo/Ba","Eo_Ba",gene_data$Cluster)
  gene_data$Cluster<-gsub(" ","",gene_data$Cluster)
  for (ct in unique(gene_data$Cluster)) {
    p<-ggplot(gene_data[gene_data$Cluster==ct],aes(x=Tissue,y=GeneExpression,col=Tissue))+
    geom_violin()+facet_wrap(vars(Guide))+ 
    ylim(c(-1,8))+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title=ct, 
        x =gene_of_interest)
    out <- dirout(paste0(basedir,gene_of_interest))
    ggsave(out(paste0(ct,"single_Cluster_expression_across_KO.pdf")),plot = p)
   
  p<-ggplot(gene_data[Guide=="NTC"][gene_data[Guide=="NTC"]$Cluster==ct],aes(x=Tissue,y=GeneExpression,col=Tissue))+
    geom_violin()+facet_wrap(vars(Guide))+ 
    ylim(c(-1,8))+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
    labs(title=ct, 
         x =gene_of_interest)
  out1 <- dirout(paste0(basedir,"NTC"))
  ggsave(out1(paste0(ct,"single_Cluster_expression_across_KO.pdf")),plot = p)
}
}

##################


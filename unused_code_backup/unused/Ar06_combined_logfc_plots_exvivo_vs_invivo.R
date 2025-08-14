source("src/00_init.R")
library(tidyselect)
library(pheatmap)
library(cowplot)
library(patchwork)
library(plotly)
library(htmlwidgets)
library(fgsea)
library(RColorBrewer)
library(clusterProfiler)
library("org.Mm.eg.db")
library(dplyr)
library(data.table)
library(msigdbr)
gseadir<-"Ar_07_ORA_from_ClusterProfiler_own_plot/"
basedir <- "Ar_06_ORA_Cluster_profiler/"
inp<-"SCRNA_33_DE_Nebula_testClustering/"
out <- dirout(basedir)
ex<-"9d"
inv<-"14d"
leuk<-"6d"
#logFC
#for (ex in c("7d","9d")){
 # for (inv in c("14d","28d")){
exvivo = fread(dirout(paste0(inp, "ex.vivo_",ex,"_useClusters"))
               ("DEG_Results_all.tsv")) %>% 
 dplyr::select(c("guide","gene_id","estimate"))#%>% (guide%in%invivo$guide)
    
invivo = fread(dirout(paste0(inp, "in.vivo_",inv,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      dplyr::select(c("guide","gene_id","estimate"))%>% filter(guide%in%exvivo$guide)
    
    leukemia=fread(dirout(paste0(inp, "leukemia_",leuk,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      dplyr::select(c("guide","gene_id","estimate"))
    

genes<-unique(exvivo$gene_id)
  combined_invivo<-merge(exvivo,invivo,by=c("guide","gene_id"))
  colnames(combined_invivo)<-c("guide","geneid","exvivo","invivo")
  
  combined_invivo<-as.data.frame(combined_invivo)
  combined_invivo["group"]<-"a:n.s"
  
  combined_invivo[which(combined_invivo['exvivo'] >= 1 & combined_invivo['invivo'] >= 1.0 ),"group"] <- "i:up"
  combined_invivo[which(combined_invivo['exvivo'] <= -1 & combined_invivo['invivo'] <= -1.0 ),"group"] <- "h:down"
  combined_invivo[which(combined_invivo['exvivo'] <= -1 & combined_invivo['invivo'] >= 1.0 ),"group"] <- "g:inv_up.ex_down"
  combined_invivo[which(combined_invivo['invivo'] <= -1 & combined_invivo['exvivo'] >= 1.0 ),"group"] <- "f:inv_down.ex_up"
  combined_invivo[which(combined_invivo['invivo'] >= -1 & combined_invivo['invivo'] <=1 & combined_invivo['exvivo'] >= 1.0 ),"group"] <- "e:inv_low.ex_up"
  combined_invivo[which(combined_invivo['invivo'] >= -1 & combined_invivo['invivo'] <=1 & combined_invivo['exvivo'] <= -1.0 ),"group"] <- "d:inv_low.ex_down"
  combined_invivo[which(combined_invivo['exvivo'] >= -1 & combined_invivo['exvivo'] <=1 & combined_invivo['invivo'] >= 1.0 ),"group"] <- "c:ex_low.in_up"
  combined_invivo[which(combined_invivo['exvivo'] >= -1 & combined_invivo['exvivo'] <=1 & combined_invivo['invivo'] <= -1.0 ),"group"] <- "b:ex_low.in_down"
  write.table(combined_invivo,out("combined_invivo_logfc.tsv"))
  head(combined_invivo)
  name<-"Ash1l"
  for (name in unique(combined_invivo$guide)){
    t<-combined_invivo[combined_invivo$guide==name,]
    head(t)
    my_col=c("gray","royalblue4","deepskyblue3","#DB7093","#556B2F","#5F9EA0","#CD5C5C","#8B6508","green")
    g<-ggplot(data=t,aes(x=exvivo,y=invivo))+
      geom_point(data=t,aes(x=exvivo,y=invivo,col=group,text=geneid))+xlim(-4,4)+ylim(-4,4)
    ggplotly(g,)
    ggsave(out(paste0(name,"exvivo_vs_invivo.png")))
  
    fig<-ggplotly(g,tooltip = "text")
    #fig <- plotly(t, x = ~exvivo, y = ~invivo, name = 'trace 0', type = 'scatter',color = ~group, text=~geneid) 
    
  
    htmlwidgets::saveWidget(fig, out(paste0(name,".html")))
  }
  combined_invivo<-combined_invivo %>% 
    mutate(entrez = mapIds(org.Mm.eg.db, keys = combined_invivo$geneid,
                           column = "ENTREZID", keytype = "SYMBOL")) %>% 
    filter(!is.na(entrez))
  head(combined_invivo)
  #############################
  #
    
  
  ggplot(combined_invivo[guide=="Copz2"], aes(, value)) +   
    geom_bar(aes(fill = variable), position = "dodge", stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  #########################################
  OrgDb<-"org.Mm.eg.db"
  ff<-data.frame()
  all_genes <- unique(combined_invivo$entrez) # directly selects the gene column
  
  combined_invivo<-as.data.table(combined_invivo)
  KO<-unique(combined_invivo$guide)[1]
  type<-unique(combined_invivo$group)[2]
  for (KO in unique(combined_invivo$guide)){
    for (type in unique(combined_invivo$group)){
      if (length(combined_invivo[combined_invivo$guide==KO & combined_invivo$group==type,]$entrez) == 0) {
        print("The data frame is empty.")
      } else {
  
        genes<-combined_invivo[combined_invivo$guide==KO & combined_invivo$group==type,]$entrez
      
      
        ego<-enrichGO(
        gene=genes,
        OrgDb,
        keyType = "ENTREZID",
        ont = "BP",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        universe=all_genes,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = 500,
        readable = FALSE,
        pool = FALSE)
        
        outdir <- dirout(paste0(basedir, KO))
        if (nrow(ego) == 0) {
          print("The data frame is empty.")
        } else {
      # Perform some operation on the non-empty data frame
      # For example, print the first few rows
      
        #plot<-dotplot(ego,orderBy="GeneRatio",title=paste0(KO,"_",type))
    
      #ggsave(outdir(paste0(KO,"_",type,".png")),h=10,w=12)
      ##################################################
      main<-
      
      t1<-ego@result
      nrow(t1[1:7])
      t1<-as.data.table(t1)
      #ff <- rbind(ff,t1[1:7], idcol = c(KO,type))
      write.tsv(t1[1:7],outdir(paste0(KO,"_",type,"_","data_table.tsv")))
        }
      }}}
  #}
    
 
 
 
 ego@result[ego@result$p.adjust<0.05,]
  
  # MSigDB R package
  library(msigdbr)
  msigdbr::msigdbr_collections() # available collections
  # Subset to Human GO-BP sets
  BP_db <- msigdbr(species = "Mus musculus", 
                   category = "C5", subcategory = "GO:BP")
  head(BP_db)
  BP<-BP_db[,c("gs_exact_source","gs_description")]
  colnames(BP)<-c("Id","description")
  BP<-as.data.frame(BP)
  t<-c(unique(BP$Id))
  
  BP<-BP[!duplicated(BP$Id),]
  
  ff<-list.files()
  GSEA_total<-data.table()
  for (KO in unique(combined_invivo$guide)){
    ff <- list.files(dirout_load("Ar_06_ORA_Cluster_profiler/")(KO), pattern=".tsv", full.names = TRUE, recursive = TRUE)
    names(ff) <- gsub("^.+\\/", "", (ff))
    ff<-lapply(ff, fread)
    GSEA_res <- rbindlist(ff, idcol = "KO")
    GSEA_total<-rbind(GSEA_total,GSEA_res)
  }
  test<-data_table.tsv
  head(GSEA_total)
  GSEA_total<-as.data.table(GSEA_total)
  GSEA_total[,type:=KO]
  GSEA_total[,type:=gsub("^.*?_","",type)]
  
  
  GSEA_total[,type:=gsub("_data_table.tsv","",type)]
  GSEA_total[,KO:=gsub("_.+$","",KO)]

  GSEA_total<-GSEA_total[!grepl('n.s', GSEA_total$type),]

  dt <- data.table(KO = rep(unique(combined_invivo$guide),each=8), type = rep(c("i:up","h:down","g:inv_up.ex_down","f:inv_down.ex_up","e:inv_low.ex_up",
                                                                         "d:inv_low.ex_down","c:ex_low.in_up"
                                                                         ,"b:ex_low.in_down"),length(unique(combined_invivo$guide))))
  
  testdf<-merge(GSEA_total,dt,by=c("KO","type"),fill=T,all.y=T)
  
  pDT<-testdf
  head(pDT)
  pDT <- hierarch.ordering(pDT, "Description", "KO", "type", TRUE)
  pDT<-pDT[-log10(p.adjust)>2,]
  pDT[,a:=as.numeric(gsub("(\\d+)/(\\d+)","\\1",GeneRatio))]
  pDT[,b:=as.numeric(gsub("(\\d+)/(\\d+)","\\2",GeneRatio))]
  pDT[,c:=as.numeric(gsub("(\\d+)/(\\d+)","\\1",BgRatio))]
  pDT[,d:=as.numeric(gsub("(\\d+)/(\\d+)","\\2",BgRatio))]
  pDT[,odds_ratio:= (a/b)/(c/d)]
  head(pDT)
  #length(unique(pDT[-log10(p.adjust)>4,]$Description))
 #pDT<-pDT[unique(pDT$KO)[10],]
  head(pDT)
  pDT<-pDT[p.adjust<=0.005]
  ggplot(pDT, aes(x=type, y=Description,size=pmin(15, -log10(p.adjust)))) +
    geom_point(aes(color=odds_ratio))+scale_fill_brewer(palette = "RdYlBu") +
    scale_size_continuous(range=c(0,3), limits = c(0,15)) +
    theme_bw(12) +
    xRot()+
    facet_grid(cols = vars(KO), scales="free") +
    
    theme(strip.text.y=element_text(angle=0),axis.text.y = element_text(size=20),
          axis.text.x = element_text(size=20)
        )
  out<-dirout(gseadir)
  ggsave(out("enrichment_trial_cut0ff_0.005.pdf"),h=20,w=20)
  
##########################################################################
## Cluster GO-BP ORA with clusterProfiler package
  all_genes <- unlist(gcSample)
  head(all_genes)
  head(gcSample)
  universe <- all_genes[Biobase::isUnique(all_genes)] # all unique genes
  head(universe)
  # List with only unique genes
  gcUnique <- lapply(gcSample, function(group_i) {
    group_i[group_i %in% universe]})
  
cp_ora<-data.table()
  head(BP_list)
  KO <-unique(combined_invivo$guide)[1]
  for (KO in unique(combined_invivo$guide)){
  cp_ora<-compareCluster(
  geneClusters = combined_invivo[combined_invivo$guide==KO,]$geneid,
  fun = "enrichGO", # ORA function to apply to each cluster
  # Arguments below are passed to enrichGO
  OrgDb = "org.Mm.eg.db", 
  keyType = "SYMBOL", 
  ont = "BP", # BP, CC, MF, or ALL for all ontologies
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.10, # do not filter by q-value
  pAdjustMethod = "BH", # p-values are adjusted within clusters
  universe = all_genes, # all genes
  minGSSize = 15, 
  maxGSSize = 500)
  #cp_ora<-rbind(cp_ora,res)
  }
  library(clusterProfiler)
  data("gcSample") # data for examples
  
head(gcSample)
# First 6 entries sorted by cluster and p-value
cp_ora@compareClusterResult %>% 
  arrange(Cluster, pvalue) %>% 
  head()


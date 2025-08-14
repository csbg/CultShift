###########
#Script that I have been working on 13-12-23
#Trying to change the naming of the enrichr lis with type only having the group info
#check from line 250
# May be I can use the ko"_"type when creating the table
#packages
source("src/00_init.R")
library(tidyselect)
library(pheatmap)
library(patchwork)
library(plotly)
library(htmlwidgets)
library(clusterProfiler)
library("org.Mm.eg.db")
library(msigdbr)
library(dplyr)

#Input and out directories

ORAdir<-"Ar_07_ORA_from_ClusterProfiler/"
inp<-"SCRNA_33_DE_Nebula_testClustering/"

#corresponding stimulation duration
ex<-"9d"
inv<-"14d"
leuk<-"6d"

# subset input data
exvivo = fread(dirout(paste0(inp, "ex.vivo_",ex,"_useClusters"))
               ("DEG_Results_all.tsv")) %>% 
  dplyr::select(c("guide","gene_id","estimate"))
    
invivo = fread(dirout(paste0(inp, "in.vivo_",inv,"_useClusters"))
               ("DEG_Results_all.tsv")) %>% 
  dplyr::select(c("guide","gene_id","estimate"))%>% 
  filter(guide %in% exvivo$guide)
    
leukemia=fread(dirout(paste0(inp, "leukemia_",leuk,"_useClusters"))
               ("DEG_Results_all.tsv")) %>% 
  dplyr::select(c("guide","gene_id","estimate"))




# combine in.vivo and ex.vivo
combined<-merge(exvivo,invivo,by=c("guide","gene_id"))
colnames(combined)<-c("guide","geneid","exvivo","invivo")

combined_invivo<-as.data.frame(combined)
combined_invivo["group"]<-"a:n.s"

combined_invivo[which(combined_invivo['exvivo'] >= 1 & 
                        combined_invivo['invivo'] >= 1.0 ),"group"] <- "i:up"
combined_invivo[which(combined_invivo['exvivo'] <= -1 & 
                        combined_invivo['invivo'] <= -1.0 ),"group"] <- "h:down"
combined_invivo[which(combined_invivo['exvivo'] <= -1 & 
                        combined_invivo['invivo'] >= 1.0 ),"group"] <- "g:invivo up & exvivo down"
combined_invivo[which(combined_invivo['invivo'] <= -1 & 
                        combined_invivo['exvivo'] >= 1.0 ),"group"] <- "f:invivo down & exvivo up"
combined_invivo[which(combined_invivo['invivo'] >= -1 & combined_invivo['invivo'] <=1 &
                        combined_invivo['exvivo'] >= 1.0 ),"group"] <- "e:invivo low & exvivo up"
combined_invivo[which(combined_invivo['invivo'] >= -1 & combined_invivo['invivo'] <=1 &
                        combined_invivo['exvivo'] <= -1.0 ),"group"] <- "d:invivo low & exvivo down"
combined_invivo[which(combined_invivo['exvivo'] >= -1 & combined_invivo['exvivo'] <=1 &
                        combined_invivo['invivo'] >= 1.0 ),"group"] <- "c:exvivo low & invivo up"
combined_invivo[which(combined_invivo['exvivo'] >= -1 & combined_invivo['exvivo'] <=1 & 
                        combined_invivo['invivo'] <= -1.0 ),"group"] <- "b:exvivo low & invivo down"
write.table(combined_invivo,out("combined_invivo_logfc.tsv"))

  name<-"Ash1l"
  for (name in unique(combined_invivo$guide)){
    t<-combined_invivo[combined_invivo$guide==name,]
    head(t)
    my_col=c("gray","royalblue4","deepskyblue3","#DB7093","#556B2F","#5F9EA0","#CD5C5C","#8B6508","green")
    g<-ggplot(data=t,aes(x=exvivo,y=invivo))+
      geom_point(data=t,aes(x=exvivo,y=invivo,col=group,text=geneid))+xlim(-4,4)+ylim(-4,4)
    g
    ggsave(out(paste0(name,"exvivo_vs_invivo.png")))
  
    fig<-ggplotly(g,tooltip = "text")
    #fig <- plotly(t, x = ~exvivo, y = ~invivo, name = 'trace 0', type = 'scatter',color = ~group, text=~geneid) 
    
  
    htmlwidgets::saveWidget(fig, out(paste0(name,".html")))
  }
  #all plots together
  ggplot(data=combined_invivo,aes(x=exvivo,y=invivo))+
    geom_point(data=combined_invivo,aes(x=exvivo,y=invivo,col=group,text=geneid))+xlim(-4,4)+ylim(-4,4)+
    facet_wrap(vars(guide),nrow=4)
  ggsave(out(paste0("exvivo_vs_invivo_all.png")))
 ?facet_wrap
  ###################
  
   combined_invivo<-combined_invivo %>% 
    mutate(entrez = mapIds(org.Mm.eg.db, keys = combined_invivo$geneid,
                           column = "ENTREZID", keytype = "SYMBOL")) %>% 
    filter(!is.na(entrez))
  head(combined_invivo)
  #############################
    
  #########################################
  OrgDb<-"org.Mm.eg.db"
  ff<-data.frame()
  all_genes <- unique(combined_invivo$entrez) # directly selects the gene column
  combined_invivo<-as.data.table(combined_invivo)
  ####################################################
  #enrichGO---------------------------
  ####################################################
  
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
        
        outdir <- dirout(paste0(ORAdir, KO))
        if (nrow(ego) == 0) {
          print("The data frame is empty.")
        } else {
          
          t1<-ego@result
          
          t1<-as.data.table(t1)
          
          write.tsv(t1,outdir(paste0(KO,"_",type,"_","data_table.tsv")))
        }
      }}}
  #}
  
  ego@result[ego@result$p.adjust<0.05,]
  
  ff<-list.files()
  ORA_total<-data.table()
  KO<-unique(combined_invivo$guide)[1]
  type<-unique(combined_invivo$group)[1]
  
  for (KO in unique(combined_invivo$guide)){
    for (type in unique(combined_invivo$group)){
      list.files(dirout_load("Ar_06_ORA_Cluster_profiler")(KO), pattern=".tsv", full.names = TRUE, recursive = TRUE)
      ff <- list.files(dirout_load("Ar_06_ORA_Cluster_profiler")(KO), pattern="MF", full.names = TRUE, recursive = TRUE)
      names(ff) <- gsub("^.+\\/", "", (ff))
      ff<-lapply(ff, fread)
      ORA_res <- rbindlist(ff, idcol = c("KO","type"))
      ORA_total<-rbind(ORA_total,ORA_res)
    }
  }
  
  head(ORA_total)
  ORA_total<-as.data.table(ORA_total)
  
  ORA_total[,type:=KO]
  ORA_total[,type:=gsub("^.*?_","",type)]
  ORA_total[,type:=gsub("_data_table.tsv","",type)]
  ORA_total[,KO:=gsub("_.+$","",KO)]
  ORA_total<-ORA_total[!grepl('n.s', ORA_total$type),]
  head(ORA_total)
  unique(ORA_total$type)
  ####
  #CHANGE THE  TYPE ACCORDINGLY
  dt <- data.table(KO = rep(unique(combined_invivo$guide),each=8), type =rep(c("i:up","h:down","g:inv_up.ex_down","f:inv_down.ex_up","e:inv_low.ex_up",
                                                                               "d:inv_low.ex_down","c:ex_low.in_up"
                                                                               ,"b:ex_low.in_down"),length(unique(combined_invivo$guide))))
  testdf<-merge(ORA_total,dt,by=c("KO","type"),fill=T,all.y=T)
  
  pDT<-testdf
  
  #pDT <- hierarch.ordering(pDT, "Description", "KO", "type", TRUE)
  head(pDT)
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
  pDT<-pDT[p.adjust<=0.05]
  ggplot(pDT, aes(x=KO, y=Description,size=pmin(15, -log10(p.adjust)))) +
    geom_point(aes(color=log2(odds_ratio)))+scale_fill_brewer(palette = "RdYlBu") +
    scale_size_continuous(range=c(0,3), limits = c(0,15)) +
    theme_bw(12) +
    xRot()+
    facet_grid(cols = vars(type), scales="free") +
    theme(strip.text.y=element_text(angle=0),axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10)
    )
  
  out_enr<-dirout(enrichment)
  ggsave(out_enr("GSEA_plot_",dbx,".pdf"),w=30,h=15)
  ##################################################################################
  
  
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
          ont = "MF",
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          universe=all_genes,
          qvalueCutoff = 1,
          minGSSize = 10,
          maxGSSize = 500,
          readable = FALSE,
          pool = FALSE)
        
        outdir <- dirout(paste0(ORAdir, KO))
        if (nrow(ego) == 0) {
          print("The data frame is empty.")
        } else {
          
          t1<-ego@result
          nrow(t1[1:7])
          t1<-as.data.table(t1)
          
          write.tsv(t1,outdir(paste0(KO,"_MF_",type,"_","data_table.tsv")))
        }
      }}}
  #}
  
  ##################################################################  
  
  enr.res.ko.list<-list()
  combined_invivo<-combined_invivo[combined_invivo$group=="g:inv_up.ex_down" | combined_invivo$group=="d:inv_low.ex_down",]
  load(PATHS$RESOURCES$Enrichr.mouse)
  
  KO<-unique(combined_invivo$guide)[1]
  type<-unique(combined_invivo$group)[2]
  for (KO in unique(combined_invivo$guide)){
    for(dbx in names(enr.terms)){
      #if (length(combined_invivo[combined_invivo$guide==KO ]$entrez) == 0) {
      #print("The data frame is empty.")
      #} else {
      # Function to perform enrichment analysis for a given group and id
      genes_down_ex<-combined_invivo[guide==KO,]%>%pull(geneid)%>%unique()
      head(genes_down_ex)
      
      enr.res_ko<-enrichr(genes_down_ex, databases = dbx)
      enr.res_ko <- bind_rows(enr.res_ko, .id="db")
      # Extract unique ids
      # Store results in the list
      enr.res.ko.list[[KO]] <- enr.res_ko
    }}
  enr.res.all <- bind_rows(enr.res.ko.list, .id="KO")
  # Combine the list into a data frame
  enr_results <- bind_rows(enr.res.ko.list)
  
  
  load(PATHS$RESOURCES$Enrichr.mouse)
  dbx<-unique(names(enr.terms))[1]
  head(unique(names(enr.terms)))
  for(ct in unique(top_table_res$id)){
    for(dbx in names(enr.terms)){
      gsea.res <- rbind(gsea.res, data.table(fgsea(
        pathways=enr.terms[[dbx]], 
        stats=with(top_table_res[id == ct], setNames(logFC, nm=genes))), 
        grp=ct,
        db=dbx))
    }
  }
  
  for(dbx in unique(gsea.res$db)){
    write.tsv(gsea.res.export[db == dbx], out("GSEA_significant_",dbx,".tsv"))
  }
  
  
  # Prepare for plotting
  for(dbx in unique(gsea.res$db)){
    pDT <- gsea.res[db == dbx]
    
    pw.display <- unique(pDT[padj < 0.05][order(NES)][, head(.SD, n=10), by=c("grp")]$pathway)
    pDT <- pDT[pathway %in% pw.display]
    pDT <- hierarch.ordering(pDT, "pathway", "grp", "NES", TRUE)
    
    ggplot(pDT, aes(x=grp, y=pathway, color=NES, size=pmin(5, -log10(padj)))) +
      geom_point() + scale_color_gradient2(low="blue", mid="white", high="red") +
      geom_point(data=pDT[padj < 0.05], shape=1, color="black") +
      scale_size_continuous(range=c(0,5), limits = c(0,5)) +
      theme_bw(12) +
      xRot() +
      facet_grid(. ~ db, space="free", scales="free") +
      theme(strip.text.y=element_text(angle=0))
    ggsave(out("GSEA_plot_",dbx,".pdf"), w=20,h=length(unique(pDT$pathway)) * 0.2 + 3, limitsize = FALSE)
  }
  
  #####################################
  #
  ## MSigDB R package
  library(msigdbr)
  msigdbr_species()  
  msigdbr::msigdbr_collections()
  look<-msigdbr::msigdbr_collections() # available collections
  # Subset to mouse
  
  BP_db <- msigdbr(species = "Mus musculus")
  unique(BP_db$gs_cat)
  
  BP_db <- msigdbr(species = "Mus musculus", 
                   category = "C2",subcategory = "CP:REACTOME")
  
  BP<-BP_db[,c("gs_description","entrez_gene")]
  head(BP)
  colnames(BP)<-c("TERM","GENE")
  BP<-as.data.frame(BP)
  head(BP)
  #BP<-BP[!duplicated(BP$Id),]
  all_genes <- unique(combined_invivo$entrez) # directly selects the gene column
  genes<-combined_invivo[combined_invivo$guide==KO & combined_invivo$group==type,]$entrez
  enr<-enricher(
    gene=genes,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe=all_genes,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    TERM2GENE=BP,
    TERM2NAME = NA
  )
  
  ??TERM2GENE
  head(enr)
  ##############################################
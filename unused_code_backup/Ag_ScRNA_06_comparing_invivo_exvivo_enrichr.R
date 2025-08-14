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
enrichment<-"Ar_07_ORA_from_enrichment/"

ORAdir<-"Ar_07_KO_DE_groups/"
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
test<-combined_invivo%>%pivot_longer(,cols=c("exvivo","invivo"),names_to="tissue",values_to = "logFC")

head(test)
dim(test)
test[test$geneid=="C1ra",]
ggplot(test[test$geneid=="C1ra",],aes(x=tissue,y=logFC))+
         geom_col()+
  facet_wrap(vars(guide))
ggsave(out("C1ra_across_KO.pdf"))

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


###############################################################
#Enrichr
##########################################################################
#all_genes <- unique(combined_invivo$entrez) # directly selects the gene column
enr.res.list<-list()

combined_invivo<-as.data.table(combined_invivo[group!="a:n.s"])
KO<-unique(combined_invivo$guide)[1]
type<-unique(combined_invivo$group)[2]
for (KO in unique(combined_invivo$guide)){
  for (type in unique(combined_invivo$group)){
    if (length(combined_invivo[combined_invivo$guide==KO & combined_invivo$group==type,]$geneid) == 0) {
      print("The data frame is empty.")
    } else {
      genes<-combined_invivo[combined_invivo$guide==KO & combined_invivo$group==type,]$geneid
      
      enr.res <- enrichr(genes,databases = c("MSigDB_Hallmark_2020", "GO_Biological_Process_2021","KEGG_2019_Mouse","WikiPathways_2019_Mouse" ))
      
      # The results will be a list, where each entry is one database. We will combine those into one long table
      enr.res <- rbindlist(enr.res, idcol = c("db"))
      ko_type<-paste0(KO,"_",type)
      # Store results in the list
      enr.res.list[[ko_type]] <- enr.res
    }
  }}


enr.res.all <-rbindlist(enr.res.list,idcol ="KO_type")

enr.res.all[,type:=KO_type]

enr.res.all[,type:=gsub("^.*?_","",type)]
enr.res.all[,type:=gsub("^.+:","",type)]
enr.res.all[,KO:=gsub("_.+$","",KO_type)]
enr.res.all[,KO:=gsub("^.+:","",KO)]
#Take the terms and plot for all--------------

databases = c("MSigDB_Hallmark_2020", "GO_Biological_Process_2021","KEGG_2019_Mouse","WikiPathways_2019_Mouse" )
unique(names(enr.terms))
dtx<-databases[1]
for (dtx in databases){
  outdir <- dirout(paste0(enrichment,dtx))
  
  plot_enr<-enr.res.all%>% filter(Odds.Ratio > 10 & 
                                    Adjusted.P.value<0.05 & 
                                    db==dtx)
  plot_enr<- enr.res.all[enr.res.all$Term %in% plot_enr$Term,]
  
  
  ggplot(plot_enr,aes(x=KO,y=substr(Term,1,30),size=log2(Odds.Ratio),
                      color=pmin(5,-log10(Adjusted.P.value))))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_point()+scale_color_gradient2(high="red", low="blue")+facet_grid(cols=vars(type))
  ggsave(outdir("enrichment_all.pdf"),w=30,h=15)
  ###############
  #specific
  types_ex<-c("invivo up & exvivo down","invivo down & exvivo up","exvivo low & invivo up","invivo low & exvivo up")
  plot_specific<-plot_enr[type %in% types_ex]
  plot_enr<- enr.res.all[enr.res.all$Term %in% plot_enr$Term,]
 
  ggplot(plot_specific,aes(x=KO,y=substr(Term,1,30),size=log2(Odds.Ratio),
                      color=pmin(5,-log10(Adjusted.P.value))))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_point()+scale_color_gradient2(high="red", low="blue")+facet_grid(cols=vars(type))
  ggsave(outdir("enrichment_specific.pdf"),w=30,h=15)
}

#####################################################################################################################

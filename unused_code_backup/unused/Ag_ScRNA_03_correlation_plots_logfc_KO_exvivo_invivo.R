source("src/00_init.R")
library(tidyselect)
library(pheatmap)


basedir<-"Ag_logfc_correlation/"
inp<-"SCRNA_33_DE_Nebula_testClustering/"

ex<-"9d"
inv<-"14d"
leuk<-"6d"
#Use clusters refers to that clusters were one of the explanatory variables while doing DE.So 
#the DE genes are the ones which are different within clusters
for (ex in c("7d","9d")){
  for (inv in c("14d","28d")){
    exvivo = fread(dirout(paste0(inp, "ex.vivo_",ex,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate","q_value"))
      exvivo<-exvivo %>% select(c("guide","gene_id","est.q_val"))
    
    invivo = fread(dirout(paste0(inp, "in.vivo_",inv,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate","q_value"))
    invivo<-invivo %>% select(c("guide","gene_id","est.q_val"))
    
    leukemia=fread(dirout(paste0(inp, "leukemia_",leuk,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate","q_value"))
    leukemia<-leukemia %>% select(c("guide","gene_id","est.q_val"))
   sign(exvivo$estimate)
       
    combined<-merge(exvivo,invivo,by=c("guide","gene_id")) %>% pivot_wider(
      names_from = guide,
      values_from=c(estimate.x,estimate.y)
    )
    
    colnames(combined)[grep("estimate.x",colnames(combined))]<-gsub("estimate.x","exvivo",colnames(combined)[grep("estimate.x",colnames(combined))])
    colnames(combined)[grep("estimate.y",colnames(combined))]<-gsub("estimate.y","invivo",colnames(combined)[grep("estimate.y",colnames(combined))])
    
    out1 <- dirout(paste0(basedir,"/exvivo_",ex,"invivo",inv))
      for (corr in c("pearson","spearman")){
        
        test<-cor(combined[grep("exvivo",colnames(combined))],y=combined[grep("invivo",colnames(combined))],
                method = corr)
        test2<-cor(combined[grep("invivo",colnames(combined))],y=combined[grep("invivo",colnames(combined))],
                   method = corr)
        
       
        pdf(out1(paste0("all_correlations_",corr,"_.pdf")))
        pheatmap(test,cluster_rows = T,cluster_cols = T,)
        pheatmap(test2,cluster_rows = T,cluster_cols = T)
        dev.off()
        stopifnot(all(gsub("^.*_","",colnames(test))==gsub("^.*_","",rownames(test))))            
        diag_correlation<-data.frame(test[row(test)==col(test)])
        colnames(diag_correlation)<-"correlation"
        rownames(diag_correlation)<-gsub("exvivo_","",rownames(test))
            
        write.tsv(data.table(data.frame(diag_correlation), keep.rownames = TRUE),
                out1(paste0("correlation_logfc_invivo_",inv,"vs_exvivo",ex,"_",corr,"_.tsv")))
      head(diag_correlation)
        ggplot(data = diag_correlation, mapping = aes(x = correlation,
         y = rownames(diag_correlation))) + geom_col(fill="deepskyblue4") + xlab(label = "Sample")+
          #scale_fill_gradient2(low = "royalblue",mid = "gray", high = "deeppink4",midpoint = ave(diag_correlation$correlation))
        
     
        ggsave(out1("pheatmap_corresponding_ko","invivo",inv,"exvivo",ex,"_",corr,"_.pdf"))
      
    }
}}


##############################

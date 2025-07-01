source("src/00_init.R")
library(tidyselect)
library(pheatmap)
library(ggpubr)
library("tidyverse")

basedir <- "SCRNA_33_DE_Nebula_testClustering/"
outB <- dirout(paste0(basedir,"/Correlations"))

#Also open the other files
#save pheatmap

for (ex in c("7d","9d")){
  for (inv in c("14d","28d")){
    exvivo_mye = fread(dirout(paste0(basedir, "ex.vivo_",ex,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
    exvivo_all = fread(paste0("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_33_DE_Nebula_testClustering/ex.vivo_",ex,
                              "_useClusters/DEG_Results_all.tsv"))
    invivo = fread(dirout(paste0(basedir, "in.vivo_",inv,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
    invivo_all = fread(paste0("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_33_DE_Nebula_testClustering/in.vivo_",inv,
                              "_useClusters/DEG_Results_all.tsv"))
    
    for (df in c("mye","all")){
      paste0("combined_",df)<-merge(paste0("exvivo_",df),paste0("invivo_",df),by=c("guide","gene_id")) %>% pivot_wider(
      names_from = guide,
      values_from=c(estimate.x,estimate.y)
    )
    colnames(combined_mye)[grep("estimate.x",colnames(combined_mye))]<-gsub("estimate.x","exvivo",colnames(combined_mye)[grep("estimate.x",colnames(combined_mye))])
    colnames(combined_mye)[grep("estimate.y",colnames(combined_mye))]<-gsub("estimate.y","invivo",colnames(combined_mye)[grep("estimate.y",colnames(combined_mye))])
    
    corrplot::corrplot(cor(x=combined_mye[grep("exvivo",colnames(combined_mye))],
                           y=combined_mye[grep("invivo",colnames(combined_mye))],method = "pearson"),diag = T)
    for (corr in c("pearson","spearman")){
      test<-cor(combined_mye[grep("exvivo",colnames(combined_mye))],y=combined_mye[grep("invivo",colnames(combined_mye))],
                method = corr)
      all_correlation<-pheatmap(test,cluster_rows = F,cluster_cols = F)
      
      temp_hm_name <- out("pheatmap.png")
    
      save_pheatmap(all_correlation, filename=temp_hm_name)
      
      out <- dirout(paste0(basedir,"exvivo_",ex,"invivo",inv))
      all_correlation
      ggsave(out("Correlations_all.png"))
      
      
      
      
      diag_correlation<-data.frame(test[row(test)==col(test)])
      colnames(diag_correlation)<-"correlation"
      rownames(diag_correlation)<-rownames(test)
      
      mine.heatmap <- ggplot(data = diag_correlation, mapping = aes(x = "correlation",
         y = rownames(diag_correlation),
         fill = correlation)) + geom_tile() + xlab(label = "Sample")+
         scale_fill_distiller(palette = "RdBu", limits = c(-1,1), na.value = "gray",
         direction = 1, labels = c(-1,-0.75, -0.5, 0.75, 1))
      
      mine.heatmap
      out <- dirout(paste0(basedir,"exvivo_",ex,"invivo",inv))
      ggsave(out("pheatmap_corresponding_ko.png"), w=15,h=15)
      
      
      
      out <- dirout(paste0(basedir,"exvivo_",ex,"invivo",inv))
      ggsave(out("pheatmap.png"), w=15,h=15)
    }}
}

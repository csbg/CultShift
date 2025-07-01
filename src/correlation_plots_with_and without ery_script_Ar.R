library(tidyselect)
library(pheatmap)

source("src/00_init.R")
basedir <- "SCRNA_33_DE_Nebula_testClustering/Output"
inp<-"SCRNA_33_DE_Nebula_testClustering/"
out <- dirout(basedir)
ex<-"9d"
inv<-"14d"
for (ex in c("7d","9d")){
  for (inv in c("14d","28d")){
    exvivo = fread(dirout(paste0(inp, "ex.vivo_",ex,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
    invivo = fread(dirout(paste0(inp, "in.vivo_",inv,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
    
    out1 <- dirout(paste0(basedir,"/exvivo_",ex,"invivo",inv))
    
    combined<-merge(exvivo,invivo,by=c("guide","gene_id")) %>% pivot_wider(
      names_from = guide,
      values_from=c(estimate.x,estimate.y)
    )
    colnames(combined)[grep("estimate.x",colnames(combined))]<-gsub("estimate.x","exvivo",colnames(combined)[grep("estimate.x",colnames(combined))])
    colnames(combined)[grep("estimate.y",colnames(combined))]<-gsub("estimate.y","invivo",colnames(combined)[grep("estimate.y",colnames(combined))])
    
    corr<-"pearson"
    for (corr in c("pearson","spearman")){
    
      test<-cor(combined[grep("exvivo",colnames(combined))],y=combined[grep("invivo",colnames(combined))],
                method = corr)
      head(test)
      png(out1("all_correlations.png"))
      all_correlation<-pheatmap(test,cluster_rows = F,cluster_cols = F)
      dev.off()
      
            
      diag_correlation<-data.frame(test[row(test)==col(test)])
      colnames(diag_correlation)<-"correlation"
      rownames(diag_correlation)<-gsub("exvivo_","",rownames(test))
            
      write.tsv(data.table(data.frame(diag_correlation), keep.rownames = TRUE),
                out1(paste0("correlation_logfc_invivo_",inv,"vs_exvivo",ex,".tsv")))
      
      mine.heatmap <- ggplot(data = diag_correlation, mapping = aes(x = "correlation",
         y = rownames(diag_correlation),
         fill = correlation)) + geom_tile() + xlab(label = "Sample")
      mine.heatmap
     
      ggsave(out1("pheatmap_corresponding_ko","invivo",inv,"exvivo",ex,".png"), w=15,h=15)
            
      
    }


   
##############################
#
############################

invivo_Nik<-fread("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_33_DE_Nebula_testClustering/ex.vivo_9d_useClusters/DEG_Results_all.tsv")
exvivo_Nik<-fread("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_33_DE_Nebula_testClustering/in.vivo_14d_useClusters/DEG_Results_all.tsv")

exvivo_Nik<-exvivo_Nik[exvivo_Nik$guide %in% invivo_Nik$guide & exvivo_Nik$gene_id %in% invivo_Nik$gene_id & exvivo_Nik$term %in% invivo_Nik$term,]
invivo_Nik<-invivo_Nik[invivo_Nik$guide %in% exvivo_Nik$guide & invivo_Nik$gene_id %in% exvivo_Nik$gene_id & invivo_Nik$term %in% exvivo_Nik$term,]


exvivo_Nik<-invivo_Nik[,c("guide","gene_id","estimate")]
invivo_Nik<-exvivo_Nik[,c("guide","gene_id","estimate")]

exvivo_Nik_wide <- exvivo_Nik %>% pivot_wider(names_from = guide,values_from = c(estimate)) %>% column_to_rownames(var="gene_id")
invivo_Nik_wide<-invivo_Nik %>% pivot_wider(names_from = guide,values_from = c(estimate)) %>% column_to_rownames(var="gene_id")


test_pearson_nik<-cor(exvivo_Nik_wide,invivo_Nik_wide,method = "pearson")
test_spearman_nik<-cor(exvivo_Nik_wide,invivo_Nik_wide,method = "spearman")

pheatmap(test_pearson_nik,cluster_rows = F,cluster_cols = F)
pheatmap(test_spearman_nik,cluster_rows = F,cluster_cols = F)

corrplot::corrplot(cor(x=exvivo_Nik_wide,y=invivo_Nik_wide,method = "pearson"),diag = T)
corrplot::corrplot(cor(x=exvivo_Nik_wide,y=invivo_Nik_wide,method = "spearman"),diag = T)

diag_pearson_nik<-data.frame(test_pearson_nik[row(test_pearson_nik)==col(test_pearson_nik)])
colnames(diag_pearson_nik)<-"pearson_corr"
rownames(diag_pearson_nik)<-rownames(test_pearson_nik)

mine.heatmap.nik <- ggplot(data = diag_pearson_nik, mapping = aes(x = "pearson_corr",
                                                                  y = rownames(diag_pearson_nik),
                                                                  fill = pearson_corr)) + geom_tile() + xlab(label = "Sample")+scale_fill_distiller(palette = "RdBu", limits = c(-1,1), na.value = "gray",
                                                                                                                                                    direction = 1, labels = c(-1,-0.75, -0.5, 0.75, 1))
mine.heatmap.nik
getwd()
out <- dirout(paste0(basedir, "compare_ko_exvivo_vs_invivo_nik"))
ggsave(out("9d_exvivo_vs_14d_invivo.png"), w=15,h=15)



plot(den, frame = FALSE, col = "blue",main = "Density plot")
out <- dirout(paste0(basedir, "compare_ko_exvivo_vs_invivo_cell_type"))
ggsave(out("9d_exvivo_vs_14d_invivo.png"), w=15,h=15)
###############################################

#density plots

###############################################

den <- density(diag_pearson$pearson_corr)

plot(den, frame = FALSE, col = "blue",main = "Density plot")
den_nik<- density(diag_pearson_nik$pearson_corr)
plot(den_nik, frame = FALSE, col = "blue",main = "Density plot")

density_compare<- cbind(diag_pearson,diag_pearson_nik$pearson_corr)
head(density_compare)
colnames(density_compare)<-c("mye","all")
density_compare<-density_compare %>% pivot_longer(cols = c(mye:all),
                                                  names_to = "cell_types",
                                                  values_to = "pearson_corr",
                                                  values_transform = list(pearson_corr=as.numeric))
head(density_compare)

#create overlaying density plots
ggplot(density_compare, aes(x=pearson_corr, fill=cell_types)) + xlim(-1,1)+
  geom_density(alpha=.25)
getwd()
out <- dirout(paste0(basedir, "TESTTcompare_ko_exvivo_vs_invivo_cell_type"))
ggsave(out("TEST9d_exvivo_vs_14d_invivo.png"), w=15,h=15)
###############

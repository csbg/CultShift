source("src/00_init.R")

library(ggpubr)
library("tidyverse")

basedir <- "SCRNA_33_DE_Nebula_testClustering/"
outB <- dirout(basedir)
ex<-"9d"

for ex in c("7d","9d"){
  for inv in c("14d","28d"){
    exvivo <- fread(dirout(paste0(basedir, "ex.vivo_",ex,"_useClusters"))("DEG_Results_all.tsv"))
    invivo <- fread(dirout(paste0(basedir, "in.vivo",inv,"_useClusters"))("DEG_Results_all.tsv"))
    colnames(invivo)<-paste0(colnames(invivo),"invivo")
    
    ######################################################
    #Filter based on common columns
    exvivo<-exvivo[exvivo$guide %in% invivo$guide & exvivo$gene_id %in% invivo$gene_id & exvivo$term %in% invivo$term,]
    invivo<-invivo[invivo$guide %in% exvivo$guide & invivo$gene_id %in% exvivo$gene_id & invivo$term %in% exvivo$term,]

    for tissue in c("exvivo","invivo"){
    wi<-tissue[,c("guide","gene_id","estimate")]%>% pivot_wider(names_from = guide,values_from = c(estimate)) %>% column_to_rownames(var="gene_id")
    
    test_pearson<-cor(wide,invivo_wide,method = "pearson")
test_spearman<-cor(exvivo_wide,invivo_wide,method = "spearman")

pheatmap(test_pearson,cluster_rows = F,cluster_cols = F)
pheatmap(test_spearman,cluster_rows = F,cluster_cols = F)


diag_pearson<-data.frame(test_pearson[row(test_pearson)==col(test_pearson)])
colnames(diag_pearson)<-"pearson_corr"
rownames(diag_pearson)<-rownames(test_pearson)


################################################################

mine.heatmap <- ggplot(data = diag_pearson, mapping = aes(x = "pearson_corr",
                                                       y = rownames(diag_pearson),
                                                       fill = pearson_corr)) + geom_tile() + xlab(label = "Sample")+scale_fill_distiller(palette = "RdBu", limits = c(-1,1), na.value = "gray",
                                                                                                                                         direction = 1, labels = c(-1,-0.75, -0.5, 0.75, 1))
  
png(mine.heatmap,""
getwd()
out <- dirout(paste0(basedir, "compare_ko_exvivo_vs_invivo"))
ggsave(out("7d_exvivo_vs_28d_invivo.png"), w=15,h=15)
corrplot::corrplot(cor(x=exvivo_wide,y=invivo_wide,method = "pearson"),diag = T)
corrplot::corrplot(cor(x=exvivo_wide,y=invivo_wide,method = "spearman"),diag = T)
#############################################################################

invivo_Nik<-fread("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_33_DE_Nebula_testClustering/ex.vivo_9d_useClusters/DEG_Results_all.tsv")
exvivo_Nik<-fread("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_33_DE_Nebula_testClustering/in.vivo_14d_useClusters/DEG_Results_all.tsv")

exvivo_Nik<-exvivo_Nik[exvivo_Nik$guide %in% invivo_Nik$guide & exvivo_Nik$gene_id %in% invivo_Nik$gene_id & exvivo_Nik$term %in% invivo_Nik$term,]
invivo_Nik<-invivo_Nik[invivo_Nik$guide %in% exvivo_Nik$guide & invivo_Nik$gene_id %in% exvivo_Nik$gene_id & invivo_Nik$term %in% exvivo_Nik$term,]


exvivo_Nik<-exvivo_Nik[,c("guide","gene_id","estimate")]
invivo_Nik<-invivo_Nik[,c("guide","gene_id","estimate")]

exvivo_Nik_wide <- exvivo_Nik %>% pivot_wider(names_from = guide,values_from = c(estimate)) %>% column_to_rownames(var="gene_id")
invivo_Nik_wide<-invivo_Nik %>% pivot_wider(names_from = guide,values_from = c(estimate)) %>% column_to_rownames(var="gene_id")


test_pearson_nik<-cor(exvivo_Nik_wide,invivo_Nik_wide,method = "pearson")
head(test_pearson_nik)
head(test_pearson)
test_spearman_nik<-cor(exvivo_Nik_wide,invivo_Nik_wide,method = "spearman")
pheatmap(test,cluster_rows = F,cluster_cols = F)
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

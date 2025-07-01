source("src/00_init.R")

library(ggpubr)
library("tidyverse")

basedir <- "SCRNA_33_DE_Nebula_testClustering/"
outB <- dirout(basedir)
for file in...invivo, do correlate with each in ex vivo
file1 <- fread(dirout(paste0(basedir, "ex.vivo_7d_useClusters"))("DEG_Results_all.tsv"))

file2<- fread(dirout(paste0(basedir, "in.vivo_14d_useClusters"))("DEG_Results_all.tsv"))
table(unique(file2$gene_id) %in% unique(file1$gene_id))
table(unique(file2$guide) %in% unique(file1$guide))
table(unique(file2$term) %in% unique(file1$term))
#Filter based on common columns
file2<-file2[file2$guide %in% file1$guide & file2$gene_id %in% file1$gene_id & file2$term %in% file1$term,]
file1<-file1[file1$guide %in% file2$guide & file1$gene_id %in% file2$gene_id & file1$term %in% file2$term,]
head(file1)


correlation <- cor(file1$estimate, file2$estimate, method = "pearson")
?cor



plot(file1$estimate, file2$estimate)

cor_ko_file1 <- file1 %>% pivot_wider(names_from = guide,values_from = estimate)
head(cor_ko_file1)
exvivo<-file1[,c("guide","gene_id","estimate")]
invivo<-file2[,c("guide","gene_id","estimate")]
exvivo_wide <- exvivo %>% pivot_wider(names_from = guide,values_from = c(estimate))
invivo_wide<-invivo %>% pivot_wider(names_from = guide,values_from = c(estimate))
exvivo_wide<-exvivo_wide %>% column_to_rownames(var="gene_id")
invivo_wide<-invivo_wide %>% column_to_rownames(var="gene_id")
test_pearson<-cor(exvivo_wide,invivo_wide,method = "pearson")
test_spearman<-cor(exvivo_wide,invivo_wide,method = "spearman")
pheatmap(test_pearson,cluster_rows = F,cluster_cols = F)
pheatmap(test_spearman,cluster_rows = F,cluster_cols = F)




#create correlation heatma
corrplot::corrplot(cor(x=exvivo_wide,y=invivo_wide,method = "pearson"),diag = T)
corrplot::corrplot(cor(x=exvivo_wide,y=invivo_wide,method = "spearman"),diag = T)
#############################################################################

file1_Nik<-fread("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_33_DE_Nebula_testClustering/ex.vivo_7d_useClusters/DEG_Results_all.tsv")
file2_Nik<-fread("/media/AGFORTELNY/PROJECTS/TfCf/Analysis/SCRNA_33_DE_Nebula_testClustering/in.vivo_14d_useClusters/DEG_Results_all.tsv")

file2_Nik<-file2_Nik[file2_Nik$guide %in% file1_Nik$guide & file2_Nik$gene_id %in% file1_Nik$gene_id & file2_Nik$term %in% file1_Nik$term,]
file1_Nik<-file1_Nik[file1_Nik$guide %in% file2_Nik$guide & file1_Nik$gene_id %in% file2_Nik$gene_id & file1_Nik$term %in% file2_Nik$term,]


exvivo_Nik<-file1_Nik[,c("guide","gene_id","estimate")]
invivo_Nik<-file2_Nik[,c("guide","gene_id","estimate")]
exvivo_Nik_wide <- exvivo_Nik %>% pivot_wider(names_from = guide,values_from = c(estimate))
invivo_Nik_wide<-invivo_Nik %>% pivot_wider(names_from = guide,values_from = c(estimate))
exvivo_Nik_wide<-exvivo_Nik_wide %>% column_to_rownames(var="gene_id")
invivo_Nik_wide<-invivo_Nik_wide %>% column_to_rownames(var="gene_id")
test_pearson<-cor(exvivo_Nik_wide,invivo_Nik_wide,method = "pearson")
test_spearman<-cor(exvivo_Nik_wide,invivo_Nik_wide,method = "spearman")
pheatmap(test_pearson,cluster_rows = F,cluster_cols = F)
pheatmap(test_spearman,cluster_rows = F,cluster_cols = F)
corrplot::corrplot(cor(x=exvivo_Nik_wide,y=invivo_Nik_wide,method = "pearson"),diag = T)
corrplot::corrplot(cor(x=exvivo_Nik_wide,y=invivo_Nik_wide,method = "spearman"),diag = T)

library(tidyselect)
library(pheatmap)
library(ggrepel)

source("src/00_init.R")
basedir <- "SCRNA_33_DE_Nebula_testClustering/Output"
inp<-"SCRNA_33_DE_Nebula_testClustering/"
out <- dirout(basedir)
ex<-"9d"
inv<-"14d"
leuk<-"6d"
#for (ex in c("7d","9d")){
#  for (inv in c("14d","28d")){
    exvivo = fread(dirout(paste0(inp, "ex.vivo_",ex,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
    head(exvivo)
    invivo = fread(dirout(paste0(inp, "in.vivo_",inv,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
    leukemia=fread(dirout(paste0(inp, "leukemia_",leuk,"_useClusters"))("DEG_Results_all.tsv")) %>% 
      select(c("guide","gene_id","estimate"))
    
    out1 <- dirout(paste0(basedir,"/invivo_",inv,"_exvivo_",ex))
    out2<- dirout(paste0(basedir,"/exvivo_",ex,"_leukemia_",leuk))
    
    combined_in_ex<-merge(exvivo,invivo,by=c("guide","gene_id")) %>% pivot_wider(
      names_from = guide,
      values_from=c(estimate.x,estimate.y))
    head(combined_in_ex)
    combined_ex_le<-merge(exvivo,leukemia,by=c("guide","gene_id")) %>% pivot_wider(
      names_from = guide,
      values_from=c(estimate.x,estimate.y))
    
    
    colnames(combined_in_ex)[grep("estimate.x",colnames(combined_in_ex))]<-gsub("estimate.x","exvivo",colnames(combined_in_ex)[grep("estimate.x",colnames(combined_in_ex))])
    colnames(combined_in_ex)[grep("estimate.y",colnames(combined_in_ex))]<-gsub("estimate.y","invivo",colnames(combined_in_ex)[grep("estimate.y",colnames(combined_in_ex))])
    
    colnames(combined_ex_le)[grep("estimate.x",colnames(combined_ex_le))]<-gsub("estimate.x","ex.vivo",colnames(combined_ex_le)[grep("estimate.x",colnames(combined_ex_le))])
    colnames(combined_ex_le)[grep("estimate.y",colnames(combined_ex_le))]<-gsub("estimate.y","leuk",colnames(combined_ex_le)[grep("estimate.y",colnames(combined_ex_le))])
    
    
    for (corr in c("pearson","spearman")){
      corr<-"pearson"
      test<-cor(combined_in_ex[grep("exvivo",colnames(combined_in_ex))],y=combined_in_ex[grep("invivo",colnames(combined_in_ex))],
                method = corr)
      head(test)
      png(out1(paste0("all_correlations_ex_in",corr,"_.png")))
      all_correlation<-pheatmap(test,cluster_rows = F,cluster_cols = F)
      dev.off()
      
            
      diag_correlation_in_ex<-data.frame(test[row(test)==col(test)])
      colnames(diag_correlation_in_ex)<-"correlation"
      rownames(diag_correlation_in_ex)<-gsub("invivo_","",colnames(test))
      diag_correlation_in_ex["gene"]<-rownames(diag_correlation_in_ex)
      
      
      diag_correlation_in_ex<-diag_correlation_in_ex%>% filter(gene %in% res_ex$gene)%>% arrange(desc(correlation))
      head(diag_correlation_in_ex)
      t<-diag_correlation_in_ex$gene
      t
      ggplot(diag_correlation_in_ex, aes(x=reorder(gene,correlation), y=correlation)) +
        geom_bar(stat='identity') + ylim(-0.1,1)+
        coord_flip()
      ggsave(out1("correlation_exvivo_invivo_logfc.png"), w=20,h=10)
      head(pDT.ntc_sum_med_diff)
      t<-reorder(diag_correlation_in_ex$gene,correlation)
      pDT.ntc_sum_med_diff<-pDT.ntc_sum_med_diff[match(t,pDT.ntc_sum_med_diff$gene)]
      head(pDT.ntc_sum_med_diff)
      ggplot(pDT.ntc_sum_med_diff, aes(x=ex.vivo, y=gene)) + 
        geom_point(size=1, shape=21, fill="white")+
        geom_segment(data=pDT.ntc_sum_med_diff, mapping=aes(x=ex.vivo, y=gene), 
                                                            xend=ex.vivo+diff, yend=gene),arrow = arrow(length = unit(0.09, "cm")), 
                     size=0.09, color="black")+
        theme_bw(12)
        #geom_text_repel(aes(label=gene,max.overlaps =40))
      
      t<-diag_correlation_in_ex$gene
      
      source("Ar02_SCRNA_50_01_02_Trajectories_comparison_only_plot.R")
      write.tsv(data.table(data.frame(diag_correlation), keep.rownames = TRUE),
                out1(paste0("correlation_logfc_leuk_",leuk,"vs_invivo",inv,corr,".tsv")))
      
      mine.heatmap <- ggplot(data = diag_correlation, mapping = aes(x = "correlation",
         y = rownames(diag_correlation),
         fill = correlation)) + geom_tile() + xlab(label = "Sample")
      mine.heatmap
     
      ggsave(out1("pheatmap_corresponding_ko","leuk",leuk,"invivo",inv,".png"), w=15,h=15)
            
      
    }
head(diag)

##################
###7
out2 <- dirout(paste0(basedir,"/exvivo_",ex,"leukemia",leuk))

combined_ex_le<-merge(exvivo,leukemia,by=c("guide","gene_id")) %>% pivot_wider(
  names_from = guide,
  values_from=c(estimate.x,estimate.y)
)
head(combined_ex_le)
colnames(combined_ex_le)[grep("estimate.x",colnames(combined_ex_le))]<-gsub("estimate.x","exvivo",colnames(combined_ex_le)[grep("estimate.x",colnames(combined_ex_le))])
colnames(combined_ex_le)[grep("estimate.y",colnames(combined_ex_le))]<-gsub("estimate.y","leuk",colnames(combined_ex_le)[grep("estimate.y",colnames(combined_ex_le))])
head(exvivo)
########


test2<-cor(combined_ex_le[grep("exvivo",colnames(combined_ex_le))],y=combined_ex_le[grep("leuk",colnames(combined_ex_le))],
          method = corr)

all_correlation1<-pheatmap(test2,cluster_rows = F,cluster_cols = F)
dev.off()


diag_correlation_2<-data.frame(test2[row(test2)==col(test2)])
colnames(diag_correlation_2)<-"correlation"
rownames(diag_correlation_2)<-gsub("exvivo_","",rownames(test2))

write.tsv(data.table(data.frame(diag_correlation_2), keep.rownames = TRUE),
          out1(paste0("correlation_logfc_leuk_",leuk,"vs_invivo",inv,corr,".tsv")))

mine.heatmap <- ggplot(data = diag_correlation_2, mapping = aes(x = "correlation",
                                                              y = rownames(diag_correlation),
                                                              fill = correlation)) + geom_tile() + xlab(label = "Sample")
mine.heatmap

ggsave(out1("pheatmap_corresponding_ko","leuk",leuk,"invivo",inv,".png"), w=15,h=15)

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
diag_correlation["gene"]<-rownames(diag_correlation)
head(diag_correlation_2)
diag_correlation<-diag_correlation[gene %in% unique(diag_correlation_2$gene)]
diag_correlation_2["gene"]<-rownames(diag_correlation_2)
compare<-merge(diag_correlation,diag_correlation_2,by="gene")
head(compare)
ggplot(compare,aes(x=correlation.x,y=correlation.y))+
  geom_point()+theme_bw(12) +
  geom_text_repel(aes(label=gene),max.overlaps = 20)
ggsave(out1("invivo_to_leuk_vs_exvivo_to_leuk.png"),w=8,h=8)


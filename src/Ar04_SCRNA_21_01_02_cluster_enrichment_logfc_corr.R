source("src/00_init.R")
require(ggrepel)


basedir_corr <- "SCRNA_33_DE_Nebula_testClustering/Output/exvivo_9dinvivo14d/"
out_corr <- dirout(basedir_corr)
outs<-dirout("SCRNA_50_02_ProjectionTrajectories/Corr_enrichment_vs_logfc")

correlation_logfc<-fread(out_corr("correlation_logfc_invivo_14dvs_exvivo9d.tsv"))
colnames(correlation_logfc)<-c("gene","correlation")
head(correlation_logfc)


ex.vivo_9d<- fread(dirout_load(paste0("SCRNA_50_02_ProjectionTrajectories/cluster.enrichments"))("Cluster_enrichments_basic_ex.vivo_withMixscape_9d.tsv"))
in.vivo_14d<- fread(dirout_load(paste0("SCRNA_50_02_ProjectionTrajectories/cluster.enrichments"))("Cluster_enrichments_basic_in.vivo_withMixscape_14d.tsv"))

in.vivo_14d<-in.vivo_14d %>%filter(Clusters %in% unique(ex.vivo_9d$Clusters)) %>% filter(gene %in% unique(ex.vivo_9d$gene)) %>% filter(Clusters !="GMP (late)")
ex.vivo_9d<-ex.vivo_9d %>%filter(Clusters %in% unique(in.vivo_14d$Clusters)) %>% filter(gene %in% unique(in.vivo_14d$gene)) %>% filter(Clusters !="GMP (late)")

in.vivo_14d_arrange<- in.vivo_14d %>% select(gene,Clusters,log2OR) %>%
  group_by(gene, Clusters)
in.vivo_14d_arrange<-as.data.table(in.vivo_14d_arrange)
in.vivo_14d<- in.vivo_14d_arrange[with(in.vivo_14d_arrange,order(gene,Clusters))]
head(in.vivo_14d)

ex.vivo_9d_arrange<- ex.vivo_9d %>% select(gene,Clusters,log2OR) %>% 
  group_by(gene, Clusters)
ex.vivo_9d_arrange<-as.data.table(ex.vivo_9d_arrange)
ex.vivo_9d<- ex.vivo_9d_arrange[with(ex.vivo_9d_arrange,order(gene,Clusters))]

head(ex.vivo_9d)
head(in.vivo_14d)
combined<-merge(ex.vivo_9d,in.vivo_14d,by=c("gene","Clusters")) %>% pivot_wider(
  names_from = gene,
  values_from=c(log2OR.x,log2OR.y)
)
head(combined)

colnames(combined)[grep("log2OR.x",colnames(combined))]<-gsub("log2OR.x","exvivo",colnames(combined)[grep("log2OR.x",colnames(combined))])
colnames(combined)[grep("log2OR.y",colnames(combined))]<-gsub("log2OR.y","invivo",colnames(combined)[grep("log2OR.y",colnames(combined))])
unique(ex.vivo_9d$Clusters
)
unique(in.vivo_14d$Clusters)

for (corr in c("pearson","spearman")){  
  test<-cor(combined[grep("exvivo",colnames(combined))],y=combined[grep("invivo",colnames(combined))],
            method = "spearman")
  head(test)
  pdf(out("Basic_pheatmap_across_all","exvivo9d","exvivo14d",corr,"basic.pdf"))  
  all_correlation<-pheatmap(test,cluster_rows = F,cluster_cols = F)
  dev.off()
  all_correlation
  
  
  
  diag_correlation<-data.frame(test[row(test)==col(test)])
  colnames(diag_correlation)<-"correlation"
  rownames(diag_correlation)<-gsub("exvivo_","",rownames(test))
  
  
  
  mine.heatmap <- ggplot(data = diag_correlation, mapping = aes(x = "correlation",
                                                                y = rownames(diag_correlation),
                                                                fill = correlation)) + geom_tile(aes(fill = correlation)) + xlab(label = "Sample")+
    ylab(label = "guide")
  mine.heatmap
  
  ggsave(out("Basic_pheatmap_corresponding_celltype","exvivo9d","exvivo14d",corr,"basic.png"), w=15,h=15)
  
  #scatter cell type enrichment vs logfc
  diag_correlation["gene"]<-rownames(diag_correlation)
  
 
  correlation_logfc<-fread(out_corr("correlation_logfc_invivo_14dvs_exvivo9d.tsv"))
  colnames(correlation_logfc)<-c("gene","correlation")
  corr_enrichment_vs_logfc<-merge(diag_correlation,correlation_logfc,by="gene")
  colnames(corr_enrichment_vs_logfc)<-c("gene","enrichment","logfc")
  
  ggplot(corr_enrichment_vs_logfc, aes(x=enrichment, y=logfc)) + 
    geom_point() + 
    theme_bw(12) +
    geom_text_repel(aes(label=gene))
  ggsave(outs(paste0("Scatter_enrichment",corr,"_vs_logfc",".pdf")),w=8,h=8)
}




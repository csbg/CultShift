
source("src/00_init.R")

base.dir <- "SCRNA_50_02_ProjectionTrajectories/"

# . Cluster enrichment analyses ---------------------------------------------
tx <- "in.vivo"
inDir <- dirout_load("SCRNA_21_01_ClusterEnrichments")
#for(tx in names(inDir.funcs)){
# out directory
out <- dirout(paste0(base.dir, "/", "cluster.enrichments_21_01/"))
out1<-dirout(paste0(base.dir, "/", "cluster_enrichment_correlation_invivo_vs_exvivo"))
# Collect enrichment scores
typex<-"broad_in.vivo_withMixscape"
gsub("Guides_Fisher_Mixscape_(.+).pdf", "\\1", list.files(inDir(""), pattern="Guides_Fisher_Mixscape_.*.pdf"))
for(typex in gsub("Guides_Fisher_Mixscape_(.+).pdf", "\\1", list.files(inDir(""), pattern="Guides_Fisher_Mixscape_.*.pdf"))){
  fish.file <- inDir("Guides_Fisher_Mixscape_",typex,".tsv")
  fish.file
  if(!file.exists(fish.file)) next
  
  fish.full <- fread(fish.file)
  
  
  fish.full[mixscape_class == "Pu.1", mixscape_class := "Spi1"]
  
  grep("S", unique(fish.full$mixscape_class), value = TRUE)
  #fish.full <- merge(fish.full, unique(SANN[,c("sample_broad", "timepoint"),with=F]), by.x="sample", by.y="sample_broad")
  timex <- "14d"
  
  for(timex in c(unique(fish.full$sample))){
    fish <- copy(fish.full)
    typex
    timex
    #fish <- fish[mixscape_class == "Rbbp4"]
    if(timex != "all") fish <- fish[sample == timex]
    
    # summarize across NTCs
    fish <- fish[, .(
      log2OR=mean(log2OR), 
      dir=length(unique(sign(log2OR[padj < 0.01]))) <= 1, 
      #dir=length(unique(sign(log2OR)))==1, 
      padj=sum(padj < 0.01),
      N=.N,
      Ncells=sum(unique(guide.cells))), by=c("sample", "Clusters", "mixscape_class")]
    fish[dir == FALSE, padj := 0]
    fish[dir == FALSE, log2OR := NA]
    
    # legacy
    fish[, gene := mixscape_class]
    
    # Summarize across guides
     fish[, gene := gsub("_.+", "", mixscape_class)]
    # fish[gene == "Pu.1", gene := "Spi1"]
    # fish <- fish[, .(
    #   log2OR=mean(log2OR, na.rm=TRUE), 
    #   dir=length(unique(sign(log2OR[!is.na(log2OR)])))==1, 
    #   padj=sum(padj), 
    #   N=sum(N)), by=c("sample", "Clusters", "gene")]
    # fish[dir == FALSE, padj := 0]
    # fish[dir == FALSE, log2OR := NA]
    
    # summarize across samples
    # fish <- fish[, .(
    #   log2OR=mean(log2OR, na.rm=TRUE), 
    #   dir=length(unique(sign(log2OR[!is.na(log2OR)])))==1, 
    #   padj=sum(padj), 
    #   N=sum(N)), by=c("Clusters", "gene")]
    # fish[dir == FALSE, padj := 0]
    # fish[dir == FALSE, log2OR := NA]
    head(fish)
    #
    head(ex.vivo_9d)
    fish<-fish[gene %in% unique(ex.vivo_9d$gene)]
    # setup for plotting
    fish[padj == 0 | is.na(log2OR), log2OR := 0]
    fish[, sig.perc := padj / N]
    fish[,log2OR_cap := pmin(abs(log2OR), 5) * sign(log2OR)]
    fish <- hierarch.ordering(dt = fish, toOrder = "gene", orderBy = "Clusters", value.var = "log2OR")
    head(fish)
    fish[, Clusters := gsub("^Gran$", "Gran.", Clusters)]
    #fish[, Clusters := cleanCelltypes(Clusters)]
    #fish <- hierarch.ordering(dt = fish, toOrder = "Clusters", orderBy = "gene", value.var = "log2OR")
    ggplot(fish, aes(x=gene, y=Clusters, size=sig.perc, color=log2OR_cap)) + 
      themeNF(rotate=TRUE) +
      scale_color_gradient2(name="log2(OR)",low="blue", midpoint = 0, high="red") +
      scale_size_continuous(name="% sign.", range = c(0,5)) +
      geom_point() +
      geom_point(shape=1, color="lightgrey") +
      xlab("Gene") + ylab("Cell type")
    ggsaveNF(
      out("Cluster_enrichments_",typex,"_", timex, "Filtered.pdf"), 
      w=length(unique(fish$gene))*0.05 + 0.5,
      h=length(unique(fish$Clusters))*0.05 + 0.5)
    write.tsv(fish, out("Cluster_enrichments_Filtered",typex,"_", timex,".tsv"))
  }
}
#}

ex.vivo_9d<- fread(dirout_load(paste0("SCRNA_50_02_ProjectionTrajectories/cluster.enrichments"))("Cluster_enrichments_broad_ex.vivo_withMixscape_9d.tsv"))
in.vivo_14d<- fread(dirout_load(paste0("SCRNA_50_02_ProjectionTrajectories/cluster.enrichments"))("Cluster_enrichments_broad_in.vivo_withMixscape_14d.tsv"))

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


combined<-merge(ex.vivo_9d,in.vivo_14d,by=c("gene","Clusters")) %>% pivot_wider(
  names_from = gene,
  values_from=c(log2OR.x,log2OR.y)
)
head(combined)

colnames(combined)[grep("log2OR.x",colnames(combined))]<-gsub("log2OR.x","exvivo",colnames(combined)[grep("log2OR.x",colnames(combined))])
colnames(combined)[grep("log2OR.y",colnames(combined))]<-gsub("log2OR.y","invivo",colnames(combined)[grep("log2OR.y",colnames(combined))])

for (corr in c("pearson","spearman")){  
  test<-cor(combined[grep("exvivo",colnames(combined))],y=combined[grep("invivo",colnames(combined))],
            method = "spearman")
  
  colnames(test)<-gsub("invivo_","",colnames(test))
  rownames(test)<-gsub("exvivo_","",rownames(test))
  head(test)
  test<- replace(test,is.na(t(test)), 0))
?pheatmap
pdf(out1("Broad_pheatmap_across_all","exvivo9d","exvivo14d",corr,".pdf")) 
  all_correlation<-pheatmap(test,cluster_rows = T,cluster_cols = T,main = "invivo_vs_exvivo")
  dev.off()
 
  
  
  
  diag_correlation<-data.frame(test[row(test)==col(test)])
  colnames(diag_correlation)<-"correlation"
  rownames(diag_correlation)<-gsub("exvivo_","",rownames(test))
  
  
  
  mine.heatmap <- ggplot(data = diag_correlation, mapping = aes(x = "correlation",
         y = rownames(diag_correlation),
         fill = correlation)) + geom_tile(aes(fill = correlation)) + xlab(label = "Sample")+
                                                     ylab(label = "guide")
      mine.heatmap
     
      ggsave(out1("Broad_pheatmap_corresponding_celltype","exvivo9d","exvivo14d",corr,".png"), w=15,h=15)
            
      
}
